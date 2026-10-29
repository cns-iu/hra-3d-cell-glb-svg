#!/usr/bin/env python3
"""
run_all_combos.py: Mass-generate HRA outputs for every organ / sex / tool
combination in HRApop that has a reference-organ GLB.

Always generates the distribution JSON (stages 10 + 20), since every other
output depends on it. Everything else is opt-in:

    --glbs   annotated GLB + node CSV (stage 30; the two come from the same
             marker-placement run, so they are always written together)
    --pngs   PNG charts (stage 50)
    --svgs   SVG charts (stage 51)
    --all    all of the above

config.yaml is never modified: the config is overridden in memory per combo
(10,000 markers). HRApop, the supertree, and the reference-organ list are
fetched once, live, and reused for every combination. Top cell types with no
supertree match fall back to OTHER_COLOR instead of crashing stage 30.

Writes under outputs/annotated organs/<organ>/ (output_name =
<stem>-<uberon>-all-as-<tool>-sphere-hra-pop):
    json/<output_name>_distribution.json
    glb/<output_name>.glb            (--glbs)
    csv/<output_name>.csv            (--glbs)
    png/{with,without}-legend/...    (--pngs)
    svg/{with,without}-legend/...    (--svgs)
(Stage 20 also leaves a scratch json/<stem>_distribution.json holding
whichever combo ran last; the per-tool file above is the real output.)

Usage (this file can live in HRA/ or HRA/scripts/):
    python run_all_combos.py --dry-run --all      # show what would run
    python run_all_combos.py --svgs               # JSON + SVGs
    python run_all_combos.py --all                # everything
    python run_all_combos.py --glbs --only heart,kidney
    python run_all_combos.py --pngs --force       # regenerate even if present
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import importlib.util
import io
import json
import os
import shutil
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Tuple

HRA_ROOT = Path(__file__).resolve().parent
if not (HRA_ROOT / "config.yaml").exists():  # placed inside HRA/scripts/
    HRA_ROOT = HRA_ROOT.parent
SCRIPTS_DIR = HRA_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import requests
import trimesh

import shared
from palette import OTHER_COLOR

NUM_NODES = 10_000
MAX_TRIANGLES = 1_000_000

VERBOSE = False
_REAL_LOAD_CONFIG = shared.load_config


# =============================================================================
# Helpers
# =============================================================================

def load_stage(name: str, cfg: dict):
    """Import a stage script with load_config() returning `cfg` (in memory)."""
    shared.load_config = lambda _path: cfg
    try:
        spec = importlib.util.spec_from_file_location(name, SCRIPTS_DIR / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        shared.load_config = _REAL_LOAD_CONFIG
    return module


def parse_organ_label(label: str) -> Tuple[str, str]:
    """'Left kidney' -> ('kidney', 'left'); 'heart' -> ('heart', '')."""
    low = label.strip()
    for side in ("Left", "Right"):
        if low.startswith(side + " "):
            return low[len(side) + 1:], side.lower()
    return low, ""


def fetch_organs(api_url: str) -> Dict[Tuple[str, str], str]:
    """(organ label lower, sex lower) -> GLB URL, from the live reference-organs API."""
    resp = requests.get(api_url, timeout=120)
    resp.raise_for_status()
    out: Dict[Tuple[str, str], str] = {}
    for o in resp.json():
        url = str(o.get("object", {}).get("file", "")).strip()
        if url:
            out[(str(o.get("label", "")).strip().lower(), str(o.get("sex", "")).strip().lower())] = url
    return out


def fetch_rows(base_cfg: dict, stage20):
    """Live HRApop rows; falls back to the cached CSV only with a loud warning."""
    try:
        rows, raw = stage20.fetch_hra_pop_rows(base_cfg["apis"]["hra_pop"])
        print(f"HRApop: {len(rows)} rows fetched LIVE")
        return rows, raw
    except Exception as exc:
        cache = HRA_ROOT.parent / "data" / "hra_pop_raw.csv"
        if not cache.exists():
            raise
        print(f"WARNING: live HRApop fetch failed ({exc}).")
        print(f"WARNING: using cached CSV {cache}, data may be out of date.")
        raw = cache.read_bytes()
        return stage20.rows_from_csv_bytes(raw), raw


def fetch_supertree(url: str):
    return shared.load_supertree(url)


# =============================================================================
# Plan
# =============================================================================

def build_plan(rows, organs, base_cfg, stage20, only: List[str]):
    combos = sorted({(r["organ"], r["sex"], r["tool"]) for r in rows})
    plan, skipped = [], []
    shape = base_cfg["markers"]["shape"]
    for label, sex, tool in combos:
        name, side = parse_organ_label(label)
        if only and name.lower() not in only and label.lower() not in only:
            continue
        url = organs.get((label.lower(), sex.lower()))
        if not url:
            skipped.append((label, sex, tool, "no reference-organ GLB"))
            continue
        glb_filename = url.rsplit("/", 1)[-1]
        stem = Path(glb_filename).stem
        uberon = stage20.extract_organ_uberon_id(rows, label, tool, sex)
        output_name = "-".join(p for p in [stem, uberon, "all-as", tool.lower(), shape, "hra-pop"] if p)
        base = Path(base_cfg["output"]["amended_folder"]) / name.lower()
        plan.append(dict(
            label=label, name=name, side=side, sex=sex, tool=tool,
            glb_filename=glb_filename, url=url, stem=stem, output_name=output_name,
            base=base,
            glb_path=base / "glb" / f"{output_name}.glb",
            csv_path=base / "csv" / f"{output_name}.csv",
            generic_json=base / "json" / f"{stem}_distribution.json",
            tool_json=base / "json" / f"{output_name}_distribution.json",
        ))
    return plan, skipped


def json_done(item) -> bool:
    """Per-tool JSON present. A JSON from an earlier run (the generic file for
    this same tool) is adopted as complete."""
    if item["tool_json"].exists():
        return True
    generic = item["generic_json"]
    if generic.exists():
        try:
            if json.loads(generic.read_text(encoding="utf-8")).get("tool", "").lower() == item["tool"].lower():
                shutil.copyfile(generic, item["tool_json"])
                return True
        except Exception:
            pass
    return False


def glb_done(item) -> bool:
    return (item["glb_path"].exists() and item["glb_path"].stat().st_size > 0
            and item["csv_path"].exists() and item["csv_path"].stat().st_size > 0)


def charts_done(item, kind: str) -> bool:
    folder = item["base"] / kind
    return all((folder / sub / f"{item['output_name']}-{sub}.{kind}").exists() for sub in ("with-legend", "without-legend"))


# =============================================================================
# Steps for one combination
# =============================================================================

def make_cfg(base_cfg: dict, item: dict) -> dict:
    cfg = copy.deepcopy(base_cfg)
    cfg["organ"].update(name=item["name"], sex=item["sex"], side=item["side"], glb_filename=item["glb_filename"])
    cfg["filters"]["tool"] = item["tool"]
    cfg["markers"]["num_nodes"] = NUM_NODES
    cfg["output"]["save_hra_pop_csv"] = False
    return cfg


def step_json(item, cfg, rows, raw, supertree_pair, organs_folder: Path) -> None:
    """Stages 10 + 20: reference GLB download and distribution JSON."""
    s10 = load_stage("10_download_organs", cfg)
    organs_folder.mkdir(parents=True, exist_ok=True)
    s10.download_glb(item["url"], organs_folder)

    s20 = load_stage("20_fetch_cell_distribution", cfg)
    s20.fetch_hra_pop_rows = lambda _url: (rows, raw)
    s20.load_supertree = lambda _path: supertree_pair
    s20.main()
    shutil.copyfile(item["generic_json"], item["tool_json"])


def use_tool_json(item) -> None:
    """Stages 30/50/51 read the generic <stem>_distribution.json; point it at this tool's data."""
    shutil.copyfile(item["tool_json"], item["generic_json"])


def step_glb(item, cfg, supertree_pair) -> int:
    """Stage 30: annotated GLB + node CSV from the same marker-placement run."""
    use_tool_json(item)
    s30 = load_stage("30_generate_markers", cfg)
    s30.load_supertree = lambda _path: supertree_pair
    s30.resolve_reference_organ_glb_url = lambda _fn, _api: item["url"]
    original_build = s30.build_color_map
    s30.build_color_map = lambda *a, **k: {
        label: (hex_color or OTHER_COLOR) for label, hex_color in original_build(*a, **k).items()
    }
    s30.main()

    scene = trimesh.load(item["glb_path"], force="scene")
    return sum(len(g.faces) for g in scene.geometry.values() if isinstance(g, trimesh.Trimesh))


def step_chart(item, cfg, supertree_pair, stage_name: str) -> None:
    """Stage 50 (PNG) or 51 (SVG) from this combo's JSON."""
    use_tool_json(item)
    stage = load_stage(stage_name, cfg)
    stage.load_supertree = lambda _path: supertree_pair
    stage.main()


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--glbs", action="store_true", help="annotated GLB + node CSV (stage 30)")
    ap.add_argument("--pngs", action="store_true", help="PNG charts (stage 50)")
    ap.add_argument("--svgs", action="store_true", help="SVG charts (stage 51)")
    ap.add_argument("--all", action="store_true", help="GLBs/CSVs, PNGs and SVGs")
    ap.add_argument("--dry-run", action="store_true", help="list what would run and exit")
    ap.add_argument("--force", action="store_true", help="regenerate even if the output exists")
    ap.add_argument("--only", default="", help="comma-separated organ names, e.g. heart,kidney")
    ap.add_argument("--verbose", action="store_true", help="show the stage output instead of hiding it")
    args = ap.parse_args()
    if args.all:
        args.glbs = args.pngs = args.svgs = True
    global VERBOSE
    VERBOSE = args.verbose
    only = [s.strip().lower() for s in args.only.split(",") if s.strip()]

    os.chdir(HRA_ROOT)
    base_cfg = _REAL_LOAD_CONFIG(HRA_ROOT / "config.yaml")
    stage20 = load_stage("20_fetch_cell_distribution", base_cfg)

    organs = fetch_organs(base_cfg["apis"]["reference_organs"])
    rows, raw = fetch_rows(base_cfg, stage20)
    plan, skipped = build_plan(rows, organs, base_cfg, stage20, only)

    # What each combination still needs.
    jobs = []
    for item in plan:
        steps = []
        if args.force or not json_done(item):
            steps.append("json")
        if args.glbs and (args.force or not glb_done(item)):
            steps.append("glb+csv")
        if args.pngs and (args.force or not charts_done(item, "png")):
            steps.append("png")
        if args.svgs and (args.force or not charts_done(item, "svg")):
            steps.append("svg")
        if steps:
            jobs.append((item, steps))

    print(f"\n{len(plan)} combinations with a GLB; {len(jobs)} need work; {len(skipped)} skipped (no GLB). "
          f"Markers per organ: {NUM_NODES:,}.")
    for label, sex, tool, why in skipped:
        print(f"  skip: {label} / {sex} / {tool}, {why}")
    if args.dry_run:
        for item, steps in jobs:
            print(f"  would run: {item['label']} / {item['sex']} / {item['tool']}: {', '.join(steps)}")
        return
    if not jobs:
        print("Nothing to do.")
        return

    supertree_pair = fetch_supertree(base_cfg["apis"]["cell_hierarchy"])
    organs_folder = Path(base_cfg["output"]["organs_folder"])

    ok, failed, heavy = 0, [], []
    for i, (item, steps) in enumerate(jobs, 1):
        tag = f"{item['label']} / {item['sex']} / {item['tool']}"
        print(f"[{i}/{len(jobs)}] {tag} ({', '.join(steps)}) ...", end=" ", flush=True)
        t0 = time.time()
        buf = io.StringIO()
        cfg = make_cfg(base_cfg, item)
        tris = None
        try:
            with (contextlib.nullcontext() if VERBOSE else contextlib.redirect_stdout(buf)):
                if "json" in steps:
                    step_json(item, cfg, rows, raw, supertree_pair, organs_folder)
                if "glb+csv" in steps:
                    tris = step_glb(item, cfg, supertree_pair)
                if "png" in steps:
                    step_chart(item, cfg, supertree_pair, "50_generate_pngs")
                if "svg" in steps:
                    step_chart(item, cfg, supertree_pair, "51_generate_svgs")
        except Exception as exc:
            print(f"FAILED ({exc})")
            failed.append((tag, "".join(traceback.format_exception_only(type(exc), exc)).strip(), buf.getvalue()[-600:]))
            continue
        note = ""
        if tris is not None:
            note = f", {tris:,} triangles" + ("  <-- OVER 1M TRIANGLES" if tris > MAX_TRIANGLES else "")
            if tris > MAX_TRIANGLES:
                heavy.append((tag, tris))
        print(f"ok{note}, {time.time() - t0:.0f}s")
        ok += 1

    print(f"\nDone: {ok} ok, {len(failed)} failed" + (f", {len(heavy)} over {MAX_TRIANGLES:,} triangles." if args.glbs else "."))
    for tag, tris in heavy:
        print(f"  heavy: {tag}, {tris:,}")
    for tag, err, tail in failed:
        print(f"  failed: {tag}, {err}")
        if tail.strip():
            print("    last output: " + tail.strip().replace("\n", "\n    "))


if __name__ == "__main__":
    main()
