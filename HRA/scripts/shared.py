#!/usr/bin/env python3
"""
shared.py — Constants and utilities shared across projects (A5, HRA, ...).

Mostly project-agnostic: no fonts, no layout opinions live here.
build_top9_granular_color_scheme() takes a palette argument so each
project supplies its own OTHER_COLOR/BG_MODE (see each project's
palette.py) while sharing the same ranking/grouping/shading algorithm.
The one exception is CT_COLORS below — the 10 base colors are meant to
be identical across every project (so "endothelial cell" is the same
color family in A5 and HRA), so that dict lives here rather than being
duplicated per project.
"""

from __future__ import annotations

import colorsys
import csv
import re
import yaml
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd


# =============================================================================
# Hex color helpers
# =============================================================================

# Fixed AS/2 supertree label -> base hex color. Every project must use
# these exact colors for these exact labels — "endothelial cell" is
# always #e31a1c everywhere in the atlas, never reassigned by whichever
# order categories happen to appear in for a given organ/sex/tool block
# (that positional-first-appearance scheme is what BASE_HUES used to do,
# and it's why the same AS/2 category could get a different color in
# different blocks). Granular cell types under a category are shaded
# variants of that category's fixed color (see _shades_for_group) and
# must never be assigned the raw base color itself.
#
# Source: CT Colors.xlsx (ColorBrewer "Paired", n=10:
# https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10). Keyed
# on the current 10-label AS/2 supertree (supertree.csv, Aug 2026
# revision) — if the supertree's AS/2 label set changes again, this dict
# needs new entries, not a re-indexed list.
CT_COLORS: Dict[str, str] = {
    "connective tissue cell": "#1f78b4",
    "endothelial cell":       "#e31a1c",
    "epithelial cell":        "#ff7f00",
    "germ line cell":         "#6a3d9a",
    "hematopoietic cell":     "#fb9a99",
    "mesenchymal cell":       "#a6cee3",
    "neural cell":            "#33a02c",
    "melanocyte":             "#cab2d6",
    "muscle cell":            "#afabab",
    "precursor cell":         "#b2df8a",
}


def hex_to_rgb(hex_color: str) -> Tuple[int, int, int]:
    h = hex_color.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def hex_to_float(hex_color: str) -> Tuple[float, float, float]:
    r, g, b = hex_to_rgb(hex_color)
    return r / 255.0, g / 255.0, b / 255.0


def _hex_to_hls(hex_color: str) -> Tuple[float, float, float]:
    r, g, b = hex_to_float(hex_color)
    return colorsys.rgb_to_hls(r, g, b)


def _hls_to_hex(h: float, l: float, s: float) -> str:
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return "#{:02x}{:02x}{:02x}".format(
        int(round(r * 255)), int(round(g * 255)), int(round(b * 255))
    )


MIN_LIGHTNESS_DELTA = 0.05  # minimum lightness gap a shade must keep from
                            # its group's unshaded base color (see below)

# Both projects share one canonical set of 9 base hues (see each project's
# palette.py) so "muscle cell" etc. means the same color family everywhere
# in the atlas. What differs per project is how those hues get shaded into
# per-cell-type variants: A5 renders on a light background and can use
# genuinely dark shades for contrast, HRA renders in VR on a near-black
# background, where a shade dark enough to read well on paper would all but
# disappear. bg_mode shifts the whole lightness range brighter for "dark"
# without changing the hue/saturation, so it's still recognizably the same
# base color family, just legible on either background.
BG_LIGHTNESS_RANGES: Dict[str, Tuple[float, float]] = {
    "light": (0.20, 0.88),
    "dark":  (0.45, 0.92),
}


def _shades_for_group(base_hex: str, n: int, bg_mode: str = "light") -> List[str]:
    """n evenly-spaced HSL-lightness shades of base_hex, darkest first,
    within the lightness floor/ceiling for bg_mode ("light" or "dark").

    Guaranteed to never return base_hex itself — a granular cell type
    sharing its supertree parent's exact color would be visually
    indistinguishable from "just show the category color", defeating the
    point of a granular-level legend. Two cases previously let this slip
    through: a group with only one top-9 member (nothing to shade against,
    so the base was returned unchanged), and a group with an odd member
    count, where symmetric lightness spacing put the middle-ranked member
    exactly back on the base lightness by construction. Both are now
    nudged at least MIN_LIGHTNESS_DELTA away from the base.
    """
    h, l0, s = _hex_to_hls(base_hex)
    floor, ceiling = BG_LIGHTNESS_RANGES.get(bg_mode, BG_LIGHTNESS_RANGES["light"])

    def _clamp(l: float) -> float:
        return min(ceiling, max(floor, l))

    def _push_off_base(l: float) -> float:
        # Keep l on whichever side of l0 it already landed on, but widen
        # the gap to at least MIN_LIGHTNESS_DELTA if it's too close. (An
        # earlier version of this pushed *toward* l0 instead of away from
        # it — l - DELTA when l >= l0 moves a lighter shade down toward
        # the base rather than further from it — which could round back
        # to base_hex's exact rgb after hex quantization.)
        if l >= l0 and l - l0 < MIN_LIGHTNESS_DELTA:
            l = l0 + MIN_LIGHTNESS_DELTA
        elif l < l0 and l0 - l < MIN_LIGHTNESS_DELTA:
            l = l0 - MIN_LIGHTNESS_DELTA
        return _clamp(l)

    if n <= 1:
        midpoint = (floor + ceiling) / 2
        l = l0 + MIN_LIGHTNESS_DELTA if l0 <= midpoint else l0 - MIN_LIGHTNESS_DELTA
        return [_hls_to_hex(h, _clamp(l), s)]

    lo, hi = max(floor, l0 - 0.28), min(ceiling, l0 + 0.28)
    if hi <= lo:
        lo, hi = floor + 0.05, ceiling - 0.05
    step = (hi - lo) / (n - 1)
    ls = [_push_off_base(lo + i * step) for i in range(n)]

    # The push-off-base nudge above can shove a shade into its neighbor
    # (rare — needs a large n bunched near a clamp boundary, e.g. every
    # top-9 slot in one supertree category). Enforce a minimum gap between
    # consecutive shades too, sweeping darkest-to-lightest.
    MIN_STEP = 0.015
    for i in range(1, len(ls)):
        if ls[i] - ls[i - 1] < MIN_STEP:
            ls[i] = min(ceiling, ls[i - 1] + MIN_STEP)

    return [_hls_to_hex(h, l, s) for l in ls]


# =============================================================================
# Top-9-granular + Other color scheme
# =============================================================================
# Rank granular (level-3) cell types organ-wide by count. The top 9 each
# belong to a level-2 parent category; each parent's color is a fixed
# lookup in `palette` (CT_COLORS) by category name — not by first-
# appearance order, so a given category is the same color in every organ/
# sex/tool block, not just within one. Members under the same parent get
# HSL-lightness shades of that color (darkest = highest count).
# Everything outside the top 9 folds into one "Other" bucket (`other_color`).
# A top-9 member whose supertree lookup finds no AS/2 category at all
# (lvl2_of() returns "") is neither shaded nor folded into Other — it's
# returned with hex="" so the caller can flag it for review instead of
# silently guessing a color or silently dropping it from the legend.

# Known spelling-variant duplicates in the raw HRApop cell_label field,
# folded into one canonical label before ranking (e.g. a typo'd duplicate
# of "Ventricular cardiomyocyte" in the source data).
CELL_LABEL_ALIASES: Dict[str, str] = {
    "ventricular cardiomycoyte": "Ventricular cardiomyocyte",
}


def canonicalize_cell_label(label: str) -> str:
    return CELL_LABEL_ALIASES.get(label.strip().lower(), label)


def build_top9_granular_color_scheme(
    organ_cell_counts: Dict[str, float],
    cell_id_map: Dict[str, str],
    supertree: Dict[str, str],
    palette: Dict[str, str],
    other_color: str,
    bg_mode: str = "light",
) -> Tuple[Dict[str, str], List[Tuple[str, str, str, float]]]:
    """
    Build the top-9-granular + Other color scheme for one organ.

    Args:
        organ_cell_counts: cell_label -> organ-wide raw cell count.
        cell_id_map:       normalized cell_label key -> compact CL ID.
        supertree:         compact CL ID -> level-2 supertree label.
        palette:           AS/2 category name -> fixed hex color (pass
                            CT_COLORS). Every project must pass the same
                            dict so a given category means the same color
                            family everywhere — use bg_mode to adapt
                            legibility instead of swapping in different
                            hues. A category with no entry in `palette`
                            falls back to `other_color`.
        other_color:       color for the aggregate "Other" bucket, and
                            the fallback for any top-9 category missing
                            from `palette`.
        bg_mode:           "light" or "dark" — see BG_LIGHTNESS_RANGES;
                            controls how far shades range from each base
                            hue so they stay legible on that background.

    Returns:
        color_map: cell_label -> hex color (unmatched top-9 labels map to "").
        legend_entries: ordered (label, hex, level2_label, count) — the
            top-9 granular cell types + one "Other" entry, sorted by
            count descending with Other pinned last. A top-9 label with
            no supertree match at all gets hex="" and level2_label="" —
            distinct from "Other", which is deliberately not one of
            these entries but the aggregate of everything ranked below
            the top 9.
    """
    merged_counts: Dict[str, float] = defaultdict(float)
    for label, count in organ_cell_counts.items():
        merged_counts[canonicalize_cell_label(label)] += count
    organ_cell_counts = dict(merged_counts)

    ranked = sorted(organ_cell_counts.items(), key=lambda kv: kv[1], reverse=True)
    top9 = ranked[:9]
    rest = ranked[9:]

    def lvl2_of(label: str) -> str:
        # cell_id_map may only have an entry under the original spelling,
        # so also try any alias that folds into this canonical label.
        keys_to_try = [normalize_cell_label_key(label)]
        for alias, canonical in CELL_LABEL_ALIASES.items():
            if canonical == label:
                keys_to_try.append(normalize_cell_label_key(alias))
        for key in keys_to_try:
            cl_id = cell_id_map.get(key, "")
            if cl_id:
                return supertree.get(cl_id, "")
        return ""

    ordered_lvl2: List[str] = []
    for label, _ in top9:
        lvl2 = lvl2_of(label)
        if lvl2 and lvl2 not in ordered_lvl2:
            ordered_lvl2.append(lvl2)

    lvl2_color = {
        lvl2: palette.get(lvl2, other_color)
        for lvl2 in ordered_lvl2
    }

    groups: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    for label, count in top9:
        groups[lvl2_of(label)].append((label, count))

    color_map: Dict[str, str] = {}
    legend_entries: List[Tuple[str, str, str, float]] = []

    for lvl2 in ordered_lvl2:
        members = sorted(groups[lvl2], key=lambda x: x[1], reverse=True)
        shades = _shades_for_group(lvl2_color[lvl2], len(members), bg_mode=bg_mode)
        for (label, count), shade in zip(members, shades):
            color_map[label] = shade
            legend_entries.append((label, shade, lvl2, count))

    # Top-9 members whose supertree lookup found no AS/2 category at all
    # (lvl2_of() == "") used to be silently dropped here: groups[""] was
    # populated but "" is never in ordered_lvl2, so the loop above never
    # visited it. Surface them instead — hex="" is the signal to the
    # caller that this needs a human decision, not an automatic color.
    for label, count in sorted(groups.get("", []), key=lambda x: x[1], reverse=True):
        color_map[label] = ""
        legend_entries.append((label, "", "", count))

    other_total = 0.0
    for label, count in rest:
        color_map[label] = other_color
        other_total += count
    legend_entries.append(("Other", other_color, "", other_total))

    # Grouping by level-2 parent above can bury a higher-count member of a
    # different parent further down, so re-sort display order by count,
    # descending, with Other pinned last. Colors/shades are unaffected.
    other_entry = legend_entries.pop()
    legend_entries.sort(key=lambda entry: entry[3], reverse=True)
    legend_entries.append(other_entry)

    return color_map, legend_entries


# =============================================================================
# Supertree loading
# =============================================================================

def load_supertree(path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Load supertree CSV (local path or Google Sheets export URL).

    Returns:
        cell_to_label: CL:XXXXXXX -> supertree label (AS/2 or AS/3 fallback)
        cell_to_id:    CL:XXXXXXX -> supertree CL ID (AS/2 or AS/3 fallback)

    Filters to azimuth/celltypist/popv (+ a few other) sources only, and
    skips the root 'cell' row.
    """
    sources = {"azimuth", "celltypist", "popv", "pan-human-azimuth", "frmatch", "sc_proteomics"}
    cell_to_label: Dict[str, str] = {}
    cell_to_id:    Dict[str, str] = {}

    df = pd.read_csv(path, encoding="utf-8-sig")

    for _, row in df.iterrows():
        source = str(row.get("CT/1 - Sources", "")).strip().lower()
        if source not in sources:
            continue

        leaf_id = None
        for i in range(12, 0, -1):
            v = str(row.get(f"AS/{i}/ID", "")).strip()
            if v and v.lower() != "nan":
                leaf_id = compact_ontology_id(v)
                break
        if not leaf_id:
            continue

        as2_label = str(row.get("AS/2/LABEL", "")).strip()
        as2_id    = str(row.get("AS/2/ID",    "")).strip()
        as3_label = str(row.get("AS/3/LABEL", "")).strip()
        as3_id    = str(row.get("AS/3/ID",    "")).strip()

        as2_label = "" if as2_label.lower() == "nan" else as2_label
        as2_id    = "" if as2_id.lower()    == "nan" else as2_id
        as3_label = "" if as3_label.lower() == "nan" else as3_label
        as3_id    = "" if as3_id.lower()    == "nan" else as3_id

        if as2_label:
            level_label = as2_label
            level_id    = compact_ontology_id(as2_id) if as2_id else f"ASCTB-TEMP-{as2_label}"
        elif as3_label:
            level_label = as3_label
            level_id    = compact_ontology_id(as3_id) if as3_id else f"ASCTB-TEMP-{as3_label}"
        else:
            continue

        if leaf_id not in cell_to_label:
            cell_to_label[leaf_id] = level_label
            cell_to_id[leaf_id]    = level_id

    return cell_to_label, cell_to_id


# =============================================================================
# Cell type display names
# =============================================================================

CELL_TYPE_FULL_NAMES: Dict[str, str] = {
    "fibroblast":                       "Fibroblast",
    "endothelial":                      "Endothelial Cell",
    "cardiomyocyte":                    "Cardiomyocyte",
    "macrophage":                       "Macrophage",
    "smooth_muscle":                    "Smooth Muscle Cell",
    "pericyte":                         "Pericyte",
    "t_cell":                           "T Cell",
    "b_cell":                           "B Cell",
    "mast_cell":                        "Mast Cell",
    "neuronal":                         "Neuronal Cell",
    "adipocyte":                        "Adipocyte",
    "epithelial":                       "Epithelial Cell",
    "proximal_tubule":                  "Proximal Tubule Cell",
    "distal_convoluted_tubule":         "Distal Convoluted Tubule Cell",
    "loop_of_henle":                    "Loop of Henle Cell",
    "collecting_duct":                  "Collecting Duct Cell",
    "podocyte":                         "Podocyte",
    "mesangial":                        "Mesangial Cell",
}

# =============================================================================
# Organ total cell counts (for legend / JSON metadata)
# =============================================================================

ORGAN_CELL_COUNTS: Dict[str, Dict] = {
    "heart": {
        "total_cells":  7_500_000_000,
        "notes":        "The human heart contains an estimated 2-3 billion cardiac muscle cells, but these account for less than a third of the total cell number in the heart.",
        "source":       "Cell Communications in the Heart",
        "doi":          "10.1161/CIRCULATIONAHA.108.847731",
    },
    "kidney": {
        "total_cells":  1_500_000_000,
        "notes":        "Estimated 1.5 billion cells per kidney.",
        "source":       "Lote, 2012",
        "doi":          "",
    },
}

# =============================================================================
# Manual AS label -> GLB node overrides
# =============================================================================

AS_LABEL_OVERRIDES: Dict[str, List[str]] = {
    "posteromedial head of posterior papillary muscle of left ventricle": [
        "VH_F_papillary_muscle_of_heart_posmed",
        "VH_M_papillary_muscle_of_heart_posmed",
    ],
}

# =============================================================================
# GLB filename resolution
# =============================================================================

def resolve_glb_filename(organ_name: str, organ_sex: str, side: str = "") -> str:
    sex_code  = "f" if organ_sex.lower() == "female" else "m"
    name_slug = organ_name.lower().strip().replace(" ", "-")
    side      = (side or "").strip().lower()
    side_code = {"left": "-l", "right": "-r"}.get(side, "")
    return f"3d-vh-{sex_code}-{name_slug}{side_code}.glb"


def resolve_organ_filter_label(organ_name: str, side: str = "") -> str:
    """HRApop CSV rows label lateralized organs (kidney, ureter, ovary,
    knee, mammary gland, ...) as 'Left <organ>' / 'Right <organ>', not
    plain 'kidney'. Row-filtering must match on this full label — an
    exact match on organ_name alone would find zero rows, and a
    substring match would silently pull in both sides' data at once."""
    side = (side or "").strip().lower()
    if side in ("left", "right"):
        return f"{side.title()} {organ_name.strip()}"
    return organ_name.strip()


# =============================================================================
# Config loading
# =============================================================================

def load_config(config_path: Path) -> dict:
    with config_path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


# =============================================================================
# Text normalisation helpers
# =============================================================================

def normalize_match_text(text: str) -> str:
    text = text.lower()
    text = re.sub(r"[^a-z0-9\s]", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def normalize_cell_label_key(label: str) -> str:
    return re.sub(r"[\s\-]+", "_", label.strip().lower())


def safe_name(text: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_]", "_", text)


# =============================================================================
# Ontology / data helpers
# =============================================================================

def compact_ontology_id(url_or_id: str) -> str:
    """
    Compact a full ontology URL to a short CURIE-style ID, e.g.
    http://purl.obolibrary.org/obo/CL_0000057 -> CL:0000057
    """
    s = (url_or_id or "").strip()
    tail = s.rstrip("/").split("/")[-1]
    match = re.match(r"^([A-Za-z]+)[_:](\d+)$", tail)
    if match:
        return f"{match.group(1).upper()}:{match.group(2)}"
    return tail


def parse_float(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def values_match_filter(value: str, filter_value: str) -> bool:
    """Exact match (case-insensitive). A substring check here would let
    e.g. tool='azimuth' also match 'pan-human-azimuth' rows, silently
    doubling every row-sum downstream."""
    if not filter_value:
        return True
    return value.strip().lower() == filter_value.strip().lower()


def auto_match_as_to_glb_nodes(
    as_labels: List[str],
    glb_node_names: List[str],
    organ_name: str,
) -> Dict[str, Optional[List[str]]]:
    organ_slug = normalize_match_text(organ_name)

    norm_to_glb: Dict[str, str] = {}
    for node in glb_node_names:
        norm = normalize_match_text(node)
        for prefix in [f"vh_f_{organ_slug}_", f"vh_m_{organ_slug}_",
                       f"vh_f_", f"vh_m_", organ_slug + "_", organ_slug]:
            if norm.startswith(prefix):
                norm = norm[len(prefix):]
                break
        norm_to_glb[norm] = node

    result: Dict[str, Optional[List[str]]] = {}

    for as_label in as_labels:
        if as_label.lower() in AS_LABEL_OVERRIDES:
            result[as_label] = AS_LABEL_OVERRIDES[as_label.lower()]
            print(f"  ✓  '{as_label}' → {result[as_label]} (manual override)")
            continue

        norm_as = normalize_match_text(as_label)
        for prefix in [organ_slug + " ", organ_slug]:
            if norm_as.startswith(prefix):
                norm_as = norm_as[len(prefix):]
                break

        matches: List[str] = []

        for norm_node, original in norm_to_glb.items():
            if norm_as == norm_node:
                matches.append(original)

        if not matches:
            for norm_node, original in norm_to_glb.items():
                if norm_as in norm_node or norm_node in norm_as:
                    matches.append(original)

        if not matches:
            as_tokens = set(norm_as.split())
            for norm_node, original in norm_to_glb.items():
                node_tokens = set(norm_node.split())
                overlap = as_tokens & node_tokens
                if len(overlap) >= max(1, len(as_tokens) // 2):
                    matches.append(original)

        if matches:
            result[as_label] = matches
            print(f"  ✓  '{as_label}' → {matches}")
        else:
            result[as_label] = None
            print(f"  ✗  '{as_label}' → no match")

    return result