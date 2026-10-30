#!/usr/bin/env python3
"""
40_generate_legend.py: Generate a companion JSON legend file for the GLB output.

Reads:  config.yaml
        outputs/<organ>/json/<stem>_distribution.json
        outputs/<organ>/csv/<stem>-<uberon>-all-as-<tool>-<shape>-hra-pop.csv
Writes: outputs/<organ>/json/<stem>-<uberon>-all-as-<tool>-<shape>-hra-pop_legend.json
"""

from __future__ import annotations

import csv
import datetime
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).parent))
from shared import (
    load_config,
    ORGAN_CELL_COUNTS,
    resolve_glb_filename,
    load_supertree,
    resolve_organ_filter_label,
    build_top9_granular_color_scheme,
)
from palette import LEVEL2_TOP9_COLORS, OTHER_COLOR, BG_MODE, apply_to_scheme

REAL_CELL_SIZE_UM = 10.0

def build_output_stem(stem: str, uberon_id: str, tool_slug: str, shape: str) -> str:
    parts = [stem]
    if uberon_id:
        parts.append(uberon_id)
    parts += ["all-as", tool_slug, shape, "hra-pop"]
    return "-".join(parts)


def load_distribution(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def load_nodes_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as f:
        return [dict(row) for row in csv.DictReader(f)]


def get_plotted_counts_by_as(rows: List[Dict[str, str]]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for row in rows:
        as_label = row.get("as_label", "").strip()
        counts[as_label] = counts.get(as_label, 0) + 1
    return counts


def main() -> None:
    config = load_config(Path(__file__).parent.parent / "config.yaml")

    organ_name  = config["organ"]["name"]
    organ_sex   = config["organ"]["sex"]
    tool        = config["filters"]["tool"]
    total_nodes_requested = config["markers"]["num_nodes"]
    sphere_radius = config["markers"]["sphere_radius"]
    shape         = config["markers"]["shape"]
    supertree_path = config["apis"]["cell_hierarchy"]

    glb_filename = (config["organ"].get("glb_filename") or "").strip()
    organ_side   = (config["organ"].get("side") or "").strip()
    if not glb_filename:
        glb_filename = resolve_glb_filename(organ_name, organ_sex, organ_side)

    stem      = Path(glb_filename).stem
    tool_slug = tool.lower()

    base_folder = Path(config["output"]["amended_folder"]) / organ_name.lower()
    json_folder = base_folder / "json"
    csv_folder  = base_folder / "csv"
    json_folder.mkdir(parents=True, exist_ok=True)

    distribution_path = json_folder / f"{stem}_distribution.json"
    if not distribution_path.exists():
        raise FileNotFoundError(f"Distribution JSON not found: {distribution_path}.")

    dist_data  = load_distribution(distribution_path)
    uberon_id  = dist_data.get("uberon_id", "")
    output_name = build_output_stem(stem, uberon_id, tool_slug, shape)

    nodes_csv_path = csv_folder  / f"{output_name}.csv"
    output_json    = json_folder / f"{output_name}_legend.json"

    if not nodes_csv_path.exists():
        raise FileNotFoundError(f"Nodes CSV not found: {nodes_csv_path}.")

    nodes_rows     = load_nodes_csv(nodes_csv_path)
    plotted_by_as  = get_plotted_counts_by_as(nodes_rows)
    organ_info = ORGAN_CELL_COUNTS.get(organ_name.lower(), {})
    organ_total_cells = organ_info.get("total_cells", 1)

    # Same inputs + palette as 30/50, so legend_entries matches theirs.
    print(f"Loading supertree from: {supertree_path}")
    supertree, _ = load_supertree(supertree_path)
    _color_map, legend_entries = build_top9_granular_color_scheme(
        organ_cell_counts=dist_data.get("organ_cell_counts", {}),
        cell_id_map=dist_data.get("cell_id_map", {}),
        supertree=supertree,
        palette=LEVEL2_TOP9_COLORS,
        other_color=OTHER_COLOR,
        bg_mode=BG_MODE,
    )
    # Cell type colors come from data/Supertree_colors_all_organs_A5.xlsx.
    _color_map, legend_entries = apply_to_scheme(
        _color_map, legend_entries,
        resolve_organ_filter_label(organ_name, (config["organ"].get("side") or "").strip()), organ_sex, tool,
    )

    total_plotted = sum(plotted_by_as.values())

    legend_json = {
        "metadata": {
            "organ": organ_name, "sex": organ_sex, "tool": tool,
            "uberon_id": uberon_id, "glb_filename": glb_filename,
            "generated_at": datetime.datetime.utcnow().isoformat() + "Z",
        },
        "organ_total_cell_count": {
            "organ":       organ_name,
            "total_cells": organ_info.get("total_cells", "unknown"),
            "notes":       organ_info.get("notes", ""),
            "source":      organ_info.get("source", ""),
            "doi":         organ_info.get("doi", ""),
        },
        "cell_type_color_key": {
            "cell type": [label for label, _hex, _lvl2, _count in legend_entries],
            "colors": {
                label: hex_color for label, hex_color, _lvl2, _count in legend_entries
            },
        },
        "anatomical_structures": [
            {
                "name":              as_label,
                "glb_node":          dist_data["as_to_glb_node"].get(as_label, ""),
                "hra_pop_cell_count": int(dist_data["as_total_counts"].get(as_label, 0)),
                "allocated_markers": dist_data["node_allocation"].get(as_label, 0),
                "plotted":           plotted_by_as.get(as_label, 0),
            }
            for as_label in dist_data["as_labels"]
        ],
        "cell_count_comparison": {
            "total_plotted":       total_plotted,
            "total_requested":     total_nodes_requested,
            "real_cells_per_marker": organ_total_cells // max(total_plotted, 1),
        },
        "size_reference": {
            "assumed_cell_size_um": REAL_CELL_SIZE_UM,
            # GLB/glTF linear units are meters per spec, and sphere_radius
            "marker_radius_um":   round(sphere_radius * 1e6, 2),
            "marker_diameter_um": round(sphere_radius * 1e6 * 2, 2),
            "marker_size_ratio":  f"{sphere_radius * 1e6 * 2:.0f}μm|{REAL_CELL_SIZE_UM:.0f}μm",
            "marker_scale_factor": round((sphere_radius * 1e6 * 2) / REAL_CELL_SIZE_UM, 1),
            "note": (
                "All cell types assumed to be 10 µm in diameter. Markers are "
                "rendered far larger than true biological scale for VR "
                "visibility, marker_scale_factor is how many times larger "
                "than a real cell each marker is (marker_diameter_um / "
                "assumed_cell_size_um)."
            ),
        },
    }

    output_json.write_text(
        json.dumps(legend_json, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Legend JSON saved to: {output_json}")


if __name__ == "__main__":
    main()