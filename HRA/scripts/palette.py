#!/usr/bin/env python3
"""
palette.py: HRA's colors for the top-9-granular + Other scheme.

The 10 fixed AS/2 category colors (LEVEL2_TOP9_COLORS) are
shared.CT_COLORS, the same fixed name->hex lookup A5 uses, so a given
supertree category (e.g. "muscle cell") means the same color, in the
same order-independent way, in both projects. "endothelial cell" is
always #e31a1c, never reassigned based on which order categories happen
to appear in for a given organ/sex/tool block. HRA previously had its
own distinct bright palette; that's gone in favor of one shared identity
across the whole atlas.

What HRA supplies locally is OTHER_COLOR and BG_MODE, which tells
shared.build_top9_granular_color_scheme() to shade those same hues for a
dark VR background (#0F1117) instead: BG_MODE="dark" shifts the whole
shading lightness range brighter (see shared.BG_LIGHTNESS_RANGES) so no
shade gets dark enough to disappear into near-black, while keeping the
exact same hue/saturation as A5's version of that category.

Imported by 30, 40, 50, 60, all four pass these into
build_top9_granular_color_scheme() so GLB markers, the CSV hex_color
column, the JSON legend, and the color spreadsheet all match.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from shared import CT_COLORS

LEVEL2_TOP9_COLORS: Dict[str, str] = CT_COLORS

# Matches the "Other" color in data/Supertree_colors_all_organs_A5.xlsx; also
# used for cell types with no supertree match.
OTHER_COLOR: str = "#b15928"

# Passed as bg_mode= into shared.build_top9_granular_color_scheme().
BG_MODE: str = "dark"


# --- Cell type colors from data/Supertree_colors_all_organs_A5.xlsx ---------
# For each organ / sex / tool the sheet lists the top 9 cell types and their hex
# colors. apply_to_scheme() swaps those into the computed color scheme; cell
# types with a blank color (no supertree match) get OTHER_COLOR.

SHEET_PATH = Path(__file__).parent.parent.parent / "data" / "Supertree_colors_all_organs_A5.xlsx"

# (organ label lower, sex lower, tool lower) -> {cell type: hex or ""}
_CACHE: Optional[Dict[Tuple[str, str, str], Dict[str, str]]] = None


def _load() -> Dict[Tuple[str, str, str], Dict[str, str]]:
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    import openpyxl  # imported here so the rest of the pipeline works without it

    if not SHEET_PATH.exists():
        raise FileNotFoundError(f"Color sheet not found: {SHEET_PATH}")
    ws = openpyxl.load_workbook(SHEET_PATH, read_only=True, data_only=True).active
    table: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    organ = sex = tool = ""
    for row in ws.iter_rows(min_row=2, values_only=True):
        row = tuple(row) + (None,) * (9 - len(row))  # read-only mode can return short rows
        # Organ / sex / tool are only filled on the first row of each block.
        organ = str(row[0]).strip() if row[0] else organ
        sex = str(row[1]).strip() if row[1] else sex
        tool = str(row[2]).strip() if row[2] else tool
        cell_type = row[6]
        if not cell_type or not organ:
            continue
        hex_color = str(row[7]).strip() if row[7] else ""
        table.setdefault((organ.lower(), sex.lower(), tool.lower()), {})[str(cell_type).strip()] = hex_color
    _CACHE = table
    return table


def sheet_colors_for(organ_label: str, sex: str, tool: str) -> Dict[str, str]:
    return _load().get((organ_label.strip().lower(), sex.strip().lower(), tool.strip().lower()), {})


def apply_to_scheme(
    color_map: Dict[str, str],
    legend_entries: List[Tuple[str, str, str, float]],
    organ_label: str,
    sex: str,
    tool: str,
) -> Tuple[Dict[str, str], List[Tuple[str, str, str, float]]]:
    """organ_label is the HRApop label, e.g. 'Left kidney' (shared.resolve_organ_filter_label)."""
    sheet = sheet_colors_for(organ_label, sex, tool)
    if not sheet:
        print(f"  Warning: no rows in the color sheet for {organ_label} / {sex} / {tool}; using computed colors.")

    def pick(label: str, computed: str) -> str:
        if label == "Other":
            return OTHER_COLOR
        if label in sheet:
            return sheet[label] or OTHER_COLOR
        return computed or OTHER_COLOR

    new_map = dict(color_map)
    for label in list(new_map):
        if label in sheet:
            new_map[label] = pick(label, new_map[label])
        elif not new_map[label]:
            new_map[label] = OTHER_COLOR
    new_entries = [(label, pick(label, hex_color), lvl2, count) for label, hex_color, lvl2, count in legend_entries]
    return new_map, new_entries
