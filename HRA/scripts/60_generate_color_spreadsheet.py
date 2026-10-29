#!/usr/bin/env python3
"""
60_generate_color_spreadsheet.py: Build/refresh the cumulative Supertree
color reference spreadsheet for every organ/sex/tool combination this
project's pipeline has processed so far.

Unlike scripts 10-50, this stage isn't scoped to the single organ named in
config.yaml, it rescans every distribution JSON already on disk under
outputs/annotated organs/*/json/*_distribution.json (one per organ/sex/tool
run through 20_fetch_cell_distribution.py) and recomputes each one's
top-9-granular + Other color scheme with THIS project's palette, so the
spreadsheet always reflects exactly what 30/40/50 actually rendered for
every organ processed to date. Running the pipeline for a new organ and
re-running this stage grows the spreadsheet by one more (Organ, Sex, Tool)
block, nothing needs to be re-entered by hand.

Reads:  outputs/annotated organs/<organ>/json/<stem>_distribution.json  (one per organ/sex/tool run through script 20)
Writes: outputs/Supertree_colors_all_organs_HRA.xlsx
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

sys.path.insert(0, str(Path(__file__).parent))
from shared import load_config, load_supertree, build_top9_granular_color_scheme
from palette import LEVEL2_TOP9_COLORS, OTHER_COLOR, BG_MODE

PROJECT_LABEL = "HRA"
BG_MODE_LABEL = "dark-background" if BG_MODE == "dark" else "light-background"


# =============================================================================
# Gather one row-group per (Organ, Sex, Tool) distribution JSON on disk
# =============================================================================

def find_distribution_files(amended_folder: Path) -> List[Path]:
    if not amended_folder.exists():
        return []
    return sorted(amended_folder.glob("*/json/*_distribution.json"))


def build_rows_for_distribution(
    dist_path: Path,
    supertree: Dict[str, str],
) -> Tuple[List[Tuple], List[Tuple]]:
    dist_data = json.loads(dist_path.read_text(encoding="utf-8"))
    organ = dist_data.get("organ", "")
    sex   = dist_data.get("sex", "")
    tool  = dist_data.get("tool", "")
    organ_cell_counts = dist_data.get("organ_cell_counts", {})
    cell_id_map       = dist_data.get("cell_id_map", {})

    if not organ_cell_counts:
        print(f"  Skipping {dist_path.name}, no organ_cell_counts in file.")
        return [], []

    _color_map, legend_entries = build_top9_granular_color_scheme(
        organ_cell_counts=organ_cell_counts,
        cell_id_map=cell_id_map,
        supertree=supertree,
        palette=LEVEL2_TOP9_COLORS,
        other_color=OTHER_COLOR,
        bg_mode=BG_MODE,
    )

    # legend_entries is (label, hex, level2_label, count), count-descending
    # with "Other" pinned last. LEVEL2_TOP9_COLORS (CT_COLORS) is now a
    # fixed name -> hex dict, so Base Color is a direct lookup, no more
    # reconstructing "which slot of a positional palette did this category
    # land in for this block", which was the actual bug behind the same
    # category getting different colors in different organ/sex/tool
    # blocks. "Other" and a top-9 member with no supertree match both
    # have level2_label == "" and so both get a blank Base Color; they're
    # told apart by Hex Color, "Other" always has OTHER_COLOR, an
    # unmatched member has hex_color == "" (see shared.py) and gets
    # listed in `unmatched` below for the review tab instead of guessing.
    rows: List[Tuple] = []
    unmatched: List[Tuple] = []
    for rank, (label, hex_color, lvl2, count) in enumerate(legend_entries, start=1):
        base_color = LEVEL2_TOP9_COLORS.get(lvl2, "") if lvl2 else ""
        rows.append((organ, sex, tool, rank, label, lvl2, base_color, hex_color, count))
        if label != "Other" and not hex_color:
            unmatched.append((organ, sex, tool, rank, label, count))
    return rows, unmatched


# =============================================================================
# xlsx output, same styling as the manually-built reference spreadsheet
# =============================================================================

def write_xlsx(
    blocks: List[List[Tuple]],
    out_path: Path,
    unmatched: Optional[List[Tuple]] = None,
) -> None:
    """blocks: one list of row-tuples per (Organ, Sex, Tool) distribution
    file, in the order they should appear. Matches the reference
    spreadsheet's own layout exactly: no title/subtitle banner, row 1 is
    blank, row 2 is the header (bold, no fill), data starts row 3. A blank
    row separates each block, and the Organ column is merged per block
    (one merge per organ/sex/tool group, not spanning across multiple
    blocks of the same organ, every block repeats its own organ name)."""
    wb = Workbook()
    ws = wb.active
    ws.title = "Sheet1"

    thin   = Side(style="thin", color="FFA8A8A8")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    HEADER_FONT = Font(name="Arial", size=11, bold=True, color="FF000000")
    BODY_FONT   = Font(name="Arial", size=11, color="FF000000")
    CENTER      = Alignment(horizontal="center", vertical="center")

    # Column order: Organ, Sex, Tool, Rank, Supertree Label, Base Color,
    # Cell Type, Hex Color, Cell Count. Row 1 is intentionally left blank.
    headers = ["Organ", "Sex", "Tool", "Rank", "Supertree label", "Base Color",
               "Cell Type", "Hex Color", "Cell Count"]
    for c, h in enumerate(headers, start=1):
        cell = ws.cell(row=2, column=c, value=h)
        cell.font = HEADER_FONT
        cell.alignment = CENTER

    r = 3
    for block_idx, block_rows in enumerate(blocks):
        if not block_rows:
            continue

        if block_idx > 0:
            r += 1  # blank separator row between every block, left unstyled

        block_start = r
        for organ, sex, tool, rank, cell_type, lvl2, base_color, hex_color, count in block_rows:
            # Column 1 (Organ) is intentionally left blank here, it's
            # filled in once per merged range below.
            values = [None, sex, tool, rank, lvl2, base_color, cell_type, hex_color, count]
            for c, v in enumerate(values, start=1):
                cell = ws.cell(row=r, column=c, value=v)
                cell.font = BODY_FONT
                cell.border = border

            ws.cell(row=r, column=4).alignment = CENTER  # Rank, centered

            if base_color:
                ws.cell(row=r, column=6).fill = PatternFill("solid", fgColor="FF" + base_color.lstrip("#"))
            if hex_color:
                ws.cell(row=r, column=8).fill = PatternFill("solid", fgColor="FF" + hex_color.lstrip("#"))

            r += 1
        block_end = r - 1

        organ_name = block_rows[0][0]
        if block_end > block_start:
            ws.merge_cells(start_row=block_start, start_column=1, end_row=block_end, end_column=1)
        top_cell = ws.cell(row=block_start, column=1, value=organ_name)
        top_cell.font = BODY_FONT
        top_cell.border = border

    widths = {"A": 19.75, "B": 9.38, "C": 20.38, "D": 9.25, "E": 18.75,
              "F": 15.13, "G": 63.0, "H": 14.75, "I": 14.63}
    for col, w in widths.items():
        ws.column_dimensions[col].width = w
    ws.freeze_panes = "A3"

    if unmatched:
        ws2 = wb.create_sheet("Unclassified — Review")
        ws2_headers = ["Organ", "Sex", "Tool", "Rank", "Cell Type", "Cell Count"]
        for c, h in enumerate(ws2_headers, start=1):
            cell = ws2.cell(row=1, column=c, value=h)
            cell.font = HEADER_FONT
            cell.alignment = CENTER
        note = ws2.cell(
            row=2, column=1,
            value=(
                "These are top-9-by-count granular cell types whose supertree lookup found no "
                "AS/2 category at all under the current supertree, not folded into \"Other\" "
                "(which is a different, deliberate bucket for everything ranked below the top 9) "
                "and not auto-colored. They need either a supertree fix or a manual color call."
            ),
        )
        note.font = Font(name="Arial", size=10, italic=True, color="FF535353")
        ws2.merge_cells(start_row=2, start_column=1, end_row=2, end_column=len(ws2_headers))
        r2 = 3
        for organ, sex, tool, rank, label, count in unmatched:
            values = [organ, sex, tool, rank, label, count]
            for c, v in enumerate(values, start=1):
                cell = ws2.cell(row=r2, column=c, value=v)
                cell.font = BODY_FONT
            ws2.cell(row=r2, column=4).alignment = CENTER
            r2 += 1
        for col, w in {"A": 19.75, "B": 9.38, "C": 20.38, "D": 9.25,
                        "E": 45, "F": 14.63}.items():
            ws2.column_dimensions[col].width = w
        ws2.freeze_panes = "A3"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_path)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    config = load_config(Path(__file__).parent.parent / "config.yaml")
    supertree_path = config["apis"]["cell_hierarchy"]
    amended_folder = Path(config["output"]["amended_folder"])
    out_path = Path("outputs") / f"Supertree_colors_all_organs_{PROJECT_LABEL}.xlsx"

    dist_files = find_distribution_files(amended_folder)
    if not dist_files:
        print(f"No distribution JSONs found under {amended_folder}/*/json/. "
              f"Run 20_fetch_cell_distribution.py for at least one organ first, skipping.")
        return

    print(f"Loading supertree from: {supertree_path}")
    supertree, _ = load_supertree(supertree_path)

    blocks: List[List[Tuple]] = []
    all_unmatched: List[Tuple] = []
    for dist_path in dist_files:
        print(f"Processing: {dist_path}")
        rows, unmatched = build_rows_for_distribution(dist_path, supertree)
        if rows:
            blocks.append(rows)
        all_unmatched.extend(unmatched)

    if not blocks:
        print("No rows produced, nothing to write.")
        return

    # Group blocks by organ (each organ's blocks get merged into one Organ
    # cell), sub-sorted by sex/tool for a stable, readable layout.
    blocks.sort(key=lambda rows: (rows[0][0].lower(), rows[0][1].lower(), rows[0][2].lower()))

    write_xlsx(blocks, out_path, unmatched=all_unmatched)
    total_rows = sum(len(b) for b in blocks)
    print(f"\nWrote {total_rows} rows across {len(blocks)} organ/sex/tool combination(s) to: {out_path}")
    if all_unmatched:
        print(f"{len(all_unmatched)} cell type(s) had no supertree match, see the "
              f"\"Unclassified — Review\" tab.")


if __name__ == "__main__":
    main()