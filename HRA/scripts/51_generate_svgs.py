#!/usr/bin/env python3
"""
51_generate_svgs.py: SVG versions of the HRA grouped horizontal bar charts.

Builds exactly the same Plotly figures as 50_generate_pngs.py (it imports that
script's chart code), so the SVGs match the PNGs: dark VR background, Nunito
body text and title. The only difference is vector output with the
Nunito font embedded as a base64 @font-face rule.

Note: Adobe Illustrator ignores embedded @font-face fonts and substitutes a
system font on import (install Nunito locally for Illustrator).
Browsers render the embedded fonts correctly.

Reads:  config.yaml
        outputs/<organ>/json/<stem>_distribution.json
        fonts/Nunito-VariableFont_wght.ttf
Writes: outputs/<organ>/svg/with-legend/<output_name>-with-legend.svg
        outputs/<organ>/svg/without-legend/<output_name>-without-legend.svg

Requires: plotly, kaleido
    pip install plotly kaleido
"""

from __future__ import annotations

import base64
import importlib.util
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from shared import (
    load_config, resolve_glb_filename, load_supertree, build_top9_granular_color_scheme,
    resolve_organ_filter_label,
)
from palette import LEVEL2_TOP9_COLORS, OTHER_COLOR, BG_MODE, apply_to_scheme

# Chart code lives in 50_generate_pngs.py; reuse it so the SVGs always match.
_spec = importlib.util.spec_from_file_location("hra_png_chart", Path(__file__).parent / "50_generate_pngs.py")
png = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(png)

FONTS_DIR = Path(__file__).parent.parent / "fonts"
EMBEDDED_FONTS = [
    # (css family, file, css font-weight, format)
    ("Nunito", FONTS_DIR / "Nunito-VariableFont_wght.ttf", "200 1000", "truetype"),
]


def font_face_css() -> str:
    faces = []
    for family, path, weight, fmt in EMBEDDED_FONTS:
        if not path.exists():
            print(f"  Warning: font not found at {path}; falling back to Arial for {family}.")
            continue
        b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        faces.append(
            f"@font-face {{font-family:'{family}';src:url(data:font/{path.suffix.lstrip('.')};base64,{b64}) "
            f"format('{fmt}');font-weight:{weight};font-style:normal;}}"
        )
    return f'<defs><style type="text/css">{"".join(faces)}</style></defs>' if faces else ""


def save_svg(fig, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        # Explicit width/height: some kaleido versions ignore fig.layout size.
        fig.write_image(str(output_path), format="svg", width=fig.layout.width, height=fig.layout.height)
        svg_text = output_path.read_text(encoding="utf-8")
        css = font_face_css()
        if css:
            svg_text = re.sub(r"(<svg[^>]*>)", lambda m: m.group(1) + css, svg_text, count=1)
        output_path.write_text(svg_text, encoding="utf-8")
        print(f"Saved SVG: {output_path}")
    except Exception as e:
        print(f"Error saving SVG (is kaleido installed?): {e}")
        print("Install with: pip install kaleido")


def main() -> None:
    config = load_config(Path(__file__).parent.parent / "config.yaml")

    organ_name = config["organ"]["name"]
    organ_sex = config["organ"]["sex"]
    tool = config["filters"]["tool"]
    shape = config["markers"]["shape"]
    supertree_path = config["apis"]["cell_hierarchy"]

    glb_filename = (config["organ"].get("glb_filename") or "").strip()
    organ_side = (config["organ"].get("side") or "").strip()
    if not glb_filename:
        glb_filename = resolve_glb_filename(organ_name, organ_sex, organ_side)

    stem = Path(glb_filename).stem
    base_folder = Path(config["output"]["amended_folder"]) / organ_name.lower()
    svg_with = base_folder / "svg" / "with-legend"
    svg_without = base_folder / "svg" / "without-legend"

    distribution_path = base_folder / "json" / f"{stem}_distribution.json"
    if not distribution_path.exists():
        raise FileNotFoundError(
            f"Distribution JSON not found: {distribution_path}. Run 20_fetch_cell_distribution.py first."
        )

    dist_data = json.loads(distribution_path.read_text(encoding="utf-8"))
    output_name = png.build_output_stem(stem, dist_data.get("uberon_id", ""), tool.lower(), shape)

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
    top9_labels = [label for label, _hex, _lvl2, _count in legend_entries if label != "Other"]

    as_labels, cell_type_counts, as_totals = png.prepare_chart_data(dist_data, top9_labels)
    wrapped_labels = [png.wrap_label(l) for l in as_labels]

    print(f"\nGenerating SVGs for: {output_name}")
    for show_legend, folder, suffix in ((True, svg_with, "with-legend"), (False, svg_without, "without-legend")):
        fig = png.make_figure(
            organ=organ_name, sex=organ_sex, tool=tool,
            as_labels=as_labels, wrapped_labels=wrapped_labels,
            cell_type_counts=cell_type_counts, as_totals=as_totals,
            legend_entries=legend_entries, show_legend=show_legend,
        )
        save_svg(fig, folder / f"{output_name}-{suffix}.svg")


if __name__ == "__main__":
    main()
