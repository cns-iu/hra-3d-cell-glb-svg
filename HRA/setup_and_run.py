#!/usr/bin/env python3
"""
setup_and_run.py: Run the full HRA pipeline (10 -> 20 -> 30 -> 40 -> 50 -> 51 -> 60)
in order, against this folder's own config.yaml.

Stage 50 writes the PNG charts and stage 51 writes matching SVGs.

Stage 60 is different from the rest: it isn't scoped to the single organ
config.yaml points at. It rescans every distribution JSON already on disk
(one per organ/sex/tool previously run through this pipeline) and rebuilds
the full cumulative Supertree color spreadsheet from scratch, so it always
reflects every organ processed so far, not just the one from this run.

Chdirs into its own directory first, so config.yaml's relative output
paths always resolve under HRA/ regardless of invocation directory.

If any stage raises, the pipeline stops there, the traceback shows
which stage failed.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"


STAGES = [
    "10_download_organs",
    "20_fetch_cell_distribution",
    "30_generate_markers",
    "40_generate_legend",
    "50_generate_pngs",
    "51_generate_svgs",
    "60_generate_color_spreadsheet",
]


def run_stage(stage_name: str) -> None:
    module_path = SCRIPTS_DIR / f"{stage_name}.py"
    print(f"\n{'=' * 70}\n{stage_name}\n{'=' * 70}")
    spec = importlib.util.spec_from_file_location(stage_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.main()


def main() -> None:
    os.chdir(PROJECT_ROOT)
    for stage in STAGES:
        run_stage(stage)
    print(f"\n{'=' * 70}\nHRA pipeline complete.\n{'=' * 70}")


if __name__ == "__main__":
    main()