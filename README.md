# HRA 3D Cell Population Pipeline

Generates annotated 3D GLB files, per-marker CSVs, distribution JSONs, and PNG/SVG bar charts for HRA organs, with cell type markers placed inside anatomical structures. Data comes from [HRApop](https://humanatlas.io), the HRA reference-organ API, and the [HRA 3D Cell Generation API](https://apps.humanatlas.io/api/v1/mesh-3d-cell-population).

Cell type colors come from `data/Supertree_colors_all_organs_A5.xlsx` (unmatched cell types use the "Other" color); charts use the dark VR background and Nunito font.

## Setup

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Usage

**One organ / sex / tool**, edit `HRA/config.yaml`, then:

```bash
cd HRA
python setup_and_run.py
```

**Every organ / sex / tool in HRApop that has a reference-organ GLB** (config.yaml is not modified; 10,000 markers per organ):

```bash
cd HRA/scripts
python run_all_combos.py --dry-run --all     # show what would run
python run_all_combos.py --all               # JSON + CSV + GLB + PNG + SVG
python run_all_combos.py --pngs --svgs       # JSON + charts only
python run_all_combos.py --glbs --only heart,kidney
python run_all_combos.py --all --force       # regenerate existing outputs
```

The distribution JSON is always generated; `--glbs`, `--pngs`, `--svgs` (or `--all`) add the rest. Markers use a 10k count to keep GLBs light for the headset (some organs, e.g. right mammary gland, exceed 1M triangles). Skin has no GLB because the cell API returns no points for its mesh.

## Output structure

```
HRA/outputs/
├── downloaded organs/            reference organ GLBs (not committed)
├── glb report.csv                triangle count for every generated GLB
└── <organ>/
    ├── glb/                      GLB with cell markers
    ├── csv/                      per-marker position, cell type, cell ID, color
    ├── json/                     distribution data and legend metadata
    ├── png/{with,without}-legend/
    └── svg/{with,without}-legend/
```

File names: `<organ stem>-<uberon>-all-as-<tool>-sphere-hra-pop.<ext>`, e.g. `3d-vh-f-heart-UBERON_0000948-all-as-azimuth-sphere-hra-pop.glb`.

## Scripts (`HRA/scripts/`)

| Script | What it does |
|---|---|
| `10_download_organs.py` | Downloads the organ GLB from the HRA reference organs API |
| `20_fetch_cell_distribution.py` | Fetches cell type distributions per anatomical structure from HRApop |
| `30_generate_markers.py` | Places 3D cell markers; exports annotated GLB + CSV |
| `40_generate_legend.py` | Writes the legend JSON |
| `50_generate_pngs.py` | Grouped horizontal bar chart PNGs (fixed 1400×1376) |
| `51_generate_svgs.py` | Matching SVGs (Nunito embedded) |
| `60_generate_color_spreadsheet.py` | Rebuilds the cumulative color spreadsheet from all distribution JSONs |
| `palette.py` | Colors, background mode, and the color-sheet loader |
| `run_all_combos.py` | Mass generation for every organ / sex / tool |

`shared scripts/shared.py` holds helpers shared with other projects.

## Sample outputs

`HRA/outputs/` contains only heart, kidney, lung and liver for now, as examples of the chart layout and colors.

## Data

`data/hra_pop_raw.csv` is a cached HRApop download (used if the live API is unavailable); `data/Supertree_colors_all_organs_A5.xlsx` is the cell type color legend.

## Supported organs

Heart, kidney, knee, large intestine, liver, lung, lymph node, mammary gland, ovary, pancreas, prostate, skin, small intestine, spleen, thymus, ureter, urinary bladder, uterus, for the tools available in HRApop for each (Azimuth, CellTypist, FRMatch, Pan-Human Azimuth, PopV, SC Proteomics), female and male where present.
