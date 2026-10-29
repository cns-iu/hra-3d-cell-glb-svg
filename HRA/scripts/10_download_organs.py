#!/usr/bin/env python3
"""
10_download_organs.py: Download organ GLB files from the HRA reference-organs API.

Reads config.yaml to determine which organ/sex to download.
Saves GLB files to the configured output folder.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from urllib.parse import urlparse

import requests

sys.path.insert(0, str(Path(__file__).parent))  # shared.py lives next to this script
from shared import load_config, safe_name, resolve_glb_filename


def fetch_reference_organs(api_url: str) -> list:
    """Fetch the full list of reference organs from the HRA API."""
    print(f"Fetching reference organs from: {api_url}")
    response = requests.get(api_url, timeout=120)
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, list):
        raise ValueError("Reference organs API did not return a list.")
    return data


def matches_organ_config(organ: dict, name_filter: str, sex_filter: str) -> bool:
    """Return True if this organ entry matches the configured organ name and sex."""
    organ_label = str(organ.get("label", "")).lower()
    organ_sex = str(organ.get("sex", "")).lower()

    name_match = name_filter.lower() in organ_label
    sex_match = sex_filter.lower() in organ_sex or sex_filter == ""

    return name_match and sex_match


def download_glb(glb_url: str, output_folder: Path) -> Path:
    """Download a single GLB file into output_folder. Skips if already present."""
    parsed = urlparse(glb_url)
    file_name = os.path.basename(parsed.path)
    file_path = output_folder / file_name

    if file_path.exists():
        print(f"  Already exists, skipping: {file_name}")
        return file_path

    print(f"  Downloading: {file_name}")
    response = requests.get(glb_url, timeout=120)
    if response.status_code == 200:
        file_path.write_bytes(response.content)
        print(f"  Saved: {file_path}")
    else:
        print(f"  Warning: Failed to download {file_name} (HTTP {response.status_code})")

    return file_path


def main() -> None:
    config = load_config(Path(__file__).parent.parent / "config.yaml")

    api_url = config["apis"]["reference_organs"]
    organ_name = config["organ"]["name"]
    organ_sex = config["organ"]["sex"]
    glb_filename = (config["organ"].get("glb_filename") or "").strip()
    output_folder = Path(config["output"]["organs_folder"])

    output_folder.mkdir(parents=True, exist_ok=True)

    organs = fetch_reference_organs(api_url)
    print(f"Total organs returned by API: {len(organs)}")

    # If a specific GLB filename is configured, download only that one
    if glb_filename:
        matched = [
            o for o in organs
            if str(o.get("object", {}).get("file", "")).endswith(glb_filename)
        ]
        if not matched:
            raise ValueError(
                f"Could not find configured glb_filename '{glb_filename}' in API response."
            )
        for organ in matched:
            glb_url = organ["object"]["file"]
            download_glb(glb_url, output_folder)
        print(f"\nDownloaded {len(matched)} GLB file(s) to: {output_folder}")
        return

    # Otherwise filter by organ name + sex
    matched = [o for o in organs if matches_organ_config(o, organ_name, organ_sex)]
    print(f"Organs matching name='{organ_name}' sex='{organ_sex}': {len(matched)}")

    if not matched:
        raise ValueError(
            f"No organs matched name='{organ_name}' and sex='{organ_sex}'. "
            "Check config.yaml or set glb_filename explicitly."
        )

    for organ in matched:
        glb_url = organ.get("object", {}).get("file", "")
        if glb_url:
            download_glb(glb_url, output_folder)

    print(f"\nDownloaded {len(matched)} GLB file(s) to: {output_folder}")


if __name__ == "__main__":
    main()
