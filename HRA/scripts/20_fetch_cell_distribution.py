#!/usr/bin/env python3
"""
20_fetch_cell_distribution.py: Fetch cell type data from HRApop API for all
anatomical structures of the configured organ, auto-match AS labels to GLB
node names, then output per-AS distributions and proportional node allocations
as JSON for use by script 30.

Reads:  config.yaml
        outputs/downloaded organs/<glb_filename>  (to extract GLB node names)
        <cell_hierarchy CSV>                       (supertree for mapping)
Writes: outputs/<organ>/json/<stem>_distribution.json
"""

from __future__ import annotations

import csv
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import trimesh

sys.path.insert(0, str(Path(__file__).parent))  # shared.py lives next to this script
from shared import (
    resolve_glb_filename,
    resolve_organ_filter_label,
    load_config,
    normalize_cell_label_key,
    compact_ontology_id,
    parse_float,
    safe_name,
    values_match_filter,
    auto_match_as_to_glb_nodes,
    load_supertree,
)


# =============================================================================
# Fetching
# =============================================================================

def fetch_url_bytes(url: str, accept: str) -> bytes:
    request = urllib.request.Request(
        url=url,
        headers={"Accept": accept, "User-Agent": "hra-cell-pipeline/1.0"},
        method="GET",
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return response.read()


def rows_from_csv_bytes(raw: bytes) -> List[Dict[str, str]]:
    decoded = raw.decode("utf-8-sig", errors="replace")
    return [dict(row) for row in csv.DictReader(decoded.splitlines())]


# grlc-hosted endpoints (like apis.hra_pop) run this same query under the
# hood. When grlc itself is down, hit the underlying SPARQL endpoint
# directly with the query grlc would have run, see the HuBMAP workaround:
# github.com/hubmapconsortium/ccf-grlc/blob/main/grlc-workaround.md, query
# source: github.com/hubmapconsortium/ccf-grlc/blob/main/hra-pop/
# cell_types_in_anatomical_structurescts_per_as.rq
# Unlike the plain GitHub CSV export, this returns organ_id (UBERON/FMA),
# matching the columns fetch_hra_pop_rows() already expects.
HRA_POP_SPARQL_QUERY = """
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX ccf: <http://purl.org/ccf/>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX UBERON: <http://purl.obolibrary.org/obo/UBERON_>
PREFIX FMA: <http://purl.org/sig/ont/fma/fma>
PREFIX HRA: <https://purl.humanatlas.io/collection/hra-api>
PREFIX HRApop: <https://purl.humanatlas.io/graph/hra-pop>

SELECT DISTINCT ?organ ?organ_id ?as ?as_label ?sex ?tool ?modality ?cell_id ?cell_label ?cell_count ?cell_percentage ?dataset_count
FROM HRA:
FROM HRApop:
WHERE {
  ?as ccf:has_cell_summary [
      ccf:cell_annotation_method ?tool ;
      ccf:cell_source_label ?as_label ;
      ccf:sex ?sex ;
      ccf:modality ?modality ;
      ccf:aggregated_summary_count ?dataset_count ;
      ccf:has_cell_summary_row [
        ccf:cell_id ?cell_id ;
        ccf:cell_label ?summary_cell_label ;
        ccf:cell_count ?cell_count ;
        ccf:percentage_of_total ?cell_percentage ;
      ]
  ] .
  FILTER(STRSTARTS(STR(?as), STR(UBERON:)) || STRSTARTS(STR(?as), STR(FMA:)))

  OPTIONAL {
    ?cell_id rdfs:label ?cell_rdfs_label .
  }
  BIND(COALESCE(?cell_rdfs_label, ?summary_cell_label) as ?cell_label)

  ?refOrganPiece ccf:representation_of ?as ;
      ccf:has_reference_organ ?refOrgan .
  ?refOrgan ccf:organ_owner_sex ?sex .
  ?refOrgan skos:prefLabel ?organ .
  ?refOrgan ccf:representation_of ?organ_id .
}
ORDER BY ?sex ?as ?tool DESC(?cell_count)
"""


def fetch_hra_pop_via_sparql(sparql_endpoint: str) -> bytes:
    query_string = urllib.parse.urlencode({"query": HRA_POP_SPARQL_QUERY})
    url = f"{sparql_endpoint}?{query_string}"
    return fetch_url_bytes(url, accept="text/csv")


def fetch_hra_pop_rows(api_url: str) -> Tuple[List[Dict[str, str]], bytes]:
    """Try CSV then JSON formats. Returns (rows, raw_bytes)."""
    errors: List[str] = []

    for url in [api_url, api_url + ".csv", api_url + "?_format=csv"]:
        try:
            raw = fetch_url_bytes(url, accept="text/csv")
            decoded = raw.decode("utf-8-sig", errors="replace")
            reader = csv.DictReader(decoded.splitlines())
            rows = [dict(row) for row in reader]
            if rows and "cell_label" in rows[0]:
                return rows, raw
        except Exception as exc:
            errors.append(f"CSV attempt failed for {url}: {exc}")

    for url in [api_url, api_url + ".json", api_url + "?_format=json"]:
        try:
            raw = fetch_url_bytes(url, accept="application/json")
            data = json.loads(raw.decode("utf-8", errors="replace"))
            if isinstance(data, list):
                return [dict(r) for r in data], raw
            if isinstance(data, dict):
                for key in ("results", "data", "rows"):
                    maybe = data.get(key)
                    if isinstance(maybe, list):
                        return [dict(r) for r in maybe], raw
        except Exception as exc:
            errors.append(f"JSON attempt failed for {url}: {exc}")

    raise RuntimeError("Could not fetch HRApop rows.\n" + "\n".join(errors))


# =============================================================================
# GLB node extraction
# =============================================================================

def get_glb_node_names(glb_path: Path) -> List[str]:
    loaded = trimesh.load(str(glb_path), force="scene")
    if isinstance(loaded, trimesh.Scene):
        return list(loaded.graph.nodes_geometry)
    return []


# =============================================================================
# UBERON / organ ID extraction
# =============================================================================

def extract_organ_uberon_id(
    rows: List[Dict[str, str]],
    organ: str,
    tool: str,
    sex: str,
) -> str:
    for row in rows:
        if not values_match_filter(row.get("organ", ""), organ):
            continue
        if not values_match_filter(row.get("tool", ""), tool):
            continue
        if not values_match_filter(row.get("sex", ""), sex):
            continue
        raw_id = str(row.get("organ_id", "")).strip()
        if raw_id:
            tail = raw_id.rstrip("/").split("/")[-1]
            return tail.replace(":", "_")
    return ""


# =============================================================================
# Per-AS distribution building
# =============================================================================

def get_all_as_labels(
    rows: List[Dict[str, str]],
    organ: str,
    tool: str,
    sex: str,
    modality: str,
) -> List[str]:
    seen = set()
    result = []
    for row in rows:
        if not values_match_filter(row.get("organ", ""), organ):
            continue
        if not values_match_filter(row.get("tool", ""), tool):
            continue
        if not values_match_filter(row.get("sex", ""), sex):
            continue
        if not values_match_filter(row.get("modality", ""), modality):
            continue
        as_label = str(row.get("as_label", "")).strip()
        if as_label and as_label not in seen:
            seen.add(as_label)
            result.append(as_label)
    return result


def filter_rows_for_as(
    rows: List[Dict[str, str]],
    organ: str,
    as_label: str,
    tool: str,
    sex: str,
    modality: str,
) -> List[Dict[str, str]]:
    return [
        row for row in rows
        if values_match_filter(row.get("organ", ""), organ)
        and values_match_filter(row.get("as_label", ""), as_label)
        and values_match_filter(row.get("tool", ""), tool)
        and values_match_filter(row.get("sex", ""), sex)
        and values_match_filter(row.get("modality", ""), modality)
    ]


def build_distribution_for_as(rows: List[Dict[str, str]]) -> Dict[str, float]:
    distribution: Dict[str, float] = {}
    for row in rows:
        label = str(row.get("cell_label", "")).strip()
        if not label:
            continue
        pct   = parse_float(row.get("cell_percentage", ""))
        count = parse_float(row.get("cell_count", ""))
        weight = pct if pct > 0 else count
        if weight > 0:
            distribution[label] = distribution.get(label, 0.0) + weight
    return distribution


def normalize_distribution(distribution: Dict[str, float]) -> Dict[str, float]:
    total = sum(distribution.values())
    if total <= 0:
        return {}
    return {k: v / total for k, v in distribution.items()}


def get_total_cell_count_for_as(rows: List[Dict[str, str]]) -> float:
    return sum(parse_float(row.get("cell_count", "")) for row in rows)


def build_organ_cell_counts(
    rows: List[Dict[str, str]],
    organ: str,
    tool: str,
    sex: str,
    modality: str,
    matched_as_labels: List[str],
) -> Dict[str, float]:
    """Organ-wide raw cell counts per cell type, summed across all matched AS."""
    organ_counts: Dict[str, float] = {}
    for as_label in matched_as_labels:
        as_rows = filter_rows_for_as(rows, organ, as_label, tool, sex, modality)
        for row in as_rows:
            label = str(row.get("cell_label", "")).strip()
            if not label:
                continue
            count = parse_float(row.get("cell_count", ""))
            if count > 0:
                organ_counts[label] = organ_counts.get(label, 0.0) + count
    return organ_counts


def build_cell_id_map(rows: List[Dict[str, str]]) -> Dict[str, str]:
    """Map normalized cell label key → compact CL ID (e.g. CL:0000057)."""
    cell_id_map: Dict[str, str] = {}
    for row in rows:
        label  = str(row.get("cell_label", "")).strip()
        raw_id = str(row.get("cell_id",    "")).strip()
        if not label or not raw_id:
            continue
        key = normalize_cell_label_key(label)
        if key not in cell_id_map:
            cell_id_map[key] = compact_ontology_id(raw_id)
    return cell_id_map


# =============================================================================
# Supertree level ranking
# =============================================================================

def build_color_order(
    organ_cell_counts: Dict[str, float],
    cell_id_map: Dict[str, str],
    supertree: Dict[str, str],
) -> List[str]:
    """
    Aggregate organ-wide cell counts per supertree level, then return levels
    ordered from most to least cells. Informational only, actual coloring
    is done by shared.build_top9_granular_color_scheme() (ranks granular
    cell types directly, not by level).
    """
    level_counts: Dict[str, float] = {}
    for cell_label, count in organ_cell_counts.items():
        cl_id = cell_id_map.get(normalize_cell_label_key(cell_label), "")
        level = supertree.get(cl_id)
        if level:
            level_counts[level] = level_counts.get(level, 0.0) + count

    return [l for l, _ in sorted(level_counts.items(), key=lambda x: x[1], reverse=True)]


# =============================================================================
# Node allocation
# =============================================================================

def allocate_nodes(
    as_total_counts: Dict[str, float],
    total_nodes: int,
) -> Dict[str, int]:
    grand_total = sum(as_total_counts.values())
    if grand_total <= 0:
        equal = total_nodes // len(as_total_counts)
        return {k: equal for k in as_total_counts}

    allocated: Dict[str, int] = {}
    remainder = total_nodes
    as_list   = list(as_total_counts.items())

    for i, (as_label, count) in enumerate(as_list):
        if i == len(as_list) - 1:
            allocated[as_label] = remainder
        else:
            n = round(total_nodes * count / grand_total)
            allocated[as_label] = n
            remainder -= n

    return allocated


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    config = load_config(Path(__file__).parent.parent / "config.yaml")

    api_url        = config["apis"]["hra_pop"]
    supertree_path = config["apis"]["cell_hierarchy"]
    organ_name     = config["organ"]["name"]
    organ_sex      = config["organ"]["sex"]
    glb_filename   = (config["organ"].get("glb_filename") or "").strip()
    tool           = config["filters"]["tool"]
    modality       = config["filters"].get("modality", "")
    total_nodes    = config["markers"]["num_nodes"]
    save_csv       = config["output"].get("save_hra_pop_csv", False)
    organs_folder  = Path(config["output"]["organs_folder"])

    organ_folder = Path(config["output"]["amended_folder"]) / organ_name.lower()
    json_folder  = organ_folder / "json"
    json_folder.mkdir(parents=True, exist_ok=True)

    organ_side = (config["organ"].get("side") or "").strip()
    if not glb_filename:
        glb_filename = resolve_glb_filename(organ_name, organ_sex, organ_side)
        print(f"Auto-resolved GLB filename: {glb_filename}")

    # HRApop CSV rows label lateralized organs as "Left kidney" / "Right
    # kidney", not plain "kidney", use this for all row-matching below so
    # config's plain organ name + side still finds the right rows without
    # accidentally pulling in both sides' data.
    organ_filter = resolve_organ_filter_label(organ_name, organ_side)
    if organ_filter != organ_name:
        print(f"Matching HRApop rows against organ='{organ_filter}' (side='{organ_side}')")

    stem     = Path(glb_filename).stem
    glb_path = organs_folder / glb_filename

    if not glb_path.exists():
        raise FileNotFoundError(
            f"GLB not found: {glb_path}. Run 10_download_organs.py first."
        )

    # Load supertree
    print(f"Loading supertree from: {supertree_path}")
    supertree, _ = load_supertree(supertree_path)
    print(f"Supertree entries loaded: {len(supertree)}")

    # Get GLB node names
    print(f"Reading GLB node names from: {glb_path}")
    glb_node_names = get_glb_node_names(glb_path)
    print(f"GLB nodes found: {len(glb_node_names)}")
    for n in sorted(glb_node_names):
        print(f"  {n}")

    # Fetch HRApop data. Fallback chain if the live grlc API is down:
    # 1) query the underlying SPARQL endpoint directly with the same query
    #    grlc runs (see HRA_POP_SPARQL_QUERY above), the HuBMAP-documented
    #    workaround for a downed grlc.io, and the only fallback that still
    #    includes organ_id/UBERON like the live API does.
    # 2) the local cache under "data/" (shared reference-data folder at
    #    the repo root), for fully offline runs.
    sparql_endpoint = config["apis"].get("hra_pop_sparql", "")
    cache_path = Path(__file__).parent.parent.parent / "data" / "hra_pop_raw.csv"
    print(f"\nFetching HRApop data from: {api_url}")
    try:
        rows, raw = fetch_hra_pop_rows(api_url)
    except Exception as api_exc:
        rows = None
        if sparql_endpoint:
            print(f"  Live fetch failed ({api_exc}); trying SPARQL endpoint: {sparql_endpoint}")
            try:
                raw = fetch_hra_pop_via_sparql(sparql_endpoint)
                rows = rows_from_csv_bytes(raw)
            except Exception as sparql_exc:
                print(f"  SPARQL fallback failed ({sparql_exc})")
        if rows is None:
            if not cache_path.exists():
                raise
            print(f"  Falling back to cached CSV: {cache_path}")
            raw = cache_path.read_bytes()
            rows = rows_from_csv_bytes(raw)
    print(f"Total rows fetched: {len(rows)}")

    if save_csv:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_bytes(raw)
        print(f"Saved raw HRApop CSV to: {cache_path}")

    # Extract UBERON ID
    uberon_id = extract_organ_uberon_id(rows, organ_filter, tool, organ_sex)
    print(f"\nOrgan UBERON ID: {uberon_id or '(not found)'}")

    # Get all AS labels
    as_labels = get_all_as_labels(rows, organ_filter, tool, organ_sex, modality)
    print(f"\nHRApop AS labels found: {len(as_labels)}")
    for label in as_labels:
        print(f"  {label}")

    if not as_labels:
        raise ValueError(
            f"No AS labels found for organ='{organ_filter}', "
            f"tool='{tool}', sex='{organ_sex}'. Check config.yaml."
        )

    # Auto-match AS labels to GLB nodes
    print(f"\nAuto-matching AS labels to GLB nodes (organ='{organ_name}'):")
    as_to_glb_node = auto_match_as_to_glb_nodes(as_labels, glb_node_names, organ_name)

    matched_as_labels = [l for l, n in as_to_glb_node.items() if n is not None]
    unmatched         = [l for l, n in as_to_glb_node.items() if n is None]

    print(f"\nMatched: {len(matched_as_labels)} | Unmatched (skipped): {len(unmatched)}")

    if not matched_as_labels:
        raise ValueError(
            "No AS labels could be matched to GLB nodes. "
            "Check organ name and GLB file."
        )

    # Build per-AS distributions
    per_as_distributions: Dict[str, Dict[str, float]] = {}
    per_as_normalized:    Dict[str, Dict[str, float]] = {}
    as_total_counts:      Dict[str, float]            = {}

    for as_label in matched_as_labels:
        as_rows    = filter_rows_for_as(rows, organ_filter, as_label, tool, organ_sex, modality)
        dist       = build_distribution_for_as(as_rows)
        total_count = get_total_cell_count_for_as(as_rows)

        per_as_distributions[as_label] = dist
        per_as_normalized[as_label]    = normalize_distribution(dist)
        as_total_counts[as_label]      = total_count

        print(f"  {as_label}: {len(dist)} cell types, total count={total_count:.0f}")

    # Organ-wide cell counts and cell ID map
    organ_cell_counts = build_organ_cell_counts(
        rows, organ_filter, tool, organ_sex, modality, matched_as_labels
    )
    cell_id_map = build_cell_id_map(rows)
    print(f"\nOrgan-wide cell type count entries: {len(organ_cell_counts)}")
    print(f"Cell ID map entries: {len(cell_id_map)}")

    # Cell type color order
    color_order = build_color_order(organ_cell_counts, cell_id_map, supertree)
    print(f"\nSupertree levels found for this organ ({len(color_order)}):")
    for i, b in enumerate(color_order):
        print(f"  [{i}] {b}")

    # Node allocation
    node_allocation = allocate_nodes(as_total_counts, total_nodes)
    print(f"\nNode allocation (total={total_nodes}):")
    for as_label, n in node_allocation.items():
        print(f"  {as_label}: {n} nodes")

    # Save distribution JSON
    output_path = json_folder / f"{stem}_distribution.json"
    output_data = {
        "stem":                stem,
        "organ":               organ_name,
        "sex":                 organ_sex,
        "tool":                tool,
        "modality":            modality,
        "uberon_id":           uberon_id,
        "cell_type_color_order": color_order,
        "as_labels":           matched_as_labels,
        "as_to_glb_node":      {k: v for k, v in as_to_glb_node.items() if v is not None},
        "as_node_counts":      {k: len(v) for k, v in as_to_glb_node.items() if v is not None},
        "unmatched_as_labels": unmatched,
        "per_as_distributions": per_as_distributions,
        "per_as_normalized":   per_as_normalized,
        "as_total_counts":     as_total_counts,
        "node_allocation":     node_allocation,
        "organ_cell_counts":   organ_cell_counts,
        "cell_id_map":         cell_id_map,
    }

    output_path.write_text(
        json.dumps(output_data, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\nDistribution saved to: {output_path}")


if __name__ == "__main__":
    main()
