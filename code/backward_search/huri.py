#!/usr/bin/env python3
"""Materialize HuRI positive evidence and optional scoped non-report evidence."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from .engine import canonical_pair, utc_now


ENSEMBL_RE = re.compile(r"^ENSG\d+(?:\.\d+)?$", re.IGNORECASE)


def _ensembl(value: Any) -> str:
    token = str(value).strip().upper()
    return token.split(".", 1)[0] if ENSEMBL_RE.fullmatch(token) else ""


def read_huri_pairs(path: Path) -> list[tuple[str, str]]:
    """Read headered or headerless HuRI TSV pairs as canonical Ensembl IDs."""
    pairs: set[tuple[str, str]] = set()
    with path.resolve().open("r", encoding="utf-8-sig") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                continue
            left, right = _ensembl(fields[0]), _ensembl(fields[1])
            if left and right and left != right:
                pairs.add(tuple(sorted((left, right))))
    if not pairs:
        raise ValueError(f"no Ensembl interaction pairs were found in {path}")
    return sorted(pairs)


def _read_scope_identifiers(path: Path) -> set[str]:
    values: set[str] = set()
    with path.resolve().open("r", encoding="utf-8-sig") as handle:
        for raw_line in handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue
            for field in re.split(r"[\t,;]", raw_line.strip()):
                ensembl = _ensembl(field)
                if ensembl:
                    values.add(ensembl)
                    break
    if not values:
        raise ValueError(f"no Ensembl identifiers were found in {path}")
    return values


def build_huri_search_evidence(
    interactions_file: Path,
    mapping_file: Path,
    graph_nodes: pd.DataFrame,
    output_directory: Path,
    *,
    positive_bayes_factor: float,
    nonreported_bayes_factor: float,
    nonreported_scope: str = "neutral",
    screened_genes_file: Path | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build sparse positive factors plus an optional node scope for negatives.

    ``nonreported_scope`` is one of:
      * ``neutral``: no negative evidence from absence;
      * ``huri_interactor_set``: weak-negative eligible only when both mapped
        genes occur somewhere in the positive HuRI network (a broad assumption);
      * ``screened_gene_file``: eligible only when both mapped genes occur in a
        caller-supplied screen/testability list.
    """
    if positive_bayes_factor <= 1:
        raise ValueError("HuRI positive Bayes factor must be greater than 1")
    if not 0 < nonreported_bayes_factor <= 1:
        raise ValueError("HuRI non-report Bayes factor must lie in (0,1]")
    allowed_modes = {"neutral", "huri_interactor_set", "screened_gene_file"}
    if nonreported_scope not in allowed_modes:
        raise ValueError(f"HuRI non-report scope must be one of {sorted(allowed_modes)}")
    if nonreported_scope == "screened_gene_file" and screened_genes_file is None:
        raise ValueError("screened_gene_file mode requires a screened genes file")

    huri_pairs = read_huri_pairs(interactions_file)
    mapping = pd.read_csv(mapping_file.resolve(), sep="\t", dtype=str).fillna("")
    required = {"symbol", "hpa_ensembl_gene_id"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"HuRI orthology mapping lacks {sorted(required - set(mapping.columns))}")
    graph_symbols = set(graph_nodes["symbol"].astype(str))
    if "node_type" in graph_nodes.columns:
        graph_symbols = set(
            graph_nodes.loc[
                graph_nodes["node_type"].astype(str).str.casefold().eq("protein"), "symbol"
            ].astype(str)
        )
    human_to_mouse: dict[str, set[str]] = {}
    mouse_to_human: dict[str, str] = {}
    for row in mapping.to_dict("records"):
        symbol = str(row["symbol"]).strip()
        human = _ensembl(row["hpa_ensembl_gene_id"])
        if symbol in graph_symbols and human:
            human_to_mouse.setdefault(human, set()).add(symbol)
            mouse_to_human[symbol] = human

    positive_rows: dict[tuple[str, str], dict[str, Any]] = {}
    huri_interactors: set[str] = set()
    for human_a, human_b in huri_pairs:
        huri_interactors.update((human_a, human_b))
        for node_a in human_to_mouse.get(human_a, set()):
            for node_b in human_to_mouse.get(human_b, set()):
                if node_a == node_b:
                    continue
                pair = canonical_pair(node_a, node_b)
                positive_rows[pair] = {
                    "node_a": pair[0], "node_b": pair[1],
                    "bayes_factor": float(positive_bayes_factor),
                    "human_ensembl_a": human_a, "human_ensembl_b": human_b,
                    "huri_status": "reported_positive",
                }

    if nonreported_scope == "neutral":
        scoped_human: set[str] = set()
        scope_note = "HuRI non-reporting is neutral."
    elif nonreported_scope == "huri_interactor_set":
        scoped_human = huri_interactors
        scope_note = (
            "Weak negatives apply only when both orthologs occur somewhere in HuRI. "
            "This is an assumption, not an explicit failed-pair table."
        )
    else:
        scoped_human = _read_scope_identifiers(screened_genes_file.resolve())
        scope_note = (
            "Weak negatives apply only when both orthologs occur in the supplied "
            "screen/testability list."
        )
    scope_symbols = sorted(
        (symbol for symbol, human in mouse_to_human.items() if human in scoped_human),
        key=lambda item: (item.casefold(), item),
    )

    output_directory = output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    factor_path = output_directory / "huri_positive_factors.tsv.gz"
    scope_path = output_directory / "huri_nonreport_scope_nodes.tsv"
    pd.DataFrame(
        positive_rows.values(),
        columns=["node_a", "node_b", "bayes_factor", "human_ensembl_a", "human_ensembl_b", "huri_status"],
    ).to_csv(factor_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame({
        "symbol": scope_symbols,
        "human_ensembl_gene_id": [mouse_to_human[symbol] for symbol in scope_symbols],
    }).to_csv(scope_path, sep="\t", index=False)
    summary = {
        "schema_version": 1,
        "created_at": utc_now(),
        "source_file": str(interactions_file.resolve()),
        "mapping_file": str(mapping_file.resolve()),
        "source_huri_pair_count": len(huri_pairs),
        "graph_protein_count": len(graph_symbols),
        "mapped_graph_protein_count": len(mouse_to_human),
        "reported_positive_graph_pair_count": len(positive_rows),
        "nonreported_scope": nonreported_scope,
        "nonreported_scope_node_count": len(scope_symbols),
        "positive_bayes_factor": float(positive_bayes_factor),
        "nonreported_bayes_factor": float(nonreported_bayes_factor),
        "positive_factor_file": str(factor_path),
        "scope_nodes_file": str(scope_path),
        "scope_note": scope_note,
        "scientific_caveat": (
            "HuRI.tsv reports detected interactions, not a complete pair-level table of "
            "failed tests. An absent pair is therefore not a demonstrated noninteraction."
        ),
    }
    summary_path = output_directory / "huri_evidence_summary.json"
    summary["summary_file"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    stream = {
        "id": "huri_binary_interaction",
        "file": str(factor_path),
        "factor_column": "bayes_factor",
        "weight": 1.0,
        "missing_bayes_factor": 1.0,
        "scope_nodes_file": str(scope_path) if nonreported_scope != "neutral" else "",
        "scope_symbol_column": "symbol",
        "scoped_missing_bayes_factor": float(nonreported_bayes_factor),
    }
    return stream, summary
