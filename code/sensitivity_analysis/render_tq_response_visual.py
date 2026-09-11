#!/usr/bin/env python3
"""Insert compact Tq-response data into the reusable visualization fragment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=(
            PROJECT_ROOT
            / "results/sensitivity_analysis/evidence_information_2026-08-24/tq_response_curves"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    frames = [
        pd.read_csv(input_dir / "node_tq_response_curves.tsv", sep="\t"),
        pd.read_csv(input_dir / "edge_tq_response_curves_seed891.tsv", sep="\t"),
    ]
    curves = pd.concat(frames, ignore_index=True)
    compact = []
    for row in curves.itertuples(index=False):
        compact.append(
            {
                "stage": row.stage,
                "stream_id": row.stream_id,
                "label": row.stream_label,
                "control": row.control_label.replace("×", "x"),
                "multiplier": round(float(row.multiplier), 8),
                "bf": round(float(row.geometric_mean_effective_bayes_factor), 8),
                "info": round(float(row.total_information_kl_bits), 8),
                "enabled": bool(row.enabled_by_default),
                "weight": round(float(row.weight), 6),
            }
        )
    template_path = SCRIPT_PATH.with_name("tq_response_visual.template.html")
    fragment = template_path.read_text(encoding="utf-8")
    payload = json.dumps(compact, separators=(",", ":"), ensure_ascii=False).replace(
        "</", "<\\/"
    )
    fragment = fragment.replace("__TQRV_DATA__", payload)
    if "__TQRV_DATA__" in fragment:
        raise RuntimeError("visualization data placeholder was not replaced")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(fragment, encoding="utf-8")
    print(f"Wrote {len(compact):,} response points to {args.output}")


if __name__ == "__main__":
    main()
