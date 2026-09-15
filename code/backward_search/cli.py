#!/usr/bin/env python3
"""Command-line interface for Version 2 backward frontier search."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .engine import (
    build_node_sequence_fasta,
    build_external_target_evidence,
    collect_alphapulldown_scores,
    import_structural_scores,
    initialize_run,
    load_run_state,
    step_run,
)
from .development import record_structural_import, render_development_report
from .legacy_screen import LEGACY_PROTOCOL_ID, import_legacy_ppi_screen
from .pipeline import (
    collect_pipeline_structural_scores,
    default_pipeline_configuration,
    initialize_end_to_end_pipeline,
    load_pipeline_state,
    refresh_pipeline_state,
    render_pipeline_trace,
    step_pipeline,
    submit_pipeline_structural_round,
)


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)
    init = commands.add_parser("init", help="Create a new resumable search run.")
    init.add_argument("--config", type=Path, required=True)
    init.add_argument("--run-dir", type=Path, required=True)
    init.add_argument(
        "--development",
        action="store_true",
        help="Pause after every auditable decision instead of every complete round.",
    )
    step = commands.add_parser(
        "step",
        help="Advance one complete round, or one decision when development mode is active.",
    )
    step.add_argument("--run-dir", type=Path, required=True)
    status = commands.add_parser("status", help="Print the saved run state.")
    status.add_argument("--run-dir", type=Path, required=True)
    imported = commands.add_parser("import", help="Import standardized structural scores.")
    imported.add_argument("--run-dir", type=Path, required=True)
    imported.add_argument("--scores", type=Path, required=True)
    collect = commands.add_parser("collect", help="Collect Biowulf AlphaPulldown ipTM outputs.")
    collect.add_argument("--run-dir", type=Path, required=True)
    trace = commands.add_parser(
        "trace",
        help="Regenerate and print the path to a development run's HTML trace.",
    )
    trace.add_argument("--run-dir", type=Path, required=True)
    fasta = commands.add_parser("build-fasta", help="Build a node-symbol FASTA from a UniProt table.")
    fasta.add_argument("--nodes", type=Path, required=True)
    fasta.add_argument("--uniprot", type=Path, required=True)
    fasta.add_argument("--output", type=Path, required=True)
    fasta.add_argument("--symbol-column", default="symbol")
    target = commands.add_parser(
        "build-target-evidence",
        help="Append an external target's HPA factors to an internal pair catalog.",
    )
    target.add_argument("--project-root", type=Path, required=True)
    target.add_argument("--target", required=True)
    target.add_argument("--base-evidence", type=Path, required=True)
    target.add_argument("--base-factor-column", default="bayes_factor")
    target.add_argument("--output", type=Path, required=True)
    legacy = commands.add_parser(
        "import-legacy-screen",
        help="Validate ppi_screen scores/structures and transfer them into the V2 pair cache.",
    )
    legacy.add_argument("--screen-dir", type=Path, required=True)
    legacy.add_argument("--cache", type=Path, required=True)
    legacy.add_argument("--audit-dir", type=Path)
    legacy.add_argument("--protocol-id", default=LEGACY_PROTOCOL_ID)
    legacy.add_argument("--allow-scores-without-structure", action="store_true")

    pipeline_config = commands.add_parser(
        "pipeline-config",
        help="Write the complete default node/edge/backward-search configuration.",
    )
    pipeline_config.add_argument("--output", type=Path, required=True)
    pipeline_init = commands.add_parser(
        "pipeline-init",
        help="Run node selection and cheap edge integration, then initialize the frontier search.",
    )
    pipeline_init.add_argument("--config", type=Path, required=True)
    pipeline_init.add_argument("--run-dir", type=Path, required=True)
    pipeline_step = commands.add_parser(
        "pipeline-step", help="Advance a complete pipeline by one visible decision stage."
    )
    pipeline_step.add_argument("--run-dir", type=Path, required=True)
    pipeline_step.add_argument(
        "--until-checkpoint",
        action="store_true",
        help="Continue until AlphaPulldown input is needed or the search completes.",
    )
    pipeline_status = commands.add_parser(
        "pipeline-status", help="Print the persisted complete-pipeline state."
    )
    pipeline_status.add_argument("--run-dir", type=Path, required=True)
    pipeline_status.add_argument("--full", action="store_true", help="Print the full scientific summaries.")
    pipeline_submit = commands.add_parser(
        "pipeline-submit", help="Submit the current AlphaPulldown round on Biowulf."
    )
    pipeline_submit.add_argument("--run-dir", type=Path, required=True)
    pipeline_collect = commands.add_parser(
        "pipeline-collect", help="Collect the current AlphaPulldown round into the pair cache."
    )
    pipeline_collect.add_argument("--run-dir", type=Path, required=True)
    pipeline_trace = commands.add_parser(
        "pipeline-trace", help="Regenerate the complete pipeline's stepwise HTML trace."
    )
    pipeline_trace.add_argument("--run-dir", type=Path, required=True)
    return root


def main(argv: list[str] | None = None) -> int:
    arguments = parser().parse_args(argv)

    def pipeline_summary(state: dict[str, object]) -> dict[str, object]:
        return {
            "status": state.get("status"),
            "run_directory": state.get("run_directory"),
            "selected_node_count": state.get("selected_node_count"),
            "search_round": state.get("search_round"),
            "next_search_stage": state.get("next_search_stage"),
            "pending_unique_structural_pairs": state.get("pending_unique_structural_pairs"),
            "solution_count": state.get("solution_count"),
            "completion_reason": state.get("completion_reason", ""),
            "development_trace": state.get("development_trace", ""),
        }
    if arguments.command == "pipeline-config":
        configuration = default_pipeline_configuration()
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(configuration, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"configuration": str(arguments.output.resolve())}, indent=2))
        return 0
    if arguments.command == "import-legacy-screen":
        print(json.dumps(import_legacy_ppi_screen(
            arguments.screen_dir,
            arguments.cache,
            protocol_id=arguments.protocol_id,
            require_ranked_structure=not arguments.allow_scores_without_structure,
            audit_directory=arguments.audit_dir,
        ), indent=2))
        return 0
    if arguments.command == "pipeline-init":
        configuration = json.loads(arguments.config.read_text(encoding="utf-8"))
        state = initialize_end_to_end_pipeline(configuration, arguments.run_dir)
        print(json.dumps(pipeline_summary(state), indent=2, default=str))
        return 0
    if arguments.command == "pipeline-step":
        state = step_pipeline(arguments.run_dir, until_checkpoint=arguments.until_checkpoint)
        print(json.dumps(pipeline_summary(state), indent=2, default=str))
        return 0
    if arguments.command == "pipeline-status":
        state = refresh_pipeline_state(arguments.run_dir)
        print(json.dumps(state if arguments.full else pipeline_summary(state), indent=2, default=str))
        return 0
    if arguments.command == "pipeline-submit":
        state = submit_pipeline_structural_round(arguments.run_dir)
        print(json.dumps(pipeline_summary(state), indent=2, default=str))
        return 0
    if arguments.command == "pipeline-collect":
        state = collect_pipeline_structural_scores(arguments.run_dir)
        print(json.dumps(pipeline_summary(state), indent=2, default=str))
        return 0
    if arguments.command == "pipeline-trace":
        report = render_pipeline_trace(arguments.run_dir)
        print(json.dumps({"development_trace": str(report)}, indent=2))
        return 0
    if arguments.command == "init":
        state = initialize_run(
            arguments.config,
            arguments.run_dir,
            development_mode=True if arguments.development else None,
        )
    elif arguments.command == "step":
        state = step_run(arguments.run_dir)
    elif arguments.command == "status":
        state = load_run_state(arguments.run_dir)
    elif arguments.command == "import":
        state = import_structural_scores(arguments.run_dir, arguments.scores)
        state = record_structural_import(arguments.run_dir.resolve(), state)
    elif arguments.command == "collect":
        score_path = collect_alphapulldown_scores(arguments.run_dir)
        state = import_structural_scores(arguments.run_dir, score_path)
        state = record_structural_import(arguments.run_dir.resolve(), state)
    elif arguments.command == "trace":
        report = render_development_report(arguments.run_dir)
        print(json.dumps({"development_trace": str(report)}, indent=2))
        return 0
    elif arguments.command == "build-fasta":
        print(json.dumps(build_node_sequence_fasta(
            arguments.nodes,
            arguments.uniprot,
            arguments.output,
            symbol_column=arguments.symbol_column,
        ), indent=2))
        return 0
    else:
        print(json.dumps(build_external_target_evidence(
            arguments.project_root,
            arguments.target,
            arguments.base_evidence,
            arguments.output,
            base_factor_column=arguments.base_factor_column,
        ), indent=2))
        return 0
    print(json.dumps({
        "status": state["status"],
        "round": state["round"],
        "beam_size": len(state.get("beam", [])),
        "solution_count": len(state.get("solutions", [])),
        "pending_unique_structural_pairs": state.get("pending_unique_structural_pairs", 0),
        "round_directory": state.get("round_directory", ""),
        "completion_reason": state.get("completion_reason", ""),
        "development_mode": bool(state.get("development_mode", False)),
        "next_development_stage": state.get("development_stage", ""),
        "development_step_count": state.get("development_step_count", 0),
        "development_trace": str(arguments.run_dir.resolve() / "development_trace.html")
        if state.get("development_mode") else "",
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
