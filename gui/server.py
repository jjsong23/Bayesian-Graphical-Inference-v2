#!/usr/bin/env python3
"""Local HTTP server for the Graphical Bayesian Inference workbench."""

from __future__ import annotations

import argparse
import json
import mimetypes
import os
import sys
import threading
import traceback
import uuid
import webbrowser
from dataclasses import dataclass, field
from datetime import datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, unquote, urlparse

from workflow_engine import (
    PROJECT_ROOT,
    WorkflowCancelled,
    default_configuration,
    load_registry,
    run_workflow,
)
from incremental_edge_cache import current_cache_counts
from evidence_inspector import inspect_edge_evidence, inspect_node_evidence


CODE_ROOT = PROJECT_ROOT / "code"
if str(CODE_ROOT) not in sys.path:
    sys.path.insert(0, str(CODE_ROOT))
from backward_search.pipeline import (  # noqa: E402
    collect_pipeline_structural_scores,
    default_pipeline_configuration,
    initialize_end_to_end_pipeline,
    refresh_pipeline_state,
    render_pipeline_trace,
    step_pipeline,
    submit_pipeline_structural_round,
)


WEB_ROOT = Path(__file__).resolve().parent / "web"


@dataclass
class Job:
    job_id: str
    status: str = "queued"
    message: str = "Queued"
    progress: float = 0.0
    configuration: dict[str, Any] = field(default_factory=dict)
    preview: dict[str, Any] | None = None
    summary: dict[str, Any] | None = None
    output_directory: str | None = None
    error: str | None = None
    cancel_event: threading.Event = field(default_factory=threading.Event, repr=False)


JOBS: dict[str, Job] = {}
JOBS_LOCK = threading.Lock()
ANALYSIS_LOCK = threading.Lock()
TERMINAL_JOB_STATUSES = {"complete", "failed", "cancelled"}
V2_RUN_ROOT = PROJECT_ROOT / "runs"


def resolve_evidence_root() -> Path:
    explicit = os.environ.get("GBI_EVIDENCE_ROOT", "").strip()
    candidates = [Path(explicit)] if explicit else []
    candidates.extend([PROJECT_ROOT, PROJECT_ROOT.parent / "graphical_bayesian_inference"])
    for candidate in candidates:
        if (candidate / "data/node_selection/mouse_signaling_nodes_liberal.tsv").is_file():
            return candidate.resolve()
    return PROJECT_ROOT


EVIDENCE_ROOT = resolve_evidence_root()


@dataclass
class V2Job:
    job_id: str
    configuration: dict[str, Any] = field(default_factory=dict)
    status: str = "queued"
    operation: str = "initialize"
    message: str = "Queued"
    progress: float = 0.0
    pipeline: dict[str, Any] | None = None
    error: str | None = None
    busy: bool = False
    cancel_event: threading.Event = field(default_factory=threading.Event, repr=False)

    @property
    def run_directory(self) -> Path:
        return V2_RUN_ROOT / self.job_id


V2_JOBS: dict[str, V2Job] = {}
V2_JOBS_LOCK = threading.Lock()


def public_job(job: Job) -> dict[str, Any]:
    return {
        "job_id": job.job_id,
        "status": job.status,
        "message": job.message,
        "progress": job.progress,
        "preview": job.preview,
        "summary": job.summary,
        "files": job.preview.get("files", []) if job.preview else [],
        "error": job.error,
        "cancel_requested": job.cancel_event.is_set(),
    }


def public_v2_job(job: V2Job) -> dict[str, Any]:
    return {
        "job_id": job.job_id,
        "status": job.status,
        "operation": job.operation,
        "message": job.message,
        "progress": job.progress,
        "pipeline": job.pipeline,
        "error": job.error,
        "busy": job.busy,
        "cancel_requested": job.cancel_event.is_set(),
        "run_directory": str(job.run_directory),
    }


def recover_v2_job(job_id: str) -> V2Job | None:
    with V2_JOBS_LOCK:
        existing = V2_JOBS.get(job_id)
    if existing is not None:
        return existing
    state_path = V2_RUN_ROOT / job_id / "pipeline_state.json"
    if not state_path.is_file():
        return None
    try:
        state = refresh_pipeline_state(V2_RUN_ROOT / job_id)
    except Exception as exc:  # noqa: BLE001 - expose persisted failure cleanly
        state = {"status": "failed", "error": str(exc)}
    job = V2Job(
        job_id=job_id,
        status=str(state.get("status", "ready")),
        operation="idle",
        message="Recovered saved pipeline",
        progress=1.0,
        pipeline=state,
        error=state.get("error"),
    )
    with V2_JOBS_LOCK:
        return V2_JOBS.setdefault(job_id, job)


def execute_v2_operation(job_id: str, operation: str) -> None:
    job = recover_v2_job(job_id)
    if job is None:
        return

    def progress(message: str, fraction: float) -> None:
        if job.cancel_event.is_set():
            raise RuntimeError("Version 2 pipeline cancelled")
        with V2_JOBS_LOCK:
            job.message = message
            job.progress = max(0.0, min(1.0, float(fraction)))

    lock_acquired = False
    try:
        while not ANALYSIS_LOCK.acquire(timeout=0.2):
            if job.cancel_event.is_set():
                raise RuntimeError("Version 2 pipeline cancelled")
            with V2_JOBS_LOCK:
                job.status = "queued"
                job.message = "Waiting for current analysis"
        lock_acquired = True
        with V2_JOBS_LOCK:
            job.busy = True
            job.operation = operation
            job.status = "running"
            job.message = {
                "initialize": "Starting the complete pipeline",
                "step": "Advancing one decision stage",
                "advance": "Running to the next structural checkpoint",
                "submit": "Submitting the current AlphaPulldown round",
                "collect": "Collecting structural results and continuing",
                "trace": "Refreshing the development trace",
            }[operation]
            job.progress = 0.0
            job.error = None
        if operation == "initialize":
            state = initialize_end_to_end_pipeline(
                job.configuration,
                job.run_directory,
                project_root=EVIDENCE_ROOT,
                progress=progress,
                cancel_requested=job.cancel_event.is_set,
            )
        elif operation == "step":
            state = step_pipeline(job.run_directory)
        elif operation == "advance":
            state = step_pipeline(job.run_directory, until_checkpoint=True)
        elif operation == "submit":
            state = submit_pipeline_structural_round(job.run_directory)
        elif operation == "collect":
            collect_pipeline_structural_scores(job.run_directory)
            state = step_pipeline(job.run_directory, until_checkpoint=True)
        elif operation == "trace":
            render_pipeline_trace(job.run_directory)
            state = refresh_pipeline_state(job.run_directory)
        else:
            raise ValueError(f"unknown Version 2 operation: {operation}")
        with V2_JOBS_LOCK:
            job.pipeline = state
            job.status = str(state.get("status", "ready"))
            job.message = (
                "Complete"
                if state.get("status") == "complete"
                else "Ready for the next action"
            )
            job.progress = 1.0
    except Exception as exc:  # noqa: BLE001 - HTTP boundary reports scientific failures
        traceback.print_exc()
        cancelled = job.cancel_event.is_set()
        with V2_JOBS_LOCK:
            job.status = "cancelled" if cancelled else "failed"
            job.message = "Cancelled" if cancelled else f"{operation.capitalize()} failed"
            job.error = None if cancelled else str(exc)
        if job.run_directory.joinpath("pipeline_state.json").is_file():
            try:
                with V2_JOBS_LOCK:
                    job.pipeline = refresh_pipeline_state(job.run_directory)
            except Exception:
                pass
    finally:
        with V2_JOBS_LOCK:
            job.busy = False
        if lock_acquired:
            ANALYSIS_LOCK.release()


def start_v2_operation(job: V2Job, operation: str) -> None:
    with V2_JOBS_LOCK:
        if job.busy:
            raise RuntimeError("this Version 2 run already has an active operation")
        job.busy = True
        job.operation = operation
        job.status = "queued"
        job.message = "Queued"
        job.cancel_event.clear()
    threading.Thread(
        target=execute_v2_operation,
        args=(job.job_id, operation),
        daemon=True,
    ).start()


def execute_job(job_id: str) -> None:
    with JOBS_LOCK:
        job = JOBS[job_id]

    def check_cancel() -> None:
        if job.cancel_event.is_set():
            raise WorkflowCancelled("Analysis cancelled by user")

    def progress(message: str, fraction: float) -> None:
        check_cancel()
        with JOBS_LOCK:
            current = JOBS[job_id]
            current.message = message
            current.progress = max(0.0, min(1.0, float(fraction)))

    lock_acquired = False
    try:
        while not ANALYSIS_LOCK.acquire(timeout=0.2):
            check_cancel()
            with JOBS_LOCK:
                current = JOBS[job_id]
                current.status = "queued"
                current.message = "Waiting for current analysis"
        lock_acquired = True
        check_cancel()
        with JOBS_LOCK:
            current = JOBS[job_id]
            current.status = "running"
            current.message = "Starting"
        result = run_workflow(
            job.configuration,
            project_root=PROJECT_ROOT,
            run_id=job_id,
            progress=progress,
            cancel_requested=job.cancel_event.is_set,
        )
        check_cancel()
        with JOBS_LOCK:
            job = JOBS[job_id]
            job.status = "complete"
            job.message = "Complete"
            job.progress = 1.0
            job.preview = result.preview
            job.summary = result.summary
            job.output_directory = str(result.output_directory)
    except WorkflowCancelled:
        with JOBS_LOCK:
            job = JOBS[job_id]
            job.status = "cancelled"
            job.message = "Cancelled"
            job.error = None
    except Exception as exc:  # noqa: BLE001 - boundary must report scientific failures
        traceback.print_exc()
        with JOBS_LOCK:
            job = JOBS[job_id]
            job.status = "failed"
            job.message = "Analysis failed"
            job.error = str(exc)
    finally:
        if lock_acquired:
            ANALYSIS_LOCK.release()


class WorkbenchHandler(BaseHTTPRequestHandler):
    server_version = "GBIWorkbench/1.0"

    def log_message(self, format: str, *args: object) -> None:
        print(f"[{self.log_date_time_string()}] {format % args}")

    def send_json(self, payload: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(payload, default=str).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def read_json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        if length <= 0 or length > 2_000_000:
            raise ValueError("request body is empty or too large")
        raw = self.rfile.read(length)
        payload = json.loads(raw.decode("utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("request body must be a JSON object")
        return payload

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        if path == "/api/health":
            self.send_json({"status": "ok"})
            return
        if path == "/api/v2/config":
            try:
                self.send_json({
                    "registry": load_registry(EVIDENCE_ROOT),
                    "defaults": default_pipeline_configuration(EVIDENCE_ROOT),
                    "run_root": str(V2_RUN_ROOT),
                    "evidence_root": str(EVIDENCE_ROOT),
                })
            except (ValueError, FileNotFoundError) as exc:
                self.send_json({"error": str(exc)}, HTTPStatus.SERVICE_UNAVAILABLE)
            return
        if path.startswith("/api/v2/runs/"):
            parts = [part for part in path.split("/") if part]
            if len(parts) == 4:
                job = recover_v2_job(parts[3])
                if job is None:
                    self.send_json({"error": "Version 2 run not found"}, HTTPStatus.NOT_FOUND)
                else:
                    if not job.busy and job.run_directory.joinpath("pipeline_state.json").is_file():
                        try:
                            job.pipeline = refresh_pipeline_state(job.run_directory)
                            job.status = str(job.pipeline.get("status", job.status))
                        except Exception as exc:  # noqa: BLE001
                            job.error = str(exc)
                    self.send_json(public_v2_job(job))
                return
            if len(parts) == 6 and parts[4] == "files":
                self.send_v2_file(parts[3], parts[5])
                return
        if path == "/api/config":
            registry = load_registry(PROJECT_ROOT)
            cache = current_cache_counts(PROJECT_ROOT)
            with JOBS_LOCK:
                active_job = next(
                    (
                        public_job(job)
                        for job in JOBS.values()
                        if job.status not in TERMINAL_JOB_STATUSES
                    ),
                    None,
                )
            self.send_json(
                {
                    "registry": registry,
                    "defaults": default_configuration(registry),
                    "active_job": active_job,
                    "project": {
                        "name": "Graphical Bayesian Inference",
                        "node_candidates": 9170,
                        "seed_catalog_nodes": 891,
                        "seed_unique_pairs": 396495,
                        "cached_incremental_pairs": cache["incremental_pairs"],
                        "cached_pair_hypotheses": (
                            cache["seed_pairs"] + cache["incremental_pairs"]
                        ),
                        "cached_incremental_nodes": cache["profiled_nodes"],
                    },
                }
            )
            return
        if path == "/api/jobs/active":
            with JOBS_LOCK:
                active_job = next(
                    (
                        public_job(job)
                        for job in JOBS.values()
                        if job.status not in TERMINAL_JOB_STATUSES
                    ),
                    None,
                )
            self.send_json({"active_job": active_job})
            return
        if path.startswith("/api/jobs/"):
            parts = [part for part in path.split("/") if part]
            if len(parts) == 3:
                job_id = parts[2]
                with JOBS_LOCK:
                    job = JOBS.get(job_id)
                    payload = public_job(job) if job else None
                if payload is None:
                    self.send_json({"error": "job not found"}, HTTPStatus.NOT_FOUND)
                else:
                    self.send_json(payload)
                return
            if len(parts) == 5 and parts[3] == "files":
                self.send_job_file(parts[2], parts[4])
                return
            if len(parts) == 5 and parts[3] == "evidence":
                self.send_evidence_inspection(parts[2], parts[4], parsed.query)
                return
        if path.startswith("/api/"):
            self.send_json({"error": "not found"}, HTTPStatus.NOT_FOUND)
            return
        self.send_static(path)

    def send_evidence_inspection(
        self,
        job_id: str,
        evidence_kind: str,
        raw_query: str,
    ) -> None:
        with JOBS_LOCK:
            job = JOBS.get(job_id)
            status = job.status if job else None
            output_directory = (
                Path(job.output_directory)
                if job and job.output_directory
                else None
            )
        if job is None:
            self.send_json({"error": "job not found"}, HTTPStatus.NOT_FOUND)
            return
        if status != "complete" or output_directory is None:
            self.send_json(
                {"error": "Evidence can be inspected only after the run is complete"},
                HTTPStatus.CONFLICT,
            )
            return
        query = parse_qs(raw_query, keep_blank_values=True)

        def one(name: str) -> str:
            value = query.get(name, [""])[0].strip()
            if len(value) > 160:
                raise ValueError(f"{name} is too long")
            return value

        try:
            if evidence_kind == "node":
                payload = inspect_node_evidence(
                    output_directory,
                    one("symbol"),
                    project_root=PROJECT_ROOT,
                )
            elif evidence_kind == "edge":
                payload = inspect_edge_evidence(
                    output_directory,
                    one("node_a"),
                    one("node_b"),
                    project_root=PROJECT_ROOT,
                )
            else:
                self.send_json({"error": "unknown evidence kind"}, HTTPStatus.NOT_FOUND)
                return
            self.send_json(payload)
        except (ValueError, FileNotFoundError) as exc:
            self.send_json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)

    def do_POST(self) -> None:  # noqa: N802
        path = unquote(urlparse(self.path).path)
        parts = [part for part in path.split("/") if part]
        if path == "/api/v2/runs":
            try:
                payload = self.read_json()
                configuration = payload.get("configuration", payload)
                if not isinstance(configuration, dict):
                    raise ValueError("configuration must be an object")
                job_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8]
                job = V2Job(job_id=job_id, configuration=configuration)
                with V2_JOBS_LOCK:
                    V2_JOBS[job_id] = job
                start_v2_operation(job, "initialize")
                self.send_json(public_v2_job(job), HTTPStatus.ACCEPTED)
            except (ValueError, json.JSONDecodeError, RuntimeError) as exc:
                self.send_json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)
            return
        if len(parts) == 5 and parts[:3] == ["api", "v2", "runs"]:
            job = recover_v2_job(parts[3])
            action = parts[4]
            if job is None:
                self.send_json({"error": "Version 2 run not found"}, HTTPStatus.NOT_FOUND)
                return
            if action == "cancel":
                with V2_JOBS_LOCK:
                    job.cancel_event.set()
                    if job.busy:
                        job.status = "cancelling"
                        job.message = "Cancellation requested"
                self.send_json(public_v2_job(job), HTTPStatus.ACCEPTED)
                return
            if action not in {"step", "advance", "submit", "collect", "trace"}:
                self.send_json({"error": "unknown Version 2 action"}, HTTPStatus.NOT_FOUND)
                return
            try:
                start_v2_operation(job, action)
                self.send_json(public_v2_job(job), HTTPStatus.ACCEPTED)
            except RuntimeError as exc:
                self.send_json({"error": str(exc)}, HTTPStatus.CONFLICT)
            return
        if len(parts) == 4 and parts[:2] == ["api", "jobs"] and parts[3] == "cancel":
            job_id = parts[2]
            with JOBS_LOCK:
                job = JOBS.get(job_id)
                if job is None:
                    payload = None
                elif job.status not in TERMINAL_JOB_STATUSES:
                    job.cancel_event.set()
                    job.status = "cancelling"
                    job.message = "Cancellation requested"
                    payload = public_job(job)
                else:
                    payload = public_job(job)
            if payload is None:
                self.send_json({"error": "job not found"}, HTTPStatus.NOT_FOUND)
            else:
                self.send_json(payload, HTTPStatus.ACCEPTED)
            return
        if path != "/api/jobs":
            self.send_json({"error": "not found"}, HTTPStatus.NOT_FOUND)
            return
        try:
            payload = self.read_json()
            configuration = payload.get("configuration", payload)
            if not isinstance(configuration, dict):
                raise ValueError("configuration must be an object")
            job_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:8]
            job = Job(job_id=job_id, configuration=configuration)
            with JOBS_LOCK:
                active = next(
                    (
                        current
                        for current in JOBS.values()
                        if current.status not in TERMINAL_JOB_STATUSES
                    ),
                    None,
                )
                if active is None:
                    JOBS[job_id] = job
                    conflict = None
                else:
                    conflict = public_job(active)
            if conflict is not None:
                self.send_json(
                    {
                        "error": (
                            "Another analysis is already active. Cancel it before "
                            "starting a new run."
                        ),
                        "active_job": conflict,
                    },
                    HTTPStatus.CONFLICT,
                )
                return
            thread = threading.Thread(target=execute_job, args=(job_id,), daemon=True)
            self.send_json(public_job(job), HTTPStatus.ACCEPTED)
            thread.start()
        except (ValueError, json.JSONDecodeError) as exc:
            self.send_json({"error": str(exc)}, HTTPStatus.BAD_REQUEST)

    def send_v2_file(self, job_id: str, requested_name: str) -> None:
        job = recover_v2_job(job_id)
        if job is None:
            self.send_json({"error": "Version 2 run not found"}, HTTPStatus.NOT_FOUND)
            return
        allowed = {
            "trace": job.run_directory / "backward_search" / "development_trace.html",
            "pipeline-state": job.run_directory / "pipeline_state.json",
            "selected-nodes": job.run_directory / "selected_nodes.tsv",
            "graph-nodes": job.run_directory / "graph_nodes.tsv",
            "configuration": job.run_directory / "submitted_configuration.json",
        }
        file_path = allowed.get(requested_name)
        if file_path is None or not file_path.is_file():
            self.send_json({"error": "file not found"}, HTTPStatus.NOT_FOUND)
            return
        mime_type = mimetypes.guess_type(file_path.name)[0] or "application/octet-stream"
        body = file_path.read_bytes()
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", f"{mime_type}; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        if requested_name != "trace":
            self.send_header("Content-Disposition", f'attachment; filename="{file_path.name}"')
        self.end_headers()
        self.wfile.write(body)

    def send_job_file(self, job_id: str, requested_name: str) -> None:
        safe_requested = Path(requested_name).name
        if safe_requested != requested_name:
            self.send_json({"error": "invalid file name"}, HTTPStatus.BAD_REQUEST)
            return
        with JOBS_LOCK:
            job = JOBS.get(job_id)
            output_directory = Path(job.output_directory) if job and job.output_directory else None
        if output_directory is None:
            self.send_json({"error": "job output is not available"}, HTTPStatus.NOT_FOUND)
            return
        file_path = (output_directory / safe_requested).resolve()
        if file_path.parent != output_directory.resolve() or not file_path.is_file():
            self.send_json({"error": "file not found"}, HTTPStatus.NOT_FOUND)
            return
        mime_type = mimetypes.guess_type(file_path.name)[0] or "application/octet-stream"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", mime_type)
        self.send_header("Content-Length", str(file_path.stat().st_size))
        self.send_header("Content-Disposition", f'attachment; filename="{file_path.name}"')
        self.end_headers()
        with file_path.open("rb") as handle:
            while block := handle.read(1024 * 1024):
                self.wfile.write(block)

    def send_static(self, requested_path: str) -> None:
        relative = requested_path.lstrip("/") or "index.html"
        file_path = (WEB_ROOT / relative).resolve()
        if WEB_ROOT.resolve() not in file_path.parents and file_path != WEB_ROOT.resolve():
            self.send_error(HTTPStatus.FORBIDDEN)
            return
        if file_path.is_dir():
            file_path = file_path / "index.html"
        if not file_path.is_file():
            file_path = WEB_ROOT / "index.html"
        body = file_path.read_bytes()
        mime_type = mimetypes.guess_type(file_path.name)[0] or "application/octet-stream"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", f"{mime_type}; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(body)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    server = ThreadingHTTPServer((args.host, args.port), WorkbenchHandler)
    url = f"http://{args.host}:{server.server_port}"
    print(f"Graphical Bayesian Inference Version 2 workbench: {url}/v2.html")
    if not args.no_browser:
        threading.Timer(0.5, lambda: webbrowser.open(f"{url}/v2.html")).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
