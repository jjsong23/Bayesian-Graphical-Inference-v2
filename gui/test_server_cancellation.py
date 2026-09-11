"""Regression tests for cooperative workbench job cancellation."""

from __future__ import annotations

import threading
import time
import unittest
from unittest.mock import patch

import server


class ServerCancellationTests(unittest.TestCase):
    def setUp(self) -> None:
        with server.JOBS_LOCK:
            server.JOBS.clear()
        self.assertFalse(server.ANALYSIS_LOCK.locked())

    def tearDown(self) -> None:
        with server.JOBS_LOCK:
            server.JOBS.clear()
        if server.ANALYSIS_LOCK.locked():
            server.ANALYSIS_LOCK.release()

    def test_running_job_cancels_at_workflow_checkpoint(self) -> None:
        started = threading.Event()

        def cancellable_workflow(*args, cancel_requested=None, **kwargs):
            started.set()
            while not cancel_requested():
                time.sleep(0.005)
            raise server.WorkflowCancelled("cancelled in test")

        job = server.Job(job_id="running-cancel", configuration={})
        with server.JOBS_LOCK:
            server.JOBS[job.job_id] = job
        with patch.object(server, "run_workflow", side_effect=cancellable_workflow):
            worker = threading.Thread(target=server.execute_job, args=(job.job_id,))
            worker.start()
            self.assertTrue(started.wait(timeout=2))
            job.cancel_event.set()
            worker.join(timeout=2)

        self.assertFalse(worker.is_alive())
        self.assertEqual(job.status, "cancelled")
        self.assertEqual(job.message, "Cancelled")
        self.assertIsNone(job.error)

    def test_waiting_job_can_cancel_before_workflow_starts(self) -> None:
        job = server.Job(job_id="queued-cancel", configuration={})
        with server.JOBS_LOCK:
            server.JOBS[job.job_id] = job
        server.ANALYSIS_LOCK.acquire()
        with patch.object(server, "run_workflow") as mocked_workflow:
            worker = threading.Thread(target=server.execute_job, args=(job.job_id,))
            worker.start()
            job.cancel_event.set()
            worker.join(timeout=2)

        self.assertFalse(worker.is_alive())
        mocked_workflow.assert_not_called()
        self.assertEqual(job.status, "cancelled")


if __name__ == "__main__":
    unittest.main()
