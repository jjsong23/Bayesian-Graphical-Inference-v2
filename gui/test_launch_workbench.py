"""Tests for the detached local-workbench launcher."""

from __future__ import annotations

import json
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from launch_workbench import is_healthy, workbench_url


class _ConfigHandler(BaseHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        return

    def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API
        payload = json.dumps({"defaults": {}, "registry": {}}).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)


class LaunchWorkbenchTests(unittest.TestCase):
    def test_health_check_accepts_workbench_configuration_payload(self) -> None:
        server = ThreadingHTTPServer(("127.0.0.1", 0), _ConfigHandler)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            url = workbench_url("127.0.0.1", server.server_port)
            self.assertTrue(is_healthy(url, timeout=2.0))
        finally:
            server.shutdown()
            server.server_close()
            thread.join(timeout=2.0)

    def test_health_check_rejects_unreachable_port(self) -> None:
        self.assertFalse(is_healthy("http://127.0.0.1:1", timeout=0.1))


if __name__ == "__main__":
    unittest.main()
