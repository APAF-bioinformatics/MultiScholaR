#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time
from datetime import datetime, timezone


DEFAULT_INTERVAL_SECONDS = 300


STATIC_REMINDER = [
    "Codex-in-chat is the only outer supervisor.",
    "Do not run an autonomous repo-level supervisor.",
    "Passive pollers may record status and reminder context only.",
    "After each bounded inner-loop checkpoint: poll, inspect, verify, decide, then relaunch manually if safe.",
    "Stop only for a genuine strategic blocker.",
    "Commit only at module boundaries and keep production commits scoped to prod files.",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def run_status(repo_root: pathlib.Path) -> dict:
    completed = subprocess.run(
        [sys.executable, "tools/refactor/stabilization-loop.py", "status"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(completed.stdout)


def write_text(path: pathlib.Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def append_text(path: pathlib.Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(content)


def render_reminder(status: dict) -> str:
    current = status.get("currentItem") or {}
    runtime = status.get("runtime") or {}
    summary = status.get("summary") or {}
    last_progress = current.get("lastProgress") or {}
    target_meta = last_progress.get("target") or {}

    lines = [
        "# Stabilization Reminder",
        "",
        f"- Updated: `{now_iso()}`",
        "",
        "## Static Rules",
    ]
    for reminder in STATIC_REMINDER:
        lines.append(f"- {reminder}")

    lines.extend(
        [
            "",
            "## Current Loop State",
            f"- Loop status: `{summary.get('status')}`",
            f"- Runtime status: `{runtime.get('status')}`",
            f"- Current item: `{current.get('id')}`",
            f"- Target: `{current.get('targetPath')}`",
            f"- Last checkpoint: `{current.get('lastCheckpoint')}`",
            f"- Last summary: `{current.get('lastSummary')}`",
        ]
    )

    if target_meta:
        lines.append(f"- Next classifier step: `{target_meta.get('next_step')}`")

    lines.extend(
        [
            "",
            "## Immediate Reminder",
            "- If the inner loop is idle after a completed checkpoint, inspect the changed files and rerun the focused gate before launching the next inner loop.",
            "- If the target family reaches a true module boundary, commit the production family files before moving on.",
            "- If a blocker appears, distinguish ordinary verification work from a real strategic blocker before stopping.",
            "",
        ]
    )
    return "\n".join(lines)


def run_poller(args) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    interval = max(5, int(args.interval_seconds))
    runtime_dir = repo_root / ".god-module-stabilization"
    reminder_path = runtime_dir / "context-reminder.md"
    state_log = runtime_dir / "poll.log"
    state_json = runtime_dir / "poll-state.json"

    while True:
        status = run_status(repo_root)
        snapshot = {
            "updatedAt": now_iso(),
            "status": status,
            "staticReminder": STATIC_REMINDER,
        }
        write_text(state_json, json.dumps(snapshot, indent=2) + "\n")
        write_text(reminder_path, render_reminder(status))
        append_text(
            state_log,
            f"=== {now_iso()} ===\n{json.dumps(status, indent=2)}\n\n",
        )
        time.sleep(interval)


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--repo-root", default=".")
    run_parser.add_argument("--interval-seconds", type=int, default=DEFAULT_INTERVAL_SECONDS)

    args = parser.parse_args()
    return run_poller(args)


if __name__ == "__main__":
    raise SystemExit(main())
