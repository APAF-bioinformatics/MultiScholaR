#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import time
from datetime import datetime, timezone


DEFAULT_INTERVAL_SECONDS = 30
DEFAULT_SETTLE_SECONDS = 5


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: pathlib.Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def append_log(path: pathlib.Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(f"[{now_iso()}] {message}\n")


def run_json(cmd: list[str], cwd: pathlib.Path, allow_failure: bool = False) -> dict:
    completed = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0 and not allow_failure:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    stdout = completed.stdout.strip()
    if not stdout:
        return {
            "ok": completed.returncode == 0,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }
    try:
        data = json.loads(stdout)
    except json.JSONDecodeError:
        data = {
            "ok": completed.returncode == 0,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }
    data.setdefault("returncode", completed.returncode)
    data.setdefault("stderr", completed.stderr)
    data.setdefault("stdout", completed.stdout)
    return data


def run_text(cmd: list[str], cwd: pathlib.Path, allow_failure: bool = False) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if completed.returncode != 0 and not allow_failure:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return completed


def loop_status(repo_root: pathlib.Path) -> dict:
    return run_json([sys.executable, "tools/refactor/stabilization-loop.py", "status"], repo_root)


def read_loop_state_from_status(status: dict) -> dict:
    state_path_value = status.get("statePath")
    if not state_path_value:
        return {}
    state_path = pathlib.Path(state_path_value)
    if not state_path.exists():
        return {}
    try:
        return json.loads(state_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def scoped_items(status: dict) -> list[dict]:
    state = read_loop_state_from_status(status)
    items = state.get("items") or []
    prefixes = state.get("scopePrefixes") or []
    if not prefixes:
        return list(items)

    def in_scope(item: dict) -> bool:
        target_path = item.get("targetPath") or ""
        work_target_path = item.get("workTargetPath") or ""
        return any(
            target_path.startswith(prefix) or work_target_path.startswith(prefix)
            for prefix in prefixes
        )

    return [item for item in items if in_scope(item)]


def scoped_open_items(status: dict) -> list[dict]:
    return [
        item
        for item in scoped_items(status)
        if item.get("status") in {"pending", "in_progress"}
    ]


def loop_run(
    repo_root: pathlib.Path,
    target: str | None = None,
    scope_prefixes: list[str] | None = None,
    *,
    clear_stop_file: bool = False,
) -> dict:
    cmd = [sys.executable, "tools/refactor/stabilization-loop.py", "run", "--iteration-limit", "1"]
    if clear_stop_file:
        cmd.append("--clear-stop-file")
    if target:
        cmd.extend(["--target", target])
    for prefix in scope_prefixes or []:
        cmd.extend(["--scope-prefix", prefix])
    return run_json(cmd, repo_root, allow_failure=True)


def loop_run_with_max(
    repo_root: pathlib.Path,
    max_iterations: int,
    target: str | None = None,
    scope_prefixes: list[str] | None = None,
    *,
    clear_stop_file: bool = False,
) -> dict:
    cmd = [
        sys.executable,
        "tools/refactor/stabilization-loop.py",
        "run",
        "--iteration-limit",
        "1",
        "--max-iterations",
        str(max_iterations),
    ]
    if clear_stop_file:
        cmd.append("--clear-stop-file")
    if target:
        cmd.extend(["--target", target])
    for prefix in scope_prefixes or []:
        cmd.extend(["--scope-prefix", prefix])
    return run_json(cmd, repo_root, allow_failure=True)


def git_status_porcelain(repo_root: pathlib.Path) -> list[str]:
    completed = run_text(["git", "status", "--short"], repo_root)
    return [line for line in completed.stdout.splitlines() if line.strip()]


def description_changed(repo_root: pathlib.Path) -> bool:
    completed = run_text(["git", "status", "--short", "--", "DESCRIPTION"], repo_root)
    return bool(completed.stdout.strip())


def is_production_relpath(path_str: str) -> bool:
    return path_str == "DESCRIPTION" or (path_str.startswith("R/") and path_str.endswith(".R"))


def to_repo_relpath(repo_root: pathlib.Path, path_str: str) -> str | None:
    path = pathlib.Path(path_str)
    if not path.is_absolute():
        rel = path.as_posix()
        return rel if is_production_relpath(rel) else None
    try:
        rel = path.relative_to(repo_root).as_posix()
    except ValueError:
        return None
    return rel if is_production_relpath(rel) else None


def dirty_production_paths(repo_root: pathlib.Path) -> list[str]:
    lines = git_status_porcelain(repo_root)
    dirty: list[str] = []
    seen: set[str] = set()
    for line in lines:
        entry = line[3:] if len(line) > 3 else ""
        if " -> " in entry:
            entry = entry.split(" -> ", 1)[1]
        rel = entry.strip()
        if rel and is_production_relpath(rel) and rel not in seen:
            dirty.append(rel)
            seen.add(rel)
    return dirty


def scoped_dirty_production_paths(repo_root: pathlib.Path, status: dict) -> list[str]:
    state = read_loop_state_from_status(status)
    prefixes = state.get("scopePrefixes") or []
    dirty = dirty_production_paths(repo_root)
    if not prefixes:
        return dirty
    scoped: list[str] = []
    for rel in dirty:
        abs_path = (repo_root / rel).resolve()
        if any(str(abs_path).startswith(prefix) for prefix in prefixes):
            scoped.append(rel)
    return scoped


def current_family_files(repo_root: pathlib.Path, item: dict) -> list[str]:
    progress = item.get("lastProgress") or {}
    family = progress.get("family") or {}
    helper_files = family.get("helper_files") or []
    last_files_changed = item.get("lastFilesChanged") or []
    files = []
    target = item.get("targetPath")
    if target:
        files.append(target)
    for helper in helper_files:
        if helper not in files:
            files.append(helper)
    for changed_file in last_files_changed:
        rel = to_repo_relpath(repo_root, changed_file)
        if not rel:
            continue
        changed_path = str((repo_root / rel).resolve())
        if changed_path not in files:
            files.append(changed_path)
    return files


def stage_commit_scope(repo_root: pathlib.Path, item: dict) -> list[str]:
    files = current_family_files(repo_root, item)
    dirty = set(dirty_production_paths(repo_root))
    paths_to_add: list[str] = []
    seen: set[str] = set()
    for path_str in files:
        rel = to_repo_relpath(repo_root, path_str)
        if not rel or rel in seen:
            continue
        candidate = repo_root / rel
        if rel in dirty or candidate.exists():
            paths_to_add.append(rel)
            seen.add(rel)
    if "DESCRIPTION" in dirty and "DESCRIPTION" not in seen:
        paths_to_add.append("DESCRIPTION")
    if not paths_to_add:
        return []
    run_text(["git", "add", "--"] + paths_to_add, repo_root)
    return paths_to_add


def has_staged_changes(repo_root: pathlib.Path) -> bool:
    completed = run_text(["git", "diff", "--cached", "--name-only"], repo_root)
    return bool(completed.stdout.strip())


def commit_message_for_item(item: dict) -> str:
    title = item.get("title") or item.get("id") or "stabilization target"
    target_path = item.get("targetPath") or ""
    target_name = pathlib.Path(target_path).stem if target_path else title
    return f"Stabilize {target_name}"


def commit_item(repo_root: pathlib.Path, item: dict, log_path: pathlib.Path) -> bool:
    staged = stage_commit_scope(repo_root, item)
    if not staged:
        append_log(log_path, f"commit skipped for {item.get('id')}: no production files to stage")
        return False
    if not has_staged_changes(repo_root):
        append_log(log_path, f"commit skipped for {item.get('id')}: staged scope produced no diff")
        return False
    message = commit_message_for_item(item)
    completed = run_text(["git", "commit", "-m", message], repo_root, allow_failure=True)
    append_log(
        log_path,
        f"git commit for {item.get('id')} returned {completed.returncode}: {message}",
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"git commit failed for {item.get('id')}:\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return True


def next_target_override(status: dict) -> str | None:
    current_item = status.get("currentItem") or {}
    item_status = current_item.get("status")
    target = current_item.get("targetPath")
    if target and item_status in {"pending", "in_progress"}:
        return target
    for item in scoped_open_items(status):
        target = item.get("targetPath")
        if target:
            return target
    return None


def should_extend_iteration_cap(status: dict) -> bool:
    summary = status.get("summary") or {}
    return summary.get("status") == "paused_max_iterations" and bool(scoped_open_items(status))


def should_commit_boundary(status: dict, supervisor_state: dict) -> bool:
    current_item = status.get("currentItem") or {}
    if current_item.get("status") != "done":
        return False
    committed = set(supervisor_state.get("committedItems") or [])
    return current_item.get("id") not in committed


def is_terminal_repo_state(repo_root: pathlib.Path, status: dict) -> bool:
    summary = status.get("summary") or {}
    return (
        summary.get("status") == "completed"
        and not scoped_open_items(status)
        and not scoped_dirty_production_paths(repo_root, status)
    )


def is_deliberate_hold_state(status: dict) -> bool:
    summary = status.get("summary") or {}
    return summary.get("status") == "stopped"


def is_blocked_hold_state(status: dict) -> bool:
    summary = status.get("summary") or {}
    return summary.get("status") == "blocked"


def is_runtime_overrun(status: dict) -> bool:
    runtime = status.get("runtime") or {}
    return runtime.get("status") == "overrun"


def supervisor_snapshot(status: dict, *, mode: str, note: str | None = None) -> dict:
    return {
        "updatedAt": now_iso(),
        "mode": mode,
        "note": note,
        "loopSummary": status.get("summary"),
        "runtime": status.get("runtime"),
        "currentItem": status.get("currentItem"),
    }


def write_supervisor_state(path: pathlib.Path, *, status: dict, mode: str, note: str | None, committed_items: list[str]) -> None:
    snapshot = supervisor_snapshot(status, mode=mode, note=note)
    snapshot["committedItems"] = committed_items
    write_json(path, snapshot)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["run"])
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--interval-seconds", type=int, default=DEFAULT_INTERVAL_SECONDS)
    parser.add_argument("--settle-seconds", type=int, default=DEFAULT_SETTLE_SECONDS)
    parser.add_argument("--state-path", default=".god-module-stabilization/supervisor-state.json")
    parser.add_argument("--log-path", default=".god-module-stabilization/supervisor.log")
    parser.add_argument("--scope-prefix", action="append", default=[])
    parser.add_argument("--clear-stop-file", action="store_true")
    args = parser.parse_args()

    repo_root = pathlib.Path(args.repo_root).resolve()
    state_path = (repo_root / args.state_path).resolve()
    log_path = (repo_root / args.log_path).resolve()

    supervisor_state = {
        "startedAt": now_iso(),
        "committedItems": [],
    }
    append_log(log_path, "supervisor started")

    while True:
        status = loop_status(repo_root)
        runtime = status.get("runtime") or {}
        summary = status.get("summary") or {}
        current_item = status.get("currentItem") or {}

        write_supervisor_state(
            state_path,
            status=status,
            mode="observe",
            note="heartbeat",
            committed_items=supervisor_state["committedItems"],
        )

        if is_terminal_repo_state(repo_root, status):
            append_log(log_path, "repo stabilization loop reports completed")
            return 0

        scoped_dirty = scoped_dirty_production_paths(repo_root, status)
        if summary.get("status") == "completed" and scoped_dirty:
            dirty_note = ", ".join(scoped_dirty[:6])
            append_log(
                log_path,
                f"loop reports completed but scoped production files remain dirty; supervisor exiting into blocked hold ({dirty_note})",
            )
            write_supervisor_state(
                state_path,
                status=status,
                mode="blocked",
                note=f"completed queue but dirty scoped production files remain: {dirty_note}",
                committed_items=supervisor_state["committedItems"],
            )
            return 0

        if runtime.get("status") == "running":
            time.sleep(max(5, args.interval_seconds))
            continue

        if is_runtime_overrun(status):
            append_log(log_path, "loop runtime exceeded phase timeout; supervisor exiting into blocked hold")
            write_supervisor_state(
                state_path,
                status=status,
                mode="blocked",
                note="loop runtime exceeded phase timeout; manual diagnosis or explicit recovery required",
                committed_items=supervisor_state["committedItems"],
            )
            return 0

        if should_commit_boundary(status, supervisor_state):
            committed = commit_item(repo_root, current_item, log_path)
            if current_item["id"] not in supervisor_state["committedItems"]:
                supervisor_state["committedItems"].append(current_item["id"])
            note = f"committed {current_item['id']}" if committed else f"acknowledged existing boundary for {current_item['id']}"
            write_supervisor_state(
                state_path,
                status=status,
                mode="commit",
                note=note,
                committed_items=supervisor_state["committedItems"],
            )
            time.sleep(max(2, args.settle_seconds))
            continue

        if is_deliberate_hold_state(status):
            if args.clear_stop_file:
                target_override = next_target_override(status)
                append_log(
                    log_path,
                    f"clearing deliberate stop hold and relaunching {target_override or 'next-open-item'}",
                )
                write_supervisor_state(
                    state_path,
                    status=status,
                    mode="launch",
                    note=f"clearing stop hold for {target_override or 'next-open-item'}",
                    committed_items=supervisor_state["committedItems"],
                )
                result = loop_run(
                    repo_root,
                    target_override,
                    args.scope_prefix,
                    clear_stop_file=True,
                )
                append_log(
                    log_path,
                    f"clear-stop loop run returned code={result.get('returncode')} ok={result.get('ok')} target={target_override or 'next-open-item'}",
                )
                if result.get("returncode") not in (0, None):
                    write_supervisor_state(
                        state_path,
                        status=loop_status(repo_root),
                        mode="blocked",
                        note=f"clear-stop loop run failed for {target_override or 'next-open-item'}",
                        committed_items=supervisor_state["committedItems"],
                    )
                    return 1
                time.sleep(max(2, args.settle_seconds))
                continue
            append_log(log_path, "loop reports deliberate stopped hold state; supervisor exiting without relaunch")
            write_supervisor_state(
                state_path,
                status=status,
                mode="hold",
                note="loop is stopped; manual clear-stop-file resume required",
                committed_items=supervisor_state["committedItems"],
            )
            return 0

        if is_blocked_hold_state(status):
            append_log(log_path, "loop reports blocked hold state; supervisor exiting without relaunch")
            write_supervisor_state(
                state_path,
                status=status,
                mode="blocked",
                note="loop is blocked; manual diagnosis or explicit recovery required",
                committed_items=supervisor_state["committedItems"],
            )
            return 0

        if should_extend_iteration_cap(status):
            next_cap = max(int(summary.get("iteration") or 0) + 25, int(summary.get("maxIterations") or 0) + 25)
            target_override = next_target_override(status)
            append_log(
                log_path,
                f"extending loop max-iterations to {next_cap} for {target_override or 'next-open-item'}",
            )
            write_supervisor_state(
                state_path,
                status=status,
                mode="launch",
                note=f"extending max-iterations to {next_cap} for {target_override or 'next-open-item'}",
                committed_items=supervisor_state["committedItems"],
            )
            result = loop_run_with_max(repo_root, next_cap, target_override, args.scope_prefix)
            append_log(
                log_path,
                f"loop run with extended max returned code={result.get('returncode')} ok={result.get('ok')} target={target_override or 'next-open-item'}",
            )
            if result.get("returncode") not in (0, None):
                write_supervisor_state(
                    state_path,
                    status=loop_status(repo_root),
                    mode="blocked",
                    note=f"loop run command failed while extending max-iterations for {target_override or 'next-open-item'}",
                    committed_items=supervisor_state["committedItems"],
                )
                return 1
            time.sleep(max(2, args.settle_seconds))
            continue

        target_override = next_target_override(status)
        append_log(
            log_path,
            f"launching next loop iteration for {target_override or 'next-open-item'} from status={summary.get('status')}",
        )
        write_supervisor_state(
            state_path,
            status=status,
            mode="launch",
            note=f"launching {target_override or 'next-open-item'}",
            committed_items=supervisor_state["committedItems"],
        )
        result = loop_run(repo_root, target_override, args.scope_prefix)
        append_log(
            log_path,
            f"loop run returned code={result.get('returncode')} ok={result.get('ok')} target={target_override or 'next-open-item'}",
        )
        if result.get("returncode") not in (0, None):
            write_supervisor_state(
                state_path,
                status=loop_status(repo_root),
                mode="blocked",
                note=f"loop run command failed for {target_override or 'next-open-item'}",
                committed_items=supervisor_state["committedItems"],
            )
            return 1
        time.sleep(max(2, args.settle_seconds))


if __name__ == "__main__":
    raise SystemExit(main())
