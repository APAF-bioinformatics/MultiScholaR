#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import pathlib
import re
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone

from stabilization_schema import OpenAISchemaCompatibilityError, validate_openai_output_schema


STATE_SCHEMA_VERSION = 2
HEARTBEAT_INTERVAL_SECONDS = 1.0
STALE_HEARTBEAT_SECONDS = 15.0
DEFAULT_EXECUTOR_TIMEOUT_MS = 2700000
DEFAULT_OVERRIDES_PATH = "tools/refactor/stabilization-target-overrides.json"
DEFAULT_HANDOVER_EXCERPT_LINES = 40
DEFAULT_HANDOVER_EXCERPT_CHARS = 2500
DEFAULT_SUMMARY_EXCERPT_CHARS = 800
HEADING_RE = re.compile(r"^###\s+(\d+)\.\s+(.+?)\s*$")
R_LINK_RE = re.compile(r"\((/[^)]+/R/[^):]+\.R):\d+\)")
HANDOVER_LINK_RE = re.compile(r"\((/[^)]+/tools/refactor/HANDOVER[^):]+\.md):\d+\)")


@dataclass
class LoopItem:
    id: str
    bucketNumber: int
    title: str
    targetPath: str
    handoverPath: str | None
    status: str
    workTargetPath: str | None = None
    focusedGateCommands: list[list[str]] | None = None
    supplementalReviewCommands: list[list[str]] | None = None
    promptHints: list[str] | None = None
    startedAt: str | None = None
    completedAt: str | None = None
    lastCheckpoint: str | None = None
    lastSummary: str | None = None
    lastProgress: dict | None = None
    lastReview: dict | None = None
    lastFilesChanged: list[str] | None = None


class ScriptExecutionError(RuntimeError):
    def __init__(
        self,
        *,
        phase: str,
        returncode: int,
        stdout_text: str,
        stderr_text: str,
        stdout_path: str,
        stderr_path: str,
        command: list[str],
        timed_out: bool = False,
    ) -> None:
        self.phase = phase
        self.returncode = returncode
        self.stdout_text = stdout_text
        self.stderr_text = stderr_text
        self.stdout_path = stdout_path
        self.stderr_path = stderr_path
        self.command = command
        self.timed_out = timed_out
        super().__init__(stderr_text or stdout_text or f"{phase} exited with code {returncode}")


class ProtocolValidationError(RuntimeError):
    def __init__(self, *, phase: str, message: str) -> None:
        self.phase = phase
        super().__init__(message)


def now_utc() -> datetime:
    return datetime.now(timezone.utc)


def now_iso() -> str:
    return now_utc().isoformat()


def parse_iso(value: str | None) -> datetime | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


def read_json(path: pathlib.Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: pathlib.Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def append_jsonl(path: pathlib.Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(value) + "\n")


def read_tail(path_str: str | None, limit: int = 4000) -> str:
    if not path_str:
        return ""
    path = pathlib.Path(path_str)
    if not path.exists():
        return ""
    text = path.read_text(encoding="utf-8", errors="replace")
    return text[-limit:]


def truncate_text(text: str | None, max_chars: int) -> str:
    if not text:
        return ""
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 3].rstrip() + "..."


def pid_is_alive(pid: int | None) -> bool:
    if not isinstance(pid, int) or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def normalize_argv_matrix(value, *, field_name: str) -> list[list[str]]:
    if not isinstance(value, list):
        raise ProtocolValidationError(phase="protocol", message=f"{field_name} must be a list")
    normalized: list[list[str]] = []
    for entry in value:
        if not isinstance(entry, list) or not entry:
            raise ProtocolValidationError(
                phase="protocol",
                message=f"{field_name} entries must be non-empty argv arrays",
            )
        argv: list[str] = []
        for token in entry:
            if not isinstance(token, str) or not token:
                raise ProtocolValidationError(
                    phase="protocol",
                    message=f"{field_name} argv entries must be non-empty strings",
                )
            argv.append(token)
        normalized.append(argv)
    return normalized


def phase_state(active_run: dict | None, phase: str | None = None) -> dict | None:
    if not active_run:
        return None
    phase_name = phase or active_run.get("phase")
    if phase_name not in {"executor", "reviewer"}:
        return None
    return active_run.get(phase_name)


def parse_backlog(backlog_path: pathlib.Path) -> list[LoopItem]:
    lines = backlog_path.read_text(encoding="utf-8").splitlines()
    items: list[LoopItem] = []
    current_match: re.Match[str] | None = None
    current_block: list[str] = []

    def flush_bucket(match: re.Match[str] | None, block_lines: list[str]) -> None:
        if match is None:
            return
        block = "\n".join(block_lines)
        file_paths = list(dict.fromkeys(R_LINK_RE.findall(block)))
        if not file_paths:
            return
        handover_matches = list(dict.fromkeys(HANDOVER_LINK_RE.findall(block)))
        lower = block.lower()
        status = "pending"
        if "completed in live `r/`" in lower or "no longer a blocker" in lower:
            status = "done"
        elif "next active target" in lower:
            status = "in_progress"

        title = match.group(2).strip()
        slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
        items.append(
            LoopItem(
                id=f"{match.group(1)}-{slug}",
                bucketNumber=int(match.group(1)),
                title=title,
                targetPath=file_paths[0],
                handoverPath=handover_matches[0] if handover_matches else None,
                workTargetPath=None,
                focusedGateCommands=None,
                supplementalReviewCommands=None,
                promptHints=None,
                status=status,
            )
        )

    for line in lines:
        heading_match = HEADING_RE.match(line)
        if heading_match:
            flush_bucket(current_match, current_block)
            current_match = heading_match
            current_block = [line]
            continue
        if current_match is not None:
            current_block.append(line)
    flush_bucket(current_match, current_block)
    return items


def summarize_state(state: dict) -> dict:
    counts = {"pending": 0, "in_progress": 0, "done": 0, "blocked": 0, "failed": 0}
    items = [item for item in state["items"] if item_in_scope(state, item)]
    for item in items:
        counts[item["status"]] = counts.get(item["status"], 0) + 1
    runtime = build_runtime_snapshot(state)
    return {
        "schemaVersion": state["schemaVersion"],
        "status": state["status"],
        "runtimeStatus": runtime["status"],
        "iteration": state["iteration"],
        "maxIterations": state["maxIterations"],
        "currentItemId": state["currentItemId"],
        "itemStatus": counts,
        "totalItems": len(items),
    }


def item_in_scope(state: dict, item: dict) -> bool:
    prefixes = state.get("scopePrefixes") or []
    if not prefixes:
        return True
    target_path = item.get("targetPath") or ""
    work_target_path = item.get("workTargetPath") or ""
    return any(target_path.startswith(prefix) or work_target_path.startswith(prefix) for prefix in prefixes)


def next_open_item(state: dict) -> dict | None:
    current_item_id = state.get("currentItemId")
    if current_item_id:
        current_item = next((item for item in state["items"] if item["id"] == current_item_id), None)
        if current_item and current_item["status"] in {"in_progress", "pending"} and item_in_scope(state, current_item):
            return current_item
    for status in ("in_progress", "pending"):
        for item in state["items"]:
            if item["status"] == status and item_in_scope(state, item):
                return item
    return None


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def default_overrides_path(repo_root: pathlib.Path) -> pathlib.Path:
    return (repo_root / DEFAULT_OVERRIDES_PATH).resolve()


def load_target_overrides(repo_root: pathlib.Path, overrides_path: pathlib.Path | None) -> dict:
    path = overrides_path or default_overrides_path(repo_root)
    if not path.exists():
        return {}
    data = read_json(path)
    if not isinstance(data, dict):
        raise ProtocolValidationError(phase="protocol", message="target overrides file must contain an object")
    items = data.get("items", {})
    if not isinstance(items, dict):
        raise ProtocolValidationError(phase="protocol", message="target overrides items must be an object")
    return items


def normalize_override_path(repo_root: pathlib.Path, value: str | None) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ProtocolValidationError(phase="protocol", message="override paths must be non-empty strings")
    return str(resolve_repo_path(repo_root, value))


def apply_item_override(item: dict, override: dict, repo_root: pathlib.Path) -> bool:
    changed = False
    handover_path = normalize_override_path(repo_root, override.get("handoverPath"))
    if handover_path is not None and item.get("handoverPath") != handover_path:
        item["handoverPath"] = handover_path
        changed = True

    work_target_path = normalize_override_path(repo_root, override.get("workTargetPath"))
    if work_target_path is not None and item.get("workTargetPath") != work_target_path:
        item["workTargetPath"] = work_target_path
        changed = True

    if "focusedGateCommands" in override:
        focused_gate_commands = normalize_argv_matrix(
            override.get("focusedGateCommands", []),
            field_name="focusedGateCommands",
        )
        if item.get("focusedGateCommands") != focused_gate_commands:
            item["focusedGateCommands"] = focused_gate_commands
            changed = True

    if "supplementalReviewCommands" in override:
        supplemental_review_commands = normalize_argv_matrix(
            override.get("supplementalReviewCommands", []),
            field_name="supplementalReviewCommands",
        )
        if item.get("supplementalReviewCommands") != supplemental_review_commands:
            item["supplementalReviewCommands"] = supplemental_review_commands
            changed = True

    if "promptHints" in override:
        prompt_hints = ensure_string_list(override.get("promptHints", []), field_name="promptHints")
        if item.get("promptHints") != prompt_hints:
            item["promptHints"] = prompt_hints
            changed = True

    return changed


def apply_repo_overrides(state: dict, repo_root: pathlib.Path) -> bool:
    overrides_path_str = state.get("overridesPath")
    overrides_path = pathlib.Path(overrides_path_str).resolve() if isinstance(overrides_path_str, str) else default_overrides_path(repo_root)
    if state.get("overridesPath") != str(overrides_path):
        state["overridesPath"] = str(overrides_path)
        changed = True
    else:
        changed = False
    overrides = load_target_overrides(repo_root, overrides_path)
    for item in state.get("items", []):
        override = overrides.get(item.get("id"))
        if isinstance(override, dict):
            changed = apply_item_override(item, override, repo_root) or changed
    return changed


def find_item_by_override(state: dict, target_override: str | None, item_id_override: str | None) -> dict | None:
    if item_id_override:
        return next((item for item in state["items"] if item["id"] == item_id_override), None)
    if target_override:
        resolved_target = str(pathlib.Path(target_override).resolve())
        return next((item for item in state["items"] if pathlib.Path(item["targetPath"]).resolve() == pathlib.Path(resolved_target)), None)
    return None


def ensure_override_item(state: dict, target_override: str | None, item_id_override: str | None) -> dict | None:
    override_item = find_item_by_override(state, target_override, item_id_override)
    if override_item is not None:
        return override_item
    if target_override is None:
        return None

    target_path = str(pathlib.Path(target_override).resolve())
    synthetic_item = {
        "id": f"manual-{slugify(pathlib.Path(target_path).stem)}",
        "bucketNumber": 0,
        "title": f"Manual Target: {pathlib.Path(target_path).name}",
        "targetPath": target_path,
        "handoverPath": None,
        "workTargetPath": None,
        "focusedGateCommands": None,
        "supplementalReviewCommands": None,
        "promptHints": None,
        "status": "pending",
        "startedAt": None,
        "completedAt": None,
        "lastCheckpoint": None,
        "lastSummary": None,
        "lastProgress": None,
        "lastReview": None,
    }
    state["items"].append(synthetic_item)
    return synthetic_item


def queue_targets(state: dict, targets: list[str], *, item_ids: list[str] | None = None) -> list[dict]:
    queued: list[dict] = []
    item_ids = item_ids or []

    for item_id in item_ids:
        item = ensure_override_item(state, None, item_id)
        if item is None:
            raise ProtocolValidationError(
                phase="protocol",
                message=f"cannot queue unknown item id: {item_id}",
            )
        if item["status"] == "done":
            item["completedAt"] = None
            item["status"] = "pending"
        elif item["status"] not in {"pending", "in_progress"}:
            item["status"] = "pending"
        queued.append(item)

    for target in targets:
        item = ensure_override_item(state, target, None)
        if item is None:
            raise ProtocolValidationError(
                phase="protocol",
                message=f"cannot queue target: {target}",
            )
        if item["status"] == "done":
            item["completedAt"] = None
            item["status"] = "pending"
        elif item["status"] not in {"pending", "in_progress"}:
            item["status"] = "pending"
        queued.append(item)

    if any(item["status"] in {"pending", "in_progress"} for item in state["items"]):
        if state.get("status") in {"completed", "blocked", "failed", "stopped"}:
            state["status"] = "pending"
        if state.get("currentItemId") and not any(
            item["id"] == state["currentItemId"] and item["status"] == "in_progress"
            for item in state["items"]
        ):
            state["currentItemId"] = None
    state["updatedAt"] = now_iso()
    return queued


def run_status_snapshot(repo_root: pathlib.Path, backlog_path: pathlib.Path, target_path: str) -> dict:
    status_script = pathlib.Path(__file__).resolve().parent / "stabilization-status.py"
    completed = subprocess.run(
        [
            sys.executable,
            str(status_script),
            "--target",
            target_path,
            "--backlog",
            str(backlog_path),
            "--json",
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(completed.stdout)


def script_command(script_path: pathlib.Path, payload_path: pathlib.Path) -> list[str]:
    if script_path.suffix == ".py":
        return [sys.executable, str(script_path), str(payload_path)]
    if script_path.suffix in {".js", ".cjs", ".mjs"}:
        return ["/usr/bin/node", str(script_path), str(payload_path)]
    return [str(script_path), str(payload_path)]


def classify_runner_error(error_text: str) -> str:
    lowered = error_text.lower()
    session_init_markers = (
        "failed to create session",
        "failed to initialize session",
        "thread/start failed",
    )
    if any(marker in lowered for marker in session_init_markers):
        return "blocked"
    network_markers = (
        "failed to lookup address information",
        "temporary failure in name resolution",
        "dns error",
        "error sending request for url",
        "stream disconnected before completion",
    )
    if any(marker in lowered for marker in network_markers):
        return "blocked"
    protocol_markers = (
        "failed to parse",
        "protocol validation failed",
        "protocol schema validation failed",
        "invalid executor result",
        "invalid reviewer result",
        "invalid_json_schema",
        "response_format",
    )
    if any(marker in lowered for marker in protocol_markers):
        return "blocked"
    return "failed"


def classify_runner_reason_code(error_text: str, *, phase: str, timed_out: bool = False) -> str:
    if timed_out:
        return f"blocked_{phase}_timeout"
    lowered = error_text.lower()
    if (
        "failed to create session" in lowered
        or "failed to initialize session" in lowered
        or "thread/start failed" in lowered
    ):
        return f"blocked_{phase}_session_init"
    if "invalid_json_schema" in lowered or "protocol schema validation failed" in lowered:
        return "blocked_protocol_schema"
    if "protocol validation failed" in lowered or "invalid executor result" in lowered or "invalid reviewer result" in lowered:
        return "blocked_protocol"
    if "failed to lookup address information" in lowered or "temporary failure in name resolution" in lowered or "dns error" in lowered:
        return "blocked_network"
    status = classify_runner_error(error_text)
    if status == "blocked":
        return "blocked_runner"
    return status


def validate_local_executor_schema(runner_script: pathlib.Path) -> None:
    default_runner = (pathlib.Path(__file__).resolve().parent / "stabilization-codex-runner.py").resolve()
    if runner_script != default_runner:
        return
    schema_path = pathlib.Path(__file__).resolve().parent / "stabilization-executor-output.schema.json"
    try:
        schema = read_json(schema_path)
    except json.JSONDecodeError as error:
        raise ProtocolValidationError(
            phase="protocol",
            message=f"protocol schema validation failed: invalid JSON in {schema_path}: {error}",
        ) from error
    except FileNotFoundError as error:
        raise ProtocolValidationError(
            phase="protocol",
            message=f"protocol schema validation failed: schema file not found: {error.filename}",
        ) from error
    try:
        validate_openai_output_schema(schema)
    except OpenAISchemaCompatibilityError as error:
        raise ProtocolValidationError(
            phase="protocol",
            message=f"protocol schema validation failed: {error}",
        ) from error


def ensure_string_list(value, *, field_name: str) -> list[str]:
    if not isinstance(value, list):
        raise ProtocolValidationError(phase="protocol", message=f"{field_name} must be a list")
    normalized: list[str] = []
    for entry in value:
        if not isinstance(entry, str):
            raise ProtocolValidationError(phase="protocol", message=f"{field_name} entries must be strings")
        normalized.append(entry)
    return normalized


def production_files_from_changed_list(changed_files: list[str]) -> list[str]:
    normalized: list[str] = []
    seen: set[str] = set()
    for path_str in changed_files:
        if path_str == "DESCRIPTION" or (path_str.startswith("R/") and path_str.endswith(".R")):
            if path_str not in seen:
                normalized.append(path_str)
                seen.add(path_str)
    return normalized


def canonicalize_replay_commands(result: dict) -> tuple[dict, bool]:
    verification = result.get("verification")
    if isinstance(verification, dict):
        replay_commands = verification.get("replayCommands", [])
        if not isinstance(replay_commands, list):
            raise ProtocolValidationError(
                phase="protocol",
                message="verification.replayCommands must be a list",
            )
        canonical_commands: list[dict] = []
        for command in replay_commands:
            if not isinstance(command, dict):
                raise ProtocolValidationError(
                    phase="protocol",
                    message="verification.replayCommands entries must be objects",
                )
            argv = command.get("argv")
            if not isinstance(argv, list) or not argv:
                raise ProtocolValidationError(
                    phase="protocol",
                    message="verification.replayCommands entries must include a non-empty argv array",
                )
            if not all(isinstance(token, str) and token for token in argv):
                raise ProtocolValidationError(
                    phase="protocol",
                    message="verification.replayCommands argv entries must be non-empty strings",
                )
            canonical_entry = {"argv": argv}
            if command.get("label") is not None:
                canonical_entry["label"] = str(command["label"])
            if command.get("phase") is not None:
                canonical_entry["phase"] = str(command["phase"])
            if command.get("replayMode") is not None:
                replay_mode = str(command["replayMode"])
                if replay_mode not in {"verify", "operate"}:
                    raise ProtocolValidationError(
                        phase="protocol",
                        message="verification.replayCommands replayMode must be verify or operate when provided",
                    )
                canonical_entry["replayMode"] = replay_mode
            canonical_commands.append(canonical_entry)
        display = verification.get("display", [])
        if not isinstance(display, list) or not all(isinstance(entry, str) for entry in display):
            raise ProtocolValidationError(
                phase="protocol",
                message="verification.display must be a list of strings when provided",
            )
        return {
            "replayCommands": canonical_commands,
            "display": display,
            "legacyFallbackUsed": bool(verification.get("legacyFallbackUsed", False)),
        }, False

    tests_run = result.get("testsRun")
    if tests_run is None:
        raise ProtocolValidationError(
            phase="protocol",
            message="executor result must include verification.replayCommands or legacy testsRun",
        )
    normalized_tests = ensure_string_list(tests_run, field_name="testsRun")
    return {
        "replayCommands": [{"replayCommand": command} for command in normalized_tests],
        "display": normalized_tests,
        "legacyFallbackUsed": True,
    }, True


def validate_executor_result(result: dict) -> dict:
    if not isinstance(result, dict):
        raise ProtocolValidationError(phase="executor", message="invalid executor result: top-level value must be an object")
    normalized = dict(result)
    status = normalized.get("status")
    if status not in {"completed", "blocked", "failed"}:
        raise ProtocolValidationError(phase="executor", message="invalid executor result: status must be completed, blocked, or failed")
    target = normalized.get("target")
    checkpoint = normalized.get("checkpoint")
    target_state = normalized.get("targetStatusAfterIteration")
    summary = normalized.get("summary")
    if not isinstance(target, str) or not target.strip():
        raise ProtocolValidationError(phase="executor", message="invalid executor result: target must be a non-empty string")
    if not isinstance(checkpoint, str) or not checkpoint.strip():
        raise ProtocolValidationError(phase="executor", message="invalid executor result: checkpoint must be a non-empty string")
    if target_state not in {"in_progress", "done", "blocked"}:
        raise ProtocolValidationError(
            phase="executor",
            message="invalid executor result: targetStatusAfterIteration must be in_progress, done, or blocked",
        )
    if not isinstance(summary, str) or not summary.strip():
        raise ProtocolValidationError(phase="executor", message="invalid executor result: summary must be a non-empty string")
    normalized["filesChanged"] = ensure_string_list(normalized.get("filesChanged", []), field_name="filesChanged")
    normalized["notes"] = ensure_string_list(normalized.get("notes", []), field_name="notes")
    verification, used_legacy = canonicalize_replay_commands(normalized)
    normalized["verification"] = verification
    if used_legacy:
        normalized["testsRun"] = verification["display"]
    elif "testsRun" in normalized and normalized["testsRun"] is not None:
        normalized["testsRun"] = ensure_string_list(normalized["testsRun"], field_name="testsRun")
    else:
        normalized["testsRun"] = verification["display"]
    if "reasonCode" in normalized and normalized["reasonCode"] is not None and not isinstance(normalized["reasonCode"], str):
        raise ProtocolValidationError(phase="executor", message="invalid executor result: reasonCode must be a string when provided")
    return normalized


def validate_reviewer_result(result: dict) -> dict:
    if not isinstance(result, dict):
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: top-level value must be an object")
    normalized = dict(result)
    status = normalized.get("status")
    if status not in {"approved", "needs_changes", "blocked"}:
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: status must be approved, needs_changes, or blocked")
    summary = normalized.get("summary")
    if not isinstance(summary, str) or not summary.strip():
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: summary must be a non-empty string")
    normalized["testsRun"] = ensure_string_list(normalized.get("testsRun", []), field_name="testsRun")
    normalized["issues"] = ensure_string_list(normalized.get("issues", []), field_name="issues")
    normalized["notes"] = ensure_string_list(normalized.get("notes", []), field_name="notes")
    if not isinstance(normalized.get("testFailures", []), list):
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: testFailures must be a list")
    verification = normalized.get("verification", {})
    if verification is None:
        verification = {}
    if not isinstance(verification, dict):
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: verification must be an object when provided")
    replayed_commands = verification.get("replayedCommands", [])
    if not isinstance(replayed_commands, list):
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: verification.replayedCommands must be a list")
    normalized["verification"] = verification
    target_state_override = normalized.get("targetStateOverride")
    if target_state_override is not None and target_state_override not in {"in_progress", "done", "blocked"}:
        raise ProtocolValidationError(
            phase="reviewer",
            message="invalid reviewer result: targetStateOverride must be in_progress, done, or blocked when provided",
        )
    if "reasonCode" in normalized and normalized["reasonCode"] is not None and not isinstance(normalized["reasonCode"], str):
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: reasonCode must be a string when provided")
    status_snapshot = normalized.get("statusSnapshot")
    if status_snapshot is not None and not isinstance(status_snapshot, dict):
        raise ProtocolValidationError(phase="reviewer", message="invalid reviewer result: statusSnapshot must be an object when provided")
    return normalized


def blocked_reason_code(executor_result: dict | None, reviewer_result: dict | None, fallback_status: str = "blocked") -> str:
    if reviewer_result:
        reviewer_reason = reviewer_result.get("reasonCode")
        if reviewer_result.get("status") != "approved":
            return str(reviewer_reason or f"{reviewer_result['status']}_reviewer")
        if reviewer_reason and reviewer_reason != "approved":
            return str(reviewer_reason)
    if executor_result and executor_result.get("status") != "completed":
        return str(executor_result.get("reasonCode") or f"{executor_result['status']}_executor")
    if executor_result and executor_result.get("targetStatusAfterIteration") == "blocked":
        return str(executor_result.get("reasonCode") or "blocked_executor_target")
    return fallback_status


def runtime_dir_for_state(state_path: pathlib.Path) -> pathlib.Path:
    return state_path.parent.resolve()


def migrate_state(state: dict, state_path: pathlib.Path) -> tuple[dict, bool]:
    changed = False
    if state.get("schemaVersion") != STATE_SCHEMA_VERSION:
        state["schemaVersion"] = STATE_SCHEMA_VERSION
        changed = True
    runtime_dir = str(runtime_dir_for_state(state_path))
    if state.get("runtimeDir") != runtime_dir:
        state["runtimeDir"] = runtime_dir
        changed = True
    if "activeRun" not in state:
        state["activeRun"] = None
        changed = True
    if "history" not in state:
        state["history"] = []
        changed = True
    if not isinstance(state.get("scopePrefixes"), list):
        state["scopePrefixes"] = []
        changed = True
    repo_root_str = state.get("repoRoot")
    if isinstance(repo_root_str, str) and repo_root_str:
        repo_root = pathlib.Path(repo_root_str).resolve()
        default_path = str(default_overrides_path(repo_root))
        if state.get("overridesPath") != default_path:
            state["overridesPath"] = default_path
            changed = True
    return state, changed


def build_runtime_snapshot(state: dict) -> dict:
    active_run = state.get("activeRun")
    if not active_run:
        return {
            "status": "idle",
            "reason": "no_active_run",
            "loopPid": None,
            "childPid": None,
            "phase": None,
            "heartbeatAt": None,
            "heartbeatAgeSeconds": None,
        }

    heartbeat_at = parse_iso(active_run.get("heartbeatAt"))
    heartbeat_age = None
    if heartbeat_at is not None:
        heartbeat_age = round((now_utc() - heartbeat_at).total_seconds(), 1)

    current_phase = active_run.get("phase")
    current_phase_state = phase_state(active_run, current_phase) or {}
    loop_pid = active_run.get("loopPid")
    child_pid = current_phase_state.get("pid")
    loop_alive = pid_is_alive(loop_pid)
    child_alive = pid_is_alive(child_pid)
    phase_started_at = parse_iso(current_phase_state.get("startedAt"))
    phase_elapsed = None
    if phase_started_at is not None:
        phase_elapsed = round((now_utc() - phase_started_at).total_seconds(), 1)
    timeout_ms = current_phase_state.get("timeoutMs")
    timeout_seconds = round(timeout_ms / 1000, 1) if isinstance(timeout_ms, int) and timeout_ms > 0 else None

    if timeout_seconds is not None and phase_elapsed is not None and phase_elapsed > timeout_seconds:
        status = "overrun"
        reason = "phase_timeout_exceeded"
    elif child_alive:
        status = "running"
        reason = "child_process_alive"
    elif heartbeat_age is not None and heartbeat_age <= STALE_HEARTBEAT_SECONDS:
        status = "running"
        reason = "heartbeat_recent_pid_unconfirmed"
    elif loop_alive and (heartbeat_age is None or heartbeat_age <= STALE_HEARTBEAT_SECONDS):
        status = "running"
        reason = "loop_heartbeat_recent"
    elif loop_alive:
        status = "stale"
        reason = "loop_alive_but_heartbeat_expired"
    else:
        status = "stale"
        reason = "no_live_processes"

    return {
        "status": status,
        "reason": reason,
        "loopPid": loop_pid,
        "childPid": child_pid,
        "phase": current_phase,
        "phaseStartedAt": current_phase_state.get("startedAt"),
        "phaseElapsedSeconds": phase_elapsed,
        "phaseTimeoutSeconds": timeout_seconds,
        "heartbeatAt": active_run.get("heartbeatAt"),
        "heartbeatAgeSeconds": heartbeat_age,
        "runId": active_run.get("runId"),
        "itemId": active_run.get("itemId"),
        "targetPath": active_run.get("targetPath"),
    }


def effective_work_target(item: dict) -> str:
    return item.get("workTargetPath") or item["targetPath"]


def display_commands(commands: list[list[str]] | None) -> list[str]:
    return [" ".join(command) for command in commands or []]


def build_supplemental_reviewer_commands(item: dict) -> list[dict]:
    commands: list[dict] = []
    seen: set[tuple[str, ...]] = set()

    def append(argv_matrix: list[list[str]] | None, *, label_prefix: str) -> None:
        if not argv_matrix:
            return
        for index, argv in enumerate(argv_matrix, start=1):
            key = tuple(argv)
            if key in seen:
                continue
            seen.add(key)
            label = label_prefix if len(argv_matrix) == 1 else f"{label_prefix} #{index}"
            commands.append(
                {
                    "argv": list(argv),
                    "label": label,
                    "phase": "post",
                    "replayMode": "verify",
                }
            )

    append(item.get("focusedGateCommands"), label_prefix="focused gate")
    append(item.get("supplementalReviewCommands"), label_prefix="supplemental review command")
    return commands


def read_handover_excerpt(handover_path: str | None) -> str:
    if not handover_path:
        return ""
    path = pathlib.Path(handover_path)
    if not path.exists():
        return ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    excerpt = "\n".join(lines[:DEFAULT_HANDOVER_EXCERPT_LINES]).strip()
    return truncate_text(excerpt, DEFAULT_HANDOVER_EXCERPT_CHARS)


def maybe_work_target_snapshot(repo_root: pathlib.Path, backlog_path: pathlib.Path, item: dict) -> dict | None:
    work_target = item.get("workTargetPath")
    if not work_target or work_target == item["targetPath"]:
        return None
    try:
        return run_status_snapshot(repo_root, backlog_path, work_target)
    except subprocess.SubprocessError:
        return None


def build_executor_prompt(
    state: dict,
    item: dict,
    state_path: pathlib.Path,
    log_path: pathlib.Path,
    repo_root: pathlib.Path,
    backlog_path: pathlib.Path,
) -> str:
    work_target = effective_work_target(item)
    work_snapshot = maybe_work_target_snapshot(repo_root, backlog_path, item)
    handover_excerpt = read_handover_excerpt(item.get("handoverPath"))
    prompt_hints = item.get("promptHints") or []
    focused_gate_commands = display_commands(item.get("focusedGateCommands"))
    last_summary = truncate_text(item.get("lastSummary"), DEFAULT_SUMMARY_EXCERPT_CHARS)

    lines = [
        "$god-module-stabilization",
        "",
        "Continue the installed god-module-stabilization skill as one loop-driven stabilization iteration.",
        "Do exactly one bounded checkpoint on the target below and then stop.",
        "",
        f"Repository root: {state['repoRoot']}",
        f"Loop state path: {state_path}",
        f"Loop log path: {log_path}",
        f"Backlog path: {state['backlogPath']}",
        f"Target bucket: {item['bucketNumber']}. {item['title']}",
        f"Bucket target file: {item['targetPath']}",
        f"Primary work file: {work_target}",
        f"Target handover: {item.get('handoverPath') or 'none'}",
        "",
        "Current checkpoint brief:",
        f"- Latest checkpoint: {item.get('lastCheckpoint') or 'none'}",
        f"- Latest reason code: {item.get('lastReasonCode') or 'none'}",
        f"- Latest summary: {last_summary or 'none'}",
    ]

    progress_line = (((item.get("lastProgress") or {}).get("progress_line")) or "").strip()
    if progress_line:
        lines.append(f"- Bucket progress: {progress_line}")
    if work_snapshot is not None:
        lines.append(f"- Work-file progress: {work_snapshot.get('progress_line', '').strip()}")
        target_snapshot = work_snapshot.get("target") or {}
        labels = ", ".join(target_snapshot.get("labels") or [])
        next_step = target_snapshot.get("next_step") or "none"
        lines.append(f"- Work-file classification: {labels or 'none'}")
        lines.append(f"- Work-file next step: {next_step}")
    if focused_gate_commands:
        lines.append("- Focused gate replay commands:")
        lines.extend([f"  - {command}" for command in focused_gate_commands])
    if prompt_hints:
        lines.append("- Control hints:")
        lines.extend([f"  - {hint}" for hint in prompt_hints])
    if handover_excerpt:
        lines.extend(
            [
                "",
                "Current handover excerpt:",
                handover_excerpt,
            ]
        )

    lines.extend(
        [
            "",
            "Required behavior:",
            "- Use $god-module-stabilization in stabilize mode.",
            "- Read the backlog and the named handover before editing, but do not reread the full historical handover trail when the compact brief above already identifies the active stop point.",
            "- Treat the bucket target as the public wrapper identity and the primary work file as the structural work surface for this iteration.",
            "- Perform exactly one bounded seam, staged wave/apply, or equivalent clean checkpoint.",
            "- Prefer structural seams in the primary work file over new review-only checkpoints when a structural stop point is already documented.",
            "- Rerun the focused gate for the target.",
            "- Update the handover and backlog if the checkpoint changes the stop point.",
            "- Stop after the checkpoint. Do not continue into a second seam.",
            "- Do not create, remove, or edit loop-control artifacts directly; leave loop-state, loop logs, and stop files to the outer loop.",
            "- Return structured JSON only.",
            "- Use a machine-readable verification object. Do not append prose to replay commands.",
            "- Only include replay-safe verification commands in verification.replayCommands when the reviewer can rerun them without mutating repo state.",
            "- Do not classify stage/apply/commit/extract operations as replay-safe verification. Either omit them from replayCommands or mark them with replayMode=\"operate\" so the reviewer will record but not rerun them.",
            '- Put replayable verification commands in verification.replayCommands as argv arrays, for example:',
            '  {"verification":{"replayCommands":[{"argv":["Rscript","tools/test_with_renv.R","tests/testthat/test-prot-04-design.R"],"label":"focused design gate","phase":"post","replayMode":"verify"}],"display":["focused design gate passed"]}}',
            "- Keep human-oriented commentary in summary/notes only, not inside replay commands.",
            "",
            "When deciding targetStatusAfterIteration:",
            '- use "in_progress" if the target still needs more stabilization work',
            '- use "done" only if this backlog target is genuinely complete',
            '- use "blocked" if safe progress is blocked',
            "",
            "Always include reasonCode. Use a short string such as completed_checkpoint, blocked_gate_failure, or null when no code is useful.",
        ]
    )
    return "\n".join(lines)


def build_initial_state(repo_root: pathlib.Path, backlog_path: pathlib.Path, state_path: pathlib.Path, log_path: pathlib.Path, stop_file: pathlib.Path, max_iterations: int) -> dict:
    items = [asdict(item) for item in parse_backlog(backlog_path) if item.targetPath]
    current = next((item for item in items if item["status"] == "in_progress"), None)
    if current is None:
        current = next((item for item in items if item["status"] == "pending"), None)
    return {
        "schemaVersion": STATE_SCHEMA_VERSION,
        "createdAt": now_iso(),
        "updatedAt": now_iso(),
        "repoRoot": str(repo_root),
        "backlogPath": str(backlog_path),
        "statePath": str(state_path),
        "logPath": str(log_path),
        "stopFile": str(stop_file),
        "overridesPath": str(default_overrides_path(repo_root)),
        "runtimeDir": str(runtime_dir_for_state(state_path)),
        "status": "in_progress",
        "iteration": 0,
        "maxIterations": max_iterations,
        "scopePrefixes": [],
        "currentItemId": current["id"] if current else None,
        "items": items,
        "history": [],
        "activeRun": None,
    }


def resolve_repo_path(repo_root: pathlib.Path, candidate: str) -> pathlib.Path:
    path = pathlib.Path(candidate)
    if path.is_absolute():
        return path.resolve()
    return (repo_root / path).resolve()


def load_or_create_state(args) -> tuple[pathlib.Path, dict]:
    repo_root = pathlib.Path(args.repo_root).resolve()
    backlog_path = resolve_repo_path(repo_root, args.backlog_path)
    state_path = resolve_repo_path(repo_root, args.state_path)
    log_path = resolve_repo_path(repo_root, args.log_path)
    stop_file = resolve_repo_path(repo_root, args.stop_file)
    if state_path.exists():
        state = read_json(state_path)
        state, changed = migrate_state(state, state_path)
        changed = apply_repo_overrides(state, repo_root) or changed
        normalized_scope_prefixes = [str(resolve_repo_path(repo_root, prefix)) for prefix in (args.scope_prefix or [])]
        if normalized_scope_prefixes and state.get("scopePrefixes") != normalized_scope_prefixes:
            state["scopePrefixes"] = normalized_scope_prefixes
            state["updatedAt"] = now_iso()
            changed = True
        if args.max_iterations > state.get("maxIterations", 0):
            state["maxIterations"] = args.max_iterations
            state["updatedAt"] = now_iso()
            changed = True
        if changed:
            write_json(state_path, state)
        return state_path, state
    state = build_initial_state(
        repo_root,
        backlog_path,
        state_path,
        log_path,
        stop_file,
        args.max_iterations,
    )
    apply_repo_overrides(state, repo_root)
    state["scopePrefixes"] = [str(resolve_repo_path(repo_root, prefix)) for prefix in (args.scope_prefix or [])]
    write_json(state_path, state)
    return state_path, state


def print_progress(progress_snapshot: dict, reviewer_status: str, checkpoint: str) -> None:
    print(
        f"{progress_snapshot['progress_line']} review={reviewer_status} checkpoint={checkpoint}",
        file=sys.stderr,
    )


def build_run_artifacts(
    state_path: pathlib.Path,
    item: dict,
    iteration_number: int,
    *,
    executor_timeout_ms: int,
    review_timeout_ms: int,
) -> dict:
    runtime_dir = runtime_dir_for_state(state_path)
    base = runtime_dir / f"{item['id']}-iter-{iteration_number:03d}"
    return {
        "executor": {
            "payloadPath": str(base.with_name(f"{base.name}-executor-payload.json")),
            "stdoutPath": str(base.with_name(f"{base.name}-executor.stdout.log")),
            "stderrPath": str(base.with_name(f"{base.name}-executor.stderr.log")),
            "pid": None,
            "startedAt": None,
            "completedAt": None,
            "returncode": None,
            "timeoutMs": executor_timeout_ms,
        },
        "reviewer": {
            "payloadPath": str(base.with_name(f"{base.name}-reviewer-payload.json")),
            "stdoutPath": str(base.with_name(f"{base.name}-reviewer.stdout.log")),
            "stderrPath": str(base.with_name(f"{base.name}-reviewer.stderr.log")),
            "pid": None,
            "startedAt": None,
            "completedAt": None,
            "returncode": None,
            "timeoutMs": review_timeout_ms,
        },
    }


def start_active_run(
    state: dict,
    state_path: pathlib.Path,
    item: dict,
    *,
    executor_timeout_ms: int,
    review_timeout_ms: int,
) -> dict:
    iteration_number = state["iteration"] + 1
    artifacts = build_run_artifacts(
        state_path,
        item,
        iteration_number,
        executor_timeout_ms=executor_timeout_ms,
        review_timeout_ms=review_timeout_ms,
    )
    active_run = {
        "runId": f"{item['id']}-iter-{iteration_number:03d}",
        "loopPid": os.getpid(),
        "itemId": item["id"],
        "targetPath": item["targetPath"],
        "workTargetPath": effective_work_target(item),
        "startedAt": now_iso(),
        "heartbeatAt": now_iso(),
        "phase": "executor",
        "executor": artifacts["executor"],
        "reviewer": artifacts["reviewer"],
    }
    state["activeRun"] = active_run
    state["updatedAt"] = now_iso()
    write_json(state_path, state)
    return active_run


def touch_active_run(state: dict, state_path: pathlib.Path) -> None:
    active_run = state.get("activeRun")
    if not active_run:
        return
    active_run["heartbeatAt"] = now_iso()
    active_run["loopPid"] = os.getpid()
    state["updatedAt"] = now_iso()
    write_json(state_path, state)


def update_phase_state(
    state: dict,
    state_path: pathlib.Path,
    *,
    phase: str,
    pid: int | None = None,
    returncode: int | None = None,
    started: bool = False,
    completed: bool = False,
) -> None:
    active_run = state.get("activeRun")
    if not active_run:
        return
    active_run["phase"] = phase
    phase_info = active_run[phase]
    if started:
        phase_info["startedAt"] = phase_info["startedAt"] or now_iso()
    if pid is not None:
        phase_info["pid"] = pid
    if returncode is not None:
        phase_info["returncode"] = returncode
    if completed:
        phase_info["completedAt"] = now_iso()
        phase_info["pid"] = None
    touch_active_run(state, state_path)


def run_script_live(
    script_path: pathlib.Path,
    payload_path: pathlib.Path,
    repo_root: pathlib.Path,
    state: dict,
    state_path: pathlib.Path,
    *,
    phase: str,
    timeout_ms: int | None = None,
) -> dict:
    command = script_command(script_path, payload_path)
    active_run = state["activeRun"]
    phase_info = active_run[phase]
    stdout_path = pathlib.Path(phase_info["stdoutPath"])
    stderr_path = pathlib.Path(phase_info["stderrPath"])
    stdout_path.parent.mkdir(parents=True, exist_ok=True)

    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open("w", encoding="utf-8") as stderr_handle:
        process = subprocess.Popen(
            command,
            cwd=repo_root,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            start_new_session=True,
        )
        update_phase_state(state, state_path, phase=phase, pid=process.pid, started=True)
        started_monotonic = time.monotonic()

        while True:
            returncode = process.poll()
            if returncode is not None:
                break
            if timeout_ms is not None and timeout_ms > 0:
                elapsed_ms = (time.monotonic() - started_monotonic) * 1000
                if elapsed_ms > timeout_ms:
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    try:
                        returncode = process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        returncode = process.wait()
                    update_phase_state(state, state_path, phase=phase, returncode=returncode, completed=True)
                    stdout_text = stdout_path.read_text(encoding="utf-8", errors="replace")
                    stderr_text = stderr_path.read_text(encoding="utf-8", errors="replace")
                    timeout_seconds = round(timeout_ms / 1000, 1)
                    raise ScriptExecutionError(
                        phase=phase,
                        returncode=returncode,
                        stdout_text=stdout_text,
                        stderr_text=f"{phase} timed out after {timeout_seconds} seconds",
                        stdout_path=str(stdout_path),
                        stderr_path=str(stderr_path),
                        command=command,
                        timed_out=True,
                    )
            touch_active_run(state, state_path)
            time.sleep(HEARTBEAT_INTERVAL_SECONDS)

    update_phase_state(state, state_path, phase=phase, returncode=returncode, completed=True)
    stdout_text = stdout_path.read_text(encoding="utf-8", errors="replace")
    stderr_text = stderr_path.read_text(encoding="utf-8", errors="replace")

    if returncode != 0:
        raise ScriptExecutionError(
            phase=phase,
            returncode=returncode,
            stdout_text=stdout_text,
            stderr_text=stderr_text,
            stdout_path=str(stdout_path),
            stderr_path=str(stderr_path),
            command=command,
        )

    try:
        return parse_structured_stdout(stdout_text)
    except json.JSONDecodeError as error:
        raise ScriptExecutionError(
            phase=phase,
            returncode=returncode,
            stdout_text=stdout_text,
            stderr_text=f"Failed to parse {phase} JSON output: {error}",
            stdout_path=str(stdout_path),
            stderr_path=str(stderr_path),
            command=command,
        ) from error


def parse_structured_stdout(stdout_text: str) -> dict:
    stripped = stdout_text.strip()
    if not stripped:
        raise json.JSONDecodeError("No JSON output found", stdout_text, 0)
    try:
        return json.loads(stripped)
    except json.JSONDecodeError as original_error:
        parsed_objects: list[dict] = []
        for line in stdout_text.splitlines():
            candidate = line.strip()
            if not candidate:
                continue
            try:
                parsed = json.loads(candidate)
            except json.JSONDecodeError:
                continue
            if isinstance(parsed, dict):
                parsed_objects.append(parsed)
        if parsed_objects:
            first = parsed_objects[0]
            if all(obj == first for obj in parsed_objects[1:]):
                return first
        raise original_error


def current_item(state: dict) -> dict | None:
    current_id = state.get("currentItemId")
    return next((item for item in state["items"] if item["id"] == current_id), None)


def latest_history_entry_for_item(state: dict, item_id: str | None) -> dict | None:
    if not item_id:
        return None
    for entry in reversed(state.get("history", [])):
        if entry.get("itemId") == item_id:
            return entry
    return None


def reconcile_legacy_error_state(state_path: pathlib.Path, state: dict) -> dict | None:
    if state.get("activeRun"):
        return None
    item = current_item(state)
    if item is None or item.get("status") != "in_progress":
        return None
    latest_entry = latest_history_entry_for_item(state, item.get("id"))
    if not latest_entry or not latest_entry.get("error"):
        return None
    if latest_entry.get("runtimeEvent") == "legacy_idle_error_reconciled":
        return None

    error_text = latest_entry["error"]
    summary = (
        "Recovered legacy ambiguous loop state from pre-liveness metadata; "
        "the latest recorded run ended with an error and no active process remains."
    )
    if error_text:
        summary = f"{summary}\n\nLast recorded error:\n{error_text[-4000:]}"

    item["status"] = "blocked"
    item["lastSummary"] = summary
    state["status"] = "blocked"
    state["updatedAt"] = now_iso()

    history_entry = {
        "iteration": state["iteration"],
        "itemId": item["id"],
        "targetPath": item["targetPath"],
        "executor": None,
        "reviewer": None,
        "progress": None,
        "runtimeEvent": "legacy_idle_error_reconciled",
        "error": summary,
        "timestamp": now_iso(),
    }
    state["history"].append(history_entry)
    append_jsonl(pathlib.Path(state["logPath"]).resolve(), history_entry)
    write_json(state_path, state)
    return {
        "event": "legacy_idle_error_reconciled",
        "summary": summary,
        "itemId": item["id"],
        "targetPath": item["targetPath"],
    }


def reconcile_stale_active_run(state_path: pathlib.Path, state: dict) -> dict | None:
    active_run = state.get("activeRun")
    if not active_run:
        return None
    runtime = build_runtime_snapshot(state)
    if runtime["status"] != "stale":
        return None

    log_path = pathlib.Path(state["logPath"]).resolve()
    item = next((entry for entry in state["items"] if entry["id"] == active_run.get("itemId")), None)
    active_phase = active_run.get("phase")
    phase_info = phase_state(active_run, active_phase) or {}
    stderr_tail = read_tail(phase_info.get("stderrPath"))
    stdout_tail = read_tail(phase_info.get("stdoutPath"))
    tail = stderr_tail or stdout_tail
    summary = (
        f"Recovered stale active run during {active_phase or 'unknown'} phase; "
        "no live loop or child process remained."
    )
    if tail:
        summary = f"{summary}\n\nLast log tail:\n{tail}"

    if item is not None:
        item["status"] = "blocked"
        item["lastSummary"] = summary
    state["status"] = "blocked"
    state["activeRun"] = None
    state["updatedAt"] = now_iso()

    history_entry = {
        "iteration": state["iteration"],
        "itemId": active_run.get("itemId"),
        "targetPath": active_run.get("targetPath"),
        "executor": None,
        "reviewer": None,
        "progress": None,
        "runtimeEvent": "stale_active_run_reconciled",
        "runtimeSnapshot": runtime,
        "error": summary,
        "timestamp": now_iso(),
    }
    state["history"].append(history_entry)
    append_jsonl(log_path, history_entry)
    write_json(state_path, state)
    return {
        "event": "stale_active_run_reconciled",
        "summary": summary,
        "runtime": runtime,
        "itemId": active_run.get("itemId"),
        "targetPath": active_run.get("targetPath"),
    }


def run_loop(args) -> dict:
    state_path, state = load_or_create_state(args)
    reconciled = reconcile_stale_active_run(state_path, state)
    if reconciled is not None:
        state = read_json(state_path)

    log_path = pathlib.Path(state["logPath"]).resolve()
    repo_root = pathlib.Path(state["repoRoot"]).resolve()
    backlog_path = pathlib.Path(state["backlogPath"]).resolve()
    stop_file = pathlib.Path(state["stopFile"]).resolve()
    iteration_limit = args.iteration_limit
    runtime = build_runtime_snapshot(state)
    if runtime["status"] == "running":
        return {
            "ok": False,
            "statePath": str(state_path),
            "error": "a stabilization loop run is already active for this state file",
            "runtime": runtime,
        }

    if args.clear_stop_file and stop_file.exists():
        stop_file.unlink()

    override_item = ensure_override_item(state, args.target, args.item_id)
    if override_item is not None:
        for item in state["items"]:
            if item["id"] != override_item["id"] and item["status"] == "in_progress":
                item["status"] = "pending"
        if override_item["status"] == "done":
            override_item["completedAt"] = None
        state["currentItemId"] = override_item["id"]
        if state["iteration"] < state["maxIterations"]:
            override_item["status"] = "in_progress"
            state["status"] = "in_progress"
        state["updatedAt"] = now_iso()
        write_json(state_path, state)

    while iteration_limit > 0:
        if stop_file.exists():
            state["status"] = "stopped"
            state["updatedAt"] = now_iso()
            write_json(state_path, state)
            return {"ok": True, "statePath": str(state_path), "summary": summarize_state(state), "reconciled": reconciled}

        if state["iteration"] >= state["maxIterations"]:
            state["status"] = "paused_max_iterations"
            state["updatedAt"] = now_iso()
            write_json(state_path, state)
            return {"ok": True, "statePath": str(state_path), "summary": summarize_state(state), "reconciled": reconciled}

        item = next_open_item(state)
        if item is None:
            state["status"] = "completed"
            state["currentItemId"] = None
            state["updatedAt"] = now_iso()
            write_json(state_path, state)
            return {"ok": True, "statePath": str(state_path), "summary": summarize_state(state), "reconciled": reconciled}

        item["status"] = "in_progress"
        item["startedAt"] = item["startedAt"] or now_iso()
        state["currentItemId"] = item["id"]
        state["status"] = "in_progress"
        start_active_run(
            state,
            state_path,
            item,
            executor_timeout_ms=args.executor_timeout_ms,
            review_timeout_ms=args.review_timeout_ms,
        )
        state = read_json(state_path)
        item = current_item(state) or item

        runner_script = pathlib.Path(args.runner_script or (pathlib.Path(__file__).resolve().parent / "stabilization-codex-runner.py")).resolve()
        reviewer_script = pathlib.Path(args.reviewer_script or (pathlib.Path(__file__).resolve().parent / "stabilization-reviewer.py")).resolve()

        try:
            validate_local_executor_schema(runner_script)
            work_target = effective_work_target(item)
            executor_payload = {
                "prompt": build_executor_prompt(state, item, state_path, log_path, repo_root, backlog_path),
                "repoRoot": str(repo_root),
                "backlogPath": str(backlog_path),
                "statePath": str(state_path),
                "logPath": str(log_path),
                "target": work_target,
                "bucketTarget": item["targetPath"],
                "workTarget": work_target,
                "handoverPath": item.get("handoverPath"),
                "item": item,
                "model": args.model,
                "sandbox": args.sandbox,
                "approvalPolicy": args.approval_policy,
                "networkAccessEnabled": args.network_access,
            }
            executor_payload_path = pathlib.Path(state["activeRun"]["executor"]["payloadPath"])
            write_json(executor_payload_path, executor_payload)
            executor_result = run_script_live(
                runner_script,
                executor_payload_path,
                repo_root,
                state,
                state_path,
                phase="executor",
                timeout_ms=args.executor_timeout_ms,
            )
            executor_result = validate_executor_result(executor_result)

            state["iteration"] += 1
            iteration_limit -= 1

            reviewer_payload = {
                "repoRoot": str(repo_root),
                "targetPath": item["targetPath"],
                "backlogPath": str(backlog_path),
                "statusScriptPath": str(pathlib.Path(__file__).resolve().parent / "stabilization-status.py"),
                "executorResult": executor_result,
                "supplementalReviewCommands": build_supplemental_reviewer_commands(item),
                "timeoutMs": args.review_timeout_ms,
            }
            reviewer_payload_path = pathlib.Path(state["activeRun"]["reviewer"]["payloadPath"])
            write_json(reviewer_payload_path, reviewer_payload)
            reviewer_result = run_script_live(
                reviewer_script,
                reviewer_payload_path,
                repo_root,
                state,
                state_path,
                phase="reviewer",
                timeout_ms=args.review_timeout_ms,
            )
            reviewer_result = validate_reviewer_result(reviewer_result)
            progress_snapshot = reviewer_result.get("statusSnapshot") or run_status_snapshot(repo_root, backlog_path, item["targetPath"])

            item["lastCheckpoint"] = executor_result["checkpoint"]
            item["lastSummary"] = executor_result["summary"]
            item["lastProgress"] = progress_snapshot
            item["lastReview"] = reviewer_result
            item["lastFilesChanged"] = production_files_from_changed_list(executor_result.get("filesChanged", []))
            item["lastReasonCode"] = blocked_reason_code(executor_result, reviewer_result, fallback_status="approved")

            effective_target_state = reviewer_result.get("targetStateOverride") or executor_result["targetStatusAfterIteration"]

            if executor_result["status"] != "completed" or reviewer_result["status"] != "approved":
                item["status"] = "blocked" if reviewer_result["status"] == "blocked" else "failed"
                state["status"] = "blocked"
            elif effective_target_state == "done":
                item["status"] = "done"
                item["completedAt"] = now_iso()
            elif effective_target_state == "blocked":
                item["status"] = "blocked"
                state["status"] = "blocked"
            else:
                item["status"] = "in_progress"
                state["status"] = "in_progress"

            history_entry = {
                "iteration": state["iteration"],
                "itemId": item["id"],
                "targetPath": item["targetPath"],
                "executor": executor_result,
                "reviewer": reviewer_result,
                "progress": progress_snapshot,
                "artifacts": state["activeRun"],
                "reasonCode": item.get("lastReasonCode"),
                "timestamp": now_iso(),
            }
            state["history"].append(history_entry)
            append_jsonl(log_path, history_entry)
            state["activeRun"] = None
            state["updatedAt"] = now_iso()
            write_json(state_path, state)
            print_progress(progress_snapshot, reviewer_result["status"], executor_result["checkpoint"])

            if state["status"] == "blocked":
                return {"ok": True, "statePath": str(state_path), "summary": summarize_state(state), "reconciled": reconciled}
        except (ScriptExecutionError, ProtocolValidationError) as error:
            state["iteration"] += 1
            iteration_limit -= 1
            if isinstance(error, ScriptExecutionError):
                error_text = error.stderr_text or error.stdout_text or str(error)
                artifacts = {
                    "stdoutPath": error.stdout_path,
                    "stderrPath": error.stderr_path,
                    "command": error.command,
                    "returncode": error.returncode,
                }
                phase = error.phase
            else:
                error_text = str(error)
                artifacts = None
                phase = error.phase
            item["status"] = "blocked" if isinstance(error, ScriptExecutionError) and error.timed_out else classify_runner_error(error_text)
            item["lastSummary"] = error_text
            if isinstance(error, ProtocolValidationError):
                item["lastReasonCode"] = classify_runner_reason_code(error_text, phase=phase, timed_out=False)
            elif isinstance(error, ScriptExecutionError):
                item["lastReasonCode"] = classify_runner_reason_code(
                    error_text,
                    phase=phase,
                    timed_out=error.timed_out,
                )
            else:
                item["lastReasonCode"] = blocked_reason_code(None, None, fallback_status=item["status"])
            state["status"] = "blocked"
            history_entry = {
                "iteration": state["iteration"],
                "itemId": item["id"],
                "targetPath": item["targetPath"],
                "executor": None,
                "reviewer": None,
                "progress": None,
                "error": error_text,
                "phase": phase,
                "artifacts": artifacts,
                "reasonCode": item.get("lastReasonCode"),
                "timestamp": now_iso(),
            }
            state["history"].append(history_entry)
            append_jsonl(log_path, history_entry)
            state["activeRun"] = None
            state["updatedAt"] = now_iso()
            write_json(state_path, state)
            return {"ok": True, "statePath": str(state_path), "summary": summarize_state(state), "reconciled": reconciled}

    state["updatedAt"] = now_iso()
    write_json(state_path, state)
    return {"ok": True, "statePath": str(state_path), "summary": summarize_state(state), "reconciled": reconciled}


def queue_result(args) -> dict:
    state_path, state = load_or_create_state(args)
    queued_items = queue_targets(state, args.target or [], item_ids=args.item_id or [])
    write_json(state_path, state)
    return {
        "ok": True,
        "statePath": str(state_path),
        "queued": [
            {
                "id": item["id"],
                "targetPath": item["targetPath"],
                "status": item["status"],
            }
            for item in queued_items
        ],
        "summary": summarize_state(state),
    }


def status_result(state_path: pathlib.Path) -> dict:
    if not state_path.exists():
        return {
            "ok": False,
            "statePath": str(state_path),
            "error": f"state file not found: {state_path}",
        }

    state = read_json(state_path)
    state, changed = migrate_state(state, state_path)
    repo_root_str = state.get("repoRoot")
    if isinstance(repo_root_str, str) and repo_root_str:
        changed = apply_repo_overrides(state, pathlib.Path(repo_root_str).resolve()) or changed
    if changed:
        write_json(state_path, state)
    reconciled = reconcile_stale_active_run(state_path, state)
    if reconciled is not None:
        state = read_json(state_path)
    legacy_reconciled = reconcile_legacy_error_state(state_path, state)
    if legacy_reconciled is not None:
        state = read_json(state_path)
        reconciled = legacy_reconciled if reconciled is None else {
            "event": "multiple_reconciliations",
            "reconciliations": [reconciled, legacy_reconciled],
        }

    return {
        "ok": True,
        "statePath": str(state_path),
        "summary": summarize_state(state),
        "currentItem": current_item(state),
        "runtime": build_runtime_snapshot(state),
        "reconciled": reconciled,
    }


def run_script_once(
    script_path: pathlib.Path,
    payload_path: pathlib.Path,
    repo_root: pathlib.Path,
    *,
    stdout_path: pathlib.Path,
    stderr_path: pathlib.Path,
) -> dict:
    command = script_command(script_path, payload_path)
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        command,
        cwd=repo_root,
        capture_output=True,
        text=True,
    )
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise ScriptExecutionError(
            phase="reviewer",
            returncode=completed.returncode,
            stdout_text=completed.stdout,
            stderr_text=completed.stderr,
            stdout_path=str(stdout_path),
            stderr_path=str(stderr_path),
            command=command,
        )
    try:
        return parse_structured_stdout(completed.stdout)
    except json.JSONDecodeError as error:
        raise ScriptExecutionError(
            phase="reviewer",
            returncode=completed.returncode,
            stdout_text=completed.stdout,
            stderr_text=f"Failed to parse reviewer JSON output: {error}",
            stdout_path=str(stdout_path),
            stderr_path=str(stderr_path),
            command=command,
        ) from error


def recover_blocked_review(state_path: pathlib.Path, state: dict, reviewer_script: pathlib.Path) -> dict | None:
    if state.get("status") != "blocked":
        return None
    item = current_item(state)
    if item is None or item.get("status") not in {"blocked", "failed"}:
        return None
    latest_entry = latest_history_entry_for_item(state, item.get("id"))
    if not latest_entry:
        return None
    latest_reviewer = latest_entry.get("reviewer") or {}
    latest_reviewer_status = latest_reviewer.get("status")
    if latest_reviewer_status not in {"blocked", "needs_changes"}:
        return None
    artifacts = latest_entry.get("artifacts") or {}
    reviewer_artifacts = artifacts.get("reviewer") or {}
    payload_path_str = reviewer_artifacts.get("payloadPath")
    if not payload_path_str:
        return None
    payload_path = pathlib.Path(payload_path_str)
    if not payload_path.exists():
        return None

    repo_root = pathlib.Path(state["repoRoot"]).resolve()
    runtime_dir = runtime_dir_for_state(state_path)
    base = runtime_dir / f"{item['id']}-review-recovery-{state['iteration']:03d}"
    stdout_path = base.with_name(f"{base.name}.stdout.log")
    stderr_path = base.with_name(f"{base.name}.stderr.log")

    reviewer_result = validate_reviewer_result(
        run_script_once(
            reviewer_script,
            payload_path,
            repo_root,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
        )
    )
    progress_snapshot = reviewer_result.get("statusSnapshot") or run_status_snapshot(
        repo_root,
        pathlib.Path(state["backlogPath"]).resolve(),
        item["targetPath"],
    )

    executor_result = latest_entry.get("executor") or {}
    item["lastReview"] = reviewer_result
    item["lastProgress"] = progress_snapshot
    item["lastFilesChanged"] = production_files_from_changed_list(executor_result.get("filesChanged", []))
    item["lastReasonCode"] = blocked_reason_code(executor_result, reviewer_result, fallback_status="blocked_review_recovery")

    recovery_event = "blocked_review_rerun"
    if reviewer_result["status"] == "approved" and executor_result.get("status") == "completed":
        target_state = reviewer_result.get("targetStateOverride") or executor_result.get("targetStatusAfterIteration", "in_progress")
        item["lastSummary"] = executor_result.get("summary", item.get("lastSummary"))
        if target_state == "done":
            item["status"] = "done"
            item["completedAt"] = item.get("completedAt") or now_iso()
        elif target_state == "blocked":
            item["status"] = "blocked"
        else:
            item["status"] = "in_progress"
        recovery_event = "blocked_review_recovered"
    else:
        item["status"] = "blocked" if reviewer_result["status"] == "blocked" else "failed"

    state["status"] = "completed" if item["status"] == "done" and next_open_item(state) is None else (
        "blocked" if item["status"] in {"blocked", "failed"} else "in_progress"
    )
    state["activeRun"] = None
    state["updatedAt"] = now_iso()

    history_entry = {
        "iteration": state["iteration"],
        "itemId": item["id"],
        "targetPath": item["targetPath"],
        "executor": executor_result,
        "reviewer": reviewer_result,
        "progress": progress_snapshot,
        "runtimeEvent": recovery_event,
        "reasonCode": item.get("lastReasonCode"),
        "artifacts": {
            "recovery": {
                "reviewerPayloadPath": str(payload_path),
                "stdoutPath": str(stdout_path),
                "stderrPath": str(stderr_path),
            }
        },
        "timestamp": now_iso(),
    }
    state["history"].append(history_entry)
    append_jsonl(pathlib.Path(state["logPath"]).resolve(), history_entry)
    write_json(state_path, state)
    return {
        "event": recovery_event,
        "itemId": item["id"],
        "targetPath": item["targetPath"],
        "status": item["status"],
        "reasonCode": item.get("lastReasonCode"),
    }


def recover_result(state_path: pathlib.Path, reviewer_script: pathlib.Path | None = None) -> dict:
    status = status_result(state_path)
    if not status["ok"]:
        return status
    if status.get("reconciled") is not None:
        status["recovered"] = True
        return status
    runtime = status["runtime"]
    if runtime["status"] == "running":
        status["ok"] = False
        status["error"] = "cannot recover while a live loop run is still active"
        return status
    state = read_json(state_path)
    reviewer_script = reviewer_script or (pathlib.Path(__file__).resolve().parent / "stabilization-reviewer.py")
    review_recovered = recover_blocked_review(state_path, state, reviewer_script.resolve())
    if review_recovered is not None:
        refreshed = status_result(state_path)
        refreshed["recovered"] = True
        refreshed["reviewRecovered"] = review_recovered
        return refreshed
    status["recovered"] = False
    status["message"] = "No stale active run needed recovery."
    return status


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--repo-root", default=".")
    run_parser.add_argument("--backlog-path", default="tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md")
    run_parser.add_argument("--state-path", default=".god-module-stabilization/loop-state.json")
    run_parser.add_argument("--log-path", default=".god-module-stabilization/loop.jsonl")
    run_parser.add_argument("--stop-file", default=".god-module-stabilization/loop.stop")
    run_parser.add_argument("--clear-stop-file", action="store_true")
    run_parser.add_argument("--runner-script")
    run_parser.add_argument("--reviewer-script")
    run_parser.add_argument("--iteration-limit", type=int, default=1)
    run_parser.add_argument("--max-iterations", type=int, default=25)
    run_parser.add_argument("--executor-timeout-ms", type=int, default=DEFAULT_EXECUTOR_TIMEOUT_MS)
    run_parser.add_argument("--review-timeout-ms", type=int, default=900000)
    run_parser.add_argument("--model")
    run_parser.add_argument("--sandbox", default="workspace-write")
    run_parser.add_argument("--approval-policy", default="never")
    run_parser.add_argument("--network-access", type=lambda value: str(value).lower() == "true", default=True)
    run_parser.add_argument("--target")
    run_parser.add_argument("--item-id")
    run_parser.add_argument("--scope-prefix", action="append", default=[])

    queue_parser = subparsers.add_parser("queue")
    queue_parser.add_argument("--repo-root", default=".")
    queue_parser.add_argument("--backlog-path", default="tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md")
    queue_parser.add_argument("--state-path", default=".god-module-stabilization/loop-state.json")
    queue_parser.add_argument("--log-path", default=".god-module-stabilization/loop.jsonl")
    queue_parser.add_argument("--stop-file", default=".god-module-stabilization/loop.stop")
    queue_parser.add_argument("--max-iterations", type=int, default=25)
    queue_parser.add_argument("--target", action="append", default=[])
    queue_parser.add_argument("--item-id", action="append", default=[])
    queue_parser.add_argument("--scope-prefix", action="append", default=[])

    status_parser = subparsers.add_parser("status")
    status_parser.add_argument("--state-path", default=".god-module-stabilization/loop-state.json")

    recover_parser = subparsers.add_parser("recover")
    recover_parser.add_argument("--state-path", default=".god-module-stabilization/loop-state.json")
    recover_parser.add_argument("--reviewer-script")

    args = parser.parse_args()
    if args.command == "run":
        print(json.dumps(run_loop(args), indent=2))
        return 0
    if args.command == "queue":
        print(json.dumps(queue_result(args), indent=2))
        return 0
    if args.command == "status":
        print(json.dumps(status_result(pathlib.Path(args.state_path).resolve()), indent=2))
        return 0

    reviewer_script = pathlib.Path(args.reviewer_script).resolve() if args.reviewer_script else None
    print(json.dumps(recover_result(pathlib.Path(args.state_path).resolve(), reviewer_script=reviewer_script), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
