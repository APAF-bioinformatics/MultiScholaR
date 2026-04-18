#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import re
import shlex
import subprocess
import sys

ANNOTATION_SUFFIX_RE = re.compile(
    r"\s+(?:['\"])?\((?:pre-checkpoint|post-checkpoint|pre-review|post-review|"
    r"pre-check|post-check|passed|failed|pass|fail|warning|warnings|skip|skips|expected)"
    r"[^)]*\)(?:['\"])?\s*$",
    re.IGNORECASE,
)
SHELL_CONTROL_TOKENS = {"|", "||", "&&", ";", ">", ">>", "<", "<<", "&"}
VERIFY_PHASES = {"post", "gate", "review", "check", "verify"}
OPERATE_PHASES = {"apply", "stage", "commit", "extract", "manifest"}
NON_REPLAY_SAFE_MARKERS = (
    "apply_wave.py",
    "stage_wave.py",
    "git commit",
    "git add",
    "verify_refactor.r",
    "extract_blocks.r",
)

STAGING_GENERATED_DIR_PARTS = ("tools", "refactor", "staging")


def run_status(status_script: pathlib.Path, repo_root: pathlib.Path, target_path: str, backlog_path: str) -> dict:
    cmd = [
        sys.executable,
        str(status_script),
        "--target",
        target_path,
        "--backlog",
        backlog_path,
        "--json",
    ]
    completed = subprocess.run(
        cmd,
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(completed.stdout)


def normalize_replay_command(command_entry) -> tuple[list[str] | None, str | None, str, dict]:
    metadata = {
        "label": None,
        "phase": None,
        "source": "legacy",
        "replayMode": None,
    }
    if isinstance(command_entry, dict):
        argv_value = command_entry.get("argv")
        metadata["label"] = command_entry.get("label")
        metadata["phase"] = command_entry.get("phase")
        metadata["replayMode"] = command_entry.get("replayMode")
        if argv_value is not None:
            metadata["source"] = "structured"
            if not isinstance(argv_value, list) or not argv_value:
                return None, "structured replay command must provide a non-empty argv array", "", metadata
            if not all(isinstance(token, str) and token.strip() for token in argv_value):
                return None, "structured replay argv entries must be non-empty strings", "", metadata
            if any(token in SHELL_CONTROL_TOKENS for token in argv_value):
                return None, "structured replay argv contains shell control operators", "", metadata
            argv = [token for token in argv_value]
            return argv, None, shlex.join(argv), metadata
        raw = str(command_entry.get("replayCommand") or command_entry.get("command") or "").strip()
        metadata["source"] = "legacy_dict"
    else:
        raw = str(command_entry or "").strip()
        metadata["source"] = "legacy_string"
    if not raw:
        return None, "empty replay command", raw, metadata

    sanitized = ANNOTATION_SUFFIX_RE.sub("", raw).strip()
    if not sanitized:
        return None, "replay command became empty after stripping annotation", raw, metadata

    try:
        argv = shlex.split(sanitized, posix=True)
    except ValueError as error:
        return None, f"unable to parse replay command: {error}", raw, metadata

    if not argv:
        return None, "empty argv after parsing replay command", raw, metadata
    if any(token in SHELL_CONTROL_TOKENS for token in argv):
        return None, "replay command uses shell control operators and is not safe to rerun directly", raw, metadata
    return argv, None, raw, metadata


def infer_replay_mode(argv: list[str], metadata: dict) -> tuple[str, str | None]:
    lowered_tokens = " ".join(argv).lower()
    phase = metadata.get("phase")
    lowered_phase = phase.lower() if isinstance(phase, str) else None

    if lowered_phase == "pre":
        return "operate", "pre-phase verification is not replay-safe after the checkpoint mutates repo state"

    marker = next((marker for marker in NON_REPLAY_SAFE_MARKERS if marker in lowered_tokens), None)
    if marker is not None:
        return "operate", f"command is not replay-safe after mutation ({marker})"

    explicit = metadata.get("replayMode")
    if isinstance(explicit, str) and explicit in {"verify", "operate"}:
        return explicit, None

    if lowered_phase is not None:
        if lowered_phase in VERIFY_PHASES:
            return "verify", None
        if lowered_phase in OPERATE_PHASES:
            return "operate", f"phase {lowered_phase} is non-replay-safe"

    return "verify", None


def extract_replay_commands(executor_result: dict, supplemental_commands: list | None = None) -> tuple[list, list[str], bool]:
    verification = executor_result.get("verification")
    commands: list = []
    notes: list[str] = []
    used_legacy = False

    if isinstance(verification, dict):
        replay_commands = verification.get("replayCommands")
        if isinstance(replay_commands, list) and replay_commands:
            commands.extend(replay_commands)
    if not commands:
        tests_run = executor_result.get("testsRun", [])
        if isinstance(tests_run, list) and tests_run:
            commands.extend(tests_run)
            notes.append("Using legacy testsRun fallback for reviewer replay.")
            used_legacy = True
    if isinstance(supplemental_commands, list) and supplemental_commands:
        commands.extend(supplemental_commands)
        notes.append(f"Added {len(supplemental_commands)} supplemental reviewer command(s) from loop state.")
    return commands, notes, used_legacy


def rerun_test_commands(repo_root: pathlib.Path, commands: list, timeout_ms: int) -> tuple[list[str], list[dict], list[str], list[dict]]:
    passed: list[str] = []
    failures: list[dict] = []
    notes: list[str] = []
    replayed_commands: list[dict] = []
    seen: set[str] = set()
    for command in commands:
        argv, normalize_error, raw_command, metadata = normalize_replay_command(command)
        if normalize_error is not None:
            failures.append(
                {
                    "command": raw_command,
                    "returncode": None,
                    "stdout": "",
                    "stderr": normalize_error,
                }
            )
            continue
        normalized_display = shlex.join(argv)
        if normalized_display in seen:
            notes.append(f"Skipped duplicate replay command: {normalized_display}")
            continue
        replay_mode, replay_reason = infer_replay_mode(argv, metadata)
        if replay_mode != "verify":
            if replay_reason:
                notes.append(f"Skipped non-replay-safe command ({replay_mode}): {normalized_display} [{replay_reason}]")
            else:
                notes.append(f"Skipped non-replay-safe command ({replay_mode}): {normalized_display}")
            continue
        seen.add(normalized_display)
        replayed_commands.append(
            {
                "argv": argv,
                "display": normalized_display,
                "label": metadata.get("label"),
                "phase": metadata.get("phase"),
                "source": metadata.get("source"),
                "replayMode": replay_mode,
            }
        )
        completed = subprocess.run(
            argv,
            cwd=repo_root,
            capture_output=True,
            text=True,
            timeout=max(1, timeout_ms) / 1000,
        )
        if completed.returncode == 0:
            passed.append(normalized_display)
            continue
        failures.append(
            {
                "command": normalized_display,
                "returncode": completed.returncode,
                "stdout": completed.stdout[-4000:],
                "stderr": completed.stderr[-4000:],
            }
        )
    return passed, failures, notes, replayed_commands


def resolve_changed_file_path(repo_root: pathlib.Path, changed_file: str) -> tuple[pathlib.Path, str]:
    path = pathlib.Path(changed_file)
    if path.is_absolute():
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(repo_root)
            return resolved, str(relative)
        except ValueError:
            return resolved, str(resolved)
    return (repo_root / path).resolve(), changed_file


def is_generated_staging_artifact(relative_display: str) -> bool:
    parts = pathlib.Path(relative_display).parts
    if len(parts) < len(STAGING_GENERATED_DIR_PARTS):
        return False
    if tuple(parts[: len(STAGING_GENERATED_DIR_PARTS)]) != STAGING_GENERATED_DIR_PARTS:
        return False
    filename = pathlib.Path(relative_display).name
    return filename.startswith("collate-")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("payload")
    args = parser.parse_args()

    payload_path = pathlib.Path(args.payload).resolve()
    payload = json.loads(payload_path.read_text(encoding="utf-8"))

    repo_root = pathlib.Path(payload["repoRoot"]).resolve()
    status_script = pathlib.Path(payload["statusScriptPath"]).resolve()
    target_path = payload["targetPath"]
    backlog_path = payload["backlogPath"]
    timeout_ms = int(payload.get("timeoutMs", 900000))
    executor_result = payload["executorResult"]
    supplemental_commands = payload.get("supplementalReviewCommands")

    replay_commands, replay_notes, used_legacy_fallback = extract_replay_commands(executor_result, supplemental_commands)
    passed_tests, test_failures, rerun_notes, replayed_commands = rerun_test_commands(repo_root, replay_commands, timeout_ms)
    status_snapshot = run_status(status_script, repo_root, target_path, backlog_path)

    issues: list[str] = []
    notes: list[str] = list(replay_notes) + list(rerun_notes)
    review_status = "approved"
    reason_code = "approved"
    target_state_override = None

    if test_failures:
        review_status = "blocked"
        reason_code = "blocked_reviewer_replay"
        issues.append("One or more reported gate commands failed when rerun by the reviewer.")

    target_state_after = executor_result.get("targetStatusAfterIteration", "in_progress")
    labels = status_snapshot["target"].get("labels", [])
    if (
        review_status == "approved"
        and target_state_after == "done"
        and "needs-seam-introduction" in labels
    ):
        target_state_override = "in_progress"
        reason_code = "approved_target_state_override"
        notes.append(
            "Fresh classifier still marks the target as needs-seam-introduction; reviewer downgraded the executor's done claim to in_progress."
        )

    changed_files = executor_result.get("filesChanged", [])
    missing_files: list[str] = []
    ignored_missing_files: list[str] = []
    for changed_file in changed_files:
        resolved_path, relative_display = resolve_changed_file_path(repo_root, changed_file)
        if resolved_path.exists():
            continue
        if is_generated_staging_artifact(relative_display):
            ignored_missing_files.append(relative_display)
            continue
        missing_files.append(relative_display)

    if ignored_missing_files:
        notes.extend(
            f"Ignored missing generated staging artifact reported by executor: {path}"
            for path in ignored_missing_files
        )

    if missing_files:
        review_status = "needs_changes"
        reason_code = "needs_changes_missing_changed_file"
        issues.append("Executor reported changed files that do not exist in the repo root.")
        notes.extend(missing_files)

    if review_status == "approved" and target_state_override == "in_progress":
        summary = (
            "Reviewer reran reported checks, accepted the checkpoint, and downgraded the target state to in_progress."
        )
    elif review_status == "approved":
        summary = "Reviewer reran reported checks and accepted the checkpoint."
    elif review_status == "blocked":
        summary = "Reviewer blocked the checkpoint because rerun verification failed."
    else:
        summary = "Reviewer found inconsistencies in the checkpoint report."

    print(json.dumps(
        {
            "status": review_status,
            "summary": summary,
            "reasonCode": reason_code,
            "testsRun": passed_tests,
            "verification": {
                "replayedCommands": replayed_commands,
                "legacyFallbackUsed": used_legacy_fallback,
            },
            "targetStateOverride": target_state_override,
            "issues": issues,
            "notes": notes,
            "testFailures": test_failures,
            "statusSnapshot": status_snapshot,
        },
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
