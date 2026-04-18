#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import pathlib
import selectors
import shutil
import subprocess
import sys
import tempfile

from stabilization_schema import OpenAISchemaCompatibilityError, validate_openai_output_schema

SCHEMA_SOURCE = pathlib.Path(__file__).resolve().parent / "stabilization-executor-output.schema.json"
SESSION_INIT_RETRY_MARKERS = (
    "failed to create session",
    "failed to initialize session",
    "thread/start failed",
)


def choose_runtime_parent(repo_root: pathlib.Path) -> pathlib.Path:
    override = os.environ.get("STABILIZATION_RUNNER_TMP_ROOT")
    candidates: list[pathlib.Path] = []
    if override:
        candidates.append(pathlib.Path(override).expanduser())
    candidates.append(repo_root / ".god-module-stabilization" / "runner-runtime")
    candidates.append(pathlib.Path("/tmp"))
    fallback = pathlib.Path(tempfile.gettempdir())
    candidates.append(fallback)

    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        try:
            candidate.mkdir(parents=True, exist_ok=True)
        except OSError:
            continue
        if os.access(candidate, os.W_OK | os.X_OK):
            return candidate
    raise RuntimeError("no writable runtime parent available for stabilization codex runner")


def real_codex_home() -> pathlib.Path:
    raw = os.environ.get("CODEX_HOME")
    if raw:
        return pathlib.Path(raw).expanduser().resolve()
    return (pathlib.Path.home() / ".codex").resolve()


def mirror_codex_home(runtime_root: pathlib.Path) -> pathlib.Path:
    source_root = real_codex_home()
    target_root = runtime_root / "codex-home"
    target_root.mkdir(parents=True, exist_ok=True)

    # Keep the installed skills/config readable while making sessions writable
    # inside the isolated runtime root so nested codex exec can initialize its
    # own thread without relying on the parent's read-only CODEX_HOME.
    link_entries = (
        ".personality_migration",
        "auth.json",
        "config.toml",
        "installation_id",
        "models_cache.json",
        "rules",
        "shell_snapshots",
        "skills",
        "version.json",
    )
    for name in link_entries:
        source = source_root / name
        target = target_root / name
        if not source.exists() or target.exists():
            continue
        try:
            os.symlink(source, target, target_is_directory=source.is_dir())
        except OSError:
            if source.is_dir():
                shutil.copytree(source, target, dirs_exist_ok=True)
            else:
                shutil.copy2(source, target)

    for name in ("sessions", "tmp", ".tmp", "cache", "log", "memories"):
        (target_root / name).mkdir(parents=True, exist_ok=True)

    return target_root


def build_child_env(runtime_root: pathlib.Path) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("CODEX_THREAD_ID", None)
    env.pop("CODEX_SANDBOX_NETWORK_DISABLED", None)
    env["TMPDIR"] = str(runtime_root)
    env["CODEX_HOME"] = str(mirror_codex_home(runtime_root))
    env.setdefault("HOME", str(pathlib.Path.home()))
    env.setdefault("XDG_RUNTIME_DIR", "/tmp")
    return env


def should_retry_session_init(stderr_text: str, returncode: int) -> bool:
    if returncode == 0:
        return False
    lowered = stderr_text.lower()
    return any(marker in lowered for marker in SESSION_INIT_RETRY_MARKERS)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("payload")
    args = parser.parse_args()

    payload_path = pathlib.Path(args.payload).resolve()
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    repo_root = pathlib.Path(payload["repoRoot"]).resolve()
    prompt = payload["prompt"]
    model = payload.get("model")
    sandbox = payload.get("sandbox", "workspace-write")
    approval_policy = payload.get("approvalPolicy", "never")
    network_access_enabled = payload.get("networkAccessEnabled", True)
    schema_source = pathlib.Path(
        os.environ.get("STABILIZATION_EXECUTOR_SCHEMA_PATH", str(SCHEMA_SOURCE))
    ).resolve()

    try:
        schema = json.loads(schema_source.read_text(encoding="utf-8"))
        validate_openai_output_schema(schema)
    except FileNotFoundError as error:
        print(f"protocol schema validation failed: schema file not found: {error.filename}", file=sys.stderr)
        return 2
    except json.JSONDecodeError as error:
        print(f"protocol schema validation failed: invalid JSON in schema: {error}", file=sys.stderr)
        return 2
    except OpenAISchemaCompatibilityError as error:
        print(f"protocol schema validation failed: {error}", file=sys.stderr)
        return 2

    runtime_root = pathlib.Path(
        tempfile.mkdtemp(prefix="godmod-codex-runner-", dir=str(choose_runtime_parent(repo_root)))
    )
    try:
        schema_path = runtime_root / "executor-schema.json"
        output_path = runtime_root / "last-message.json"
        schema_path.write_text(json.dumps(schema, indent=2), encoding="utf-8")

        cmd = [
            "codex",
            "exec",
            "--cd",
            str(repo_root),
            "--skip-git-repo-check",
            "--sandbox",
            sandbox,
            "--config",
            f"approval_policy={json.dumps(approval_policy)}",
            "--config",
            f"sandbox_workspace_write.network_access={'true' if network_access_enabled else 'false'}",
            "--output-schema",
            str(schema_path),
            "--output-last-message",
            str(output_path),
            "-",
        ]
        if model:
            cmd.extend(["--model", model])

        attempts = 2
        last_error: subprocess.CalledProcessError | None = None
        for attempt in range(1, attempts + 1):
            process = subprocess.Popen(
                cmd,
                cwd=repo_root,
                env=build_child_env(runtime_root),
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
            )
            assert process.stdin is not None
            assert process.stdout is not None
            assert process.stderr is not None
            process.stdin.write(prompt)
            process.stdin.close()

            selector = selectors.DefaultSelector()
            selector.register(process.stdout, selectors.EVENT_READ, data="stdout")
            selector.register(process.stderr, selectors.EVENT_READ, data="stderr")
            captured_stdout: list[str] = []
            captured_stderr: list[str] = []

            while selector.get_map():
                for key, _ in selector.select():
                    stream = key.fileobj
                    chunk = stream.readline()
                    if chunk == "":
                        selector.unregister(stream)
                        continue
                    if key.data == "stdout":
                        captured_stdout.append(chunk)
                        print(chunk, file=sys.stderr, end="")
                    else:
                        captured_stderr.append(chunk)
                        print(chunk, file=sys.stderr, end="")

            returncode = process.wait()
            stderr_text = "".join(captured_stderr)
            stdout_text = "".join(captured_stdout)

            if returncode == 0:
                print(output_path.read_text(encoding="utf-8"))
                return 0

            if returncode < 0:
                signal_name = f"signal {-returncode}"
                last_error = subprocess.CalledProcessError(
                    returncode,
                    cmd,
                    output=stdout_text,
                    stderr=f"codex exec terminated by {signal_name}\n{stderr_text}",
                )
                break

            last_error = subprocess.CalledProcessError(
                returncode,
                cmd,
                output=stdout_text,
                stderr=stderr_text,
            )
            if attempt < attempts and should_retry_session_init(stderr_text, returncode):
                print("codex runner retrying after session initialization failure", file=sys.stderr)
                continue
            break

        assert last_error is not None
        raise last_error
    finally:
        shutil.rmtree(runtime_root, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
