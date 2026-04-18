#!/usr/bin/env python3
from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile


def run(cmd: list[str]) -> None:
    print("==> " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def run_text(cmd: list[str]) -> str:
    completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return completed.stdout


def manifest_paths(manifest: pathlib.Path, emit_collate: str | None) -> tuple[list[str], list[str], list[str]]:
    r_expr = """
suppressPackageStartupMessages(library(yaml))
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
manifest <- yaml::read_yaml(commandArgs(trailingOnly = TRUE)[1])
sources <- unique(vapply(manifest$entries, function(entry) entry$source %||% "", character(1)))
sources <- sources[nzchar(sources)]
targets <- unique(vapply(Filter(function(entry) !is.null(entry$target), manifest$entries), function(entry) entry$target %||% "", character(1)))
targets <- targets[nzchar(targets)]
cat("[sources]\\n")
if (length(sources)) cat(paste(sources, collapse = "\\n"), "\\n", sep = "")
cat("[targets]\\n")
if (length(targets)) cat(paste(targets, collapse = "\\n"), "\\n", sep = "")
"""
    output = run_text(["Rscript", "-e", r_expr, str(manifest)])
    mode = None
    sources: list[str] = []
    targets: list[str] = []
    for line in output.splitlines():
        if line == "[sources]":
            mode = "sources"
            continue
        if line == "[targets]":
            mode = "targets"
            continue
        if not line.strip():
            continue
        if mode == "sources":
            sources.append(line.strip())
        elif mode == "targets":
            targets.append(line.strip())
    extras = [emit_collate] if emit_collate else []
    return sources, targets, extras


def top_level_symbols(path: pathlib.Path) -> list[str]:
    if not path.exists():
        return []
    r_expr = """
normalize_call_name <- function(expr) {
  if (is.null(expr) || !is.call(expr)) return(NULL)
  head <- expr[[1]]
  if (is.symbol(head)) return(as.character(head))
  NULL
}

extract_top_level_symbol <- function(expr) {
  call_name <- normalize_call_name(expr)
  if (is.call(expr) && length(expr) >= 3 && !is.null(call_name) && call_name %in% c("<-", "=")) {
    lhs <- expr[[2]]
    if (is.symbol(lhs)) return(as.character(lhs))
  }
  NULL
}

exprs <- parse(file = commandArgs(trailingOnly = TRUE)[1], keep.source = TRUE)
symbols <- Filter(Negate(is.null), lapply(exprs, extract_top_level_symbol))
if (length(symbols)) cat(paste(unlist(symbols), collapse = "\\n"))
"""
    output = run_text(["Rscript", "-e", r_expr, str(path)])
    return [line.strip() for line in output.splitlines() if line.strip()]


def backup_paths(repo_root: pathlib.Path, relative_paths: list[str]) -> tuple[pathlib.Path, dict[str, bool]]:
    backup_root = pathlib.Path(tempfile.mkdtemp(prefix="godmod-apply-backup-"))
    existed: dict[str, bool] = {}
    for relative in relative_paths:
        src = repo_root / relative
        existed[relative] = src.exists()
        if src.exists():
            dst = backup_root / relative
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
    return backup_root, existed


def restore_paths(repo_root: pathlib.Path, backup_root: pathlib.Path, existed: dict[str, bool]) -> None:
    for relative, did_exist in existed.items():
        dst = repo_root / relative
        backup = backup_root / relative
        if did_exist:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(backup, dst)
        elif dst.exists():
            dst.unlink()


def assert_target_preservation(repo_root: pathlib.Path, target_symbols_before: dict[str, list[str]]) -> None:
    failures: list[str] = []
    for relative, before_symbols in target_symbols_before.items():
        if not before_symbols:
            continue
        after_symbols = set(top_level_symbols(repo_root / relative))
        missing = [symbol for symbol in before_symbols if symbol not in after_symbols]
        if missing:
            failures.append(f"{relative} lost previously present top-level functions: {', '.join(missing)}")
    if failures:
        raise RuntimeError("\n".join(failures))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--emit-collate")
    parser.add_argument("--skip-check", action="store_true")
    args = parser.parse_args()

    manifest = pathlib.Path(args.manifest)
    if not manifest.exists():
        raise SystemExit(f"missing manifest: {manifest}")

    repo_root = pathlib.Path.cwd()
    sources, targets, extras = manifest_paths(manifest, args.emit_collate)
    touched_paths = sorted(set(sources + targets + extras))
    backup_root, existed = backup_paths(repo_root, touched_paths)
    target_symbols_before = {
        relative: top_level_symbols(repo_root / relative)
        for relative in targets
        if (repo_root / relative).exists()
    }

    try:
        run(["Rscript", "tools/refactor/verify_refactor.R", "--manifest", args.manifest])

        cmd = [
            "Rscript",
            "tools/refactor/extract_blocks.R",
            "--manifest",
            args.manifest,
            "--write-targets",
            "--preserve-existing-targets",
            "--rewrite-sources",
        ]
        if args.emit_collate:
            cmd.extend(["--emit-collate", args.emit_collate])
        run(cmd)

        assert_target_preservation(repo_root, target_symbols_before)

        if not args.skip_check:
            run(["Rscript", "tools/refactor/check_wave_apply.R", "--manifest", args.manifest])
    except Exception as exc:
        restore_paths(repo_root, backup_root, existed)
        print(
            "apply_wave failed and repo files were restored from backup. Confirm the manifest was reviewed, "
            "still matches current sources, and preserves previously applied helper targets.",
            file=sys.stderr,
        )
        if isinstance(exc, subprocess.CalledProcessError):
            return exc.returncode
        print(str(exc), file=sys.stderr)
        return 1

    print("wave apply complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
