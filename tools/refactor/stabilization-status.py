#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import re
from dataclasses import asdict, dataclass


FUNC_START_RE = re.compile(r"^([A-Za-z0-9_.]+)\s*(<-|=)\s*function\s*\(")
HEADING_RE = re.compile(r"^###\s+(\d+)\.\s+(.+?)\s*$")
ABS_R_LINK_RE = re.compile(r"\((/[^)]+/R/[^):]+\.R):\d+\)")


@dataclass
class FileMetrics:
    path: str
    lines: int
    top_level_functions: int
    max_function_lines: int
    module_servers: int
    observers: int
    renderers: int
    labels: list[str]
    next_step: str


@dataclass
class BacklogBucket:
    bucket_number: int
    title: str
    status: str
    primary_target: str | None
    file_paths: list[str]
    is_proteomics: bool


def read_lines(path: pathlib.Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def function_spans(lines: list[str]) -> list[int]:
    spans: list[int] = []
    index = 0
    while index < len(lines):
        if not FUNC_START_RE.match(lines[index].strip()):
            index += 1
            continue
        depth = lines[index].count("{") - lines[index].count("}")
        start = index
        index += 1
        while index < len(lines):
            depth += lines[index].count("{") - lines[index].count("}")
            if depth <= 0:
                spans.append(index - start + 1)
                index += 1
                break
            index += 1
        else:
            spans.append(len(lines) - start)
    return spans


def classify_file(path: pathlib.Path) -> FileMetrics:
    lines = read_lines(path)
    spans = function_spans(lines)
    text = "\n".join(lines)
    module_servers = len(re.findall(r"\bmoduleServer\s*\(", text))
    observers = len(re.findall(r"\bobserveEvent\s*\(", text))
    renderers = len(re.findall(r"\brender(UI|DT|Text|Plot|Table)\s*\(", text))

    labels: list[str] = []
    if module_servers and (observers >= 3 or renderers >= 3):
        labels.extend(["high-risk-wrapper", "needs-seam-introduction"])
    elif spans and max(spans, default=0) >= 250:
        labels.append("review")
    if spans and len(spans) >= 8 and module_servers == 0:
        labels.append("direct-extraction-ready")
    if not labels:
        labels.append("direct-extraction-ready" if module_servers == 0 else "review")

    if "needs-seam-introduction" in labels:
        next_step = "Introduce one top-level helper seam in the same file and rerun the focused gate."
    elif "direct-extraction-ready" in labels:
        next_step = "Draft a manifest, verify it, and stage the wave without rewriting live sources."
    else:
        next_step = "Add focused characterization before structural edits."

    return FileMetrics(
        path=str(path),
        lines=len(lines),
        top_level_functions=len(spans),
        max_function_lines=max(spans, default=0),
        module_servers=module_servers,
        observers=observers,
        renderers=renderers,
        labels=labels,
        next_step=next_step,
    )


def helper_family_files(target: pathlib.Path) -> list[pathlib.Path]:
    pattern = f"{target.stem}_*.R"
    return sorted(
        path for path in target.parent.glob(pattern)
        if path.resolve() != target.resolve()
    )


def parse_backlog(backlog_path: pathlib.Path) -> list[BacklogBucket]:
    lines = read_lines(backlog_path)
    buckets: list[BacklogBucket] = []
    current_match: re.Match[str] | None = None
    current_block: list[str] = []

    def flush_bucket(match: re.Match[str] | None, block_lines: list[str]) -> None:
        if match is None:
            return
        block_text = "\n".join(block_lines)
        title = match.group(2).strip()

        status = "pending"
        lower_block = block_text.lower()
        if "completed in live `r/`" in lower_block or "no longer a blocker" in lower_block:
            status = "done"
        elif "next active target" in lower_block:
            status = "in_progress"

        file_paths = list(dict.fromkeys(ABS_R_LINK_RE.findall(block_text)))
        primary_target = file_paths[0] if file_paths else None
        buckets.append(
            BacklogBucket(
                bucket_number=int(match.group(1)),
                title=title,
                status=status,
                primary_target=primary_target,
                file_paths=file_paths,
                is_proteomics="proteomics" in title.lower(),
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
    return buckets


def summarize_buckets(buckets: list[BacklogBucket], proteomics_only: bool = False) -> dict[str, int]:
    selected = [bucket for bucket in buckets if (bucket.is_proteomics or not proteomics_only)]
    if proteomics_only:
        selected = [bucket for bucket in buckets if bucket.is_proteomics]

    counts = {"total": len(selected), "done": 0, "in_progress": 0, "pending": 0, "blocked": 0}
    for bucket in selected:
        counts[bucket.status] = counts.get(bucket.status, 0) + 1
    return counts


def pct(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return round((numerator / denominator) * 100, 1)


def build_progress_line(target: pathlib.Path, status: dict) -> str:
    return (
        "GODMOD_PROGRESS "
        f"target={target.as_posix()} "
        f"target_pct={status['progress']['target_refactored_pct_estimate']} "
        f"repo_pct={status['progress']['repo_completed_pct_estimate']} "
        f"proteomics_pct={status['progress']['proteomics_completed_pct_estimate']} "
        f"refactored_functions={status['family']['refactored_top_level_functions']} "
        f"legacy_functions={status['family']['legacy_top_level_functions']}"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument(
        "--backlog",
        default="tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md",
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    target = pathlib.Path(args.target).resolve()
    backlog_path = pathlib.Path(args.backlog).resolve()

    target_metrics = classify_file(target)
    helper_files = helper_family_files(target)
    helper_metrics = [classify_file(path) for path in helper_files]
    refactored_functions = sum(item.top_level_functions for item in helper_metrics)
    refactored_lines = sum(item.lines for item in helper_metrics)

    buckets = parse_backlog(backlog_path)
    repo_counts = summarize_buckets(buckets, proteomics_only=False)
    proteomics_counts = summarize_buckets(buckets, proteomics_only=True)

    family_total_functions = target_metrics.top_level_functions + refactored_functions
    family_total_lines = target_metrics.lines + refactored_lines

    status = {
        "target": asdict(target_metrics),
        "family": {
            "helper_files": [item.path for item in helper_metrics],
            "helper_file_count": len(helper_metrics),
            "refactored_top_level_functions": refactored_functions,
            "refactored_lines": refactored_lines,
            "legacy_top_level_functions": target_metrics.top_level_functions,
            "legacy_lines": target_metrics.lines,
            "family_total_top_level_functions": family_total_functions,
            "family_total_lines": family_total_lines,
        },
        "backlog": {
            "path": str(backlog_path),
            "repo": repo_counts,
            "proteomics": proteomics_counts,
            "active_bucket": next(
                (
                    {
                        "bucket_number": bucket.bucket_number,
                        "title": bucket.title,
                        "primary_target": bucket.primary_target,
                    }
                    for bucket in buckets
                    if bucket.status == "in_progress"
                ),
                None,
            ),
        },
        "progress": {
            "target_refactored_pct_estimate": pct(refactored_functions, family_total_functions),
            "target_legacy_pct_estimate": pct(target_metrics.top_level_functions, family_total_functions),
            "repo_completed_pct_estimate": pct(repo_counts["done"], repo_counts["total"]),
            "proteomics_completed_pct_estimate": pct(proteomics_counts["done"], proteomics_counts["total"]),
        },
        "estimation_basis": {
            "target": "top-level function distribution across the target file and sibling helper files",
            "backlog": "bucket completion counts derived from backlog current-state markers",
        },
    }
    status["progress_line"] = build_progress_line(target, status)

    if args.json:
        print(json.dumps(status, indent=2))
    else:
        print(f"target: {target_metrics.path}")
        print(f"legacy_top_level_functions: {target_metrics.top_level_functions}")
        print(f"refactored_top_level_functions: {refactored_functions}")
        print(f"target_refactored_pct_estimate: {status['progress']['target_refactored_pct_estimate']}")
        print(f"repo_completed_pct_estimate: {status['progress']['repo_completed_pct_estimate']}")
        print(f"proteomics_completed_pct_estimate: {status['progress']['proteomics_completed_pct_estimate']}")
        print(status["progress_line"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
