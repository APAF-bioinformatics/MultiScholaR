#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import io
import json
import pathlib
import re
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone


SCHEMA_VERSION = 5
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
DEFAULT_AUDIT_DIR = ".refactor-fidelity-audit"
DEFAULT_EVENTS_PATH = "events.jsonl"
DEFAULT_DB_PATH = "audit.db"
DEFAULT_REPORTS_DIR = "reports"
DEFAULT_SUMMARY_JSON = "latest-summary.json"
DEFAULT_SUMMARY_MD = "latest-summary.md"
DEFAULT_EXCEPTIONS_JSON = "latest-exceptions.json"
DEFAULT_EXCEPTIONS_MD = "latest-exceptions.md"
DEFAULT_CLOSEOUT_JSON = "latest-closeout.json"
DEFAULT_CLOSEOUT_MD = "latest-closeout.md"
DEFAULT_BEHAVIOR_CASE_GLOB = "behavior-cases*.json"
DEFAULT_TEST_FILE_GLOB = "test-*.R"

SURFACE_OPEN_DIFFS = {
    "missing_definition",
    "extra_definition",
    "duplicate_definition",
    "definition_drift",
    "missing_export",
    "extra_export",
    "missing_file",
    "extra_file",
    "collate_drift",
}

DIFF_SEVERITY = {
    "missing_definition": "high",
    "extra_definition": "medium",
    "duplicate_definition": "high",
    "definition_drift": "medium",
    "moved_definition": "info",
    "missing_export": "high",
    "extra_export": "medium",
    "missing_file": "medium",
    "extra_file": "low",
    "collate_drift": "medium",
}

TEST_FAMILY_SPECS = (
    {
        "family": "module_contract",
        "suffix": "module-contracts",
        "surface": "demonolithed_wrappers_and_modules",
        "certainty_target": "T5_contract_replayed",
    },
    {
        "family": "characterization",
        "suffix": "characterization",
        "surface": "cross_module_and_workflow_segments",
        "certainty_target": "T5_contract_replayed",
    },
    {
        "family": "compat",
        "suffix": "compat",
        "surface": "cross_module_and_workflow_segments",
        "certainty_target": "T5_contract_replayed",
    },
    {
        "family": "golden_master",
        "suffix": "golden-master",
        "surface": "cross_module_and_workflow_segments",
        "certainty_target": "T6_end_to_end_verified",
    },
)

DEFAULT_TEST_FAMILY = {
    "family": "general_unit",
    "surface": "unspecified_general",
    "certainty_target": "T4_behavior_replayed",
}

MODULE_CONTRACT_SCENARIOS = (
    "initialization",
    "happy_path",
    "missing_prerequisite",
    "invalid_input",
    "state_restore",
    "reset_cleanup",
    "export_report_side_effects",
    "error_notification",
)

SCENARIO_PATTERNS = (
    ("export_report_side_effects", (r"\bdownload\b", r"\bexport\b", r"\breport\b", r"\barchive\b", r"\bgithub push\b", r"\bcopy\b")),
    ("state_restore", (r"\brevert", r"\brestore", r"\brestored\b", r"\bstate restore\b")),
    ("reset_cleanup", (r"\breset\b", r"\bcleanup\b", r"\bclear", r"\bclose", r"\bremove", r"\bclean")),
    ("invalid_input", (r"\binvalid\b", r"\bunsupported\b", r"\bmissing inputs?\b", r"\bmissing .* column\b", r"\berrors clearly\b")),
    ("missing_prerequisite", (r"\bwithout\b", r"\bmissing\b", r"\babsent\b", r"\bno prior\b", r"\bno valid\b", r"\bfalls back\b", r"\bskip", r"\bunavailable\b", r"\bbefore analysis\b")),
    ("error_notification", (r"\berror\b", r"\bwarn", r"\bnotification\b", r"\bfailure\b")),
    ("initialization", (r"\binitialize", r"\binitial", r"\bbootstrap\b", r"\bdefault\b", r"\brender", r"\bseed", r"\bpublic create contract\b")),
    ("happy_path", (r"\breturn", r"\bstore", r"\bapply", r"\bdelegate", r"\bupdate", r"\bformat", r"\bhydrate", r"\bsave", r"\bcreate", r"\bbuild", r"\bcomplete")),
)

TERMINAL_CURATION_STATUSES = {"auto_curated", "accepted", "resolved", "dismissed"}
CURATION_STATUS_CHOICES = (
    "candidate",
    "triaged",
    "auto_curated",
    "accepted",
    "resolved",
    "dismissed",
)

AUTO_CURATION_RULES = {
    "whitespace_only_drift": {
        "curation_status": "auto_curated",
        "disposition": "equivalent_whitespace_only",
        "note": "Auto-curated because only whitespace normalization changed.",
        "rule_name": "auto:whitespace_only_equivalent",
    },
    "comment_only_drift": {
        "curation_status": "auto_curated",
        "disposition": "equivalent_comment_only",
        "note": "Auto-curated because only comments or roxygen changed.",
        "rule_name": "auto:comment_only_equivalent",
    },
    "target_resolution_asymmetry": {
        "curation_status": "auto_curated",
        "disposition": "equivalent_target_resolution_asymmetry",
        "note": "Auto-curated because content matched exactly after fallback target resolution.",
        "rule_name": "auto:target_resolution_asymmetry",
    },
}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: pathlib.Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def append_jsonl(path: pathlib.Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(value) + "\n")


def write_text(path: pathlib.Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def connect_db(path: pathlib.Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def table_columns(conn: sqlite3.Connection, table_name: str) -> set[str]:
    rows = conn.execute(f"PRAGMA table_info({table_name})").fetchall()
    return {str(row["name"]) for row in rows}


def ensure_column(conn: sqlite3.Connection, table_name: str, column_name: str, definition: str) -> None:
    if column_name in table_columns(conn, table_name):
        return
    conn.execute(f"ALTER TABLE {table_name} ADD COLUMN {column_name} {definition}")


def create_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS audit_runs (
            run_id TEXT PRIMARY KEY,
            started_at TEXT NOT NULL,
            completed_at TEXT,
            repo_root TEXT NOT NULL,
            baseline_ref TEXT,
            target_ref TEXT,
            mode TEXT NOT NULL,
            status TEXT NOT NULL,
            summary_json TEXT
        );

        CREATE TABLE IF NOT EXISTS inventory_entities (
            entity_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            side TEXT NOT NULL,
            entity_key TEXT NOT NULL,
            entity_kind TEXT NOT NULL,
            entity_name TEXT NOT NULL,
            signature_key TEXT,
            file_path TEXT NOT NULL,
            line_start INTEGER,
            line_end INTEGER,
            exported INTEGER NOT NULL,
            collate_index INTEGER,
            hash_raw TEXT,
            hash_normalized TEXT,
            hash_ast TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE INDEX IF NOT EXISTS idx_inventory_entities_run_side
            ON inventory_entities(run_id, side);

        CREATE INDEX IF NOT EXISTS idx_inventory_entities_key
            ON inventory_entities(entity_key);

        CREATE TABLE IF NOT EXISTS inventory_files (
            file_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            side TEXT NOT NULL,
            file_path TEXT NOT NULL,
            collate_index INTEGER,
            present INTEGER NOT NULL,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE INDEX IF NOT EXISTS idx_inventory_files_run_side
            ON inventory_files(run_id, side);

        CREATE TABLE IF NOT EXISTS inventory_exports (
            export_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            side TEXT NOT NULL,
            export_kind TEXT NOT NULL,
            export_name TEXT NOT NULL,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE INDEX IF NOT EXISTS idx_inventory_exports_run_side
            ON inventory_exports(run_id, side);

        CREATE TABLE IF NOT EXISTS inventory_diffs (
            diff_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            entity_kind TEXT NOT NULL,
            entity_key TEXT NOT NULL,
            diff_class TEXT NOT NULL,
            baseline_entity_id TEXT,
            target_entity_id TEXT,
            severity TEXT,
            status TEXT NOT NULL,
            evidence_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE TABLE IF NOT EXISTS manifest_entries (
            manifest_entry_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            manifest_path TEXT NOT NULL,
            entry_id TEXT NOT NULL,
            selector_kind TEXT NOT NULL,
            selector_payload_json TEXT NOT NULL,
            source_path TEXT NOT NULL,
            target_path TEXT,
            group_name TEXT,
            action TEXT NOT NULL,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE TABLE IF NOT EXISTS manifest_comparisons (
            comparison_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            manifest_entry_id TEXT NOT NULL,
            baseline_source_resolver TEXT,
            target_resolver TEXT,
            raw_exact INTEGER,
            normalized_text_exact INTEGER,
            ast_exact INTEGER,
            hash_raw_baseline TEXT,
            hash_raw_target TEXT,
            hash_normalized_baseline TEXT,
            hash_normalized_target TEXT,
            hash_ast_baseline TEXT,
            hash_ast_target TEXT,
            status TEXT NOT NULL,
            evidence_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(manifest_entry_id) REFERENCES manifest_entries(manifest_entry_id)
        );

        CREATE TABLE IF NOT EXISTS behavior_cases (
            behavior_case_id TEXT PRIMARY KEY,
            family TEXT NOT NULL,
            entity_key TEXT NOT NULL,
            fixture_kind TEXT NOT NULL,
            fixture_path TEXT,
            normalizer_key TEXT,
            risk_level TEXT,
            catalog_path TEXT,
            case_json TEXT
        );

        CREATE INDEX IF NOT EXISTS idx_behavior_cases_entity_key
            ON behavior_cases(entity_key);

        CREATE TABLE IF NOT EXISTS behavior_results (
            behavior_result_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            behavior_case_id TEXT NOT NULL,
            baseline_status TEXT,
            target_status TEXT,
            normalized_equal INTEGER,
            warning_equal INTEGER,
            message_equal INTEGER,
            error_equal INTEGER,
            artifact_equal INTEGER,
            result_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(behavior_case_id) REFERENCES behavior_cases(behavior_case_id)
        );

        CREATE INDEX IF NOT EXISTS idx_behavior_results_run
            ON behavior_results(run_id);

        CREATE TABLE IF NOT EXISTS testing_matrix_entries (
            testing_matrix_entry_id TEXT PRIMARY KEY,
            surface TEXT NOT NULL,
            primary_evidence_json TEXT NOT NULL,
            fixture_strategy_json TEXT,
            required_scenarios_json TEXT,
            acceptance_json TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS contract_test_files (
            contract_test_file_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            file_path TEXT NOT NULL,
            family TEXT NOT NULL,
            module_family TEXT NOT NULL,
            surface TEXT NOT NULL,
            certainty_target TEXT NOT NULL,
            case_count INTEGER NOT NULL,
            metadata_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE INDEX IF NOT EXISTS idx_contract_test_files_run
            ON contract_test_files(run_id);

        CREATE TABLE IF NOT EXISTS contract_test_cases (
            contract_test_case_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            contract_test_file_id TEXT NOT NULL,
            test_name TEXT NOT NULL,
            entity_hint TEXT,
            primary_scenario TEXT,
            scenario_tags_json TEXT NOT NULL,
            certainty_target TEXT NOT NULL,
            metadata_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(contract_test_file_id) REFERENCES contract_test_files(contract_test_file_id)
        );

        CREATE INDEX IF NOT EXISTS idx_contract_test_cases_run
            ON contract_test_cases(run_id);

        CREATE TABLE IF NOT EXISTS contract_results (
            contract_result_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            contract_test_file_id TEXT NOT NULL,
            execution_mode TEXT NOT NULL,
            status TEXT NOT NULL,
            test_count INTEGER,
            failure_count INTEGER,
            skip_count INTEGER,
            result_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(contract_test_file_id) REFERENCES contract_test_files(contract_test_file_id)
        );

        CREATE INDEX IF NOT EXISTS idx_contract_results_run
            ON contract_results(run_id);

        CREATE TABLE IF NOT EXISTS exceptions (
            exception_id TEXT PRIMARY KEY,
            exception_key TEXT,
            run_id TEXT NOT NULL,
            audit_layer TEXT NOT NULL,
            entity_key TEXT NOT NULL,
            severity TEXT NOT NULL,
            exception_type TEXT NOT NULL,
            reason TEXT NOT NULL,
            baseline_ref TEXT,
            target_ref TEXT,
            evidence_json TEXT,
            curation_status TEXT NOT NULL,
            disposition TEXT,
            owner TEXT,
            review_note TEXT,
            updated_at TEXT,
            curation_rule TEXT,
            resolved_at TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        """
    )
    ensure_column(conn, "inventory_diffs", "evidence_json", "TEXT")
    ensure_column(conn, "manifest_comparisons", "hash_normalized_baseline", "TEXT")
    ensure_column(conn, "manifest_comparisons", "hash_normalized_target", "TEXT")
    ensure_column(conn, "manifest_comparisons", "evidence_json", "TEXT")
    ensure_column(conn, "behavior_cases", "catalog_path", "TEXT")
    ensure_column(conn, "behavior_cases", "case_json", "TEXT")
    ensure_column(conn, "behavior_results", "message_equal", "INTEGER")
    ensure_column(conn, "contract_test_files", "certainty_target", "TEXT")
    ensure_column(conn, "exceptions", "exception_key", "TEXT")
    ensure_column(conn, "exceptions", "baseline_ref", "TEXT")
    ensure_column(conn, "exceptions", "target_ref", "TEXT")
    ensure_column(conn, "exceptions", "review_note", "TEXT")
    ensure_column(conn, "exceptions", "updated_at", "TEXT")
    ensure_column(conn, "exceptions", "curation_rule", "TEXT")
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_exceptions_key ON exceptions(exception_key)"
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_exceptions_status ON exceptions(curation_status)"
    )
    conn.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
    conn.commit()


def resolve_paths(repo_root: pathlib.Path, audit_dir: str) -> dict[str, pathlib.Path]:
    audit_root = repo_root / audit_dir
    reports_dir = audit_root / DEFAULT_REPORTS_DIR
    return {
        "audit_root": audit_root,
        "db_path": audit_root / DEFAULT_DB_PATH,
        "events_path": audit_root / DEFAULT_EVENTS_PATH,
        "reports_dir": reports_dir,
        "summary_json": reports_dir / DEFAULT_SUMMARY_JSON,
        "summary_md": reports_dir / DEFAULT_SUMMARY_MD,
        "exceptions_json": reports_dir / DEFAULT_EXCEPTIONS_JSON,
        "exceptions_md": reports_dir / DEFAULT_EXCEPTIONS_MD,
        "closeout_json": reports_dir / DEFAULT_CLOSEOUT_JSON,
        "closeout_md": reports_dir / DEFAULT_CLOSEOUT_MD,
    }


def bootstrap_store(repo_root: pathlib.Path, audit_dir: str, *, emit_event: bool) -> dict[str, pathlib.Path]:
    paths = resolve_paths(repo_root, audit_dir)
    paths["audit_root"].mkdir(parents=True, exist_ok=True)
    paths["reports_dir"].mkdir(parents=True, exist_ok=True)
    conn = connect_db(paths["db_path"])
    try:
        create_schema(conn)
    finally:
        conn.close()

    if not paths["events_path"].exists():
        paths["events_path"].write_text("", encoding="utf-8")

    if not paths["exceptions_json"].exists():
        write_json(paths["exceptions_json"], {"exceptions": []})
    if not paths["exceptions_md"].exists():
        write_text(paths["exceptions_md"], "# Exceptions\n\nNo exceptions recorded.\n")
    if not paths["summary_json"].exists():
        write_json(
            paths["summary_json"],
            {
                "status": "bootstrap_ready",
                "schema_version": SCHEMA_VERSION,
                "repo_root": str(repo_root),
            },
        )
    if not paths["summary_md"].exists():
        write_text(
            paths["summary_md"],
            "# Refactor Fidelity Audit\n\nBootstrap complete. No runs recorded yet.\n",
        )
    if not paths["closeout_json"].exists():
        write_json(
            paths["closeout_json"],
            {
                "status": "bootstrap_ready",
                "mode": "closeout",
                "repo_root": str(repo_root),
            },
        )
    if not paths["closeout_md"].exists():
        write_text(
            paths["closeout_md"],
            "# Refactor Fidelity Closeout\n\nBootstrap complete. No closeout run recorded yet.\n",
        )

    if emit_event:
        append_jsonl(
            paths["events_path"],
            {
                "timestamp": now_iso(),
                "phase": "bootstrap",
                "status": "ready",
                "repo_root": str(repo_root),
                "artifacts": {
                    "audit_root": str(paths["audit_root"]),
                    "db_path": str(paths["db_path"]),
                    "events_path": str(paths["events_path"]),
                    "reports_dir": str(paths["reports_dir"]),
                },
            },
        )

    return paths


def ensure_list(value) -> list:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def render_inventory_summary_markdown(summary: dict) -> str:
    counts = summary["counts"]
    exports = counts["exports"]
    entity_counts = counts.get("entities", {})
    lines = [
        "# Refactor Fidelity Audit",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Side: `{summary['side']}`",
        f"- Target ref: `{summary['target_ref']}`",
        f"- Status: `{summary['status']}`",
        f"- R files: `{counts['r_files']}`",
        f"- Collate entries: `{counts['collate']}`",
        f"- Exported functions: `{exports['functions']}`",
        f"- Exported methods: `{exports['methods']}`",
        f"- Exported classes: `{exports['classes']}`",
        f"- Parse failures: `{counts['parse_failures']}`",
    ]
    if entity_counts:
        lines.append("- Entity counts:")
        for kind in sorted(entity_counts):
            lines.append(f"  - `{kind}`: `{entity_counts[kind]}`")
    if summary["parse_failures"]:
        lines.extend(["", "## Parse Failures", ""])
        for failure in summary["parse_failures"]:
            lines.append(f"- {failure}")
    return "\n".join(lines) + "\n"


def render_surface_summary_markdown(summary: dict) -> str:
    diff_counts = summary["diff_counts"]
    severity_counts = summary["severity_counts"]
    lines = [
        "# Refactor Fidelity Audit",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Baseline: `{summary['baseline']['label']}`",
        f"- Target: `{summary['target']['label']}`",
        f"- Status: `{summary['status']}`",
        f"- Total diffs: `{summary['total_diffs']}`",
        f"- Open diffs: `{summary['open_diff_count']}`",
        f"- Exception candidates: `{summary['exception_count']}`",
        f"- Baseline parse failures: `{len(summary['baseline']['parse_failures'])}`",
        f"- Target parse failures: `{len(summary['target']['parse_failures'])}`",
    ]
    if diff_counts:
        lines.extend(["", "## Drift Counts", ""])
        for diff_class in sorted(diff_counts):
            lines.append(f"- `{diff_class}`: `{diff_counts[diff_class]}`")
    if severity_counts:
        lines.extend(["", "## Severity Counts", ""])
        for severity in sorted(severity_counts):
            lines.append(f"- `{severity}`: `{severity_counts[severity]}`")
    if summary["baseline"]["parse_failures"]:
        lines.extend(["", "## Baseline Parse Failures", ""])
        for failure in summary["baseline"]["parse_failures"]:
            lines.append(f"- {failure}")
    if summary["target"]["parse_failures"]:
        lines.extend(["", "## Target Parse Failures", ""])
        for failure in summary["target"]["parse_failures"]:
            lines.append(f"- {failure}")
    if summary["sample_diffs"]:
        lines.extend(["", "## Sample Findings", ""])
        for diff in summary["sample_diffs"]:
            lines.append(
                f"- [{diff['severity']}] `{diff['diff_class']}` `{diff['entity_key']}`"
            )
    return "\n".join(lines) + "\n"


def render_manifest_summary_markdown(summary: dict) -> str:
    status_counts = summary["status_counts"]
    resolver_counts = summary["resolver_counts"]
    lines = [
        "# Refactor Fidelity Audit",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Baseline: `{summary['baseline']}`",
        f"- Target: `{summary['target']}`",
        f"- Status: `{summary['status']}`",
        f"- Manifest files: `{summary['manifest_count']}`",
        f"- Manifest entries: `{summary['entry_count']}`",
        f"- Exceptions emitted: `{summary['exception_count']}`",
    ]
    if status_counts:
        lines.extend(["", "## Comparison Status Counts", ""])
        for status_name in sorted(status_counts):
            lines.append(f"- `{status_name}`: `{status_counts[status_name]}`")
    if resolver_counts:
        lines.extend(["", "## Target Resolver Counts", ""])
        for resolver_name in sorted(resolver_counts):
            lines.append(f"- `{resolver_name}`: `{resolver_counts[resolver_name]}`")
    if summary["sample_comparisons"]:
        lines.extend(["", "## Sample Comparisons", ""])
        for comparison in summary["sample_comparisons"]:
            lines.append(
                f"- `{comparison['status']}` `{comparison['manifest_path']}::{comparison['entry_id']}` via `{comparison['target_resolver']}`"
            )
    return "\n".join(lines) + "\n"


def render_behavior_summary_markdown(summary: dict) -> str:
    status_counts = summary["status_counts"]
    family_counts = summary["family_counts"]
    lines = [
        "# Refactor Fidelity Audit",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Baseline: `{summary['baseline']}`",
        f"- Target: `{summary['target']}`",
        f"- Status: `{summary['status']}`",
        f"- Case files: `{summary['case_file_count']}`",
        f"- Cases: `{summary['case_count']}`",
        f"- Exceptions emitted: `{summary['exception_count']}`",
    ]
    if status_counts:
        lines.extend(["", "## Replay Status Counts", ""])
        for status_name in sorted(status_counts):
            lines.append(f"- `{status_name}`: `{status_counts[status_name]}`")
    if family_counts:
        lines.extend(["", "## Family Counts", ""])
        for family_name in sorted(family_counts):
            lines.append(f"- `{family_name}`: `{family_counts[family_name]}`")
    if summary["sample_results"]:
        lines.extend(["", "## Sample Results", ""])
        for result in summary["sample_results"]:
            lines.append(
                f"- `{result['comparison_status']}` `{result['case_id']}` `{result['entity_key']}`"
            )
    return "\n".join(lines) + "\n"


def render_contract_summary_markdown(summary: dict) -> str:
    lines = [
        "# Refactor Fidelity Audit",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Target: `{summary['target']}`",
        f"- Status: `{summary['status']}`",
        f"- Matrix surfaces: `{summary['matrix_surface_count']}`",
        f"- Test files cataloged: `{summary['test_file_count']}`",
        f"- Test cases cataloged: `{summary['test_case_count']}`",
        f"- Executed files: `{summary['executed_file_count']}`",
        f"- Exceptions emitted: `{summary['exception_count']}`",
    ]
    for title, key in (
        ("Family Counts", "family_counts"),
        ("Module Family Counts", "module_family_counts"),
        ("Surface Counts", "surface_counts"),
        ("Scenario Counts", "scenario_counts"),
        ("Execution Status Counts", "execution_status_counts"),
    ):
        values = summary.get(key, {})
        if not values:
            continue
        lines.extend(["", f"## {title}", ""])
        for name in sorted(values):
            lines.append(f"- `{name}`: `{values[name]}`")
    if summary.get("sample_files"):
        lines.extend(["", "## Sample Files", ""])
        for item in summary["sample_files"]:
            lines.append(
                f"- `{item['family']}` `{item['module_family']}` `{item['file_path']}`"
            )
    return "\n".join(lines) + "\n"


def render_closeout_summary_markdown(summary: dict) -> str:
    readiness = summary["readiness"]
    lines = [
        "# Refactor Fidelity Closeout",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Target: `{summary['target']}`",
        f"- Pinned baseline: `{summary['baseline']}`",
        f"- Mainline: `{summary['mainline']}`",
        f"- Status: `{summary['status']}`",
        f"- Coverage-campaign ready: `{str(readiness['coverage_campaign_ready']).lower()}`",
        f"- Proven parity: `{str(summary['proven_parity']).lower()}`",
        f"- Open exceptions: `{readiness['open_exception_count']}`",
        f"- High-severity open exceptions: `{readiness['high_open_exception_count']}`",
        f"- Auto-curated exceptions: `{summary['auto_curated_count']}`",
    ]

    gate_order = (
        "pinned_baseline_explicit",
        "surface_gate",
        "manifest_gate",
        "behavior_gate",
        "contract_gate",
        "high_severity_gate",
    )
    gate_labels = {
        "pinned_baseline_explicit": "Pinned Baseline Explicit",
        "surface_gate": "Surface Gate",
        "manifest_gate": "Manifest Gate",
        "behavior_gate": "Behavior Gate",
        "contract_gate": "Contract Gate",
        "high_severity_gate": "High Severity Gate",
    }
    lines.extend(["", "## Gates", ""])
    for key in gate_order:
        lines.append(f"- `{gate_labels[key]}`: `{str(readiness[key]).lower()}`")

    blockers = ensure_list(summary.get("blockers"))
    if blockers:
        lines.extend(["", "## Blockers", ""])
        for blocker in blockers:
            lines.append(f"- {blocker}")

    lines.extend(["", "## Component Runs", ""])
    for name in sorted(summary["component_runs"]):
        lines.append(f"- `{name}`: `{summary['component_runs'][name]}`")

    lines.extend(["", "## Component Statuses", ""])
    for name in ("baseline_surface", "mainline_surface", "manifest", "behavior", "contracts"):
        component = summary["components"][name]
        lines.append(f"- `{name}`: `{component['status']}`")

    return "\n".join(lines) + "\n"


def exception_is_open(record: dict) -> bool:
    return record.get("curation_status") not in TERMINAL_CURATION_STATUSES


def render_exception_report_markdown(report: dict) -> str:
    summary = report["summary"]
    exceptions = report["exceptions"]
    lines = [
        "# Exceptions",
        "",
        f"- Repo root: `{report['repo_root']}`",
        f"- Total exceptions: `{summary['total_count']}`",
        f"- Open exceptions: `{summary['open_count']}`",
        f"- High-severity open exceptions: `{summary['high_open_count']}`",
        f"- Decision ready: `{str(summary['decision_ready']).lower()}`",
    ]

    for title, key in (
        ("Status Counts", "status_counts"),
        ("Severity Counts", "severity_counts"),
        ("Layer Counts", "layer_counts"),
        ("Type Counts", "type_counts"),
        ("Disposition Counts", "disposition_counts"),
    ):
        values = summary.get(key, {})
        if not values:
            continue
        lines.extend(["", f"## {title}", ""])
        for name in sorted(values):
            lines.append(f"- `{name}`: `{values[name]}`")

    open_high = [record for record in exceptions if record["severity"] == "high" and record["is_open"]]
    if open_high:
        lines.extend(["", "## High-Severity Open Exceptions", ""])
        for record in open_high[:25]:
            lines.append(
                f"- `{record['audit_layer']}` `{record['exception_type']}` `{record['entity_key']}` "
                f"[{record['curation_status']}] ({record['baseline_ref']} -> {record['target_ref']}): {record['reason']}"
            )

    unresolved = [record for record in exceptions if record["is_open"]]
    if unresolved:
        lines.extend(["", "## Open Exceptions", ""])
        for record in unresolved[:50]:
            suggestion = record.get("suggested_disposition")
            suffix = f" Suggested: `{suggestion}`." if suggestion else ""
            lines.append(
                f"- [{record['severity']}] `{record['audit_layer']}` `{record['exception_type']}` "
                f"`{record['entity_key']}` [{record['curation_status']}]"
                f"{suffix}"
            )
    elif exceptions:
        lines.extend(["", "## Open Exceptions", "", "- None."])
    else:
        lines.extend(["", "No exceptions recorded."])

    return "\n".join(lines) + "\n"


def render_exceptions_markdown(exceptions: list[dict]) -> str:
    normalized = []
    for record in exceptions:
        item = dict(record)
        item.setdefault("is_open", exception_is_open(item))
        item.setdefault("baseline_ref", item.get("baseline_ref"))
        item.setdefault("target_ref", item.get("target_ref"))
        item.setdefault("suggested_disposition", suggested_disposition_for_exception(item.get("exception_type", "")))
        normalized.append(item)
    report = build_exception_report(repo_root=".", records=normalized, filters={})
    return render_exception_report_markdown(report)


def load_json_payload(text: str) -> dict:
    stripped = text.strip()
    if not stripped:
        raise RuntimeError("inventory snapshot emitted no JSON")

    try:
        return json.loads(stripped)
    except json.JSONDecodeError:
        pass

    decoder = json.JSONDecoder()
    for index in range(len(stripped) - 1, -1, -1):
        if stripped[index] not in "[{":
            continue
        try:
            payload, end = decoder.raw_decode(stripped[index:])
        except json.JSONDecodeError:
            continue
        if stripped[index + end :].strip():
            continue
        return payload

    raise RuntimeError("inventory snapshot emitted non-JSON output")


def run_inventory_snapshot(repo_root: pathlib.Path) -> dict:
    script_path = SCRIPT_DIR / "fidelity_inventory_snapshot.R"
    command = ["Rscript", "--vanilla", str(script_path), "--repo-root", str(repo_root)]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "inventory snapshot failed"
        raise RuntimeError(message)
    return load_json_payload(completed.stdout)


def run_manifest_compare(manifest_path: pathlib.Path, baseline_root: pathlib.Path, target_root: pathlib.Path) -> dict:
    script_path = SCRIPT_DIR / "fidelity_manifest_compare.R"
    command = [
        "Rscript",
        "--vanilla",
        str(script_path),
        "--manifest",
        str(manifest_path),
        "--baseline-root",
        str(baseline_root),
        "--target-root",
        str(target_root),
    ]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "manifest comparison failed"
        raise RuntimeError(message)
    return load_json_payload(completed.stdout)


def relative_label(path: pathlib.Path, repo_root: pathlib.Path) -> str:
    try:
        return str(path.resolve().relative_to(repo_root.resolve()))
    except ValueError:
        return str(path.resolve())


def list_manifest_paths(repo_root: pathlib.Path, explicit_paths: list[str] | None) -> list[pathlib.Path]:
    if explicit_paths:
        manifests = []
        for value in explicit_paths:
            manifest_path = pathlib.Path(value)
            if not manifest_path.is_absolute():
                manifest_path = repo_root / manifest_path
            manifests.append(manifest_path.resolve())
        return manifests

    return sorted((repo_root / "tools" / "refactor").glob("manifest-*.yml"))


def list_behavior_case_paths(repo_root: pathlib.Path, explicit_paths: list[str] | None) -> list[pathlib.Path]:
    if explicit_paths:
        case_paths = []
        for value in explicit_paths:
            case_path = pathlib.Path(value)
            if not case_path.is_absolute():
                case_path = repo_root / case_path
            case_paths.append(case_path.resolve())
        return case_paths

    return sorted((repo_root / "tools" / "refactor").glob(DEFAULT_BEHAVIOR_CASE_GLOB))


def stable_json(value) -> str:
    return json.dumps(value, sort_keys=True)


def make_exception_key(audit_layer: str, exception_type: str, entity_key: str) -> str:
    return f"{audit_layer}::{exception_type}::{entity_key}"


def exception_auto_rule(exception_type: str) -> dict | None:
    return AUTO_CURATION_RULES.get(exception_type)


def suggested_disposition_for_exception(exception_type: str) -> str | None:
    rule = exception_auto_rule(exception_type)
    if rule:
        return str(rule["disposition"])

    suggestions = {
        "manual_merge_expected": "expected_manual_merge",
        "missing_target_block": "needs_manifest_lineage_review",
        "selector_ambiguous": "unsupported_audit_shape",
        "semantic_mismatch": "fix_required",
        "behavior_runner_failure": "unsupported_audit_shape",
        "behavior_status_mismatch": "fix_required",
        "behavior_error_mismatch": "fix_required",
        "behavior_output_mismatch": "fix_required",
        "behavior_warning_mismatch": "review_warning_contract",
        "behavior_message_mismatch": "review_message_contract",
        "behavior_stdout_mismatch": "review_stdout_contract",
        "surface_missing_export": "intentional_api_change",
        "surface_extra_export": "intentional_api_change",
        "surface_collate_drift": "intentional_collate_change",
        "surface_definition_drift": "fix_required",
        "surface_missing_definition": "fix_required",
        "surface_extra_definition": "fix_required",
        "surface_duplicate_definition": "fix_required",
        "surface_missing_file": "fix_required",
        "surface_extra_file": "fix_required",
        "contract_test_failure": "fix_required",
        "contract_runner_failure": "unsupported_audit_shape",
    }
    return suggestions.get(exception_type)


def load_behavior_cases(
    repo_root: pathlib.Path, explicit_paths: list[str] | None
) -> tuple[list[pathlib.Path], list[dict]]:
    case_paths = list_behavior_case_paths(repo_root, explicit_paths)
    if not case_paths:
        raise RuntimeError("No behavior case files found for behavior audit")

    cases: list[dict] = []
    seen_case_ids: dict[str, str] = {}
    for case_path in case_paths:
        payload = json.loads(case_path.read_text(encoding="utf-8"))
        raw_cases = ensure_list(payload.get("cases"))
        if not raw_cases:
            raise RuntimeError(f"Behavior case file contains no cases: {case_path}")
        for raw_case in raw_cases:
            if not isinstance(raw_case, dict):
                raise RuntimeError(f"Behavior case must be an object in {case_path}")
            case = dict(raw_case)
            case_id = case.get("id")
            if not case_id:
                raise RuntimeError(f"Behavior case is missing id in {case_path}")
            if case_id in seen_case_ids:
                raise RuntimeError(
                    f"Behavior case id {case_id} is duplicated in "
                    f"{seen_case_ids[case_id]} and {case_path}"
                )
            seen_case_ids[case_id] = str(case_path)
            case["family"] = case.get("family") or "pure_helper"
            case["entity_key"] = case.get("entity_key") or f"behavior::{case_id}"
            case["fixture_kind"] = case.get("fixture_kind") or "inline_setup"
            case["normalizer_key"] = case.get("normalizer_key") or (
                "inline_expr" if case.get("normalize_expr") else "identity"
            )
            case["risk_level"] = case.get("risk_level") or "medium"
            case["catalog_path"] = relative_label(case_path, repo_root)
            case["_catalog_file"] = case_path
            case["compare_warnings"] = bool(case.get("compare_warnings", True))
            case["compare_messages"] = bool(case.get("compare_messages", False))
            case["compare_output"] = bool(case.get("compare_output", False))
            cases.append(case)

    return case_paths, cases


def run_behavior_case(case_path: pathlib.Path, case_id: str, repo_root: pathlib.Path, side: str) -> dict:
    script_path = SCRIPT_DIR / "fidelity_behavior_replay.R"
    command = [
        "Rscript",
        "--vanilla",
        str(script_path),
        "--case-file",
        str(case_path),
        "--case-id",
        case_id,
        "--repo-root",
        str(repo_root),
        "--side",
        side,
    ]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "behavior replay failed"
        return {
            "case_id": case_id,
            "side": side,
            "status": "runner_error",
            "raw_result_json": None,
            "normalized_result_json": None,
            "warnings": [],
            "messages": [],
            "output": [],
            "error": {
                "message": message,
                "class": ["behavior_runner_error"],
            },
        }
    try:
        payload = load_json_payload(completed.stdout)
    except RuntimeError as exc:
        payload = {
            "case_id": case_id,
            "side": side,
            "status": "runner_error",
            "raw_result_json": None,
            "normalized_result_json": None,
            "warnings": [],
            "messages": [],
            "output": [],
            "error": {
                "message": str(exc),
                "class": ["behavior_runner_output_error"],
            },
        }
    payload["warnings"] = [str(value) for value in ensure_list(payload.get("warnings"))]
    payload["messages"] = [str(value) for value in ensure_list(payload.get("messages"))]
    payload["output"] = [str(value) for value in ensure_list(payload.get("output"))]
    payload["status"] = str(payload.get("status") or "runner_error")
    return payload


def compare_behavior_payloads(case: dict, baseline_payload: dict, target_payload: dict) -> dict:
    compare_warnings = bool(case.get("compare_warnings", True))
    compare_messages = bool(case.get("compare_messages", False))
    compare_output = bool(case.get("compare_output", False))

    baseline_status = baseline_payload.get("status") or "runner_error"
    target_status = target_payload.get("status") or "runner_error"
    artifact_equal = (
        baseline_status == "ok"
        and target_status == "ok"
        and baseline_payload.get("raw_result_json") == target_payload.get("raw_result_json")
    )
    normalized_equal = (
        baseline_status == "ok"
        and target_status == "ok"
        and baseline_payload.get("normalized_result_json") == target_payload.get("normalized_result_json")
    )
    warning_equal = (
        baseline_payload.get("warnings", []) == target_payload.get("warnings", [])
        if compare_warnings
        else True
    )
    message_equal = (
        baseline_payload.get("messages", []) == target_payload.get("messages", [])
        if compare_messages
        else True
    )
    output_equal = (
        baseline_payload.get("output", []) == target_payload.get("output", [])
        if compare_output
        else True
    )
    error_equal = stable_json(baseline_payload.get("error")) == stable_json(target_payload.get("error"))

    mismatch_reasons: list[str] = []
    if "runner_error" in {baseline_status, target_status}:
        mismatch_reasons.append("runner_failure")
        comparison_status = "runner_failure"
    elif baseline_status != target_status:
        mismatch_reasons.append("status_mismatch")
        comparison_status = "mismatch"
    elif baseline_status == "error":
        if not error_equal:
            mismatch_reasons.append("error_payload_mismatch")
        if compare_warnings and not warning_equal:
            mismatch_reasons.append("warning_mismatch")
        if compare_messages and not message_equal:
            mismatch_reasons.append("message_mismatch")
        if compare_output and not output_equal:
            mismatch_reasons.append("stdout_mismatch")
        comparison_status = "error_match" if not mismatch_reasons else "mismatch"
    else:
        if not normalized_equal:
            mismatch_reasons.append("normalized_result_mismatch")
        if compare_warnings and not warning_equal:
            mismatch_reasons.append("warning_mismatch")
        if compare_messages and not message_equal:
            mismatch_reasons.append("message_mismatch")
        if compare_output and not output_equal:
            mismatch_reasons.append("stdout_mismatch")

        if mismatch_reasons:
            comparison_status = "mismatch"
        elif artifact_equal:
            comparison_status = "exact_match"
        else:
            comparison_status = "normalized_match"

    parity_tier = {
        "exact_match": "T4_behavior_replayed_exact",
        "normalized_match": "T4_behavior_replayed_normalized",
        "error_match": "T4_behavior_replayed_error",
        "mismatch": "T3_behavior_mismatch",
        "runner_failure": "T3_behavior_runner_failure",
    }[comparison_status]

    return {
        "behavior_case_id": case["id"],
        "case_id": case["id"],
        "catalog_path": case["catalog_path"],
        "family": case["family"],
        "entity_key": case["entity_key"],
        "risk_level": case["risk_level"],
        "comparison_status": comparison_status,
        "baseline_status": baseline_status,
        "target_status": target_status,
        "normalized_equal": normalized_equal,
        "warning_equal": warning_equal,
        "message_equal": message_equal,
        "error_equal": error_equal,
        "artifact_equal": artifact_equal,
        "result_payload": {
            "case_id": case["id"],
            "catalog_path": case["catalog_path"],
            "family": case["family"],
            "entity_key": case["entity_key"],
            "risk_level": case["risk_level"],
            "compare_flags": {
                "warnings": compare_warnings,
                "messages": compare_messages,
                "output": compare_output,
            },
            "comparison_status": comparison_status,
            "parity_tier": parity_tier,
            "artifact_equal": artifact_equal,
            "normalized_equal": normalized_equal,
            "warning_equal": warning_equal,
            "message_equal": message_equal,
            "error_equal": error_equal,
            "output_equal": output_equal,
            "mismatch_reasons": mismatch_reasons,
            "baseline": baseline_payload,
            "target": target_payload,
        },
    }


def behavior_exception_type(result: dict) -> str:
    reasons = result["result_payload"]["mismatch_reasons"]
    for reason_name, exception_type in (
        ("runner_failure", "behavior_runner_failure"),
        ("status_mismatch", "behavior_status_mismatch"),
        ("error_payload_mismatch", "behavior_error_mismatch"),
        ("normalized_result_mismatch", "behavior_output_mismatch"),
        ("warning_mismatch", "behavior_warning_mismatch"),
        ("message_mismatch", "behavior_message_mismatch"),
        ("stdout_mismatch", "behavior_stdout_mismatch"),
    ):
        if reason_name in reasons:
            return exception_type
    return "behavior_review_required"


def behavior_exception_severity(exception_type: str, risk_level: str) -> str:
    severity = {
        "behavior_runner_failure": "high",
        "behavior_status_mismatch": "high",
        "behavior_error_mismatch": "high",
        "behavior_output_mismatch": "high",
        "behavior_warning_mismatch": "medium",
        "behavior_message_mismatch": "low",
        "behavior_stdout_mismatch": "low",
        "behavior_review_required": "medium",
    }[exception_type]
    risk_rank = {"low": 0, "medium": 1, "high": 2, "critical": 3}
    severity_rank = {"low": 0, "medium": 1, "high": 2}
    if risk_rank.get(risk_level, 1) >= 2 and severity_rank[severity] < 2:
        return "high"
    return severity


def behavior_exception_reason(result: dict, exception_type: str) -> str:
    baseline_status = result["baseline_status"]
    target_status = result["target_status"]
    if exception_type == "behavior_runner_failure":
        return "Behavior replay runner failed before a comparable result was produced."
    if exception_type == "behavior_status_mismatch":
        return f"Baseline replay ended in `{baseline_status}` while target replay ended in `{target_status}`."
    if exception_type == "behavior_error_mismatch":
        return "Both sides errored, but the captured error payloads differ."
    if exception_type == "behavior_output_mismatch":
        return "Normalized behavior outputs differ after applying the case normalizer."
    if exception_type == "behavior_warning_mismatch":
        return "Behavior replay warnings differ between baseline and target."
    if exception_type == "behavior_message_mismatch":
        return "Behavior replay messages differ between baseline and target."
    if exception_type == "behavior_stdout_mismatch":
        return "Behavior replay stdout differs between baseline and target."
    return "Behavior replay requires manual review."


def build_behavior_exception(
    run_id: str,
    result: dict,
    behavior_result_id: str,
    *,
    baseline_ref: str,
    target_ref: str,
) -> dict | None:
    if result["comparison_status"] in {"exact_match", "normalized_match", "error_match"}:
        return None

    exception_type = behavior_exception_type(result)
    entity_key = result["entity_key"]
    return {
        "exception_id": f"{run_id}:behavior_exception:{result['case_id']}",
        "exception_key": make_exception_key("behavior_replay", exception_type, entity_key),
        "run_id": run_id,
        "audit_layer": "behavior_replay",
        "entity_key": entity_key,
        "severity": behavior_exception_severity(exception_type, result["risk_level"]),
        "exception_type": exception_type,
        "reason": behavior_exception_reason(result, exception_type),
        "baseline_ref": baseline_ref,
        "target_ref": target_ref,
        "evidence_json": stable_json(
            {
                "behavior_result_id": behavior_result_id,
                "case_id": result["case_id"],
                "catalog_path": result["catalog_path"],
                "comparison_status": result["comparison_status"],
                "mismatch_reasons": result["result_payload"]["mismatch_reasons"],
                "baseline_status": result["baseline_status"],
                "target_status": result["target_status"],
                "result": result["result_payload"],
            }
        ),
        "curation_status": "candidate",
        "disposition": None,
        "owner": None,
        "review_note": None,
        "updated_at": now_iso(),
        "curation_rule": None,
        "resolved_at": None,
    }


def list_test_file_paths(target_root: pathlib.Path, explicit_paths: list[str] | None) -> list[pathlib.Path]:
    if explicit_paths:
        resolved_paths = []
        for value in explicit_paths:
            test_path = pathlib.Path(value)
            if not test_path.is_absolute():
                test_path = target_root / test_path
            resolved_paths.append(test_path.resolve())
        return sorted(resolved_paths)

    test_root = target_root / "tests" / "testthat"
    return sorted(test_root.glob(DEFAULT_TEST_FILE_GLOB))


def classify_test_file(file_path: pathlib.Path) -> dict:
    stem = file_path.stem
    suffix_token = stem.removeprefix("test-")
    spec = DEFAULT_TEST_FAMILY
    matched_suffix = ""
    for candidate in TEST_FAMILY_SPECS:
        if suffix_token.endswith(candidate["suffix"]):
            spec = candidate
            matched_suffix = candidate["suffix"]
            break

    module_token = suffix_token
    if matched_suffix:
        module_token = module_token[: -len(matched_suffix)].rstrip("-")
    module_token = re.sub(r"^(?:prot|pept|metab|lipid)-\d+[a-z]?[-_]?", "", module_token)
    module_token = module_token or suffix_token
    module_family = re.sub(r"[^a-z0-9]+", "_", module_token.lower()).strip("_") or "general"

    return {
        "family": spec["family"],
        "surface": spec["surface"],
        "certainty_target": spec["certainty_target"],
        "module_family": module_family,
    }


TEST_NAME_PATTERN = re.compile(
    r"test_that\(\s*(?P<literal>\"(?:[^\"\\]|\\.)*\"|'(?:[^'\\]|\\.)*')",
    re.MULTILINE | re.DOTALL,
)

ENTITY_HINT_PATTERN = re.compile(r"\b(?:mod_[a-z0-9_]+|[A-Za-z][A-Za-z0-9_.]*[A-Z_][A-Za-z0-9_.]*)\b")


def extract_test_cases(file_path: pathlib.Path) -> list[dict]:
    text = file_path.read_text(encoding="utf-8")
    cases: list[dict] = []
    for match in TEST_NAME_PATTERN.finditer(text):
        test_name = ast.literal_eval(match.group("literal"))
        lowered = test_name.lower()
        scenario_tags: list[str] = []
        for scenario_name, patterns in SCENARIO_PATTERNS:
            if any(re.search(pattern, lowered) for pattern in patterns):
                scenario_tags.append(scenario_name)
        if not scenario_tags:
            scenario_tags = ["happy_path"]

        primary_scenario = scenario_tags[0]
        entity_hint_match = ENTITY_HINT_PATTERN.search(test_name)
        entity_hint = entity_hint_match.group(0) if entity_hint_match else None
        cases.append(
            {
                "test_name": test_name,
                "entity_hint": entity_hint,
                "primary_scenario": primary_scenario,
                "scenario_tags": scenario_tags,
            }
        )
    return cases


def load_testing_matrix_entries() -> list[dict]:
    plan_path = SCRIPT_DIR / "fidelity-audit-plan.json"
    if not plan_path.exists():
        return []
    payload = json.loads(plan_path.read_text(encoding="utf-8"))
    return [entry for entry in ensure_list(payload.get("testing_matrix")) if isinstance(entry, dict)]


def run_contract_test_file(
    target_root: pathlib.Path,
    file_path: pathlib.Path,
    *,
    load_package: bool,
) -> dict:
    script_path = SCRIPT_DIR / "fidelity_contract_runner.R"
    command = [
        "Rscript",
        "--vanilla",
        str(script_path),
        "--repo-root",
        str(target_root),
        "--test-file",
        str(file_path),
    ]
    if load_package:
        command.append("--load-package")
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "contract runner failed"
        return {
            "file_path": str(file_path),
            "status": "runner_error",
            "test_count": 0,
            "failure_count": 0,
            "skip_count": 0,
            "load_package": load_package,
            "tests": [],
            "error": {
                "message": message,
                "class": ["contract_runner_error"],
            },
        }
    try:
        payload = load_json_payload(completed.stdout)
    except RuntimeError as exc:
        payload = {
            "file_path": str(file_path),
            "status": "runner_error",
            "test_count": 0,
            "failure_count": 0,
            "skip_count": 0,
            "load_package": load_package,
            "tests": [],
            "error": {
                "message": str(exc),
                "class": ["contract_runner_output_error"],
            },
        }
    payload["tests"] = [item for item in ensure_list(payload.get("tests")) if isinstance(item, dict)]
    return payload


def catalog_contract_files(target_root: pathlib.Path, explicit_paths: list[str] | None) -> list[dict]:
    test_paths = list_test_file_paths(target_root, explicit_paths)
    if not test_paths:
        raise RuntimeError("No testthat files found for contract audit")

    contract_files: list[dict] = []
    for test_path in test_paths:
        classification = classify_test_file(test_path)
        cases = extract_test_cases(test_path)
        file_record = {
            "file_path": relative_label(test_path, target_root),
            "absolute_path": test_path,
            "family": classification["family"],
            "module_family": classification["module_family"],
            "surface": classification["surface"],
            "certainty_target": classification["certainty_target"],
            "required_scenarios": list(MODULE_CONTRACT_SCENARIOS)
            if classification["family"] == "module_contract"
            else [],
            "cases": cases,
        }
        contract_files.append(file_record)

    return contract_files


def run_git(repo_root: pathlib.Path, args: list[str], *, text: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", "-C", str(repo_root), *args],
        capture_output=True,
        text=text,
        check=False,
    )


def resolve_git_ref_label(repo_root: pathlib.Path, ref: str) -> str:
    completed = run_git(repo_root, ["rev-parse", ref])
    if completed.returncode != 0:
        message = completed.stderr.strip() or f"Unable to resolve git ref {ref}"
        raise RuntimeError(message)
    sha = completed.stdout.strip()
    return f"{ref}@{sha}"


def materialize_git_ref(repo_root: pathlib.Path, ref: str, archive_paths: list[str] | None = None) -> dict:
    git_args = ["archive", ref]
    if archive_paths:
        git_args.extend(archive_paths)
    archive = run_git(repo_root, git_args, text=False)
    if archive.returncode != 0:
        message = archive.stderr.decode("utf-8", errors="replace").strip() or f"Unable to archive git ref {ref}"
        raise RuntimeError(message)

    tmpdir = tempfile.TemporaryDirectory(prefix="fidelity-ref-")
    try:
        with tarfile.open(fileobj=io.BytesIO(archive.stdout), mode="r:*") as tar_handle:
            tar_handle.extractall(tmpdir.name)
    except Exception:
        tmpdir.cleanup()
        raise

    return {
        "path": pathlib.Path(tmpdir.name),
        "label": resolve_git_ref_label(repo_root, ref),
        "cleanup": tmpdir.cleanup,
    }


def resolve_side_source(
    repo_root: pathlib.Path,
    *,
    side: str,
    ref: str | None,
    path: str | None,
    default_ref: str,
    allow_worktree: bool,
    archive_paths: list[str] | None = None,
) -> dict:
    if path and ref:
        raise RuntimeError(f"{side} path and ref are mutually exclusive")

    if path:
        resolved = pathlib.Path(path).resolve()
        return {
            "path": resolved,
            "label": str(resolved),
            "cleanup": lambda: None,
        }

    chosen_ref = ref or default_ref
    if allow_worktree and chosen_ref == "WORKTREE":
        return {
            "path": repo_root,
            "label": "WORKTREE",
            "cleanup": lambda: None,
        }

    return materialize_git_ref(repo_root, chosen_ref, archive_paths=archive_paths)


def normalize_snapshot_exports(snapshot: dict) -> dict[str, list[str]]:
    exports = snapshot.get("exports", {})
    return {
        "functions": [str(value) for value in ensure_list(exports.get("functions"))],
        "methods": [str(value) for value in ensure_list(exports.get("methods"))],
        "classes": [str(value) for value in ensure_list(exports.get("classes"))],
    }


def insert_audit_run(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    repo_root: pathlib.Path,
    baseline_ref: str | None,
    target_ref: str | None,
    mode: str,
    status: str,
    summary: dict,
) -> None:
    timestamp = now_iso()
    conn.execute(
        """
        INSERT INTO audit_runs (
            run_id, started_at, completed_at, repo_root, baseline_ref, target_ref, mode, status, summary_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            run_id,
            timestamp,
            timestamp if status != "running" else None,
            str(repo_root),
            baseline_ref,
            target_ref,
            mode,
            status,
            json.dumps(summary, sort_keys=True),
        ),
    )


def update_audit_run(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    status: str,
    summary: dict,
) -> None:
    conn.execute(
        """
        UPDATE audit_runs
        SET completed_at = ?, status = ?, summary_json = ?
        WHERE run_id = ?
        """,
        (now_iso(), status, json.dumps(summary, sort_keys=True), run_id),
    )


def insert_inventory_snapshot(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    side: str,
    snapshot: dict,
) -> dict:
    refs = {
        "entities": defaultdict(list),
        "files": {},
        "exports": {},
    }

    for index, entity in enumerate(snapshot["entities"], start=1):
        entity_id = f"{run_id}:{side}:entity:{index}"
        conn.execute(
            """
            INSERT INTO inventory_entities (
                entity_id, run_id, side, entity_key, entity_kind, entity_name, signature_key,
                file_path, line_start, line_end, exported, collate_index, hash_raw,
                hash_normalized, hash_ast
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                entity_id,
                run_id,
                side,
                entity["entity_key"],
                entity["entity_kind"],
                entity["entity_name"],
                entity.get("signature_key"),
                entity["file_path"],
                entity.get("line_start"),
                entity.get("line_end"),
                1 if entity.get("exported") else 0,
                entity.get("collate_index"),
                entity.get("hash_raw"),
                entity.get("hash_normalized"),
                entity.get("hash_ast"),
            ),
        )
        refs["entities"][entity["entity_key"]].append(entity_id)

    collate_lookup = {
        entry["file_path"]: entry.get("collate_index")
        for entry in ensure_list(snapshot.get("collate"))
    }
    for index, file_path in enumerate(ensure_list(snapshot.get("r_files")), start=1):
        file_id = f"{run_id}:{side}:file:{index}"
        conn.execute(
            """
            INSERT INTO inventory_files (
                file_id, run_id, side, file_path, collate_index, present
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                file_id,
                run_id,
                side,
                file_path,
                collate_lookup.get(file_path),
                1,
            ),
        )
        refs["files"][file_path] = file_id

    exports = normalize_snapshot_exports(snapshot)
    for export_kind, names in exports.items():
        for index, export_name in enumerate(names, start=1):
            export_id = f"{run_id}:{side}:export:{export_kind}:{index}"
            conn.execute(
                """
                INSERT INTO inventory_exports (
                    export_id, run_id, side, export_kind, export_name
                ) VALUES (?, ?, ?, ?, ?)
                """,
                (
                    export_id,
                    run_id,
                    side,
                    export_kind,
                    export_name,
                ),
            )
            refs["exports"][f"{export_kind}::{export_name}"] = export_id

    return refs


def make_diff(
    *,
    entity_kind: str,
    entity_key: str,
    diff_class: str,
    baseline_id: str | None,
    target_id: str | None,
    evidence: dict,
) -> dict:
    return {
        "entity_kind": entity_kind,
        "entity_key": entity_key,
        "diff_class": diff_class,
        "baseline_entity_id": baseline_id,
        "target_entity_id": target_id,
        "severity": DIFF_SEVERITY[diff_class],
        "status": "open" if diff_class in SURFACE_OPEN_DIFFS else "observed",
        "evidence": evidence,
    }


def entity_lookup(snapshot: dict) -> dict[str, list[dict]]:
    lookup: dict[str, list[dict]] = defaultdict(list)
    for entity in snapshot["entities"]:
        lookup[entity["entity_key"]].append(entity)
    return lookup


def file_lookup(snapshot: dict) -> dict[str, int | None]:
    return {
        entry["file_path"]: entry.get("collate_index")
        for entry in ensure_list(snapshot.get("collate"))
    }


def compare_entity_surfaces(
    baseline_snapshot: dict,
    target_snapshot: dict,
    baseline_refs: dict,
    target_refs: dict,
) -> list[dict]:
    diffs: list[dict] = []
    baseline_entities = entity_lookup(baseline_snapshot)
    target_entities = entity_lookup(target_snapshot)

    for entity_key in sorted(set(baseline_entities) | set(target_entities)):
        baseline_entries = baseline_entities.get(entity_key, [])
        target_entries = target_entities.get(entity_key, [])
        sample = (baseline_entries or target_entries)[0]
        entity_kind = sample["entity_kind"]
        baseline_ids = baseline_refs["entities"].get(entity_key, [])
        target_ids = target_refs["entities"].get(entity_key, [])

        if len(baseline_entries) > 1 or len(target_entries) > 1:
            diffs.append(
                make_diff(
                    entity_kind=entity_kind,
                    entity_key=entity_key,
                    diff_class="duplicate_definition",
                    baseline_id=baseline_ids[0] if baseline_ids else None,
                    target_id=target_ids[0] if target_ids else None,
                    evidence={
                        "baseline_count": len(baseline_entries),
                        "target_count": len(target_entries),
                        "baseline_files": sorted({entry["file_path"] for entry in baseline_entries}),
                        "target_files": sorted({entry["file_path"] for entry in target_entries}),
                    },
                )
            )

        if not baseline_entries:
            diffs.append(
                make_diff(
                    entity_kind=entity_kind,
                    entity_key=entity_key,
                    diff_class="extra_definition",
                    baseline_id=None,
                    target_id=target_ids[0] if target_ids else None,
                    evidence={
                        "target_count": len(target_entries),
                        "target_files": sorted({entry["file_path"] for entry in target_entries}),
                    },
                )
            )
            continue

        if not target_entries:
            diffs.append(
                make_diff(
                    entity_kind=entity_kind,
                    entity_key=entity_key,
                    diff_class="missing_definition",
                    baseline_id=baseline_ids[0] if baseline_ids else None,
                    target_id=None,
                    evidence={
                        "baseline_count": len(baseline_entries),
                        "baseline_files": sorted({entry["file_path"] for entry in baseline_entries}),
                    },
                )
            )
            continue

        if len(baseline_entries) == 1 and len(target_entries) == 1:
            baseline_entry = baseline_entries[0]
            target_entry = target_entries[0]

            if baseline_entry["file_path"] != target_entry["file_path"]:
                diffs.append(
                    make_diff(
                        entity_kind=entity_kind,
                        entity_key=entity_key,
                        diff_class="moved_definition",
                        baseline_id=baseline_ids[0] if baseline_ids else None,
                        target_id=target_ids[0] if target_ids else None,
                        evidence={
                            "baseline_file": baseline_entry["file_path"],
                            "target_file": target_entry["file_path"],
                        },
                    )
                )

            if baseline_entry.get("hash_ast") != target_entry.get("hash_ast"):
                diffs.append(
                    make_diff(
                        entity_kind=entity_kind,
                        entity_key=entity_key,
                        diff_class="definition_drift",
                        baseline_id=baseline_ids[0] if baseline_ids else None,
                        target_id=target_ids[0] if target_ids else None,
                        evidence={
                            "baseline_file": baseline_entry["file_path"],
                            "target_file": target_entry["file_path"],
                            "baseline_hash_ast": baseline_entry.get("hash_ast"),
                            "target_hash_ast": target_entry.get("hash_ast"),
                        },
                    )
                )

    return diffs


def compare_file_surfaces(
    baseline_snapshot: dict,
    target_snapshot: dict,
    baseline_refs: dict,
    target_refs: dict,
) -> list[dict]:
    diffs: list[dict] = []
    baseline_files = file_lookup(baseline_snapshot)
    target_files = file_lookup(target_snapshot)

    for file_path in sorted(set(baseline_files) | set(target_files)):
        baseline_id = baseline_refs["files"].get(file_path)
        target_id = target_refs["files"].get(file_path)
        if file_path not in baseline_files:
            diffs.append(
                make_diff(
                    entity_kind="r_file",
                    entity_key=file_path,
                    diff_class="extra_file",
                    baseline_id=None,
                    target_id=target_id,
                    evidence={"target_collate_index": target_files[file_path]},
                )
            )
            continue
        if file_path not in target_files:
            diffs.append(
                make_diff(
                    entity_kind="r_file",
                    entity_key=file_path,
                    diff_class="missing_file",
                    baseline_id=baseline_id,
                    target_id=None,
                    evidence={"baseline_collate_index": baseline_files[file_path]},
                )
            )
            continue
        if baseline_files[file_path] != target_files[file_path]:
            diffs.append(
                make_diff(
                    entity_kind="r_file",
                    entity_key=file_path,
                    diff_class="collate_drift",
                    baseline_id=baseline_id,
                    target_id=target_id,
                    evidence={
                        "baseline_collate_index": baseline_files[file_path],
                        "target_collate_index": target_files[file_path],
                    },
                )
            )

    return diffs


def compare_export_surfaces(
    baseline_snapshot: dict,
    target_snapshot: dict,
    baseline_refs: dict,
    target_refs: dict,
) -> list[dict]:
    diffs: list[dict] = []
    baseline_exports = normalize_snapshot_exports(baseline_snapshot)
    target_exports = normalize_snapshot_exports(target_snapshot)

    for export_kind in ("functions", "methods", "classes"):
        baseline_set = set(baseline_exports[export_kind])
        target_set = set(target_exports[export_kind])
        entity_kind = f"export_{export_kind[:-1]}"

        for export_name in sorted(baseline_set - target_set):
            export_key = f"{export_kind}::{export_name}"
            diffs.append(
                make_diff(
                    entity_kind=entity_kind,
                    entity_key=export_key,
                    diff_class="missing_export",
                    baseline_id=baseline_refs["exports"].get(export_key),
                    target_id=None,
                    evidence={"export_kind": export_kind, "export_name": export_name},
                )
            )

        for export_name in sorted(target_set - baseline_set):
            export_key = f"{export_kind}::{export_name}"
            diffs.append(
                make_diff(
                    entity_kind=entity_kind,
                    entity_key=export_key,
                    diff_class="extra_export",
                    baseline_id=None,
                    target_id=target_refs["exports"].get(export_key),
                    evidence={"export_kind": export_kind, "export_name": export_name},
                )
            )

    return diffs


def sort_diffs(diffs: list[dict]) -> list[dict]:
    severity_order = {"high": 0, "medium": 1, "low": 2, "info": 3}
    return sorted(
        diffs,
        key=lambda diff: (
            severity_order.get(diff["severity"], 99),
            diff["diff_class"],
            diff["entity_key"],
        ),
    )


def insert_inventory_diffs(conn: sqlite3.Connection, *, run_id: str, diffs: list[dict]) -> None:
    for index, diff in enumerate(diffs, start=1):
        conn.execute(
            """
            INSERT INTO inventory_diffs (
                diff_id, run_id, entity_kind, entity_key, diff_class, baseline_entity_id,
                target_entity_id, severity, status, evidence_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                f"{run_id}:diff:{index}",
                run_id,
                diff["entity_kind"],
                diff["entity_key"],
                diff["diff_class"],
                diff.get("baseline_entity_id"),
                diff.get("target_entity_id"),
                diff["severity"],
                diff["status"],
                json.dumps(diff.get("evidence", {}), sort_keys=True),
            ),
        )


def build_surface_exception(run_id: str, diff: dict, diff_id: str, *, baseline_ref: str, target_ref: str) -> dict | None:
    if diff["status"] != "open":
        return None

    exception_type = f"surface_{diff['diff_class']}"
    entity_key = diff["entity_key"]
    return {
        "exception_id": f"{run_id}:surface_exception:{diff['diff_class']}:{entity_key}",
        "exception_key": make_exception_key("surface_inventory", exception_type, entity_key),
        "run_id": run_id,
        "audit_layer": "surface_inventory",
        "entity_key": entity_key,
        "severity": diff["severity"],
        "exception_type": exception_type,
        "reason": f"Surface audit reported `{diff['diff_class']}` for `{entity_key}`.",
        "baseline_ref": baseline_ref,
        "target_ref": target_ref,
        "evidence_json": stable_json(
            {
                "diff_id": diff_id,
                "diff_class": diff["diff_class"],
                "entity_kind": diff["entity_kind"],
                "evidence": diff.get("evidence", {}),
            }
        ),
        "curation_status": "candidate",
        "disposition": None,
        "owner": None,
        "review_note": None,
        "updated_at": now_iso(),
        "curation_rule": None,
        "resolved_at": None,
    }


def insert_surface_exceptions(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    diffs: list[dict],
    baseline_ref: str,
    target_ref: str,
) -> list[dict]:
    exceptions: list[dict] = []
    for index, diff in enumerate(diffs, start=1):
        record = build_surface_exception(
            run_id,
            diff,
            f"{run_id}:diff:{index}",
            baseline_ref=baseline_ref,
            target_ref=target_ref,
        )
        if record is not None:
            exceptions.append(record)
    insert_exceptions(conn, exceptions)
    return exceptions


def normalize_manifest_comparisons(payload: dict) -> list[dict]:
    comparisons = ensure_list(payload.get("comparisons"))
    normalized: list[dict] = []
    for comparison in comparisons:
        record = dict(comparison)
        selector_payload = record.get("selector_payload")
        if selector_payload is None:
            selector_payload = {}
        record["selector_payload"] = selector_payload
        normalized.append(record)
    return normalized


def manifest_exception_severity(exception_type: str) -> str:
    return {
        "missing_target_block": "high",
        "selector_ambiguous": "high",
        "semantic_mismatch": "high",
        "manual_merge_expected": "medium",
        "ast_only_drift": "medium",
        "target_resolution_asymmetry": "low",
        "whitespace_only_drift": "low",
        "comment_only_drift": "low",
    }.get(exception_type, "medium")


def manifest_exception_reason(comparison: dict) -> str:
    status = comparison["status"]
    target_resolver = comparison.get("target_resolver") or "unresolved"
    note = comparison.get("note")
    if status == "manual_merge_expected":
        return "Manifest entry is marked manual_merge and requires curated reconciliation."
    if status == "content_missing":
        return f"Target block could not be resolved via `{target_resolver}`."
    if status == "selector_ambiguous":
        return "Target selector resolution was ambiguous."
    if status == "normalized_only":
        return f"Target block matched only after whitespace normalization via `{target_resolver}`."
    if status == "ast_only":
        return f"Target block matched only at AST level via `{target_resolver}`."
    if status == "semantic_drift_suspected":
        return "Target block resolved, but raw, normalized-text, and AST hashes all differ."
    if status == "raw_exact" and target_resolver != "selector":
        return f"Target block was exact, but only after `{target_resolver}` fallback resolution."
    if note:
        return note
    return f"Manifest comparison ended in `{status}`."


def build_manifest_exception(
    run_id: str,
    comparison: dict,
    *,
    baseline_ref: str,
    target_ref: str,
) -> dict | None:
    if comparison["status"] == "raw_exact" and (comparison.get("target_resolver") or "selector") == "selector":
        return None

    if comparison["status"] == "raw_exact":
        exception_type = "target_resolution_asymmetry"
    else:
        exception_type = comparison.get("exception_candidate") or "manifest_review_required"
    entity_key = f"manifest::{comparison['manifest_path']}::{comparison['entry_id']}"
    return {
        "exception_id": f"{run_id}:exception:{comparison['manifest_path']}:{comparison['entry_id']}",
        "exception_key": make_exception_key("manifest_fidelity", exception_type, entity_key),
        "run_id": run_id,
        "audit_layer": "manifest_fidelity",
        "entity_key": entity_key,
        "severity": manifest_exception_severity(exception_type),
        "exception_type": exception_type,
        "reason": manifest_exception_reason(comparison),
        "baseline_ref": baseline_ref,
        "target_ref": target_ref,
        "evidence_json": json.dumps(
            {
                "status": comparison["status"],
                "target_resolver": comparison.get("target_resolver"),
                "source_path": comparison.get("source_path"),
                "target_path": comparison.get("target_path"),
                "note": comparison.get("note"),
                "hash_raw_baseline": comparison.get("hash_raw_baseline"),
                "hash_raw_target": comparison.get("hash_raw_target"),
                "hash_normalized_baseline": comparison.get("hash_normalized_baseline"),
                "hash_normalized_target": comparison.get("hash_normalized_target"),
                "hash_ast_baseline": comparison.get("hash_ast_baseline"),
                "hash_ast_target": comparison.get("hash_ast_target"),
            },
            sort_keys=True,
        ),
        "curation_status": "candidate",
        "disposition": None,
        "owner": None,
        "review_note": None,
        "updated_at": now_iso(),
        "curation_rule": None,
        "resolved_at": None,
    }


def insert_exceptions(conn: sqlite3.Connection, exceptions: list[dict]) -> None:
    for record in exceptions:
        conn.execute(
            """
            INSERT INTO exceptions (
                exception_id, exception_key, run_id, audit_layer, entity_key, severity, exception_type,
                reason, baseline_ref, target_ref, evidence_json, curation_status, disposition,
                owner, review_note, updated_at, curation_rule, resolved_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                record["exception_id"],
                record.get("exception_key"),
                record["run_id"],
                record["audit_layer"],
                record["entity_key"],
                record["severity"],
                record["exception_type"],
                record["reason"],
                record.get("baseline_ref"),
                record.get("target_ref"),
                record["evidence_json"],
                record["curation_status"],
                record["disposition"],
                record["owner"],
                record.get("review_note"),
                record.get("updated_at"),
                record.get("curation_rule"),
                record["resolved_at"],
            ),
        )


def insert_manifest_results(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    manifest_comparisons: list[dict],
    baseline_ref: str,
    target_ref: str,
) -> list[dict]:
    exceptions: list[dict] = []
    for index, comparison in enumerate(manifest_comparisons, start=1):
        manifest_entry_id = f"{run_id}:manifest_entry:{index}"
        comparison_id = f"{run_id}:manifest_comparison:{index}"
        conn.execute(
            """
            INSERT INTO manifest_entries (
                manifest_entry_id, run_id, manifest_path, entry_id, selector_kind, selector_payload_json,
                source_path, target_path, group_name, action
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                manifest_entry_id,
                run_id,
                comparison["manifest_path"],
                comparison["entry_id"],
                comparison["selector_kind"],
                json.dumps(comparison["selector_payload"], sort_keys=True),
                comparison["source_path"],
                comparison.get("target_path"),
                comparison.get("group_name"),
                comparison["action"],
            ),
        )
        conn.execute(
            """
            INSERT INTO manifest_comparisons (
                comparison_id, run_id, manifest_entry_id, baseline_source_resolver, target_resolver,
                raw_exact, normalized_text_exact, ast_exact, hash_raw_baseline, hash_raw_target,
                hash_normalized_baseline, hash_normalized_target, hash_ast_baseline, hash_ast_target,
                status, evidence_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                comparison_id,
                run_id,
                manifest_entry_id,
                comparison.get("baseline_source_resolver"),
                comparison.get("target_resolver"),
                1 if comparison.get("raw_exact") else 0,
                1 if comparison.get("normalized_text_exact") else 0,
                1 if comparison.get("ast_exact") else 0,
                comparison.get("hash_raw_baseline"),
                comparison.get("hash_raw_target"),
                comparison.get("hash_normalized_baseline"),
                comparison.get("hash_normalized_target"),
                comparison.get("hash_ast_baseline"),
                comparison.get("hash_ast_target"),
                comparison["status"],
                json.dumps(
                    {
                        "note": comparison.get("note"),
                        "source_start_line": comparison.get("source_start_line"),
                        "source_end_line": comparison.get("source_end_line"),
                        "target_start_line": comparison.get("target_start_line"),
                        "target_end_line": comparison.get("target_end_line"),
                    },
                    sort_keys=True,
                ),
            ),
        )

        exception_record = build_manifest_exception(
            run_id,
            comparison,
            baseline_ref=baseline_ref,
            target_ref=target_ref,
        )
        if exception_record is not None:
            exceptions.append(exception_record)

    insert_exceptions(conn, exceptions)
    return exceptions


def upsert_behavior_cases(conn: sqlite3.Connection, cases: list[dict]) -> None:
    for case in cases:
        case_payload = {key: value for key, value in case.items() if not key.startswith("_")}
        conn.execute(
            """
            INSERT INTO behavior_cases (
                behavior_case_id, family, entity_key, fixture_kind, fixture_path,
                normalizer_key, risk_level, catalog_path, case_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(behavior_case_id) DO UPDATE SET
                family = excluded.family,
                entity_key = excluded.entity_key,
                fixture_kind = excluded.fixture_kind,
                fixture_path = excluded.fixture_path,
                normalizer_key = excluded.normalizer_key,
                risk_level = excluded.risk_level,
                catalog_path = excluded.catalog_path,
                case_json = excluded.case_json
            """,
            (
                case["id"],
                case["family"],
                case["entity_key"],
                case["fixture_kind"],
                case.get("fixture_path"),
                case["normalizer_key"],
                case["risk_level"],
                case["catalog_path"],
                stable_json(case_payload),
            ),
        )


def insert_behavior_results(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    behavior_results: list[dict],
    baseline_ref: str,
    target_ref: str,
) -> list[dict]:
    exceptions: list[dict] = []
    for index, result in enumerate(behavior_results, start=1):
        behavior_result_id = f"{run_id}:behavior_result:{index}"
        conn.execute(
            """
            INSERT INTO behavior_results (
                behavior_result_id, run_id, behavior_case_id, baseline_status, target_status,
                normalized_equal, warning_equal, message_equal, error_equal, artifact_equal, result_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                behavior_result_id,
                run_id,
                result["behavior_case_id"],
                result["baseline_status"],
                result["target_status"],
                1 if result["normalized_equal"] else 0,
                1 if result["warning_equal"] else 0,
                1 if result["message_equal"] else 0,
                1 if result["error_equal"] else 0,
                1 if result["artifact_equal"] else 0,
                stable_json(result["result_payload"]),
            ),
        )
        exception_record = build_behavior_exception(
            run_id,
            result,
            behavior_result_id,
            baseline_ref=baseline_ref,
            target_ref=target_ref,
        )
        if exception_record is not None:
            exceptions.append(exception_record)

    insert_exceptions(conn, exceptions)
    return exceptions


def upsert_testing_matrix_entries(conn: sqlite3.Connection, matrix_entries: list[dict]) -> None:
    for entry in matrix_entries:
        conn.execute(
            """
            INSERT INTO testing_matrix_entries (
                testing_matrix_entry_id, surface, primary_evidence_json, fixture_strategy_json,
                required_scenarios_json, acceptance_json
            ) VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT(testing_matrix_entry_id) DO UPDATE SET
                surface = excluded.surface,
                primary_evidence_json = excluded.primary_evidence_json,
                fixture_strategy_json = excluded.fixture_strategy_json,
                required_scenarios_json = excluded.required_scenarios_json,
                acceptance_json = excluded.acceptance_json
            """,
            (
                entry["surface"],
                entry["surface"],
                stable_json(ensure_list(entry.get("primary_evidence"))),
                stable_json(ensure_list(entry.get("fixture_strategy"))),
                stable_json(ensure_list(entry.get("required_scenarios"))),
                stable_json(entry.get("acceptance") or {}),
            ),
        )


def insert_contract_catalog(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    repo_root: pathlib.Path,
    contract_files: list[dict],
) -> dict[str, str]:
    file_ids: dict[str, str] = {}
    case_index = 1
    for file_index, file_record in enumerate(contract_files, start=1):
        file_id = f"{run_id}:contract_file:{file_index}"
        file_ids[file_record["file_path"]] = file_id
        conn.execute(
            """
            INSERT INTO contract_test_files (
                contract_test_file_id, run_id, file_path, family, module_family, surface,
                certainty_target, case_count, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                file_id,
                run_id,
                file_record["file_path"],
                file_record["family"],
                file_record["module_family"],
                file_record["surface"],
                file_record["certainty_target"],
                len(file_record["cases"]),
                stable_json(
                    {
                        "relative_path": file_record["file_path"],
                        "required_scenarios": file_record.get("required_scenarios", []),
                    }
                ),
            ),
        )
        for case in file_record["cases"]:
            conn.execute(
                """
                INSERT INTO contract_test_cases (
                    contract_test_case_id, run_id, contract_test_file_id, test_name, entity_hint,
                    primary_scenario, scenario_tags_json, certainty_target, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:contract_case:{case_index}",
                    run_id,
                    file_id,
                    case["test_name"],
                    case.get("entity_hint"),
                    case["primary_scenario"],
                    stable_json(case["scenario_tags"]),
                    file_record["certainty_target"],
                    stable_json({"module_family": file_record["module_family"]}),
                ),
            )
            case_index += 1
    return file_ids


def insert_contract_results(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    file_ids: dict[str, str],
    contract_results: list[dict],
) -> None:
    for index, result in enumerate(contract_results, start=1):
        file_id = file_ids[result["file_path"]]
        conn.execute(
            """
            INSERT INTO contract_results (
                contract_result_id, run_id, contract_test_file_id, execution_mode, status,
                test_count, failure_count, skip_count, result_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                f"{run_id}:contract_result:{index}",
                run_id,
                file_id,
                result["execution_mode"],
                result["status"],
                result.get("test_count"),
                result.get("failure_count"),
                result.get("skip_count"),
                stable_json(result["result_payload"]),
            ),
        )


def build_contract_exception(
    run_id: str,
    result: dict,
    *,
    target_ref: str,
    family: str,
) -> dict | None:
    status = result["status"]
    if status not in {"failed", "runner_error"}:
        return None

    exception_type = "contract_runner_failure" if status == "runner_error" else "contract_test_failure"
    severity = "high" if status == "runner_error" or family in {"module_contract", "compat", "golden_master"} else "medium"
    entity_key = f"contract::{result['file_path']}"
    return {
        "exception_id": f"{run_id}:contract_exception:{result['file_path']}",
        "exception_key": make_exception_key("contract_replay", exception_type, entity_key),
        "run_id": run_id,
        "audit_layer": "contract_replay",
        "entity_key": entity_key,
        "severity": severity,
        "exception_type": exception_type,
        "reason": f"Contract execution ended in `{status}` for `{result['file_path']}`.",
        "baseline_ref": None,
        "target_ref": target_ref,
        "evidence_json": stable_json(
            {
                "file_path": result["file_path"],
                "family": family,
                "execution_mode": result["execution_mode"],
                "test_count": result.get("test_count"),
                "failure_count": result.get("failure_count"),
                "skip_count": result.get("skip_count"),
                "result_payload": result["result_payload"],
            }
        ),
        "curation_status": "candidate",
        "disposition": None,
        "owner": None,
        "review_note": None,
        "updated_at": now_iso(),
        "curation_rule": None,
        "resolved_at": None,
    }


def insert_contract_exceptions(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    contract_files: list[dict],
    contract_results: list[dict],
    target_ref: str,
) -> list[dict]:
    family_by_file = {record["file_path"]: record["family"] for record in contract_files}
    exceptions: list[dict] = []
    for result in contract_results:
        record = build_contract_exception(
            run_id,
            result,
            target_ref=target_ref,
            family=family_by_file.get(result["file_path"], "general_unit"),
        )
        if record is not None:
            exceptions.append(record)
    insert_exceptions(conn, exceptions)
    return exceptions


def build_inventory_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    side: str,
    target_ref: str,
    snapshot: dict,
) -> dict:
    return {
        "run_id": run_id,
        "mode": "inventory",
        "repo_root": str(repo_root),
        "side": side,
        "target_ref": target_ref,
        "status": "completed_with_parse_failures" if snapshot["parse_failures"] else "completed",
        "counts": snapshot["counts"],
        "parse_failures": snapshot["parse_failures"],
    }


def build_surface_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    baseline_label: str,
    target_label: str,
    baseline_snapshot: dict,
    target_snapshot: dict,
    diffs: list[dict],
    exceptions: list[dict],
) -> dict:
    diff_counts = Counter(diff["diff_class"] for diff in diffs)
    severity_counts = Counter(diff["severity"] for diff in diffs)
    status = (
        "completed_with_parse_failures"
        if baseline_snapshot["parse_failures"] or target_snapshot["parse_failures"]
        else "completed"
    )
    return {
        "run_id": run_id,
        "mode": "surface",
        "repo_root": str(repo_root),
        "baseline": {
            "label": baseline_label,
            "counts": baseline_snapshot["counts"],
            "parse_failures": baseline_snapshot["parse_failures"],
        },
        "target": {
            "label": target_label,
            "counts": target_snapshot["counts"],
            "parse_failures": target_snapshot["parse_failures"],
        },
        "status": status,
        "total_diffs": len(diffs),
        "open_diff_count": sum(1 for diff in diffs if diff["status"] == "open"),
        "exception_count": len(exceptions),
        "diff_counts": dict(sorted(diff_counts.items())),
        "severity_counts": dict(sorted(severity_counts.items())),
        "sample_diffs": [
            {
                "diff_class": diff["diff_class"],
                "entity_key": diff["entity_key"],
                "severity": diff["severity"],
            }
            for diff in diffs[:25]
        ],
    }


def build_manifest_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    baseline_label: str,
    target_label: str,
    manifest_paths: list[str],
    comparisons: list[dict],
    exceptions: list[dict],
) -> dict:
    status_counts = Counter(comparison["status"] for comparison in comparisons)
    resolver_counts = Counter(
        comparison.get("target_resolver") or "null"
        for comparison in comparisons
    )
    return {
        "run_id": run_id,
        "mode": "manifest",
        "repo_root": str(repo_root),
        "baseline": baseline_label,
        "target": target_label,
        "status": "completed",
        "manifest_count": len(manifest_paths),
        "entry_count": len(comparisons),
        "exception_count": len(exceptions),
        "status_counts": dict(sorted(status_counts.items())),
        "resolver_counts": dict(sorted(resolver_counts.items())),
        "sample_comparisons": [
            {
                "manifest_path": comparison["manifest_path"],
                "entry_id": comparison["entry_id"],
                "status": comparison["status"],
                "target_resolver": comparison.get("target_resolver") or "null",
            }
            for comparison in comparisons[:25]
        ],
    }


def build_behavior_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    baseline_label: str,
    target_label: str,
    case_paths: list[pathlib.Path],
    behavior_results: list[dict],
    exceptions: list[dict],
) -> dict:
    status_counts = Counter(result["comparison_status"] for result in behavior_results)
    family_counts = Counter(result["family"] for result in behavior_results)
    return {
        "run_id": run_id,
        "mode": "behavior",
        "repo_root": str(repo_root),
        "baseline": baseline_label,
        "target": target_label,
        "status": "completed",
        "case_file_count": len(case_paths),
        "case_count": len(behavior_results),
        "exception_count": len(exceptions),
        "status_counts": dict(sorted(status_counts.items())),
        "family_counts": dict(sorted(family_counts.items())),
        "sample_results": [
            {
                "case_id": result["case_id"],
                "entity_key": result["entity_key"],
                "comparison_status": result["comparison_status"],
                "family": result["family"],
            }
            for result in behavior_results[:25]
        ],
    }


def build_contract_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    target_label: str,
    matrix_entries: list[dict],
    contract_files: list[dict],
    contract_results: list[dict],
    exceptions: list[dict],
) -> dict:
    family_counts = Counter(file_record["family"] for file_record in contract_files)
    module_family_counts = Counter(file_record["module_family"] for file_record in contract_files)
    surface_counts = Counter(file_record["surface"] for file_record in contract_files)
    scenario_counts: Counter[str] = Counter()
    for file_record in contract_files:
        for case in file_record["cases"]:
            for scenario_tag in case["scenario_tags"]:
                scenario_counts[scenario_tag] += 1

    execution_status_counts = Counter(result["status"] for result in contract_results)
    status = "completed"
    if any(result["status"] in {"failed", "runner_error"} for result in contract_results):
        status = "completed_with_failures"

    return {
        "run_id": run_id,
        "mode": "contracts",
        "repo_root": str(repo_root),
        "target": target_label,
        "status": status,
        "matrix_surface_count": len(matrix_entries),
        "test_file_count": len(contract_files),
        "test_case_count": sum(len(file_record["cases"]) for file_record in contract_files),
        "executed_file_count": len(contract_results),
        "exception_count": len(exceptions),
        "family_counts": dict(sorted(family_counts.items())),
        "module_family_counts": dict(sorted(module_family_counts.items())),
        "surface_counts": dict(sorted(surface_counts.items())),
        "scenario_counts": dict(sorted(scenario_counts.items())),
        "execution_status_counts": dict(sorted(execution_status_counts.items())),
        "sample_files": [
            {
                "file_path": file_record["file_path"],
                "family": file_record["family"],
                "module_family": file_record["module_family"],
            }
            for file_record in contract_files[:25]
        ],
    }


def parse_json_field(value):
    if value in (None, ""):
        return None
    if isinstance(value, (dict, list)):
        return value
    return json.loads(value)


def fetch_exception_records(
    conn: sqlite3.Connection,
    *,
    run_id: str | None = None,
    run_ids: list[str] | None = None,
    exception_ids: list[str] | None = None,
    exception_keys: list[str] | None = None,
    curation_statuses: list[str] | None = None,
    open_only: bool = False,
) -> list[dict]:
    clauses = []
    params: list[str] = []
    if run_id:
        clauses.append("e.run_id = ?")
        params.append(run_id)
    if run_ids:
        placeholders = ", ".join("?" for _ in run_ids)
        clauses.append(f"e.run_id IN ({placeholders})")
        params.extend(run_ids)
    if exception_ids:
        placeholders = ", ".join("?" for _ in exception_ids)
        clauses.append(f"e.exception_id IN ({placeholders})")
        params.extend(exception_ids)
    if exception_keys:
        placeholders = ", ".join("?" for _ in exception_keys)
        clauses.append(f"e.exception_key IN ({placeholders})")
        params.extend(exception_keys)
    if curation_statuses:
        placeholders = ", ".join("?" for _ in curation_statuses)
        clauses.append(f"e.curation_status IN ({placeholders})")
        params.extend(curation_statuses)
    if open_only:
        placeholders = ", ".join("?" for _ in TERMINAL_CURATION_STATUSES)
        clauses.append(f"e.curation_status NOT IN ({placeholders})")
        params.extend(sorted(TERMINAL_CURATION_STATUSES))

    where_sql = f"WHERE {' AND '.join(clauses)}" if clauses else ""
    rows = conn.execute(
        f"""
        SELECT
            e.exception_id,
            e.exception_key,
            e.run_id,
            e.audit_layer,
            e.entity_key,
            e.severity,
            e.exception_type,
            e.reason,
            e.baseline_ref,
            e.target_ref,
            e.evidence_json,
            e.curation_status,
            e.disposition,
            e.owner,
            e.review_note,
            e.updated_at,
            e.curation_rule,
            e.resolved_at,
            r.mode AS run_mode
        FROM exceptions e
        LEFT JOIN audit_runs r ON r.run_id = e.run_id
        {where_sql}
        ORDER BY
            CASE e.severity WHEN 'high' THEN 0 WHEN 'medium' THEN 1 WHEN 'low' THEN 2 ELSE 3 END,
            CASE e.curation_status WHEN 'candidate' THEN 0 WHEN 'triaged' THEN 1 ELSE 2 END,
            e.audit_layer,
            e.exception_type,
            e.entity_key
        """,
        params,
    ).fetchall()

    records: list[dict] = []
    for row in rows:
        record = dict(row)
        record["evidence"] = parse_json_field(record.pop("evidence_json"))
        record["is_open"] = exception_is_open(record)
        record["suggested_disposition"] = suggested_disposition_for_exception(record["exception_type"])
        records.append(record)
    return records


def build_exception_report(repo_root: str | pathlib.Path, records: list[dict], filters: dict) -> dict:
    status_counts = Counter(record["curation_status"] for record in records)
    severity_counts = Counter(record["severity"] for record in records)
    layer_counts = Counter(record["audit_layer"] for record in records)
    type_counts = Counter(record["exception_type"] for record in records)
    disposition_counts = Counter(
        record["disposition"] or "unset" for record in records if record.get("disposition")
    )
    open_records = [record for record in records if record["is_open"]]
    open_status_counts = Counter(record["curation_status"] for record in open_records)
    high_open_count = sum(1 for record in open_records if record["severity"] == "high")

    return {
        "generated_at": now_iso(),
        "repo_root": str(repo_root),
        "filters": filters,
        "summary": {
            "total_count": len(records),
            "open_count": len(open_records),
            "high_open_count": high_open_count,
            "decision_ready": high_open_count == 0,
            "status_counts": dict(sorted(status_counts.items())),
            "open_status_counts": dict(sorted(open_status_counts.items())),
            "severity_counts": dict(sorted(severity_counts.items())),
            "layer_counts": dict(sorted(layer_counts.items())),
            "type_counts": dict(sorted(type_counts.items())),
            "disposition_counts": dict(sorted(disposition_counts.items())),
            "auto_curated_count": sum(1 for record in records if record["curation_status"] == "auto_curated"),
        },
        "exceptions": records,
    }


def materialize_exception_reports(
    conn: sqlite3.Connection,
    *,
    repo_root: pathlib.Path,
    paths: dict[str, pathlib.Path],
    run_id: str | None = None,
    run_ids: list[str] | None = None,
    curation_statuses: list[str] | None = None,
    open_only: bool = False,
) -> dict:
    filters = {
        "run_id": run_id,
        "run_ids": run_ids or [],
        "curation_statuses": curation_statuses or [],
        "open_only": open_only,
    }
    report = build_exception_report(
        repo_root=repo_root,
        records=fetch_exception_records(
            conn,
            run_id=run_id,
            run_ids=run_ids,
            curation_statuses=curation_statuses,
            open_only=open_only,
        ),
        filters=filters,
    )
    write_json(paths["exceptions_json"], report)
    write_text(paths["exceptions_md"], render_exception_report_markdown(report))
    return report


def update_exception_record(
    conn: sqlite3.Connection,
    *,
    exception_id: str,
    curation_status: str | None = None,
    disposition: str | None = None,
    owner: str | None = None,
    review_note: str | None = None,
    curation_rule: str | None = None,
) -> None:
    current = conn.execute(
        """
        SELECT curation_status, disposition, owner, review_note
        FROM exceptions
        WHERE exception_id = ?
        """,
        (exception_id,),
    ).fetchone()
    if current is None:
        raise RuntimeError(f"Unknown exception id: {exception_id}")

    next_status = curation_status or str(current["curation_status"])
    next_disposition = disposition if disposition is not None else current["disposition"]
    next_owner = owner if owner is not None else current["owner"]
    next_note = review_note if review_note is not None else current["review_note"]
    resolved_at = now_iso() if next_status in TERMINAL_CURATION_STATUSES else None

    conn.execute(
        """
        UPDATE exceptions
        SET curation_status = ?, disposition = ?, owner = ?, review_note = ?,
            updated_at = ?, curation_rule = ?, resolved_at = ?
        WHERE exception_id = ?
        """,
        (
            next_status,
            next_disposition,
            next_owner,
            next_note,
            now_iso(),
            curation_rule,
            resolved_at,
            exception_id,
        ),
    )


def auto_curate_exceptions(
    conn: sqlite3.Connection,
    *,
    run_id: str | None = None,
) -> list[dict]:
    records = fetch_exception_records(conn, run_id=run_id, curation_statuses=["candidate"])
    updates: list[dict] = []
    for record in records:
        rule = exception_auto_rule(record["exception_type"])
        if not rule:
            continue
        update_exception_record(
            conn,
            exception_id=record["exception_id"],
            curation_status=str(rule["curation_status"]),
            disposition=str(rule["disposition"]),
            review_note=str(rule["note"]),
            curation_rule=str(rule["rule_name"]),
        )
        updates.append(
            {
                "exception_id": record["exception_id"],
                "exception_key": record.get("exception_key"),
                "rule_name": rule["rule_name"],
                "curation_status": rule["curation_status"],
                "disposition": rule["disposition"],
            }
        )
    return updates


def curate_exceptions(
    conn: sqlite3.Connection,
    *,
    exception_ids: list[str] | None,
    exception_keys: list[str] | None,
    curation_status: str | None,
    disposition: str | None,
    owner: str | None,
    review_note: str | None,
) -> list[str]:
    if not exception_ids and not exception_keys:
        raise RuntimeError("curate requires --exception-id or --exception-key")

    records = fetch_exception_records(
        conn,
        exception_ids=exception_ids,
        exception_keys=exception_keys,
    )
    if not records:
        raise RuntimeError("No exceptions matched the requested selection")

    updated_ids: list[str] = []
    for record in records:
        update_exception_record(
            conn,
            exception_id=record["exception_id"],
            curation_status=curation_status,
            disposition=disposition,
            owner=owner,
            review_note=review_note,
            curation_rule=None,
        )
        updated_ids.append(str(record["exception_id"]))
    return updated_ids


def run_audit_subcommand(arguments: list[str]) -> dict:
    command = [sys.executable, str(SCRIPT_DIR / "fidelity_audit.py"), *arguments]
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "audit subcommand failed"
        raise RuntimeError(message)
    return load_json_payload(completed.stdout)


def build_closeout_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    baseline_label: str,
    mainline_label: str,
    target_label: str,
    component_runs: dict[str, str],
    components: dict[str, dict],
    exceptions_report: dict,
    auto_curated_updates: list[dict],
    pinned_baseline_explicit: bool,
    contracts_executed: bool,
) -> dict:
    records = ensure_list(exceptions_report.get("exceptions"))
    open_records = [record for record in records if record.get("is_open")]
    open_by_layer = Counter(record["audit_layer"] for record in open_records)
    type_counts = Counter(record["exception_type"] for record in records)

    surface_gate = open_by_layer.get("surface_inventory", 0) == 0
    manifest_gate = open_by_layer.get("manifest_fidelity", 0) == 0
    behavior_gate = open_by_layer.get("behavior_replay", 0) == 0
    contract_gate = contracts_executed and open_by_layer.get("contract_replay", 0) == 0
    high_severity_gate = exceptions_report["summary"]["high_open_count"] == 0

    blockers: list[str] = []
    if not pinned_baseline_explicit:
        blockers.append("Pinned baseline was not provided explicitly for closeout.")
    if not surface_gate:
        blockers.append("Surface inventory still has open structural exceptions.")
    if not manifest_gate:
        blockers.append("Manifest fidelity still has open exceptions.")
    if not behavior_gate:
        blockers.append("Behavior replay still has open exceptions.")
    if not contracts_executed:
        blockers.append("Contract execution was not requested during closeout.")
    elif not contract_gate:
        blockers.append("Contract execution still has open exceptions.")
    if not high_severity_gate:
        blockers.append("High-severity open exceptions remain unresolved.")

    coverage_campaign_ready = (
        pinned_baseline_explicit
        and surface_gate
        and manifest_gate
        and behavior_gate
        and contract_gate
        and high_severity_gate
    )

    return {
        "run_id": run_id,
        "mode": "closeout",
        "repo_root": str(repo_root),
        "baseline": baseline_label,
        "mainline": mainline_label,
        "target": target_label,
        "status": "ready_for_coverage_campaign" if coverage_campaign_ready else "blocked",
        "proven_parity": coverage_campaign_ready,
        "auto_curated_count": len(auto_curated_updates),
        "component_runs": component_runs,
        "components": components,
        "readiness": {
            "pinned_baseline_explicit": pinned_baseline_explicit,
            "surface_gate": surface_gate,
            "manifest_gate": manifest_gate,
            "behavior_gate": behavior_gate,
            "contract_gate": contract_gate,
            "high_severity_gate": high_severity_gate,
            "coverage_campaign_ready": coverage_campaign_ready,
            "open_exception_count": exceptions_report["summary"]["open_count"],
            "high_open_exception_count": exceptions_report["summary"]["high_open_count"],
            "open_exception_layer_counts": dict(sorted(open_by_layer.items())),
            "contracts_executed": contracts_executed,
        },
        "exception_counts": {
            "total": exceptions_report["summary"]["total_count"],
            "open": exceptions_report["summary"]["open_count"],
            "high_open": exceptions_report["summary"]["high_open_count"],
            "type_counts": dict(sorted(type_counts.items())),
        },
        "blockers": blockers,
        "sample_open_exceptions": [
            {
                "audit_layer": record["audit_layer"],
                "exception_type": record["exception_type"],
                "entity_key": record["entity_key"],
                "severity": record["severity"],
                "suggested_disposition": record.get("suggested_disposition"),
            }
            for record in open_records[:25]
        ],
    }


def bootstrap_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=True)
    print(
        json.dumps(
            {
                "status": "bootstrap_ready",
                "schema_version": SCHEMA_VERSION,
                "repo_root": str(repo_root),
                "artifacts": {key: str(value) for key, value in paths.items()},
            }
        )
    )
    return 0


def inventory_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"inventory-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "inventory",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "side": args.side,
            "target_ref": args.target_ref,
        },
    )
    snapshot = run_inventory_snapshot(repo_root)
    summary = build_inventory_summary(
        run_id=run_id,
        repo_root=repo_root,
        side=args.side,
        target_ref=args.target_ref,
        snapshot=snapshot,
    )

    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=None,
            target_ref=args.target_ref,
            mode="inventory",
            status=summary["status"],
            summary=summary,
        )
        insert_inventory_snapshot(conn, run_id=run_id, side=args.side, snapshot=snapshot)
        conn.commit()
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_inventory_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "inventory",
            "status": summary["status"],
            "run_id": run_id,
            "side": args.side,
            "target_ref": args.target_ref,
            "entity_count": len(snapshot["entities"]),
            "parse_failures": len(snapshot["parse_failures"]),
        },
    )
    print(json.dumps(summary))
    return 0


def surface_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"surface-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    baseline_source = resolve_side_source(
        repo_root,
        side="baseline",
        ref=args.baseline_ref,
        path=args.baseline_path,
        default_ref="main",
        allow_worktree=False,
        archive_paths=["DESCRIPTION", "NAMESPACE", "R"],
    )
    target_source = resolve_side_source(
        repo_root,
        side="target",
        ref=args.target_ref,
        path=args.target_path,
        default_ref="WORKTREE",
        allow_worktree=True,
        archive_paths=["DESCRIPTION", "NAMESPACE", "R"],
    )

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "surface",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "baseline": baseline_source["label"],
            "target": target_source["label"],
        },
    )

    try:
        baseline_snapshot = run_inventory_snapshot(baseline_source["path"])
        target_snapshot = run_inventory_snapshot(target_source["path"])
    finally:
        baseline_source["cleanup"]()
        target_source["cleanup"]()

    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=baseline_source["label"],
            target_ref=target_source["label"],
            mode="surface",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "surface",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        baseline_refs = insert_inventory_snapshot(conn, run_id=run_id, side="baseline", snapshot=baseline_snapshot)
        target_refs = insert_inventory_snapshot(conn, run_id=run_id, side="target", snapshot=target_snapshot)

        diffs = sort_diffs(
            compare_file_surfaces(baseline_snapshot, target_snapshot, baseline_refs, target_refs)
            + compare_export_surfaces(baseline_snapshot, target_snapshot, baseline_refs, target_refs)
            + compare_entity_surfaces(baseline_snapshot, target_snapshot, baseline_refs, target_refs)
        )
        insert_inventory_diffs(conn, run_id=run_id, diffs=diffs)
        exceptions = insert_surface_exceptions(
            conn,
            run_id=run_id,
            diffs=diffs,
            baseline_ref=baseline_source["label"],
            target_ref=target_source["label"],
        )

        summary = build_surface_summary(
            run_id=run_id,
            repo_root=repo_root,
            baseline_label=baseline_source["label"],
            target_label=target_source["label"],
            baseline_snapshot=baseline_snapshot,
            target_snapshot=target_snapshot,
            diffs=diffs,
            exceptions=exceptions,
        )
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
        materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_surface_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "surface",
            "status": summary["status"],
            "run_id": run_id,
            "baseline": summary["baseline"]["label"],
            "target": summary["target"]["label"],
            "total_diffs": summary["total_diffs"],
            "open_diff_count": summary["open_diff_count"],
        },
    )
    print(json.dumps(summary))
    return 0


def manifest_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"manifest-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    baseline_source = resolve_side_source(
        repo_root,
        side="baseline",
        ref=args.baseline_ref,
        path=args.baseline_path,
        default_ref="main",
        allow_worktree=False,
        archive_paths=["DESCRIPTION", "NAMESPACE", "R"],
    )
    target_source = resolve_side_source(
        repo_root,
        side="target",
        ref=args.target_ref,
        path=args.target_path,
        default_ref="WORKTREE",
        allow_worktree=True,
        archive_paths=["DESCRIPTION", "NAMESPACE", "R"],
    )
    manifest_paths = list_manifest_paths(repo_root, args.manifest_path)
    if not manifest_paths:
        raise RuntimeError("No manifest files found for manifest audit")

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "manifest",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "baseline": baseline_source["label"],
            "target": target_source["label"],
            "manifest_count": len(manifest_paths),
        },
    )

    comparisons: list[dict] = []
    try:
        for manifest_path in manifest_paths:
            payload = run_manifest_compare(manifest_path, baseline_source["path"], target_source["path"])
            for comparison in normalize_manifest_comparisons(payload):
                comparison["manifest_path"] = relative_label(manifest_path, repo_root)
                comparisons.append(comparison)
    finally:
        baseline_source["cleanup"]()
        target_source["cleanup"]()

    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=baseline_source["label"],
            target_ref=target_source["label"],
            mode="manifest",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "manifest",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        exceptions = insert_manifest_results(
            conn,
            run_id=run_id,
            manifest_comparisons=comparisons,
            baseline_ref=baseline_source["label"],
            target_ref=target_source["label"],
        )
        summary = build_manifest_summary(
            run_id=run_id,
            repo_root=repo_root,
            baseline_label=baseline_source["label"],
            target_label=target_source["label"],
            manifest_paths=[relative_label(path, repo_root) for path in manifest_paths],
            comparisons=comparisons,
            exceptions=exceptions,
        )
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
        materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_manifest_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "manifest",
            "status": summary["status"],
            "run_id": run_id,
            "baseline": summary["baseline"],
            "target": summary["target"],
            "manifest_count": summary["manifest_count"],
            "entry_count": summary["entry_count"],
            "exception_count": summary["exception_count"],
        },
    )
    print(json.dumps(summary))
    return 0


def behavior_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"behavior-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    baseline_source = resolve_side_source(
        repo_root,
        side="baseline",
        ref=args.baseline_ref,
        path=args.baseline_path,
        default_ref="main",
        allow_worktree=False,
        archive_paths=["DESCRIPTION", "NAMESPACE", "R"],
    )
    target_source = resolve_side_source(
        repo_root,
        side="target",
        ref=args.target_ref,
        path=args.target_path,
        default_ref="WORKTREE",
        allow_worktree=True,
        archive_paths=["DESCRIPTION", "NAMESPACE", "R"],
    )
    case_paths, cases = load_behavior_cases(repo_root, args.case_path)

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "behavior",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "baseline": baseline_source["label"],
            "target": target_source["label"],
            "case_file_count": len(case_paths),
            "case_count": len(cases),
        },
    )

    behavior_results: list[dict] = []
    try:
        for case in cases:
            baseline_payload = run_behavior_case(
                case["_catalog_file"], case["id"], baseline_source["path"], "baseline"
            )
            target_payload = run_behavior_case(
                case["_catalog_file"], case["id"], target_source["path"], "target"
            )
            behavior_results.append(compare_behavior_payloads(case, baseline_payload, target_payload))
    finally:
        baseline_source["cleanup"]()
        target_source["cleanup"]()

    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=baseline_source["label"],
            target_ref=target_source["label"],
            mode="behavior",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "behavior",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        upsert_behavior_cases(conn, cases)
        exceptions = insert_behavior_results(
            conn,
            run_id=run_id,
            behavior_results=behavior_results,
            baseline_ref=baseline_source["label"],
            target_ref=target_source["label"],
        )
        summary = build_behavior_summary(
            run_id=run_id,
            repo_root=repo_root,
            baseline_label=baseline_source["label"],
            target_label=target_source["label"],
            case_paths=case_paths,
            behavior_results=behavior_results,
            exceptions=exceptions,
        )
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
        materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_behavior_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "behavior",
            "status": summary["status"],
            "run_id": run_id,
            "baseline": summary["baseline"],
            "target": summary["target"],
            "case_file_count": summary["case_file_count"],
            "case_count": summary["case_count"],
            "exception_count": summary["exception_count"],
        },
    )
    print(json.dumps(summary))
    return 0


def contracts_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"contracts-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    target_source = resolve_side_source(
        repo_root,
        side="target",
        ref=args.target_ref,
        path=args.target_path,
        default_ref="WORKTREE",
        allow_worktree=True,
        archive_paths=None,
    )

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "contracts",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "target": target_source["label"],
            "execute": bool(args.execute),
        },
    )

    try:
        matrix_entries = load_testing_matrix_entries()
        contract_files = catalog_contract_files(target_source["path"], args.test_file)
        contract_results: list[dict] = []
        if args.execute:
            for file_record in contract_files:
                execution = run_contract_test_file(
                    target_source["path"],
                    file_record["absolute_path"],
                    load_package=bool(args.load_package),
                )
                contract_results.append(
                    {
                        "file_path": file_record["file_path"],
                        "execution_mode": "load_package" if args.load_package else "test_file_only",
                        "status": execution["status"],
                        "test_count": execution.get("test_count"),
                        "failure_count": execution.get("failure_count"),
                        "skip_count": execution.get("skip_count"),
                        "result_payload": execution,
                    }
                )
    finally:
        target_source["cleanup"]()

    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=None,
            target_ref=target_source["label"],
            mode="contracts",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "contracts",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        upsert_testing_matrix_entries(conn, matrix_entries)
        file_ids = insert_contract_catalog(conn, run_id=run_id, repo_root=repo_root, contract_files=contract_files)
        insert_contract_results(conn, run_id=run_id, file_ids=file_ids, contract_results=contract_results)
        exceptions = insert_contract_exceptions(
            conn,
            run_id=run_id,
            contract_files=contract_files,
            contract_results=contract_results,
            target_ref=target_source["label"],
        )
        summary = build_contract_summary(
            run_id=run_id,
            repo_root=repo_root,
            target_label=target_source["label"],
            matrix_entries=matrix_entries,
            contract_files=contract_files,
            contract_results=contract_results,
            exceptions=exceptions,
        )
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
        materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_contract_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "contracts",
            "status": summary["status"],
            "run_id": run_id,
            "target": summary["target"],
            "test_file_count": summary["test_file_count"],
            "test_case_count": summary["test_case_count"],
            "executed_file_count": summary["executed_file_count"],
        },
    )
    print(json.dumps(summary))
    return 0


def exceptions_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    conn = connect_db(paths["db_path"])
    try:
        auto_updates = auto_curate_exceptions(conn, run_id=args.run_id) if args.auto_curate else []
        if auto_updates:
            conn.commit()
        report = materialize_exception_reports(
            conn,
            repo_root=repo_root,
            paths=paths,
            run_id=args.run_id,
            curation_statuses=args.curation_status,
            open_only=bool(args.open_only),
        )
    finally:
        conn.close()

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "exceptions",
            "status": "reported",
            "repo_root": str(repo_root),
            "run_id": args.run_id,
            "open_only": bool(args.open_only),
            "auto_curated_count": len(auto_updates),
            "open_count": report["summary"]["open_count"],
            "high_open_count": report["summary"]["high_open_count"],
        },
    )
    print(json.dumps(report))
    return 0


def curate_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    conn = connect_db(paths["db_path"])
    try:
        updated_ids = curate_exceptions(
            conn,
            exception_ids=args.exception_id,
            exception_keys=args.exception_key,
            curation_status=args.curation_status,
            disposition=args.disposition,
            owner=args.owner,
            review_note=args.note,
        )
        conn.commit()
        report = materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "exceptions",
            "status": "curated",
            "repo_root": str(repo_root),
            "updated_count": len(updated_ids),
            "curation_status": args.curation_status,
            "disposition": args.disposition,
            "owner": args.owner,
        },
    )
    print(
        json.dumps(
            {
                "status": "curated",
                "updated_count": len(updated_ids),
                "updated_ids": updated_ids,
                "open_count": report["summary"]["open_count"],
                "high_open_count": report["summary"]["high_open_count"],
            }
        )
    )
    return 0


def closeout_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"closeout-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    def source_key(path_value: str | None, ref_value: str | None, default_ref: str) -> tuple[str, str]:
        if path_value:
            return ("path", str(pathlib.Path(path_value).resolve()))
        return ("ref", ref_value or default_ref)

    baseline_key = source_key(args.baseline_path, args.baseline_ref, "main")
    mainline_key = source_key(args.mainline_path, args.mainline_ref, "main")
    baseline_explicit = bool(args.baseline_path or args.baseline_ref)

    component_runs = {
        "baseline_surface": f"{run_id}:surface:baseline",
        "mainline_surface": f"{run_id}:surface:mainline",
        "manifest": f"{run_id}:manifest",
        "behavior": f"{run_id}:behavior",
        "contracts": f"{run_id}:contracts",
    }

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "closeout",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "baseline": baseline_key[1],
            "mainline": mainline_key[1],
            "target_ref": args.target_ref or "WORKTREE",
            "contracts_execute": bool(args.contracts_execute),
        },
    )

    surface_baseline_args = [
        "surface",
        "--repo-root",
        str(repo_root),
        "--audit-dir",
        args.audit_dir,
        "--run-id",
        component_runs["baseline_surface"],
    ]
    if args.baseline_path:
        surface_baseline_args.extend(["--baseline-path", args.baseline_path])
    else:
        surface_baseline_args.extend(["--baseline-ref", args.baseline_ref or "main"])
    if args.target_path:
        surface_baseline_args.extend(["--target-path", args.target_path])
    else:
        surface_baseline_args.extend(["--target-ref", args.target_ref or "WORKTREE"])

    baseline_surface_summary = run_audit_subcommand(surface_baseline_args)

    if mainline_key == baseline_key:
        mainline_surface_summary = baseline_surface_summary
        component_runs["mainline_surface"] = component_runs["baseline_surface"]
    else:
        surface_mainline_args = [
            "surface",
            "--repo-root",
            str(repo_root),
            "--audit-dir",
            args.audit_dir,
            "--run-id",
            component_runs["mainline_surface"],
        ]
        if args.mainline_path:
            surface_mainline_args.extend(["--baseline-path", args.mainline_path])
        else:
            surface_mainline_args.extend(["--baseline-ref", args.mainline_ref or "main"])
        if args.target_path:
            surface_mainline_args.extend(["--target-path", args.target_path])
        else:
            surface_mainline_args.extend(["--target-ref", args.target_ref or "WORKTREE"])
        mainline_surface_summary = run_audit_subcommand(surface_mainline_args)

    manifest_args = [
        "manifest",
        "--repo-root",
        str(repo_root),
        "--audit-dir",
        args.audit_dir,
        "--run-id",
        component_runs["manifest"],
    ]
    if args.baseline_path:
        manifest_args.extend(["--baseline-path", args.baseline_path])
    else:
        manifest_args.extend(["--baseline-ref", args.baseline_ref or "main"])
    if args.target_path:
        manifest_args.extend(["--target-path", args.target_path])
    else:
        manifest_args.extend(["--target-ref", args.target_ref or "WORKTREE"])
    for manifest_path in ensure_list(args.manifest_path):
        manifest_args.extend(["--manifest-path", manifest_path])
    manifest_summary = run_audit_subcommand(manifest_args)

    behavior_args = [
        "behavior",
        "--repo-root",
        str(repo_root),
        "--audit-dir",
        args.audit_dir,
        "--run-id",
        component_runs["behavior"],
    ]
    if args.baseline_path:
        behavior_args.extend(["--baseline-path", args.baseline_path])
    else:
        behavior_args.extend(["--baseline-ref", args.baseline_ref or "main"])
    if args.target_path:
        behavior_args.extend(["--target-path", args.target_path])
    else:
        behavior_args.extend(["--target-ref", args.target_ref or "WORKTREE"])
    for case_path in ensure_list(args.case_path):
        behavior_args.extend(["--case-path", case_path])
    behavior_summary = run_audit_subcommand(behavior_args)

    contracts_args = [
        "contracts",
        "--repo-root",
        str(repo_root),
        "--audit-dir",
        args.audit_dir,
        "--run-id",
        component_runs["contracts"],
    ]
    if args.target_path:
        contracts_args.extend(["--target-path", args.target_path])
    else:
        contracts_args.extend(["--target-ref", args.target_ref or "WORKTREE"])
    for test_file in ensure_list(args.test_file):
        contracts_args.extend(["--test-file", test_file])
    if args.contracts_execute:
        contracts_args.append("--execute")
    if args.contracts_load_package:
        contracts_args.append("--load-package")
    contracts_summary = run_audit_subcommand(contracts_args)

    component_run_ids = sorted(set(component_runs.values()))
    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=baseline_surface_summary["baseline"]["label"],
            target_ref=baseline_surface_summary["target"]["label"],
            mode="closeout",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "closeout",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        auto_curated_updates: list[dict] = []
        for component_run_id in component_run_ids:
            auto_curated_updates.extend(auto_curate_exceptions(conn, run_id=component_run_id))
        conn.commit()

        full_exception_report = build_exception_report(
            repo_root=repo_root,
            records=fetch_exception_records(conn, run_ids=component_run_ids),
            filters={"run_ids": component_run_ids, "open_only": False},
        )
        unresolved_exception_report = build_exception_report(
            repo_root=repo_root,
            records=fetch_exception_records(conn, run_ids=component_run_ids, open_only=True),
            filters={"run_ids": component_run_ids, "open_only": True},
        )

        closeout_summary = build_closeout_summary(
            run_id=run_id,
            repo_root=repo_root,
            baseline_label=baseline_surface_summary["baseline"]["label"],
            mainline_label=mainline_surface_summary["baseline"]["label"],
            target_label=baseline_surface_summary["target"]["label"],
            component_runs=component_runs,
            components={
                "baseline_surface": {
                    "run_id": component_runs["baseline_surface"],
                    "status": baseline_surface_summary["status"],
                    "total_diffs": baseline_surface_summary["total_diffs"],
                    "open_diff_count": baseline_surface_summary["open_diff_count"],
                    "exception_count": baseline_surface_summary["exception_count"],
                },
                "mainline_surface": {
                    "run_id": component_runs["mainline_surface"],
                    "status": mainline_surface_summary["status"],
                    "total_diffs": mainline_surface_summary["total_diffs"],
                    "open_diff_count": mainline_surface_summary["open_diff_count"],
                    "exception_count": mainline_surface_summary["exception_count"],
                },
                "manifest": {
                    "run_id": component_runs["manifest"],
                    "status": manifest_summary["status"],
                    "entry_count": manifest_summary["entry_count"],
                    "exception_count": manifest_summary["exception_count"],
                },
                "behavior": {
                    "run_id": component_runs["behavior"],
                    "status": behavior_summary["status"],
                    "case_count": behavior_summary["case_count"],
                    "exception_count": behavior_summary["exception_count"],
                },
                "contracts": {
                    "run_id": component_runs["contracts"],
                    "status": contracts_summary["status"],
                    "test_file_count": contracts_summary["test_file_count"],
                    "executed_file_count": contracts_summary["executed_file_count"],
                    "exception_count": contracts_summary["exception_count"],
                },
            },
            exceptions_report=full_exception_report,
            auto_curated_updates=auto_curated_updates,
            pinned_baseline_explicit=baseline_explicit,
            contracts_executed=bool(args.contracts_execute),
        )
        update_audit_run(conn, run_id=run_id, status=closeout_summary["status"], summary=closeout_summary)
        conn.commit()
    finally:
        conn.close()

    write_json(paths["summary_json"], closeout_summary)
    write_text(paths["summary_md"], render_closeout_summary_markdown(closeout_summary))
    write_json(paths["closeout_json"], closeout_summary)
    write_text(paths["closeout_md"], render_closeout_summary_markdown(closeout_summary))
    write_json(paths["exceptions_json"], unresolved_exception_report)
    write_text(paths["exceptions_md"], render_exception_report_markdown(unresolved_exception_report))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "closeout",
            "status": closeout_summary["status"],
            "run_id": run_id,
            "component_runs": component_runs,
            "open_exception_count": closeout_summary["readiness"]["open_exception_count"],
            "high_open_exception_count": closeout_summary["readiness"]["high_open_exception_count"],
            "coverage_campaign_ready": closeout_summary["readiness"]["coverage_campaign_ready"],
        },
    )
    print(json.dumps(closeout_summary))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Refactor fidelity audit tooling")
    subparsers = parser.add_subparsers(dest="command", required=True)

    bootstrap_parser = subparsers.add_parser("bootstrap", help="Create the local audit store")
    bootstrap_parser.add_argument("--repo-root", default=".")
    bootstrap_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    bootstrap_parser.set_defaults(func=bootstrap_command)

    inventory_parser = subparsers.add_parser("inventory", help="Capture an inventory snapshot into the audit store")
    inventory_parser.add_argument("--repo-root", default=".")
    inventory_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    inventory_parser.add_argument("--side", default="target", choices=["baseline", "target"])
    inventory_parser.add_argument("--target-ref", default="HEAD")
    inventory_parser.add_argument("--run-id")
    inventory_parser.set_defaults(func=inventory_command)

    surface_parser = subparsers.add_parser("surface", help="Compare baseline and target inventory surfaces")
    surface_parser.add_argument("--repo-root", default=".")
    surface_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    surface_parser.add_argument("--baseline-ref")
    surface_parser.add_argument("--target-ref")
    surface_parser.add_argument("--baseline-path")
    surface_parser.add_argument("--target-path")
    surface_parser.add_argument("--run-id")
    surface_parser.set_defaults(func=surface_command)

    manifest_parser = subparsers.add_parser("manifest", help="Compare manifest-managed extracted blocks")
    manifest_parser.add_argument("--repo-root", default=".")
    manifest_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    manifest_parser.add_argument("--baseline-ref")
    manifest_parser.add_argument("--target-ref")
    manifest_parser.add_argument("--baseline-path")
    manifest_parser.add_argument("--target-path")
    manifest_parser.add_argument("--manifest-path", action="append")
    manifest_parser.add_argument("--run-id")
    manifest_parser.set_defaults(func=manifest_command)

    behavior_parser = subparsers.add_parser("behavior", help="Replay behavior cases against baseline and target")
    behavior_parser.add_argument("--repo-root", default=".")
    behavior_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    behavior_parser.add_argument("--baseline-ref")
    behavior_parser.add_argument("--target-ref")
    behavior_parser.add_argument("--baseline-path")
    behavior_parser.add_argument("--target-path")
    behavior_parser.add_argument("--case-path", action="append")
    behavior_parser.add_argument("--run-id")
    behavior_parser.set_defaults(func=behavior_command)

    contracts_parser = subparsers.add_parser("contracts", help="Catalog and optionally execute testthat contract families")
    contracts_parser.add_argument("--repo-root", default=".")
    contracts_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    contracts_parser.add_argument("--target-ref")
    contracts_parser.add_argument("--target-path")
    contracts_parser.add_argument("--test-file", action="append")
    contracts_parser.add_argument("--execute", action="store_true")
    contracts_parser.add_argument("--load-package", action="store_true")
    contracts_parser.add_argument("--run-id")
    contracts_parser.set_defaults(func=contracts_command)

    exceptions_parser = subparsers.add_parser("exceptions", help="Materialize the current exception ledger and optional auto-curation")
    exceptions_parser.add_argument("--repo-root", default=".")
    exceptions_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    exceptions_parser.add_argument("--run-id")
    exceptions_parser.add_argument("--curation-status", action="append", choices=CURATION_STATUS_CHOICES)
    exceptions_parser.add_argument("--open-only", action="store_true")
    exceptions_parser.add_argument("--auto-curate", action="store_true")
    exceptions_parser.set_defaults(func=exceptions_command)

    curate_parser = subparsers.add_parser("curate", help="Apply manual curation metadata to exception rows")
    curate_parser.add_argument("--repo-root", default=".")
    curate_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    curate_parser.add_argument("--exception-id", action="append")
    curate_parser.add_argument("--exception-key", action="append")
    curate_parser.add_argument("--curation-status", required=True, choices=CURATION_STATUS_CHOICES)
    curate_parser.add_argument("--disposition")
    curate_parser.add_argument("--owner")
    curate_parser.add_argument("--note")
    curate_parser.set_defaults(func=curate_command)

    closeout_parser = subparsers.add_parser("closeout", help="Run the integrated closeout audit against baseline and mainline")
    closeout_parser.add_argument("--repo-root", default=".")
    closeout_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    closeout_parser.add_argument("--baseline-ref")
    closeout_parser.add_argument("--baseline-path")
    closeout_parser.add_argument("--mainline-ref")
    closeout_parser.add_argument("--mainline-path")
    closeout_parser.add_argument("--target-ref")
    closeout_parser.add_argument("--target-path")
    closeout_parser.add_argument("--manifest-path", action="append")
    closeout_parser.add_argument("--case-path", action="append")
    closeout_parser.add_argument("--test-file", action="append")
    closeout_parser.add_argument("--contracts-execute", action="store_true")
    closeout_parser.add_argument("--contracts-load-package", action="store_true")
    closeout_parser.add_argument("--run-id")
    closeout_parser.set_defaults(func=closeout_command)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
