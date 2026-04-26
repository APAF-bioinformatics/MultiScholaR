#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import io
import json
import os
import pathlib
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import uuid
from collections import Counter, defaultdict
from datetime import datetime, timezone


SCHEMA_VERSION = 8
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
DEFAULT_COVERAGE_JSON = "latest-coverage.json"
DEFAULT_COVERAGE_MD = "latest-coverage.md"
DEFAULT_COVERAGE_COMPARE_JSON = "latest-coverage-compare.json"
DEFAULT_COVERAGE_COMPARE_MD = "latest-coverage-compare.md"
DEFAULT_COVERAGE_EVIDENCE_JSON = "latest-coverage-evidence.json"
DEFAULT_COVERAGE_EVIDENCE_MD = "latest-coverage-evidence.md"
DEFAULT_BUNDLES_JSON = "latest-bundles.json"
DEFAULT_BUNDLES_MD = "latest-bundles.md"
DEFAULT_BEHAVIOR_CASE_GLOB = "behavior-cases*.json"
DEFAULT_CURATION_GLOB = "fidelity-curations*.json"
DEFAULT_TEST_FILE_GLOB = "test-*.R"

SURFACE_OPEN_DIFFS = {
    "missing_definition",
    "duplicate_definition",
    "definition_drift",
    "missing_export",
    "extra_export",
    "missing_file",
    "collate_drift",
}

DIFF_SEVERITY = {
    "missing_definition": "high",
    "extra_definition": "low",
    "duplicate_definition": "high",
    "definition_drift": "medium",
    "moved_definition": "info",
    "missing_export": "high",
    "extra_export": "medium",
    "missing_file": "medium",
    "extra_file": "info",
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
    {
        "family": "compat",
        "suffix": "shared",
        "surface": "cross_module_and_workflow_segments",
        "certainty_target": "T5_contract_replayed",
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
    "source_lineage_gap_target_resolved": {
        "curation_status": "auto_curated",
        "disposition": "equivalent_target_resolved_without_baseline_ancestor",
        "note": "Auto-curated because the target block resolved cleanly even though the manifest source no longer maps directly onto the pinned baseline.",
        "rule_name": "auto:source_lineage_gap_target_resolved",
    },
}

REDUNDANT_SURFACE_MANIFEST_RULE = {
    "curation_status": "auto_curated",
    "disposition": "covered_by_manifest_semantic_review",
    "note": "Auto-curated because the same entity is already blocked by an open manifest semantic mismatch in this closeout batch.",
    "rule_name": "auto:surface_definition_drift_manifest_overlap",
}

BUNDLE_TYPE_PRIORITY = {
    "lineage_family": 6,
    "wrapper_entrypoint": 5,
    "manual_merge_canonical": 4,
    "duplicate_canonical": 4,
    "s4_method_surface": 3,
    "helper_surface": 2,
    "normalized_exact_surface": 1,
}

DISPOSITION_EVIDENCE_DEPTH = {
    "equivalent_target_resolved_without_baseline_ancestor": "T2_structural",
    "equivalent_target_resolution_asymmetry": "T2_structural",
    "equivalent_whitespace_only": "T3_normalized_exact",
    "equivalent_comment_only": "T3_ast_exact",
    "covered_by_manifest_semantic_review": "T4_manifest_covered",
    "helper_equivalent_under_tests": "T4_behavior_replayed",
    "method_equivalent_under_tests": "T4_behavior_replayed",
    "manual_merge_canonicalized": "T4_behavior_replayed",
    "canonical_duplicate_retired": "T4_behavior_replayed",
    "entrypoint_shell_equivalent_under_tests": "T5_contract_replayed",
}

EVIDENCE_DEPTH_PRIORITY = {
    "T2_structural": 1,
    "T3_normalized_exact": 2,
    "T3_ast_exact": 2,
    "T4_manifest_covered": 3,
    "T4_behavior_replayed": 4,
    "T5_contract_replayed": 5,
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

        CREATE TABLE IF NOT EXISTS coverage_runs (
            coverage_run_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            side TEXT NOT NULL,
            source_label TEXT NOT NULL,
            bundle_manifest_path TEXT,
            status TEXT NOT NULL,
            tool_status TEXT NOT NULL,
            test_file_count INTEGER NOT NULL,
            bundle_count INTEGER NOT NULL,
            file_count INTEGER NOT NULL,
            line_percent REAL,
            lines_total INTEGER,
            lines_covered INTEGER,
            summary_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE INDEX IF NOT EXISTS idx_coverage_runs_run
            ON coverage_runs(run_id);

        CREATE TABLE IF NOT EXISTS coverage_files (
            coverage_file_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            coverage_run_id TEXT NOT NULL,
            file_path TEXT NOT NULL,
            line_percent REAL,
            lines_total INTEGER,
            lines_covered INTEGER,
            covered_lines_json TEXT,
            uncovered_lines_json TEXT,
            result_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(coverage_run_id) REFERENCES coverage_runs(coverage_run_id)
        );

        CREATE INDEX IF NOT EXISTS idx_coverage_files_run
            ON coverage_files(run_id);

        CREATE TABLE IF NOT EXISTS coverage_bundles (
            coverage_bundle_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            bundle_id TEXT NOT NULL,
            bundle_type TEXT NOT NULL,
            source_file_path TEXT,
            primary_entity_key TEXT,
            baseline_ref TEXT,
            target_ref TEXT,
            linked_exception_count INTEGER,
            scenario_class TEXT,
            differential_replay_status TEXT,
            coverage_gate_status TEXT,
            review_note TEXT,
            priority_rank INTEGER,
            priority_score REAL,
            metadata_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id)
        );

        CREATE INDEX IF NOT EXISTS idx_coverage_bundles_run
            ON coverage_bundles(run_id);

        CREATE TABLE IF NOT EXISTS coverage_bundle_members (
            coverage_bundle_member_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            coverage_bundle_id TEXT NOT NULL,
            entity_key TEXT,
            file_path TEXT NOT NULL,
            line_start INTEGER,
            line_end INTEGER,
            metadata_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(coverage_bundle_id) REFERENCES coverage_bundles(coverage_bundle_id)
        );

        CREATE INDEX IF NOT EXISTS idx_coverage_bundle_members_run
            ON coverage_bundle_members(run_id);

        CREATE TABLE IF NOT EXISTS coverage_test_links (
            coverage_test_link_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            coverage_bundle_id TEXT,
            test_file_path TEXT NOT NULL,
            source TEXT NOT NULL,
            metadata_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(coverage_bundle_id) REFERENCES coverage_bundles(coverage_bundle_id)
        );

        CREATE INDEX IF NOT EXISTS idx_coverage_test_links_run
            ON coverage_test_links(run_id);

        CREATE TABLE IF NOT EXISTS coverage_bundle_results (
            coverage_bundle_result_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            coverage_bundle_id TEXT NOT NULL,
            status TEXT NOT NULL,
            line_percent REAL,
            lines_total INTEGER,
            lines_covered INTEGER,
            result_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(coverage_bundle_id) REFERENCES coverage_bundles(coverage_bundle_id)
        );

        CREATE INDEX IF NOT EXISTS idx_coverage_bundle_results_run
            ON coverage_bundle_results(run_id);

        CREATE TABLE IF NOT EXISTS coverage_exception_links (
            coverage_exception_link_id TEXT PRIMARY KEY,
            run_id TEXT NOT NULL,
            coverage_bundle_id TEXT NOT NULL,
            exception_id TEXT NOT NULL,
            exception_key TEXT,
            audit_layer TEXT NOT NULL,
            entity_key TEXT NOT NULL,
            curation_status TEXT NOT NULL,
            disposition TEXT,
            metadata_json TEXT,
            FOREIGN KEY(run_id) REFERENCES audit_runs(run_id),
            FOREIGN KEY(coverage_bundle_id) REFERENCES coverage_bundles(coverage_bundle_id),
            FOREIGN KEY(exception_id) REFERENCES exceptions(exception_id)
        );

        CREATE INDEX IF NOT EXISTS idx_coverage_exception_links_run
            ON coverage_exception_links(run_id);

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
    ensure_column(conn, "coverage_bundles", "primary_entity_key", "TEXT")
    ensure_column(conn, "coverage_bundles", "baseline_ref", "TEXT")
    ensure_column(conn, "coverage_bundles", "target_ref", "TEXT")
    ensure_column(conn, "coverage_bundles", "linked_exception_count", "INTEGER")
    ensure_column(conn, "coverage_bundles", "scenario_class", "TEXT")
    ensure_column(conn, "coverage_bundles", "differential_replay_status", "TEXT")
    ensure_column(conn, "coverage_bundles", "coverage_gate_status", "TEXT")
    ensure_column(conn, "coverage_bundles", "review_note", "TEXT")
    ensure_column(conn, "coverage_bundles", "priority_rank", "INTEGER")
    ensure_column(conn, "coverage_bundles", "priority_score", "REAL")
    ensure_column(conn, "coverage_bundles", "side", "TEXT")
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
    conn.execute(
        "CREATE INDEX IF NOT EXISTS idx_coverage_exception_links_exception_key ON coverage_exception_links(exception_key)"
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
        "coverage_json": reports_dir / DEFAULT_COVERAGE_JSON,
        "coverage_md": reports_dir / DEFAULT_COVERAGE_MD,
        "coverage_compare_json": reports_dir / DEFAULT_COVERAGE_COMPARE_JSON,
        "coverage_compare_md": reports_dir / DEFAULT_COVERAGE_COMPARE_MD,
        "coverage_evidence_json": reports_dir / DEFAULT_COVERAGE_EVIDENCE_JSON,
        "coverage_evidence_md": reports_dir / DEFAULT_COVERAGE_EVIDENCE_MD,
        "bundles_json": reports_dir / DEFAULT_BUNDLES_JSON,
        "bundles_md": reports_dir / DEFAULT_BUNDLES_MD,
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
    if not paths["coverage_json"].exists():
        write_json(
            paths["coverage_json"],
            {
                "status": "bootstrap_ready",
                "mode": "coverage",
                "repo_root": str(repo_root),
            },
        )
    if not paths["coverage_md"].exists():
        write_text(
            paths["coverage_md"],
            "# Refactor Fidelity Coverage\n\nBootstrap complete. No coverage run recorded yet.\n",
        )
    if not paths["coverage_compare_json"].exists():
        write_json(
            paths["coverage_compare_json"],
            {
                "status": "bootstrap_ready",
                "mode": "coverage_compare",
                "repo_root": str(repo_root),
            },
        )
    if not paths["coverage_compare_md"].exists():
        write_text(
            paths["coverage_compare_md"],
            "# Refactor Fidelity Coverage Compare\n\nBootstrap complete. No comparative coverage run recorded yet.\n",
        )
    if not paths["coverage_evidence_json"].exists():
        write_json(
            paths["coverage_evidence_json"],
            {
                "status": "bootstrap_ready",
                "mode": "coverage_evidence",
                "repo_root": str(repo_root),
            },
        )
    if not paths["coverage_evidence_md"].exists():
        write_text(
            paths["coverage_evidence_md"],
            "# Refactor Fidelity Coverage Evidence\n\nBootstrap complete. No coverage evidence run recorded yet.\n",
        )
    if not paths["bundles_json"].exists():
        write_json(
            paths["bundles_json"],
            {
                "status": "bootstrap_ready",
                "mode": "bundle_map",
                "repo_root": str(repo_root),
            },
        )
    if not paths["bundles_md"].exists():
        write_text(
            paths["bundles_md"],
            "# Refactor Fidelity Bundles\n\nBootstrap complete. No bundle map recorded yet.\n",
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


def render_coverage_summary_markdown(summary: dict) -> str:
    lines = [
        "# Refactor Fidelity Coverage",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Side: `{summary['side']}`",
        f"- Source: `{summary['source']}`",
        f"- Status: `{summary['status']}`",
        f"- Tool status: `{summary['tool_status']}`",
        f"- Test files: `{summary['test_file_count']}`",
        f"- Executed tests: `{summary.get('test_count', 0)}`",
        f"- Failing tests: `{summary.get('failure_count', 0)}`",
        f"- Skipped tests: `{summary.get('skip_count', 0)}`",
        f"- Bundles: `{summary['bundle_count']}`",
        f"- Covered files: `{summary['file_count']}`",
        f"- Lines covered: `{summary['lines_covered']}` / `{summary['lines_total']}`",
        f"- Line coverage: `{summary['line_percent_display']}`",
    ]
    if summary.get("bundle_manifest_path"):
        lines.append(f"- Bundle manifest: `{summary['bundle_manifest_path']}`")
    status_counts = summary.get("bundle_status_counts", {})
    if status_counts:
        lines.extend(["", "## Bundle Status Counts", ""])
        for status_name in sorted(status_counts):
            lines.append(f"- `{status_name}`: `{status_counts[status_name]}`")
    if summary.get("sample_bundle_results"):
        lines.extend(["", "## Sample Bundle Results", ""])
        for item in summary["sample_bundle_results"]:
            lines.append(
                f"- `{item['status']}` `{item['bundle_id']}` coverage `{item['line_percent_display']}`"
            )
    if summary.get("sample_files"):
        lines.extend(["", "## Sample File Coverage", ""])
        for item in summary["sample_files"]:
            lines.append(
                f"- `{item['file_path']}` `{item['lines_covered']}/{item['lines_total']}` `{item['line_percent_display']}`"
            )
    if summary.get("notes"):
        lines.extend(["", "## Notes", ""])
        for note in summary["notes"]:
            lines.append(f"- {note}")
    return "\n".join(lines) + "\n"


def render_coverage_compare_summary_markdown(summary: dict) -> str:
    lines = [
        "# Refactor Fidelity Coverage Compare",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Test source: `{summary['test_source']}`",
        f"- Status: `{summary['status']}`",
        f"- Test files: `{summary['test_file_count']}`",
        f"- Bundles: `{summary['bundle_count']}`",
    ]
    if summary.get("bundle_manifest_path"):
        lines.append(f"- Bundle manifest: `{summary['bundle_manifest_path']}`")
    for title, side_key in (("Baseline", "baseline"), ("Target", "target")):
        side = summary[side_key]
        lines.extend(
            [
                "",
                f"## {title}",
                "",
                f"- Source: `{side['source']}`",
                f"- Status: `{side['status']}`",
                f"- Tool status: `{side['tool_status']}`",
                f"- Executed tests: `{side.get('test_count', 0)}`",
                f"- Failing tests: `{side.get('failure_count', 0)}`",
                f"- Skipped tests: `{side.get('skip_count', 0)}`",
                f"- Covered files: `{side['file_count']}`",
                f"- Lines covered: `{side['lines_covered']}` / `{side['lines_total']}`",
                f"- Line coverage: `{side['line_percent_display']}`",
            ]
        )
        if side.get("failed_tests"):
            lines.append(f"- Failed tests: `{', '.join(side['failed_tests'][:10])}`")
    for title, key in (("Bundle Status Counts", "bundle_status_counts"), ("Delta Counts", "delta_counts")):
        values = summary.get(key, {})
        if not values:
            continue
        lines.extend(["", f"## {title}", ""])
        for name in sorted(values):
            lines.append(f"- `{name}`: `{values[name]}`")
    if summary.get("sample_bundle_results"):
        lines.extend(["", "## Sample Bundle Results", ""])
        for item in summary["sample_bundle_results"]:
            lines.append(
                f"- `{item['status']}` `{item['bundle_id']}` baseline `{format_percent(item.get('baseline_line_percent'))}` target `{format_percent(item.get('target_line_percent'))}` delta `{format_percent(item.get('line_percent_delta'))}` `{item['delta_category']}`"
            )
    if summary.get("notes"):
        lines.extend(["", "## Notes", ""])
        for note in summary["notes"]:
            lines.append(f"- {note}")
    return "\n".join(lines) + "\n"


def render_coverage_evidence_markdown(summary: dict) -> str:
    package = summary["package_coverage"]
    lines = [
        "# Refactor Fidelity Coverage Evidence",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Status: `{summary['status']}`",
        f"- Bundle manifest: `{summary['bundle_manifest_path']}`",
        f"- Coverage compare run: `{summary['coverage_compare_run_id']}`",
        f"- Bundles: `{summary['bundle_count']}`",
        f"- Justified by tests: `{summary['justified_bundle_count']}`",
        f"- Low-coverage bundles: `{summary['low_coverage_bundle_count']}`",
        f"- Comparison-gap bundles: `{summary.get('comparison_gap_bundle_count', 0)}`",
        f"- Regressed bundles: `{summary['regression_bundle_count']}`",
        f"- Exceptions emitted: `{summary['exception_count']}`",
        "",
        "## Package Coverage",
        "",
        f"- Baseline: `{package['baseline_line_percent_display']}`",
        f"- Target: `{package['target_line_percent_display']}`",
        f"- Delta: `{package['line_percent_delta_display']}`",
        f"- Target threshold: `{package['threshold_display']}`",
        f"- Package gate: `{package['gate_status']}`",
    ]
    for title, key in (
        ("Gate Counts", "gate_counts"),
        ("Delta Counts", "delta_counts"),
        ("Threshold Counts", "threshold_counts"),
    ):
        values = summary.get(key, {})
        if not values:
            continue
        lines.extend(["", f"## {title}", ""])
        for name in sorted(values):
            lines.append(f"- `{name}`: `{values[name]}`")
    if summary.get("below_target_bundles"):
        lines.extend(["", "## Below Target Bundles", ""])
        for item in summary["below_target_bundles"][:25]:
            lines.append(
                f"- `{item['bundle_id']}` `{item['bundle_type']}` target `{item['target_line_percent_display']}` threshold `{item['threshold_display']}` tests `{item['shared_test_file_count']}`"
            )
    if summary.get("comparison_gap_bundles"):
        lines.extend(["", "## Comparison Gap Bundles", ""])
        for item in summary["comparison_gap_bundles"][:25]:
            lines.append(
                f"- `{item['bundle_id']}` `{item['bundle_type']}` baseline `{item['baseline_line_percent_display']}` target `{item['target_line_percent_display']}` threshold `{item['threshold_display']}` tests `{item['shared_test_file_count']}`"
            )
    if summary.get("regressed_bundles"):
        lines.extend(["", "## Regressed Bundles", ""])
        for item in summary["regressed_bundles"][:25]:
            lines.append(
                f"- `{item['bundle_id']}` baseline `{item['baseline_line_percent_display']}` target `{item['target_line_percent_display']}` delta `{item['line_percent_delta_display']}`"
            )
    if summary.get("sample_bundles"):
        lines.extend(["", "## Sample Bundle Evidence", ""])
        for item in summary["sample_bundles"][:25]:
            lines.append(
                f"- `{item['coverage_gate_status']}` `{item['bundle_id']}` target `{item['target_line_percent_display']}` threshold `{item['threshold_display']}` tests `{item['shared_test_file_count']}`"
            )
    if summary.get("notes"):
        lines.extend(["", "## Notes", ""])
        for note in summary["notes"]:
            lines.append(f"- {note}")
    return "\n".join(lines) + "\n"


def render_bundle_map_summary_markdown(summary: dict) -> str:
    lines = [
        "# Refactor Fidelity Bundles",
        "",
        f"Run: `{summary['run_id']}`",
        "",
        f"- Mode: `{summary['mode']}`",
        f"- Repo root: `{summary['repo_root']}`",
        f"- Closeout run: `{summary['closeout_run_id']}`",
        f"- Baseline: `{summary['baseline']}`",
        f"- Target: `{summary['target']}`",
        f"- Status: `{summary['status']}`",
        f"- Bundles: `{summary['bundle_count']}`",
        f"- Linked exceptions: `{summary['linked_exception_count']}`",
        f"- Priority families: `{summary['priority_family_count']}`",
        f"- Manifest path: `{summary['bundle_manifest_path']}`",
    ]
    for title, key in (
        ("Bundle Type Counts", "bundle_type_counts"),
        ("Disposition Counts", "disposition_counts"),
        ("Scenario Class Counts", "scenario_class_counts"),
        ("Evidence Depth Counts", "evidence_depth_counts"),
    ):
        values = summary.get(key, {})
        if not values:
            continue
        lines.extend(["", f"## {title}", ""])
        for name in sorted(values):
            lines.append(f"- `{name}`: `{values[name]}`")
    if summary.get("priority_bundles"):
        lines.extend(["", "## Priority Bundles", ""])
        for item in summary["priority_bundles"][:25]:
            lines.append(
                f"- `#{item['priority_rank']}` `{item['bundle_id']}` `{item['bundle_type']}` exceptions `{item['linked_exception_count']}` tests `{item['shared_test_file_count']}`"
            )
    if summary.get("priority_families"):
        lines.extend(["", "## Priority Families", ""])
        for item in summary["priority_families"][:25]:
            lines.append(
                f"- `{item['family_path']}` bundles `{item['bundle_count']}` exceptions `{item['linked_exception_count']}`"
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
        f"- Catalog-curated exceptions: `{summary.get('catalog_curated_count', 0)}`",
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
        raise RuntimeError("audit helper emitted no JSON")

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

    raise RuntimeError("audit helper emitted non-JSON output")


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


def list_curation_paths(repo_root: pathlib.Path) -> list[pathlib.Path]:
    return sorted((repo_root / "tools" / "refactor").glob(DEFAULT_CURATION_GLOB))


def stable_json(value) -> str:
    return json.dumps(value, sort_keys=True)


def make_exception_key(audit_layer: str, exception_type: str, entity_key: str) -> str:
    return f"{audit_layer}::{exception_type}::{entity_key}"


def exception_auto_rule(exception_type: str) -> dict | None:
    return AUTO_CURATION_RULES.get(exception_type)


def manifest_entity_key(manifest_path: str, entry_id: str) -> str:
    return f"manifest::{manifest_path}::{entry_id}"


def parse_manifest_exception_entity_key(entity_key: str) -> tuple[str, str] | None:
    if not entity_key.startswith("manifest::"):
        return None
    payload = entity_key[len("manifest::") :]
    if "::" not in payload:
        return None
    manifest_path, entry_id = payload.rsplit("::", 1)
    if not manifest_path or not entry_id:
        return None
    return manifest_path, entry_id


def suggested_disposition_for_exception(exception_type: str) -> str | None:
    rule = exception_auto_rule(exception_type)
    if rule:
        return str(rule["disposition"])

    suggestions = {
        "manual_merge_expected": "expected_manual_merge",
        "missing_target_block": "needs_manifest_lineage_review",
        "selector_ambiguous": "unsupported_audit_shape",
        "semantic_mismatch": "fix_required",
        "coverage_target_below_threshold": "add_tests_or_expand_coverage",
        "coverage_target_unmeasured": "fix_bundle_mapping_or_add_tests",
        "coverage_comparison_gap": "add_shared_compat_or_characterization_tests",
        "package_coverage_below_threshold": "expand_package_coverage",
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


def dedupe_paths(paths: list[str]) -> list[str]:
    seen: set[str] = set()
    unique: list[str] = []
    for path in paths:
        if not path or path in seen:
            continue
        seen.add(path)
        unique.append(path)
    return unique


def discover_project_library_paths(project_root: pathlib.Path) -> list[str]:
    library_root = project_root / "renv" / "library"
    if not library_root.exists():
        return []
    return dedupe_paths(
        [
            str(candidate)
            for candidate in sorted(library_root.glob("*/*/*"))
            if candidate.is_dir()
        ]
    )


def merge_env_paths(preferred: list[str], *existing_values: str | None) -> str:
    merged = list(preferred)
    for value in existing_values:
        if not value:
            continue
        merged.extend(segment for segment in value.split(os.pathsep) if segment)
    return os.pathsep.join(dedupe_paths(merged))


def run_contract_test_file(
    target_root: pathlib.Path,
    file_path: pathlib.Path,
    *,
    load_package: bool,
    library_paths: list[str] | None = None,
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
    env = os.environ.copy()
    if library_paths:
        merged_paths = merge_env_paths(
            library_paths,
            env.get("R_LIBS"),
            env.get("R_LIBS_USER"),
            env.get("R_LIBS_SITE"),
        )
        env["R_LIBS"] = merged_paths
        env["R_LIBS_USER"] = merged_paths
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        env=env,
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


def shared_collate_sequence(snapshot: dict, allowed_files: set[str]) -> list[str]:
    return [
        entry["file_path"]
        for entry in ensure_list(snapshot.get("collate"))
        if entry["file_path"] in allowed_files
    ]


def entity_fingerprint(entry: dict) -> str | None:
    return entry.get("hash_ast") or entry.get("hash_normalized") or entry.get("hash_raw")


def duplicate_shape_semantically_preserved(
    baseline_entries: list[dict],
    target_entries: list[dict],
) -> bool:
    if not baseline_entries or not target_entries:
        return False
    if len(target_entries) > len(baseline_entries):
        return False

    baseline_fingerprints = {
        fingerprint
        for entry in baseline_entries
        if (fingerprint := entity_fingerprint(entry))
    }
    target_fingerprints = {
        fingerprint
        for entry in target_entries
        if (fingerprint := entity_fingerprint(entry))
    }

    if not baseline_fingerprints or not target_fingerprints:
        return False

    return baseline_fingerprints == target_fingerprints


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
        baseline_files = sorted({entry["file_path"] for entry in baseline_entries})
        target_files = sorted({entry["file_path"] for entry in target_entries})

        duplicate_shape_changed = len(baseline_entries) != len(target_entries)
        if (len(baseline_entries) > 1 or len(target_entries) > 1) and duplicate_shape_changed:
            if not duplicate_shape_semantically_preserved(baseline_entries, target_entries):
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
                            "baseline_files": baseline_files,
                            "target_files": target_files,
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
                        "target_files": target_files,
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
                        "baseline_files": baseline_files,
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

    shared_files = set(baseline_files) & set(target_files)
    baseline_sequence = shared_collate_sequence(baseline_snapshot, shared_files)
    target_sequence = shared_collate_sequence(target_snapshot, shared_files)
    if baseline_sequence != target_sequence:
        baseline_positions = {file_path: index for index, file_path in enumerate(baseline_sequence, start=1)}
        target_positions = {file_path: index for index, file_path in enumerate(target_sequence, start=1)}
        for file_path in sorted(shared_files):
            if baseline_positions[file_path] == target_positions[file_path]:
                continue
            diffs.append(
                make_diff(
                    entity_kind="r_file",
                    entity_key=file_path,
                    diff_class="collate_drift",
                    baseline_id=baseline_refs["files"].get(file_path),
                    target_id=target_refs["files"].get(file_path),
                    evidence={
                        "baseline_collate_index": baseline_files[file_path],
                        "target_collate_index": target_files[file_path],
                        "baseline_shared_index": baseline_positions[file_path],
                        "target_shared_index": target_positions[file_path],
                        "baseline_shared_sequence": baseline_sequence,
                        "target_shared_sequence": target_sequence,
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
        "source_lineage_gap": "high",
        "source_lineage_gap_target_resolved": "low",
        "source_selector_ambiguous": "high",
        "source_resolution_failure": "high",
        "manual_merge_expected": "medium",
        "ast_only_drift": "medium",
        "target_resolution_asymmetry": "low",
        "whitespace_only_drift": "low",
        "comment_only_drift": "low",
    }.get(exception_type, "medium")


def manifest_exception_reason(comparison: dict) -> str:
    status = comparison["status"]
    baseline_resolver = comparison.get("baseline_source_resolver") or "unresolved"
    target_resolver = comparison.get("target_resolver") or "unresolved"
    note = comparison.get("note")
    if status == "source_missing_target_resolved":
        return (
            f"Baseline/source block could not be resolved via `{baseline_resolver}`, "
            f"but the target block resolved via `{target_resolver}`."
            if not note
            else (
                f"Baseline/source block could not be resolved via `{baseline_resolver}`, "
                f"but the target block resolved via `{target_resolver}`: {note}"
            )
        )
    if status == "source_missing":
        return (
            f"Baseline/source block could not be resolved via `{baseline_resolver}`."
            if not note
            else f"Baseline/source block could not be resolved via `{baseline_resolver}`: {note}"
        )
    if status == "source_ambiguous":
        return (
            f"Baseline/source selector was ambiguous via `{baseline_resolver}`."
            if not note
            else f"Baseline/source selector was ambiguous via `{baseline_resolver}`: {note}"
        )
    if status == "source_unresolved":
        return (
            f"Baseline/source resolution failed via `{baseline_resolver}`."
            if not note
            else f"Baseline/source resolution failed via `{baseline_resolver}`: {note}"
        )
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
                "baseline_source_resolver": comparison.get("baseline_source_resolver"),
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


def format_percent(value: float | int | None) -> str:
    if value is None:
        return "n/a"
    return f"{float(value):.1f}%"


def coerce_int(value) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def coerce_float(value) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def normalize_int_list(values) -> list[int]:
    normalized: list[int] = []
    for value in ensure_list(values):
        coerced = coerce_int(value)
        if coerced is None:
            continue
        normalized.append(coerced)
    return sorted(set(normalized))


def resolve_repo_path(repo_root: pathlib.Path, value: str) -> pathlib.Path:
    candidate = pathlib.Path(value)
    if not candidate.is_absolute():
        candidate = repo_root / candidate
    return candidate.resolve()


def normalize_coverage_bundle_members(
    bundle_id: str,
    raw_members,
    *,
    source_file_path: str | None,
) -> list[dict]:
    members: list[dict] = []
    for member_index, raw_member in enumerate(ensure_list(raw_members), start=1):
        if not isinstance(raw_member, dict):
            continue
        file_path = raw_member.get("file_path") or raw_member.get("source_file") or source_file_path
        if not file_path:
            continue
        members.append(
            {
                "member_id": str(raw_member.get("member_id") or f"{bundle_id}:member:{member_index}"),
                "entity_key": raw_member.get("entity_key"),
                "file_path": str(file_path),
                "line_start": coerce_int(raw_member.get("line_start")),
                "line_end": coerce_int(raw_member.get("line_end")),
                "metadata": {
                    key: value
                    for key, value in raw_member.items()
                    if key not in {"member_id", "entity_key", "file_path", "source_file", "line_start", "line_end"}
                },
            }
        )
    return members


def load_coverage_bundle_manifest(
    repo_root: pathlib.Path,
    manifest_path: str | None,
) -> tuple[str | None, list[dict]]:
    if not manifest_path:
        manifest_path = str(pathlib.Path(".refactor-fidelity-audit") / "reports" / DEFAULT_BUNDLES_JSON)

    bundle_path = resolve_repo_path(repo_root, manifest_path)
    payload = json.loads(bundle_path.read_text(encoding="utf-8"))
    raw_bundles = ensure_list(payload.get("bundles"))
    bundles: list[dict] = []
    for index, raw_bundle in enumerate(raw_bundles, start=1):
        if not isinstance(raw_bundle, dict):
            continue
        bundle_id = str(raw_bundle.get("bundle_id") or raw_bundle.get("id") or f"bundle_{index}")
        bundle_type = str(raw_bundle.get("bundle_type") or "unspecified")
        source_file_path = raw_bundle.get("source_file")
        members_payload = ensure_list(raw_bundle.get("members"))
        if not members_payload and source_file_path:
            members_payload = [{"file_path": source_file_path}]
        members = normalize_coverage_bundle_members(
            bundle_id,
            members_payload,
            source_file_path=str(source_file_path) if source_file_path else None,
        )
        baseline_members = normalize_coverage_bundle_members(
            f"{bundle_id}:baseline",
            raw_bundle.get("baseline_members"),
            source_file_path=None,
        )
        target_members = normalize_coverage_bundle_members(
            f"{bundle_id}:target",
            raw_bundle.get("target_members"),
            source_file_path=None,
        )
        shared_test_files = sorted(
            {
                str(value)
                for value in ensure_list(raw_bundle.get("shared_test_files"))
                if value
            }
        )
        test_files = sorted(
            {
                str(value)
                for value in ensure_list(raw_bundle.get("test_files"))
                if value
            }
        )

        bundles.append(
            {
                "bundle_id": bundle_id,
                "bundle_type": bundle_type,
                "source_file_path": str(source_file_path) if source_file_path else None,
                "primary_entity_key": raw_bundle.get("primary_entity_key"),
                "baseline_ref": raw_bundle.get("baseline_ref"),
                "target_ref": raw_bundle.get("target_ref"),
                "baseline_paths": [str(value) for value in ensure_list(raw_bundle.get("baseline_paths")) if value],
                "target_paths": [str(value) for value in ensure_list(raw_bundle.get("target_paths")) if value],
                "linked_exception_count": coerce_int(raw_bundle.get("linked_exception_count")) or 0,
                "scenario_class": raw_bundle.get("scenario_class"),
                "differential_replay_status": raw_bundle.get("differential_replay_status"),
                "coverage_gate_status": raw_bundle.get("coverage_gate_status"),
                "review_note": raw_bundle.get("review_note"),
                "priority_rank": coerce_int(raw_bundle.get("priority_rank")),
                "priority_score": coerce_float(raw_bundle.get("priority_score")),
                "evidence_depth": raw_bundle.get("evidence_depth"),
                "test_files": test_files,
                "shared_test_files": shared_test_files,
                "members": members,
                "baseline_members": baseline_members,
                "target_members": target_members,
                "metadata": {
                    key: value
                    for key, value in raw_bundle.items()
                    if key not in {
                        "bundle_id",
                        "id",
                        "bundle_type",
                        "source_file",
                        "test_files",
                        "shared_test_files",
                        "members",
                        "baseline_members",
                        "target_members",
                    }
                },
            }
        )
    return relative_label(bundle_path, repo_root), bundles


def collect_bundle_test_files(bundles: list[dict]) -> list[str]:
    selected: list[str] = []
    for bundle in bundles:
        selected.extend(str(value) for value in ensure_list(bundle.get("test_files")) if value)
        selected.extend(str(value) for value in ensure_list(bundle.get("shared_test_files")) if value)
    return dedupe_paths(selected)


def discover_shared_compare_test_files(test_source_root: pathlib.Path) -> list[str]:
    testthat_root = test_source_root / "tests" / "testthat"
    if not testthat_root.exists():
        return []

    selected: list[str] = []
    for test_path in sorted(testthat_root.glob("test-*.R")):
        try:
            file_text = test_path.read_text(encoding="utf-8")
        except OSError:
            continue
        if "fidelity-coverage-compare: shared" not in file_text.lower():
            continue
        selected.append(relative_label(test_path, test_source_root))
    return selected


def coverage_bundle_members_for_side(bundle: dict, side: str) -> list[dict]:
    side_members = ensure_list(bundle.get(f"{side}_members"))
    if side_members:
        return side_members
    return ensure_list(bundle.get("members"))


def prepare_coverage_bundles_for_side(bundles: list[dict], side: str) -> list[dict]:
    prepared: list[dict] = []
    for bundle in bundles:
        prepared.append(
            {
                **bundle,
                "members": coverage_bundle_members_for_side(bundle, side),
                "member_source": (
                    f"{side}_members"
                    if ensure_list(bundle.get(f"{side}_members"))
                    else "members"
                ),
            }
        )
    return prepared


def normalize_coverage_file_record(record: dict, project_root: pathlib.Path) -> dict | None:
    file_path = record.get("file_path") or record.get("filename")
    if not file_path:
        return None

    candidate = pathlib.Path(str(file_path))
    normalized_path = (
        relative_label(candidate, project_root)
        if candidate.is_absolute()
        else str(candidate).replace("\\", "/")
    )
    covered_lines = normalize_int_list(record.get("covered_lines"))
    uncovered_lines = normalize_int_list(record.get("uncovered_lines"))
    observed_lines = sorted(set(covered_lines) | set(uncovered_lines))

    lines_total = coerce_int(record.get("lines_total"))
    if lines_total is None:
        lines_total = len(observed_lines)
    lines_covered = coerce_int(record.get("lines_covered"))
    if lines_covered is None:
        lines_covered = len(covered_lines)
    line_percent = coerce_float(record.get("line_percent"))
    if line_percent is None and lines_total:
        line_percent = (lines_covered / lines_total) * 100.0

    return {
        "file_path": normalized_path,
        "line_percent": line_percent,
        "lines_total": lines_total or 0,
        "lines_covered": lines_covered or 0,
        "covered_lines": covered_lines,
        "uncovered_lines": uncovered_lines,
        "result_payload": {key: value for key, value in record.items()},
    }


def run_coverage_capture(
    project_root: pathlib.Path,
    *,
    test_files: list[pathlib.Path],
    library_paths: list[str] | None = None,
    env_overrides: dict[str, str] | None = None,
) -> dict:
    script_path = SCRIPT_DIR / "fidelity_coverage_capture.R"
    command = [
        "Rscript",
        "--vanilla",
        str(script_path),
        "--project-root",
        str(project_root),
    ]
    for test_file in test_files:
        command.extend(["--test-file", str(test_file)])

    env = os.environ.copy()
    if library_paths:
        merged_paths = merge_env_paths(
            library_paths,
            env.get("R_LIBS"),
            env.get("R_LIBS_USER"),
            env.get("R_LIBS_SITE"),
        )
        env["R_LIBS"] = merged_paths
        env["R_LIBS_USER"] = merged_paths
    if env_overrides:
        env.update({key: value for key, value in env_overrides.items() if value is not None})

    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    if completed.returncode != 0:
        message = completed.stderr.strip() or completed.stdout.strip() or "coverage runner failed"
        return {
            "status": "runner_error",
            "tool_status": "runner_error",
            "project_root": str(project_root),
            "test_files": [relative_label(test_file, project_root) for test_file in test_files],
            "files": [],
            "line_percent": None,
            "lines_total": 0,
            "lines_covered": 0,
            "notes": [message],
            "error": {"message": message, "class": ["coverage_runner_error"]},
        }

    try:
        payload = load_json_payload(completed.stdout)
    except RuntimeError as exc:
        return {
            "status": "runner_error",
            "tool_status": "runner_output_error",
            "project_root": str(project_root),
            "test_files": [relative_label(test_file, project_root) for test_file in test_files],
            "files": [],
            "line_percent": None,
            "lines_total": 0,
            "lines_covered": 0,
            "notes": [str(exc)],
            "error": {"message": str(exc), "class": ["coverage_runner_output_error"]},
        }

    normalized_files = [
        item
        for item in (
            normalize_coverage_file_record(record, project_root)
            for record in ensure_list(payload.get("files"))
            if isinstance(record, dict)
        )
        if item is not None
    ]
    payload["files"] = normalized_files
    payload["status"] = str(payload.get("status") or "completed")
    payload["tool_status"] = str(
        payload.get("tool_status")
        or ("ok" if payload["status"] == "completed" else payload["status"])
    )
    payload["project_root"] = str(project_root)
    payload["test_files"] = [
        relative_label(pathlib.Path(str(value)), project_root)
        if pathlib.Path(str(value)).is_absolute()
        else str(value).replace("\\", "/")
        for value in ensure_list(payload.get("test_files"))
        if value
    ]
    payload["notes"] = [str(note) for note in ensure_list(payload.get("notes")) if str(note).strip()]
    payload["test_count"] = coerce_int(payload.get("test_count"))
    payload["failure_count"] = coerce_int(payload.get("failure_count"))
    payload["skip_count"] = coerce_int(payload.get("skip_count"))
    payload["failed_tests"] = [str(value) for value in ensure_list(payload.get("failed_tests")) if value]
    payload["lines_total"] = coerce_int(payload.get("lines_total"))
    payload["lines_covered"] = coerce_int(payload.get("lines_covered"))
    payload["line_percent"] = coerce_float(payload.get("line_percent"))
    if payload["lines_total"] is None:
        payload["lines_total"] = sum(record["lines_total"] for record in normalized_files)
    if payload["lines_covered"] is None:
        payload["lines_covered"] = sum(record["lines_covered"] for record in normalized_files)
    if payload["line_percent"] is None and payload["lines_total"]:
        payload["line_percent"] = (payload["lines_covered"] / payload["lines_total"]) * 100.0
    return payload


def materialize_coverage_workspace(
    source_root: pathlib.Path,
    *,
    test_source_root: pathlib.Path,
) -> tempfile.TemporaryDirectory:
    tmpdir = tempfile.TemporaryDirectory(prefix="fidelity-coverage-")
    workspace_root = pathlib.Path(tmpdir.name)
    copy_names = ("DESCRIPTION", "NAMESPACE", "R", "inst", "data", "tests")
    for name in copy_names:
        source_path = source_root / name
        target_path = workspace_root / name
        if not source_path.exists():
            continue
        if source_path.is_dir():
            shutil.copytree(source_path, target_path, dirs_exist_ok=True)
        else:
            target_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_path, target_path)

    if test_source_root != source_root:
        tests_source = test_source_root / "tests"
        if tests_source.exists():
            shutil.copytree(tests_source, workspace_root / "tests", dirs_exist_ok=True)

    return tmpdir


def coverage_lines_for_member(file_record: dict, member: dict) -> tuple[set[int], set[int], bool]:
    covered_lines = set(file_record.get("covered_lines", []))
    uncovered_lines = set(file_record.get("uncovered_lines", []))
    observed_lines = covered_lines | uncovered_lines
    line_start = member.get("line_start")
    line_end = member.get("line_end")

    if line_start is not None or line_end is not None:
        if not observed_lines:
            return set(), set(), True
        start = line_start if line_start is not None else min(observed_lines)
        end = line_end if line_end is not None else max(observed_lines)
        total_lines = {line for line in observed_lines if start <= line <= end}
        return total_lines, covered_lines & total_lines, False

    if observed_lines:
        return observed_lines, covered_lines, False

    lines_total = int(file_record.get("lines_total") or 0)
    lines_covered = int(file_record.get("lines_covered") or 0)
    if lines_total <= 0:
        return set(), set(), False

    synthetic_total = set(range(1, lines_total + 1))
    synthetic_covered = set(range(1, min(lines_covered, lines_total) + 1))
    return synthetic_total, synthetic_covered, False


def build_coverage_bundle_results(bundles: list[dict], coverage_files: list[dict]) -> list[dict]:
    files_by_path = {record["file_path"]: record for record in coverage_files}
    results: list[dict] = []
    for bundle in bundles:
        total_pairs: set[tuple[str, int]] = set()
        covered_pairs: set[tuple[str, int]] = set()
        missing_files: set[str] = set()
        line_detail_missing = False

        for member in bundle.get("members", []):
            file_path = member["file_path"]
            file_record = files_by_path.get(file_path)
            if file_record is None:
                missing_files.add(file_path)
                continue
            member_total, member_covered, member_line_gap = coverage_lines_for_member(file_record, member)
            line_detail_missing = line_detail_missing or member_line_gap
            total_pairs.update((file_path, line) for line in member_total)
            covered_pairs.update((file_path, line) for line in member_covered)

        if total_pairs:
            status = "completed"
        elif line_detail_missing:
            status = "line_detail_missing"
        elif missing_files:
            status = "not_observed"
        elif not bundle.get("members"):
            status = "no_members"
        else:
            status = "not_observed"

        lines_total = len(total_pairs)
        lines_covered = len(covered_pairs)
        line_percent = ((lines_covered / lines_total) * 100.0) if lines_total else None
        results.append(
            {
                "bundle_id": bundle["bundle_id"],
                "status": status,
                "line_percent": line_percent,
                "lines_total": lines_total,
                "lines_covered": lines_covered,
                "line_percent_display": format_percent(line_percent),
                "result_payload": {
                    "missing_files": sorted(missing_files),
                    "member_count": len(bundle.get("members", [])),
                    "line_detail_missing": line_detail_missing,
                    "observed_file_count": len(
                        {file_path for file_path, _ in total_pairs}
                    ),
                },
            }
        )
    return results


def insert_coverage_results(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    side: str,
    scope: str,
    source_label: str,
    bundle_manifest_path: str | None,
    coverage_payload: dict,
    bundles: list[dict],
    bundle_results: list[dict],
    selected_test_files: list[str],
    command_test_files: list[str],
) -> None:
    coverage_run_id = f"{run_id}:coverage:{scope}"
    conn.execute(
        """
        INSERT INTO coverage_runs (
            coverage_run_id, run_id, side, source_label, bundle_manifest_path, status, tool_status,
            test_file_count, bundle_count, file_count, line_percent, lines_total, lines_covered, summary_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            coverage_run_id,
            run_id,
            side,
            source_label,
            bundle_manifest_path,
            coverage_payload["status"],
            coverage_payload["tool_status"],
            len(selected_test_files),
            len(bundles),
            len(coverage_payload["files"]),
            coverage_payload.get("line_percent"),
            coverage_payload.get("lines_total"),
            coverage_payload.get("lines_covered"),
            stable_json(coverage_payload),
        ),
    )

    for index, file_record in enumerate(coverage_payload["files"], start=1):
        conn.execute(
            """
            INSERT INTO coverage_files (
                coverage_file_id, run_id, coverage_run_id, file_path, line_percent, lines_total,
                lines_covered, covered_lines_json, uncovered_lines_json, result_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                f"{run_id}:coverage_file:{scope}:{index}",
                run_id,
                coverage_run_id,
                file_record["file_path"],
                file_record.get("line_percent"),
                file_record.get("lines_total"),
                file_record.get("lines_covered"),
                stable_json(file_record.get("covered_lines", [])),
                stable_json(file_record.get("uncovered_lines", [])),
                stable_json(file_record.get("result_payload", {})),
            ),
        )

    bundle_ids: dict[str, str] = {}
    for index, bundle in enumerate(bundles, start=1):
        coverage_bundle_id = f"{run_id}:coverage_bundle:{scope}:{index}"
        bundle_ids[bundle["bundle_id"]] = coverage_bundle_id
        conn.execute(
            """
            INSERT INTO coverage_bundles (
                coverage_bundle_id, run_id, bundle_id, bundle_type, source_file_path, side, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                coverage_bundle_id,
                run_id,
                bundle["bundle_id"],
                bundle["bundle_type"],
                bundle.get("source_file_path"),
                side,
                stable_json(
                    {
                        "test_files": bundle.get("test_files", []),
                        "shared_test_files": bundle.get("shared_test_files", []),
                        "member_source": bundle.get("member_source"),
                        "metadata": bundle.get("metadata", {}),
                    }
                ),
            ),
        )
        for member_index, member in enumerate(bundle.get("members", []), start=1):
            conn.execute(
                """
                INSERT INTO coverage_bundle_members (
                    coverage_bundle_member_id, run_id, coverage_bundle_id, entity_key,
                    file_path, line_start, line_end, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:coverage_bundle_member:{scope}:{index}:{member_index}",
                    run_id,
                    coverage_bundle_id,
                    member.get("entity_key"),
                    member["file_path"],
                    member.get("line_start"),
                    member.get("line_end"),
                    stable_json(member.get("metadata", {})),
                ),
            )

    link_index = 1
    for test_file_path in command_test_files:
        conn.execute(
            """
                INSERT INTO coverage_test_links (
                    coverage_test_link_id, run_id, coverage_bundle_id, test_file_path, source, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:coverage_test_link:{scope}:{link_index}",
                    run_id,
                    None,
                    test_file_path,
                "command",
                stable_json({}),
            ),
        )
        link_index += 1

    for bundle in bundles:
        coverage_bundle_id = bundle_ids[bundle["bundle_id"]]
        bundle_test_files = dedupe_paths(
            [str(value) for value in ensure_list(bundle.get("test_files")) if value] +
            [str(value) for value in ensure_list(bundle.get("shared_test_files")) if value]
        )
        for test_file_path in bundle_test_files:
            conn.execute(
                """
                INSERT INTO coverage_test_links (
                    coverage_test_link_id, run_id, coverage_bundle_id, test_file_path, source, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:coverage_test_link:{scope}:{link_index}",
                    run_id,
                    coverage_bundle_id,
                    test_file_path,
                    "bundle_manifest",
                    stable_json({}),
                ),
            )
            link_index += 1

    for index, result in enumerate(bundle_results, start=1):
        conn.execute(
            """
                INSERT INTO coverage_bundle_results (
                    coverage_bundle_result_id, run_id, coverage_bundle_id, status,
                    line_percent, lines_total, lines_covered, result_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                f"{run_id}:coverage_bundle_result:{scope}:{index}",
                run_id,
                bundle_ids[result["bundle_id"]],
                result["status"],
                result.get("line_percent"),
                result.get("lines_total"),
                result.get("lines_covered"),
                stable_json(result.get("result_payload", {})),
            ),
        )


def insert_coverage_compare_results(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    bundles: list[dict],
    compare_bundle_results: list[dict],
) -> None:
    bundle_ids: dict[str, str] = {}
    for index, bundle in enumerate(bundles, start=1):
        coverage_bundle_id = f"{run_id}:coverage_bundle:comparison:{index}"
        bundle_ids[bundle["bundle_id"]] = coverage_bundle_id
        conn.execute(
            """
            INSERT INTO coverage_bundles (
                coverage_bundle_id, run_id, bundle_id, bundle_type, source_file_path, side, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                coverage_bundle_id,
                run_id,
                bundle["bundle_id"],
                bundle["bundle_type"],
                bundle.get("source_file_path"),
                "comparison",
                stable_json(
                    {
                        "test_files": bundle.get("test_files", []),
                        "shared_test_files": bundle.get("shared_test_files", []),
                        "metadata": bundle.get("metadata", {}),
                        "baseline_members": bundle.get("baseline_members", []),
                        "target_members": bundle.get("target_members", []),
                    }
                ),
            ),
        )

    for index, result in enumerate(compare_bundle_results, start=1):
        conn.execute(
            """
            INSERT INTO coverage_bundle_results (
                coverage_bundle_result_id, run_id, coverage_bundle_id, status,
                line_percent, lines_total, lines_covered, result_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                f"{run_id}:coverage_bundle_result:comparison:{index}",
                run_id,
                bundle_ids[result["bundle_id"]],
                result["status"],
                result.get("target_line_percent"),
                coerce_int((result.get("result_payload") or {}).get("target", {}).get("lines_total")),
                coerce_int((result.get("result_payload") or {}).get("target", {}).get("lines_covered")),
                stable_json(result.get("result_payload", {})),
            ),
        )


def build_coverage_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    side: str,
    source_label: str,
    bundle_manifest_path: str | None,
    selected_test_files: list[str],
    coverage_payload: dict,
    bundle_results: list[dict],
) -> dict:
    bundle_status_counts = Counter(result["status"] for result in bundle_results)
    summary = {
        "run_id": run_id,
        "mode": "coverage",
        "repo_root": str(repo_root),
        "side": side,
        "source": source_label,
        "status": coverage_payload["status"],
        "tool_status": coverage_payload["tool_status"],
        "bundle_manifest_path": bundle_manifest_path,
        "test_file_count": len(selected_test_files),
        "bundle_count": len(bundle_results),
        "file_count": len(coverage_payload["files"]),
        "test_count": coverage_payload.get("test_count"),
        "failure_count": coverage_payload.get("failure_count"),
        "skip_count": coverage_payload.get("skip_count"),
        "failed_tests": coverage_payload.get("failed_tests", []),
        "line_percent": coverage_payload.get("line_percent"),
        "line_percent_display": format_percent(coverage_payload.get("line_percent")),
        "lines_total": int(coverage_payload.get("lines_total") or 0),
        "lines_covered": int(coverage_payload.get("lines_covered") or 0),
        "bundle_status_counts": dict(sorted(bundle_status_counts.items())),
        "sample_bundle_results": [
            {
                "bundle_id": result["bundle_id"],
                "status": result["status"],
                "line_percent": result.get("line_percent"),
                "line_percent_display": result["line_percent_display"],
                "lines_total": result.get("lines_total"),
                "lines_covered": result.get("lines_covered"),
            }
            for result in bundle_results[:25]
        ],
        "sample_files": [
            {
                "file_path": record["file_path"],
                "line_percent": record.get("line_percent"),
                "line_percent_display": format_percent(record.get("line_percent")),
                "lines_total": record.get("lines_total"),
                "lines_covered": record.get("lines_covered"),
            }
            for record in coverage_payload["files"][:25]
        ],
        "notes": coverage_payload.get("notes", []),
    }
    return summary


def categorize_coverage_delta(value: float | None) -> str:
    if value is None:
        return "unavailable"
    if value > 0.05:
        return "improved"
    if value < -0.05:
        return "regressed"
    return "flat"


def build_coverage_compare_bundle_results(
    bundles: list[dict],
    *,
    baseline_bundle_results: list[dict],
    target_bundle_results: list[dict],
    baseline_capture_status: str = "completed",
    target_capture_status: str = "completed",
) -> list[dict]:
    baseline_index = {result["bundle_id"]: result for result in baseline_bundle_results}
    target_index = {result["bundle_id"]: result for result in target_bundle_results}
    results: list[dict] = []
    for bundle in bundles:
        bundle_id = bundle["bundle_id"]
        baseline_result = baseline_index.get(bundle_id)
        target_result = target_index.get(bundle_id)
        baseline_status = str((baseline_result or {}).get("status") or "missing")
        target_status = str((target_result or {}).get("status") or "missing")
        if baseline_capture_status != "completed":
            baseline_status = baseline_capture_status
        if target_capture_status != "completed":
            target_status = target_capture_status
        baseline_line_percent = coerce_float((baseline_result or {}).get("line_percent"))
        target_line_percent = coerce_float((target_result or {}).get("line_percent"))
        line_percent_delta = None
        if baseline_line_percent is not None and target_line_percent is not None:
            line_percent_delta = target_line_percent - baseline_line_percent

        if baseline_status == "completed" and target_status == "completed":
            status = "completed"
        elif baseline_status == "tool_unavailable" or target_status == "tool_unavailable":
            status = "tool_unavailable"
        elif baseline_status != "completed" and target_status != "completed":
            status = "both_incomplete"
        elif baseline_status != "completed":
            status = "baseline_incomplete"
        else:
            status = "target_incomplete"

        results.append(
            {
                "bundle_id": bundle_id,
                "status": status,
                "line_percent": target_line_percent,
                "lines_total": coerce_int((target_result or {}).get("lines_total")),
                "lines_covered": coerce_int((target_result or {}).get("lines_covered")),
                "line_percent_display": format_percent(target_line_percent),
                "baseline_status": baseline_status,
                "target_status": target_status,
                "baseline_line_percent": baseline_line_percent,
                "target_line_percent": target_line_percent,
                "line_percent_delta": line_percent_delta,
                "delta_category": categorize_coverage_delta(line_percent_delta),
                "result_payload": {
                    "bundle_type": bundle.get("bundle_type"),
                    "baseline": {
                        "status": baseline_status,
                        "line_percent": baseline_line_percent,
                        "lines_total": coerce_int((baseline_result or {}).get("lines_total")),
                        "lines_covered": coerce_int((baseline_result or {}).get("lines_covered")),
                        "payload": (baseline_result or {}).get("result_payload", {}),
                    },
                    "target": {
                        "status": target_status,
                        "line_percent": target_line_percent,
                        "lines_total": coerce_int((target_result or {}).get("lines_total")),
                        "lines_covered": coerce_int((target_result or {}).get("lines_covered")),
                        "payload": (target_result or {}).get("result_payload", {}),
                    },
                    "line_percent_delta": line_percent_delta,
                    "delta_category": categorize_coverage_delta(line_percent_delta),
                    "shared_test_files": bundle.get("shared_test_files", []),
                    "test_files": bundle.get("test_files", []),
                    "baseline_member_count": len(bundle.get("baseline_members", [])),
                    "target_member_count": len(bundle.get("target_members", [])),
                },
            }
        )
    return results


def build_coverage_compare_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    baseline_source_label: str,
    target_source_label: str,
    test_source_label: str,
    bundle_manifest_path: str | None,
    selected_test_files: list[str],
    baseline_payload: dict,
    target_payload: dict,
    compare_bundle_results: list[dict],
) -> dict:
    status_counts = Counter(result["status"] for result in compare_bundle_results)
    delta_counts = Counter(result["delta_category"] for result in compare_bundle_results)
    summary_status = "completed"
    if baseline_payload["status"] != "completed" or target_payload["status"] != "completed":
        if baseline_payload["status"] == "tool_unavailable" or target_payload["status"] == "tool_unavailable":
            summary_status = "tool_unavailable"
        else:
            summary_status = "completed_with_incomplete_sides"
    elif any(result["status"] != "completed" for result in compare_bundle_results):
        summary_status = "completed_with_incomplete_bundles"

    summary = {
        "run_id": run_id,
        "mode": "coverage_compare",
        "repo_root": str(repo_root),
        "baseline": {
            "source": baseline_source_label,
            "status": baseline_payload["status"],
            "tool_status": baseline_payload["tool_status"],
            "test_count": baseline_payload.get("test_count"),
            "failure_count": baseline_payload.get("failure_count"),
            "skip_count": baseline_payload.get("skip_count"),
            "failed_tests": baseline_payload.get("failed_tests", []),
            "file_count": len(baseline_payload["files"]),
            "line_percent": baseline_payload.get("line_percent"),
            "line_percent_display": format_percent(baseline_payload.get("line_percent")),
            "lines_total": int(baseline_payload.get("lines_total") or 0),
            "lines_covered": int(baseline_payload.get("lines_covered") or 0),
        },
        "target": {
            "source": target_source_label,
            "status": target_payload["status"],
            "tool_status": target_payload["tool_status"],
            "test_count": target_payload.get("test_count"),
            "failure_count": target_payload.get("failure_count"),
            "skip_count": target_payload.get("skip_count"),
            "failed_tests": target_payload.get("failed_tests", []),
            "file_count": len(target_payload["files"]),
            "line_percent": target_payload.get("line_percent"),
            "line_percent_display": format_percent(target_payload.get("line_percent")),
            "lines_total": int(target_payload.get("lines_total") or 0),
            "lines_covered": int(target_payload.get("lines_covered") or 0),
        },
        "test_source": test_source_label,
        "status": summary_status,
        "bundle_manifest_path": bundle_manifest_path,
        "test_file_count": len(selected_test_files),
        "bundle_count": len(compare_bundle_results),
        "bundle_status_counts": dict(sorted(status_counts.items())),
        "delta_counts": dict(sorted(delta_counts.items())),
        "bundle_results": compare_bundle_results,
        "sample_bundle_results": [
            {
                "bundle_id": result["bundle_id"],
                "status": result["status"],
                "baseline_status": result["baseline_status"],
                "target_status": result["target_status"],
                "baseline_line_percent": result.get("baseline_line_percent"),
                "target_line_percent": result.get("target_line_percent"),
                "line_percent_delta": result.get("line_percent_delta"),
                "delta_category": result["delta_category"],
            }
            for result in compare_bundle_results[:25]
        ],
        "selected_test_files": selected_test_files,
        "notes": dedupe_paths(
            [str(note) for note in ensure_list(baseline_payload.get("notes")) if str(note).strip()] +
            [str(note) for note in ensure_list(target_payload.get("notes")) if str(note).strip()]
        ),
    }
    return summary


def select_coverage_bundles(
    bundles: list[dict],
    *,
    bundle_ids: list[str] | None = None,
    limit: int | None = None,
) -> list[dict]:
    selected = bundles
    if bundle_ids:
        allowed = set(bundle_ids)
        selected = [bundle for bundle in selected if bundle["bundle_id"] in allowed]
    if limit is not None and limit >= 0:
        selected = selected[:limit]
    return selected


def build_coverage_evidence_exception(
    run_id: str,
    bundle: dict,
) -> dict | None:
    gate_status = str(bundle.get("coverage_gate_status") or "")
    if gate_status == "below_target":
        exception_type = "coverage_target_below_threshold"
        severity = "medium"
        reason = (
            f"Target coverage for `{bundle['bundle_id']}` is `{format_percent(bundle.get('target_line_percent'))}`, "
            f"below the required `{format_percent(bundle.get('threshold_pct'))}`."
        )
    elif gate_status == "target_unmeasured":
        exception_type = "coverage_target_unmeasured"
        severity = "high"
        reason = (
            f"Target coverage for `{bundle['bundle_id']}` could not be measured from the selected shared test subset."
        )
    elif gate_status == "comparison_gap":
        exception_type = "coverage_comparison_gap"
        severity = "medium"
        reason = (
            f"Selected shared tests cover target `{bundle['bundle_id']}` at `{format_percent(bundle.get('target_line_percent'))}`, "
            f"but baseline coverage remains `{format_percent(bundle.get('baseline_line_percent'))}`. "
            "Add shared compat or characterization tests that exercise both sides."
        )
    else:
        return None

    entity_key = f"coverage_bundle::{bundle['bundle_id']}"
    return {
        "exception_id": f"{run_id}:coverage_exception:{bundle['bundle_id']}",
        "exception_key": make_exception_key("coverage_evidence", exception_type, entity_key),
        "run_id": run_id,
        "audit_layer": "coverage_evidence",
        "entity_key": entity_key,
        "severity": severity,
        "exception_type": exception_type,
        "reason": reason,
        "baseline_ref": bundle.get("baseline_ref"),
        "target_ref": bundle.get("target_ref"),
        "evidence_json": stable_json(
            {
                "bundle_id": bundle["bundle_id"],
                "bundle_type": bundle.get("bundle_type"),
                "target_line_percent": bundle.get("target_line_percent"),
                "baseline_line_percent": bundle.get("baseline_line_percent"),
                "line_percent_delta": bundle.get("line_percent_delta"),
                "threshold_pct": bundle.get("threshold_pct"),
                "shared_test_files": bundle.get("shared_test_files", []),
                "review_note": bundle.get("review_note"),
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


def build_package_coverage_exception(
    run_id: str,
    *,
    baseline_ref: str | None,
    target_ref: str | None,
    target_line_percent: float | None,
    threshold_pct: float,
    enforce_gate: bool = True,
) -> dict | None:
    if not enforce_gate:
        return None
    if target_line_percent is None or target_line_percent + 1e-9 >= threshold_pct:
        return None
    return {
        "exception_id": f"{run_id}:coverage_exception:package",
        "exception_key": make_exception_key("coverage_evidence", "package_coverage_below_threshold", "coverage_package::target"),
        "run_id": run_id,
        "audit_layer": "coverage_evidence",
        "entity_key": "coverage_package::target",
        "severity": "medium",
        "exception_type": "package_coverage_below_threshold",
        "reason": (
            f"Target package coverage is `{format_percent(target_line_percent)}`, below the required `{format_percent(threshold_pct)}`."
        ),
        "baseline_ref": baseline_ref,
        "target_ref": target_ref,
        "evidence_json": stable_json(
            {
                "target_line_percent": target_line_percent,
                "threshold_pct": threshold_pct,
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


def insert_coverage_evidence_results(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    bundles: list[dict],
    exceptions: list[dict],
) -> None:
    bundle_ids: dict[str, str] = {}
    for index, bundle in enumerate(bundles, start=1):
        coverage_bundle_id = f"{run_id}:coverage_bundle:evidence:{index}"
        bundle_ids[bundle["bundle_id"]] = coverage_bundle_id
        conn.execute(
            """
            INSERT INTO coverage_bundles (
                coverage_bundle_id, run_id, bundle_id, bundle_type, source_file_path,
                primary_entity_key, baseline_ref, target_ref, linked_exception_count,
                scenario_class, differential_replay_status, coverage_gate_status,
                review_note, priority_rank, priority_score, side, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                coverage_bundle_id,
                run_id,
                bundle["bundle_id"],
                bundle.get("bundle_type"),
                bundle.get("source_file_path"),
                bundle.get("primary_entity_key"),
                bundle.get("baseline_ref"),
                bundle.get("target_ref"),
                bundle.get("linked_exception_count"),
                bundle.get("scenario_class"),
                bundle.get("differential_replay_status"),
                bundle.get("coverage_gate_status"),
                bundle.get("review_note"),
                bundle.get("priority_rank"),
                bundle.get("priority_score"),
                "evidence",
                stable_json(
                    {
                        "baseline_paths": bundle.get("baseline_paths", []),
                        "target_paths": bundle.get("target_paths", []),
                        "shared_test_files": bundle.get("shared_test_files", []),
                        "threshold_pct": bundle.get("threshold_pct"),
                        "baseline_status": bundle.get("baseline_status"),
                        "target_status": bundle.get("target_status"),
                        "baseline_line_percent": bundle.get("baseline_line_percent"),
                        "target_line_percent": bundle.get("target_line_percent"),
                        "line_percent_delta": bundle.get("line_percent_delta"),
                        "delta_category": bundle.get("delta_category"),
                        "regressed_vs_baseline": bundle.get("regressed_vs_baseline"),
                        "evidence_depth": bundle.get("evidence_depth"),
                    }
                ),
            ),
        )
        for member_index, member in enumerate(ensure_list(bundle.get("baseline_members")), start=1):
            conn.execute(
                """
                INSERT INTO coverage_bundle_members (
                    coverage_bundle_member_id, run_id, coverage_bundle_id, entity_key,
                    file_path, line_start, line_end, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:coverage_bundle_member:evidence:baseline:{index}:{member_index}",
                    run_id,
                    coverage_bundle_id,
                    member.get("entity_key"),
                    member["file_path"],
                    member.get("line_start"),
                    member.get("line_end"),
                    stable_json({"side": "baseline"}),
                ),
            )
        for member_index, member in enumerate(ensure_list(bundle.get("target_members")), start=1):
            conn.execute(
                """
                INSERT INTO coverage_bundle_members (
                    coverage_bundle_member_id, run_id, coverage_bundle_id, entity_key,
                    file_path, line_start, line_end, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:coverage_bundle_member:evidence:target:{index}:{member_index}",
                    run_id,
                    coverage_bundle_id,
                    member.get("entity_key"),
                    member["file_path"],
                    member.get("line_start"),
                    member.get("line_end"),
                    stable_json({"side": "target"}),
                ),
            )
        for link_index, test_file_path in enumerate(bundle.get("shared_test_files", []), start=1):
            conn.execute(
                """
                INSERT INTO coverage_test_links (
                    coverage_test_link_id, run_id, coverage_bundle_id, test_file_path, source, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:coverage_test_link:evidence:{index}:{link_index}",
                    run_id,
                    coverage_bundle_id,
                    test_file_path,
                    "coverage_evidence",
                    stable_json({}),
                ),
            )
        conn.execute(
            """
            INSERT INTO coverage_bundle_results (
                coverage_bundle_result_id, run_id, coverage_bundle_id, status,
                line_percent, lines_total, lines_covered, result_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                f"{run_id}:coverage_bundle_result:evidence:{index}",
                run_id,
                coverage_bundle_id,
                bundle.get("coverage_gate_status"),
                bundle.get("target_line_percent"),
                coerce_int((bundle.get("bundle_result_payload") or {}).get("target", {}).get("lines_total")),
                coerce_int((bundle.get("bundle_result_payload") or {}).get("target", {}).get("lines_covered")),
                stable_json(
                    {
                        "threshold_pct": bundle.get("threshold_pct"),
                        "baseline_status": bundle.get("baseline_status"),
                        "target_status": bundle.get("target_status"),
                        "baseline_line_percent": bundle.get("baseline_line_percent"),
                        "target_line_percent": bundle.get("target_line_percent"),
                        "line_percent_delta": bundle.get("line_percent_delta"),
                        "delta_category": bundle.get("delta_category"),
                        "regressed_vs_baseline": bundle.get("regressed_vs_baseline"),
                        "justified_by_tests": bundle.get("justified_by_tests"),
                        "bundle_result_payload": bundle.get("bundle_result_payload", {}),
                    }
                ),
            ),
        )

    insert_exceptions(conn, exceptions)
    link_index = 1
    exceptions_by_entity_key = {record["entity_key"]: record for record in exceptions}
    for bundle in bundles:
        exception = exceptions_by_entity_key.get(f"coverage_bundle::{bundle['bundle_id']}")
        if exception is None:
            continue
        conn.execute(
            """
            INSERT INTO coverage_exception_links (
                coverage_exception_link_id, run_id, coverage_bundle_id, exception_id,
                exception_key, audit_layer, entity_key, curation_status, disposition, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                f"{run_id}:coverage_exception_link:evidence:{link_index}",
                run_id,
                bundle_ids[bundle["bundle_id"]],
                exception["exception_id"],
                exception.get("exception_key"),
                exception["audit_layer"],
                exception["entity_key"],
                exception["curation_status"],
                exception.get("disposition"),
                stable_json({"severity": exception["severity"], "exception_type": exception["exception_type"]}),
            ),
        )
        link_index += 1


def build_coverage_evidence_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    bundle_manifest_path: str | None,
    coverage_compare_summary: dict,
    acceptance_policy: dict,
    evaluated_bundles: list[dict],
    exceptions: list[dict],
    package_gate_applies: bool,
) -> dict:
    gate_counts = Counter(bundle["coverage_gate_status"] for bundle in evaluated_bundles)
    delta_counts = Counter(bundle["delta_category"] for bundle in evaluated_bundles)
    threshold_counts = Counter(
        format_percent(bundle["threshold_pct"]) for bundle in evaluated_bundles if bundle.get("threshold_pct") is not None
    )
    baseline_line_percent = coerce_float((coverage_compare_summary.get("baseline") or {}).get("line_percent"))
    target_line_percent = coerce_float((coverage_compare_summary.get("target") or {}).get("line_percent"))
    package_delta = None
    if baseline_line_percent is not None and target_line_percent is not None:
        package_delta = target_line_percent - baseline_line_percent
    package_threshold = float(acceptance_policy["package_line_coverage_target_pct"])
    if not package_gate_applies:
        package_gate = "not_applicable_subset"
    elif (coverage_compare_summary.get("target") or {}).get("status") == "tool_unavailable":
        package_gate = "tool_unavailable"
    elif target_line_percent is None:
        package_gate = "target_unmeasured"
    elif target_line_percent + 1e-9 < package_threshold:
        package_gate = "below_target"
    else:
        package_gate = "pass"

    low_coverage_bundles = [
        bundle for bundle in evaluated_bundles
        if bundle["coverage_gate_status"] in {"below_target", "target_unmeasured"}
    ]
    comparison_gap_bundles = [
        bundle for bundle in evaluated_bundles
        if bundle["coverage_gate_status"] == "comparison_gap"
    ]
    regressed_bundles = [bundle for bundle in evaluated_bundles if bundle["regressed_vs_baseline"]]
    status = "completed"
    if (coverage_compare_summary.get("target") or {}).get("status") == "tool_unavailable":
        status = "tool_unavailable"
    elif low_coverage_bundles or comparison_gap_bundles or package_gate == "below_target":
        status = "completed_with_gaps"

    def bundle_brief(bundle: dict) -> dict:
        return {
            "bundle_id": bundle["bundle_id"],
            "bundle_type": bundle.get("bundle_type"),
            "coverage_gate_status": bundle.get("coverage_gate_status"),
            "threshold_pct": bundle.get("threshold_pct"),
            "threshold_display": format_percent(bundle.get("threshold_pct")),
            "baseline_line_percent": bundle.get("baseline_line_percent"),
            "baseline_line_percent_display": format_percent(bundle.get("baseline_line_percent")),
            "target_line_percent": bundle.get("target_line_percent"),
            "target_line_percent_display": format_percent(bundle.get("target_line_percent")),
            "line_percent_delta": bundle.get("line_percent_delta"),
            "line_percent_delta_display": format_percent(bundle.get("line_percent_delta")),
            "shared_test_file_count": len(bundle.get("shared_test_files", [])),
            "shared_test_files": bundle.get("shared_test_files", []),
            "priority_rank": bundle.get("priority_rank"),
            "review_note": bundle.get("review_note"),
        }

    return {
        "run_id": run_id,
        "mode": "coverage_evidence",
        "repo_root": str(repo_root),
        "status": status,
        "bundle_manifest_path": bundle_manifest_path,
        "coverage_compare_run_id": coverage_compare_summary.get("run_id"),
        "bundle_count": len(evaluated_bundles),
        "justified_bundle_count": sum(1 for bundle in evaluated_bundles if bundle["justified_by_tests"]),
        "low_coverage_bundle_count": len(low_coverage_bundles),
        "comparison_gap_bundle_count": len(comparison_gap_bundles),
        "regression_bundle_count": len(regressed_bundles),
        "exception_count": len(exceptions),
        "gate_counts": dict(sorted(gate_counts.items())),
        "delta_counts": dict(sorted(delta_counts.items())),
        "threshold_counts": dict(sorted(threshold_counts.items())),
        "acceptance": acceptance_policy,
        "package_coverage": {
            "baseline_line_percent": baseline_line_percent,
            "baseline_line_percent_display": format_percent(baseline_line_percent),
            "target_line_percent": target_line_percent,
            "target_line_percent_display": format_percent(target_line_percent),
            "line_percent_delta": package_delta,
            "line_percent_delta_display": format_percent(package_delta),
            "threshold_pct": package_threshold,
            "threshold_display": format_percent(package_threshold),
            "gate_status": package_gate,
            "applies": package_gate_applies,
        },
        "below_target_bundles": [bundle_brief(bundle) for bundle in low_coverage_bundles[:50]],
        "comparison_gap_bundles": [bundle_brief(bundle) for bundle in comparison_gap_bundles[:50]],
        "regressed_bundles": [bundle_brief(bundle) for bundle in regressed_bundles[:50]],
        "sample_bundles": [bundle_brief(bundle) for bundle in evaluated_bundles[:50]],
        "bundles": evaluated_bundles,
        "notes": dedupe_paths(
            [str(note) for note in ensure_list(coverage_compare_summary.get("notes")) if str(note).strip()] +
            [str(acceptance_policy.get("wrapper_policy") or "").strip()]
        ),
    }


def parse_json_field(value):
    if value in (None, ""):
        return None
    if isinstance(value, (dict, list)):
        return value
    return json.loads(value)


def load_run_summary_by_mode(
    repo_root: pathlib.Path,
    *,
    mode: str,
    run_id: str | None = None,
    audit_dir: str = DEFAULT_AUDIT_DIR,
) -> dict:
    paths = resolve_paths(repo_root, audit_dir)
    if run_id:
        conn = connect_db(paths["db_path"])
        try:
            row = conn.execute(
                "SELECT summary_json FROM audit_runs WHERE run_id = ? AND mode = ?",
                (run_id, mode),
            ).fetchone()
        finally:
            conn.close()
        if row is None:
            raise RuntimeError(f"Unknown {mode} run id: {run_id}")
        return parse_json_field(row["summary_json"]) or {}

    summary_path_by_mode = {
        "closeout": paths["closeout_json"],
        "coverage_compare": paths["coverage_compare_json"],
        "coverage_evidence": paths["coverage_evidence_json"],
    }
    summary_path = summary_path_by_mode.get(mode)
    if summary_path is None:
        raise RuntimeError(f"Unsupported summary mode: {mode}")
    if not summary_path.exists():
        raise RuntimeError(f"{mode} summary does not exist yet")
    return json.loads(summary_path.read_text(encoding="utf-8"))


def load_closeout_summary(
    repo_root: pathlib.Path,
    closeout_run_id: str | None = None,
    *,
    audit_dir: str = DEFAULT_AUDIT_DIR,
) -> dict:
    return load_run_summary_by_mode(
        repo_root,
        mode="closeout",
        run_id=closeout_run_id,
        audit_dir=audit_dir,
    )


def load_coverage_compare_summary(
    repo_root: pathlib.Path,
    coverage_compare_run_id: str | None = None,
    *,
    audit_dir: str = DEFAULT_AUDIT_DIR,
    compare_path: str | None = None,
) -> dict:
    if compare_path:
        compare_summary_path = resolve_repo_path(repo_root, compare_path)
        return json.loads(compare_summary_path.read_text(encoding="utf-8"))
    return load_run_summary_by_mode(
        repo_root,
        mode="coverage_compare",
        run_id=coverage_compare_run_id,
        audit_dir=audit_dir,
    )


def load_coverage_acceptance_policy(repo_root: pathlib.Path) -> dict:
    plan_path = repo_root / "tools" / "refactor" / "fidelity-audit-plan.json"
    if not plan_path.exists():
        return {
            "package_line_coverage_target_pct": 80,
            "auto_curated_bundle_target_pct": 80,
            "helper_bundle_target_pct": 90,
            "preserved_or_s4_bundle_target_pct": 85,
            "wrapper_policy": "Wrapper/module bundles still require contract and scenario completeness; line coverage alone is not sufficient.",
            "comparison_policy": "Baseline and target coverage must be measured using the same test subsets and recorded side by side.",
        }
    payload = json.loads(plan_path.read_text(encoding="utf-8"))
    coverage_program = payload.get("coverage_backed_auto_curation") or {}
    acceptance = coverage_program.get("acceptance") or {}
    return {
        "package_line_coverage_target_pct": coerce_float(acceptance.get("package_line_coverage_target_pct")) or 80.0,
        "auto_curated_bundle_target_pct": coerce_float(acceptance.get("auto_curated_bundle_target_pct")) or 80.0,
        "helper_bundle_target_pct": coerce_float(acceptance.get("helper_bundle_target_pct")) or 90.0,
        "preserved_or_s4_bundle_target_pct": coerce_float(acceptance.get("preserved_or_s4_bundle_target_pct")) or 85.0,
        "wrapper_policy": str(
            acceptance.get("wrapper_policy")
            or "Wrapper/module bundles still require contract and scenario completeness; line coverage alone is not sufficient."
        ),
        "comparison_policy": str(
            acceptance.get("comparison_policy")
            or "Baseline and target coverage must be measured using the same test subsets and recorded side by side."
        ),
    }


def coverage_threshold_for_bundle(bundle: dict, policy: dict) -> float:
    bundle_type = str(bundle.get("bundle_type") or "")
    if bundle_type in {"helper_surface", "normalized_exact_surface"}:
        return float(policy["helper_bundle_target_pct"])
    if bundle_type == "s4_method_surface":
        return float(policy["preserved_or_s4_bundle_target_pct"])
    return float(policy["auto_curated_bundle_target_pct"])


def coverage_gate_status_for_bundle(
    *,
    baseline_capture_status: str,
    baseline_bundle_status: str,
    baseline_line_percent: float | None,
    target_capture_status: str,
    target_bundle_status: str,
    target_line_percent: float | None,
    threshold_pct: float,
    justified_by_tests: bool,
) -> str:
    if target_capture_status == "tool_unavailable":
        return "tool_unavailable"
    if (
        justified_by_tests
        and baseline_capture_status == "test_failures"
        and target_capture_status == "completed"
        and target_bundle_status == "completed"
        and target_line_percent is not None
        and target_line_percent > 0.0 + 1e-9
    ):
        return "comparison_gap"
    if target_capture_status != "completed":
        return "target_unmeasured"
    if target_bundle_status != "completed" or target_line_percent is None:
        return "target_unmeasured"
    if target_line_percent + 1e-9 < threshold_pct:
        return "below_target"
    if (
        justified_by_tests
        and baseline_capture_status == "completed"
        and baseline_bundle_status == "completed"
        and baseline_line_percent is not None
        and baseline_line_percent <= 0.0 + 1e-9
        and target_line_percent > 0.0 + 1e-9
    ):
        return "comparison_gap"
    return "pass"


def build_coverage_evidence_bundles(
    bundles: list[dict],
    *,
    coverage_compare_summary: dict,
    acceptance_policy: dict,
) -> list[dict]:
    compare_results_by_id = {
        result["bundle_id"]: result
        for result in ensure_list(coverage_compare_summary.get("bundle_results"))
        if isinstance(result, dict) and result.get("bundle_id")
    }
    evaluated: list[dict] = []
    target_capture_status = str((coverage_compare_summary.get("target") or {}).get("status") or "missing")
    baseline_capture_status = str((coverage_compare_summary.get("baseline") or {}).get("status") or "missing")
    for bundle in bundles:
        compare_result = compare_results_by_id.get(bundle["bundle_id"], {})
        threshold_pct = coverage_threshold_for_bundle(bundle, acceptance_policy)
        baseline_line_percent = coerce_float(compare_result.get("baseline_line_percent"))
        target_line_percent = coerce_float(compare_result.get("target_line_percent"))
        line_percent_delta = coerce_float(compare_result.get("line_percent_delta"))
        delta_category = str(compare_result.get("delta_category") or categorize_coverage_delta(line_percent_delta))
        baseline_bundle_status = str(compare_result.get("baseline_status") or baseline_capture_status)
        target_bundle_status = str(compare_result.get("target_status") or target_capture_status)
        shared_test_files = dedupe_paths(
            [str(value) for value in ensure_list(bundle.get("shared_test_files")) if value] +
            [str(value) for value in ensure_list(bundle.get("test_files")) if value]
        )
        gate_status = coverage_gate_status_for_bundle(
            baseline_capture_status=baseline_capture_status,
            baseline_bundle_status=baseline_bundle_status,
            baseline_line_percent=baseline_line_percent,
            target_capture_status=target_capture_status,
            target_bundle_status=target_bundle_status,
            target_line_percent=target_line_percent,
            threshold_pct=threshold_pct,
            justified_by_tests=bool(shared_test_files),
        )
        regressed_vs_baseline = (
            baseline_line_percent is not None
            and target_line_percent is not None
            and target_line_percent + 0.05 < baseline_line_percent
        )
        evaluated.append(
            {
                "bundle_id": bundle["bundle_id"],
                "bundle_type": bundle.get("bundle_type"),
                "primary_entity_key": bundle.get("primary_entity_key"),
                "source_file_path": bundle.get("source_file_path"),
                "baseline_ref": bundle.get("baseline_ref"),
                "target_ref": bundle.get("target_ref"),
                "baseline_paths": bundle.get("baseline_paths", []),
                "target_paths": bundle.get("target_paths", []),
                "baseline_members": ensure_list(bundle.get("baseline_members")),
                "target_members": ensure_list(bundle.get("target_members")),
                "priority_rank": bundle.get("priority_rank"),
                "priority_score": bundle.get("priority_score"),
                "linked_exception_count": bundle.get("linked_exception_count", 0),
                "scenario_class": bundle.get("scenario_class"),
                "differential_replay_status": bundle.get("differential_replay_status"),
                "review_note": bundle.get("review_note"),
                "evidence_depth": bundle.get("evidence_depth"),
                "shared_test_files": shared_test_files,
                "justified_by_tests": bool(shared_test_files),
                "threshold_pct": threshold_pct,
                "coverage_gate_status": gate_status,
                "baseline_status": baseline_bundle_status,
                "target_status": target_bundle_status,
                "baseline_line_percent": baseline_line_percent,
                "target_line_percent": target_line_percent,
                "line_percent_delta": line_percent_delta,
                "delta_category": delta_category,
                "regressed_vs_baseline": regressed_vs_baseline,
                "bundle_result_payload": compare_result.get("result_payload", {}),
            }
        )

    evaluated_by_id = {item["bundle_id"]: item for item in evaluated}
    for item in evaluated:
        bundle_id = str(item.get("bundle_id") or "")
        if (
            item.get("coverage_gate_status") != "target_unmeasured"
            or not bundle_id.startswith("surface::manifest::")
            or not bundle_id.endswith("_legacy")
        ):
            continue
        active_bundle_id = f"{bundle_id[:-len('_legacy')]}_active"
        active_item = evaluated_by_id.get(active_bundle_id)
        if active_item is None or active_item.get("coverage_gate_status") != "pass":
            continue
        if item.get("bundle_type") != active_item.get("bundle_type"):
            continue
        if item.get("shared_test_files") != active_item.get("shared_test_files"):
            continue

        item["coverage_gate_status"] = "pass"
        item["baseline_status"] = active_item.get("baseline_status")
        item["target_status"] = active_item.get("target_status")
        item["baseline_line_percent"] = active_item.get("baseline_line_percent")
        item["target_line_percent"] = active_item.get("target_line_percent")
        item["line_percent_delta"] = active_item.get("line_percent_delta")
        item["delta_category"] = str(active_item.get("delta_category") or item.get("delta_category") or "flat")
        item["regressed_vs_baseline"] = bool(active_item.get("regressed_vs_baseline"))
        note = str(item.get("review_note") or "").strip()
        active_note = (
            f" Coverage accepted via active sibling `{active_bundle_id}` reaching "
            f"{format_percent(active_item.get('target_line_percent'))} under the same shared suite."
        )
        item["review_note"] = f"{note}{active_note}" if note else active_note.strip()
        payload = dict(item.get("bundle_result_payload") or {})
        payload["coverage_proxy_bundle_id"] = active_bundle_id
        item["bundle_result_payload"] = payload

    evaluated.sort(
        key=lambda item: (
            0 if item["coverage_gate_status"] in {"below_target", "target_unmeasured", "comparison_gap"} else 1,
            0 if item["regressed_vs_baseline"] else 1,
            item["priority_rank"] if item.get("priority_rank") is not None else 999999,
            item["bundle_id"],
        )
    )
    return evaluated


def split_identifier_tokens(value: str) -> list[str]:
    normalized = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", value)
    normalized = re.sub(r"[^A-Za-z0-9]+", "_", normalized)
    tokens = [token.lower() for token in normalized.split("_") if token]
    return dedupe_paths(tokens)


def extract_surface_entity_name(entity_key: str | None) -> str | None:
    if not entity_key:
        return None
    parts = str(entity_key).split("::")
    if len(parts) < 2:
        return None
    if parts[0] == "setMethod":
        return "::".join(parts[1:])
    return parts[1]


def selector_kind_to_bundle_type(selector_surface_key: str | None) -> str:
    if not selector_surface_key:
        return "normalized_exact_surface"
    if selector_surface_key.startswith("setMethod::"):
        return "s4_method_surface"
    return "helper_surface"


def bundle_type_from_disposition(disposition: str | None, selector_surface_key: str | None) -> str:
    disposition = str(disposition or "")
    if disposition == "equivalent_target_resolved_without_baseline_ancestor":
        return "lineage_family"
    if disposition == "entrypoint_shell_equivalent_under_tests":
        return "wrapper_entrypoint"
    if disposition == "helper_equivalent_under_tests":
        return "helper_surface"
    if disposition == "method_equivalent_under_tests":
        return "s4_method_surface"
    if disposition == "canonical_duplicate_retired":
        return "duplicate_canonical"
    if disposition == "manual_merge_canonicalized":
        return "manual_merge_canonical"
    if disposition in {"equivalent_target_resolution_asymmetry", "equivalent_whitespace_only", "equivalent_comment_only"}:
        return selector_kind_to_bundle_type(selector_surface_key)
    if disposition == "covered_by_manifest_semantic_review":
        return selector_kind_to_bundle_type(selector_surface_key)
    return selector_kind_to_bundle_type(selector_surface_key)


def scenario_class_for_bundle_type(bundle_type: str) -> str:
    return {
        "lineage_family": "legacy_lineage_family",
        "wrapper_entrypoint": "wrapper_contract",
        "manual_merge_canonical": "manual_merge",
        "duplicate_canonical": "compatibility_duplicate",
        "s4_method_surface": "s4_method",
        "helper_surface": "helper_unit",
        "normalized_exact_surface": "normalized_exact",
    }.get(bundle_type, "unspecified")


def differential_replay_status_for_disposition(disposition: str | None) -> str:
    disposition = str(disposition or "")
    if disposition == "entrypoint_shell_equivalent_under_tests":
        return "contract_curated"
    if disposition in {"helper_equivalent_under_tests", "method_equivalent_under_tests"}:
        return "behavior_curated"
    if disposition in {"manual_merge_canonicalized", "canonical_duplicate_retired"}:
        return "canonicalized_under_tests"
    if disposition == "covered_by_manifest_semantic_review":
        return "manifest_semantic_overlap"
    if disposition in {
        "equivalent_target_resolved_without_baseline_ancestor",
        "equivalent_target_resolution_asymmetry",
        "equivalent_whitespace_only",
        "equivalent_comment_only",
    }:
        return "structural_equivalence_only"
    return "unknown"


def normalize_exception_file_paths(record: dict) -> tuple[list[str], list[str]]:
    evidence = record.get("evidence") or {}
    source_paths: list[str] = []
    target_paths: list[str] = []
    if record.get("audit_layer") == "manifest_fidelity":
        if evidence.get("source_path"):
            source_paths.append(str(evidence["source_path"]))
        if evidence.get("target_path"):
            target_paths.append(str(evidence["target_path"]))
    elif record.get("audit_layer") == "surface_inventory":
        nested = evidence.get("evidence") or {}
        for key in ("baseline_file",):
            if nested.get(key):
                source_paths.append(str(nested[key]))
        for key in ("target_file",):
            if nested.get(key):
                target_paths.append(str(nested[key]))
        for key in ("baseline_files",):
            source_paths.extend(str(value) for value in ensure_list(nested.get(key)) if value)
        for key in ("target_files",):
            target_paths.extend(str(value) for value in ensure_list(nested.get(key)) if value)
    return sorted(set(source_paths)), sorted(set(target_paths))


def infer_module_family(text: str) -> str | None:
    lowered = text.lower()
    for family in (
        "import",
        "norm",
        "qc",
        "summary",
        "design_builder",
        "design",
        "da",
        "enrich",
        "annotation",
        "rollup",
        "plotting",
        "filemgmt",
    ):
        family_token = family.replace("_", "")
        if family_token in lowered.replace("_", ""):
            return family
    return None


def infer_lineage_cluster_name(*values: str | None) -> str:
    tokens: set[str] = set()
    for value in values:
        if value:
            tokens.update(split_identifier_tokens(str(value)))

    if not tokens:
        return "misc"
    if "observer" in tokens and "setup" in tokens:
        return "observer_setup"
    if "observer" in tokens and ("register" in tokens or "registration" in tokens):
        return "observer_register"
    if "observer" in tokens and tokens.intersection(
        {"preflight", "handoff", "finalize", "failure", "completion", "working", "shell", "run"}
    ):
        return "observer_runtime"
    if ("run" in tokens and "analysis" in tokens and "body" in tokens) or tokens.intersection(
        {"error", "failure", "notify", "message", "report", "log"}
    ):
        return "observer_runtime"
    if "renderer" in tokens or "render" in tokens:
        return "output_renderer"
    if ("output" in tokens or "results" in tokens or "summary" in tokens) and "setup" in tokens:
        return "output_setup"
    if ("output" in tokens or "results" in tokens or "summary" in tokens) and (
        "register" in tokens or "registration" in tokens
    ):
        return "output_register"
    if "formatter" in tokens or (
        "text" in tokens and tokens.intersection({"status", "summary", "analysis", "contrast", "contrasts"})
    ):
        return "formatter"
    if tokens.intersection({"download", "archive", "filename", "zip"}):
        return "download_io"
    if tokens.intersection(
        {"builder", "preparer", "propagator", "payload", "capturer", "resolver", "persister", "finalizer", "completer", "updater"}
    ):
        return "builder_resolver"
    if tokens.intersection({"reactive", "supported", "organisms", "analysis", "method", "taxon", "species", "contrast"}):
        return "reactive_state"
    if "setup" in tokens:
        return "bootstrap_setup"
    if tokens.intersection({"plot", "plotly"}):
        return "plotting"
    return "misc"


def build_range_member(
    *,
    entity_key: str | None,
    file_path: str | None,
    line_start: int | None,
    line_end: int | None,
    side: str,
    selector_kind: str | None = None,
    comparison_status: str | None = None,
) -> dict | None:
    if not file_path:
        return None
    if line_start is None and line_end is None:
        return None
    return {
        "entity_key": entity_key,
        "file_path": str(file_path),
        "line_start": line_start,
        "line_end": line_end,
        "metadata": {
            "side": side,
            "selector_kind": selector_kind,
            "comparison_status": comparison_status,
        },
    }


def append_unique_member(target: list[dict], member: dict | None) -> None:
    if member is None:
        return
    signature = (
        member.get("entity_key"),
        member.get("file_path"),
        member.get("line_start"),
        member.get("line_end"),
    )
    for existing in target:
        existing_signature = (
            existing.get("entity_key"),
            existing.get("file_path"),
            existing.get("line_start"),
            existing.get("line_end"),
        )
        if existing_signature == signature:
            return
    target.append(member)


def build_manifest_selector_index(conn: sqlite3.Connection, run_ids: list[str]) -> dict[tuple[str, str, str], dict]:
    if not run_ids:
        return {}
    rows = conn.execute(
        f"""
        SELECT
            me.run_id,
            me.manifest_path,
            me.entry_id,
            me.selector_kind,
            me.selector_payload_json,
            me.source_path,
            me.target_path,
            mc.status AS comparison_status,
            mc.evidence_json
        FROM manifest_entries me
        LEFT JOIN manifest_comparisons mc
          ON mc.manifest_entry_id = me.manifest_entry_id
         AND mc.run_id = me.run_id
        WHERE me.run_id IN ({", ".join("?" for _ in run_ids)})
        """,
        run_ids,
    ).fetchall()
    index: dict[tuple[str, str, str], dict] = {}
    for row in rows:
        selector_payload = parse_json_field(row["selector_payload_json"]) or {}
        evidence = parse_json_field(row["evidence_json"]) or {}
        surface_entity_key = selector_payload_to_surface_entity_key(str(row["selector_kind"]), selector_payload)
        index[(str(row["run_id"]), str(row["manifest_path"]), str(row["entry_id"]))] = {
            "selector_kind": str(row["selector_kind"]),
            "selector_payload": selector_payload,
            "surface_entity_key": surface_entity_key,
            "source_path": str(row["source_path"]),
            "target_path": row["target_path"],
            "comparison_status": str(row["comparison_status"] or ""),
            "source_start_line": coerce_int(evidence.get("source_start_line")),
            "source_end_line": coerce_int(evidence.get("source_end_line")),
            "target_start_line": coerce_int(evidence.get("target_start_line")),
            "target_end_line": coerce_int(evidence.get("target_end_line")),
        }
    return index


def build_inventory_entity_index(
    conn: sqlite3.Connection,
    run_ids: list[str],
) -> dict[tuple[str, str, str], dict]:
    if not run_ids:
        return {}
    rows = conn.execute(
        f"""
        SELECT
            run_id,
            side,
            entity_key,
            entity_kind,
            entity_name,
            signature_key,
            file_path,
            line_start,
            line_end,
            collate_index
        FROM inventory_entities
        WHERE run_id IN ({", ".join("?" for _ in run_ids)})
        """,
        run_ids,
    ).fetchall()
    index: dict[tuple[str, str, str], dict] = {}
    by_entity: dict[tuple[str, str], dict[str, dict]] = defaultdict(dict)

    def upsert_candidate(
        *,
        side: str,
        entity_key: str,
        file_path: str,
        candidate: dict,
    ) -> None:
        key = (side, entity_key, file_path)
        existing = index.get(key)
        if existing is None:
            index[key] = candidate
            by_entity[(side, entity_key)][file_path] = candidate
            return
        existing_start = coerce_int(existing.get("line_start"))
        existing_end = coerce_int(existing.get("line_end"))
        candidate_start = candidate.get("line_start")
        candidate_end = candidate.get("line_end")
        existing_span = (
            (existing_end - existing_start)
            if existing_start is not None and existing_end is not None
            else None
        )
        candidate_span = (
            (candidate_end - candidate_start)
            if candidate_start is not None and candidate_end is not None
            else None
        )
        if existing_span is None or (
            candidate_span is not None and candidate_span < existing_span
        ):
            index[key] = candidate
            by_entity[(side, entity_key)][file_path] = candidate
        elif file_path not in by_entity[(side, entity_key)]:
            by_entity[(side, entity_key)][file_path] = existing

    for row in rows:
        entity_key = str(row["entity_key"] or "")
        file_path = str(row["file_path"] or "")
        if not entity_key or not file_path:
            continue
        side = str(row["side"] or "")
        entity_kind = str(row["entity_kind"] or "")
        entity_name = str(row["entity_name"] or "")
        candidate = {
            "side": side,
            "entity_key": entity_key,
            "entity_kind": entity_kind,
            "entity_name": entity_name,
            "signature_key": row["signature_key"],
            "file_path": file_path,
            "line_start": coerce_int(row["line_start"]),
            "line_end": coerce_int(row["line_end"]),
            "collate_index": coerce_int(row["collate_index"]),
        }
        upsert_candidate(
            side=side,
            entity_key=entity_key,
            file_path=file_path,
            candidate=candidate,
        )
        if entity_kind and entity_name:
            upsert_candidate(
                side=side,
                entity_key=f"selector_lookup::{entity_kind}::{entity_name}",
                file_path=file_path,
                candidate=candidate,
            )
    for (side, entity_key), path_candidates in by_entity.items():
        if len(path_candidates) == 1:
            index[(side, entity_key, "")] = next(iter(path_candidates.values()))
        else:
            active_candidate = max(
                path_candidates.values(),
                key=lambda candidate: (
                    candidate.get("collate_index") is not None,
                    candidate.get("collate_index") or -1,
                    candidate.get("line_start") or -1,
                ),
            )
            index[(side, entity_key, "<active>")] = active_candidate
    return index


def build_contract_search_index(target_root: pathlib.Path) -> list[dict]:
    records = []
    for record in catalog_contract_files(target_root, None):
        file_path = target_root / record["file_path"]
        shared_compare = False
        try:
            file_text = file_path.read_text(encoding="utf-8")
            shared_compare = "fidelity-coverage-compare: shared" in file_text.lower()
        except OSError:
            file_text = ""
            shared_compare = False
        case_text = " ".join(
            [
                case["test_name"]
                + " "
                + str(case.get("entity_hint") or "")
                + " "
                + " ".join(case.get("scenario_tags", []))
                for case in record["cases"]
            ]
        ).lower()
        entity_hints = {
            str(case.get("entity_hint"))
            for case in record["cases"]
            if case.get("entity_hint")
        }
        entity_hint_tokens = {
            token.lower()
            for hint in entity_hints
            for token in split_identifier_tokens(hint)
            if token
        }
        path_tokens = split_identifier_tokens(pathlib.Path(record["file_path"]).stem)
        source_tokens = {
            token.lower()
            for token in split_identifier_tokens(file_text)
            if len(token) >= 3
        } if shared_compare else set()
        records.append(
            {
                "file_path": record["file_path"],
                "family": record["family"],
                "module_family": record["module_family"],
                "surface": record["surface"],
                "shared_compare": shared_compare,
                "entity_hints": {hint.lower() for hint in entity_hints},
                "entity_hint_tokens": entity_hint_tokens,
                "path_tokens": {token.lower() for token in path_tokens if token},
                "source_tokens": source_tokens,
                "search_text": case_text + " " + " ".join(path_tokens),
            }
        )
    return records


def gather_bundle_tokens(bundle: dict) -> set[str]:
    tokens: set[str] = set()
    for value in ensure_list(bundle.get("entity_tokens")):
        if value:
            tokens.add(str(value).lower())
    for path_value in ensure_list(bundle.get("baseline_paths")) + ensure_list(bundle.get("target_paths")):
        path_tokens = split_identifier_tokens(pathlib.Path(str(path_value)).stem)
        tokens.update(token for token in path_tokens if len(token) >= 3)
    return tokens


GENERIC_BUNDLE_TOKENS = {
    "mod",
    "module",
    "modules",
    "symbol",
    "surface",
    "lineage",
    "manifest",
    "method",
    "methods",
    "func",
    "function",
    "functions",
    "object",
    "objects",
    "norm",
    "server",
    "ui",
    "helper",
    "helpers",
    "render",
    "renderer",
    "renderers",
    "output",
    "outputs",
    "observer",
    "observers",
    "runtime",
    "register",
    "registration",
    "setup",
    "shell",
    "body",
    "public",
    "wrapper",
    "entry",
    "entrypoint",
    "misc",
    "bootstrap",
    "state",
    "reactive",
    "builder",
    "manifest",
    "refactor",
    "tools",
    "fidelity",
    "coverage",
    "compare",
    "shared",
    "drift",
    "whitespace",
    "only",
    "yml",
}


OMICS_BUNDLE_TOKENS = {
    "prot",
    "protein",
    "proteins",
    "proteomics",
    "pept",
    "peptide",
    "peptides",
    "metab",
    "metabolite",
    "metabolites",
    "metabolomics",
    "lipid",
    "lipids",
    "lipidomics",
}


OMICS_TOKEN_SYNONYM_GROUPS = (
    {"prot", "protein", "proteins", "proteomics"},
    {"pept", "peptide", "peptides"},
    {"metab", "metabolite", "metabolites", "metabolomics"},
    {"lipid", "lipids", "lipidomics"},
)


def expand_omics_token_synonyms(tokens: set[str]) -> set[str]:
    expanded = set(tokens)
    for synonym_group in OMICS_TOKEN_SYNONYM_GROUPS:
        if expanded & synonym_group:
            expanded.update(synonym_group)
    return expanded


def gather_specific_bundle_tokens(bundle: dict) -> set[str]:
    tokens = gather_bundle_tokens(bundle)
    tokens.update(
        token.lower()
        for token in split_identifier_tokens(str(bundle.get("entity_name") or ""))
        if token
    )
    tokens.update(
        token.lower()
        for token in split_identifier_tokens(str(bundle.get("primary_entity_key") or ""))
        if token
    )
    return {
        token
        for token in tokens
        if len(token) >= 4 and token not in GENERIC_BUNDLE_TOKENS
    }


def bundle_test_overlap_tokens(bundle: dict, test_record: dict) -> set[str]:
    def token_values(value) -> list:
        if isinstance(value, (set, list, tuple)):
            return list(value)
        return ensure_list(value)

    token_space = expand_omics_token_synonyms({
        str(token).lower()
        for token in token_values(test_record.get("path_tokens"))
    } | {
        str(token).lower()
        for token in token_values(test_record.get("entity_hint_tokens"))
    } | {
        str(token).lower()
        for token in token_values(test_record.get("source_tokens"))
    })
    bundle_tokens = expand_omics_token_synonyms(gather_specific_bundle_tokens(bundle))
    return {
        token
        for token in bundle_tokens
        if token in token_space
    }


def count_bundle_test_token_overlap(bundle: dict, test_record: dict) -> int:
    return len(bundle_test_overlap_tokens(bundle, test_record))


def score_test_record_for_bundle(bundle: dict, test_record: dict) -> int:
    score = 0
    entity_name = str(bundle.get("entity_name") or "").lower()
    primary_entity_key = str(bundle.get("primary_entity_key") or "").lower()
    module_family = bundle.get("module_family")
    bundle_type = str(bundle.get("bundle_type") or "")
    tokens = gather_bundle_tokens(bundle)
    overlap_count = count_bundle_test_token_overlap(bundle, test_record)

    if entity_name and entity_name.lower() in test_record["entity_hints"]:
        score += 24
    if primary_entity_key and primary_entity_key.lower() in test_record["search_text"]:
        score += 12
    for token in sorted(tokens):
        if len(token) < 4:
            continue
        if token in test_record["entity_hints"]:
            score += 10
        elif token in test_record["search_text"]:
            score += 4
    if module_family and module_family == test_record["module_family"]:
        score += 8
    if bundle_type == "wrapper_entrypoint" and test_record["family"] == "module_contract":
        score += 4
    if bundle_type in {"helper_surface", "normalized_exact_surface"} and test_record["family"] in {"general_unit", "characterization"}:
        score += 2
    if bundle_type == "s4_method_surface" and test_record["family"] in {"characterization", "compat", "general_unit"}:
        score += 2
    if overlap_count:
        score += overlap_count * 6
    return score


def build_shared_test_candidates(bundle: dict, contract_index: list[dict]) -> list[dict]:
    entity_name = str(bundle.get("entity_name") or "").lower()
    primary_entity_key = str(bundle.get("primary_entity_key") or "").lower()
    bundle_type = str(bundle.get("bundle_type") or "")
    candidates: list[dict] = []
    for test_record in contract_index:
        if not test_record.get("shared_compare"):
            continue
        family = str(test_record.get("family") or "")
        exact_entity = bool(entity_name and entity_name in test_record["entity_hints"])
        exact_primary = bool(primary_entity_key and primary_entity_key in test_record["search_text"])
        overlap_tokens = bundle_test_overlap_tokens(bundle, test_record)
        overlap_count = len(overlap_tokens)
        bundle_omics_tokens = gather_specific_bundle_tokens(bundle) & OMICS_BUNDLE_TOKENS
        non_omics_overlap_count = len(
            [token for token in overlap_tokens if token not in OMICS_BUNDLE_TOKENS]
        )
        if not exact_entity and not exact_primary and overlap_count == 0:
            continue
        if (
            bundle_omics_tokens
            and not exact_entity
            and not exact_primary
            and not (overlap_tokens & bundle_omics_tokens)
        ):
            continue
        if bundle_type == "wrapper_entrypoint":
            if family not in {"characterization", "compat", "module_contract"}:
                continue
            if not exact_entity and not exact_primary and (overlap_count < 2 or non_omics_overlap_count < 1):
                continue
        elif bundle_type == "lineage_family":
            if family not in {"characterization", "compat", "module_contract", "general_unit"}:
                continue
            if not exact_entity and not exact_primary and (overlap_count < 2 or non_omics_overlap_count < 1):
                continue
        elif bundle_type in {"helper_surface", "normalized_exact_surface", "s4_method_surface", "manual_merge_canonical", "duplicate_canonical"}:
            if family not in {"characterization", "compat", "general_unit"}:
                continue
            if not exact_entity and not exact_primary and (overlap_count < 2 or non_omics_overlap_count < 1):
                continue
        elif not exact_entity and not exact_primary and (overlap_count < 2 or non_omics_overlap_count < 1):
            continue
        candidates.append(
            {
                "file_path": str(test_record["file_path"]),
                "family": family,
                "score": score_test_record_for_bundle(bundle, test_record),
                "token_overlap": overlap_count,
                "exact_entity": exact_entity,
                "exact_primary": exact_primary,
            }
        )
    return candidates


def select_shared_tests_for_bundle(bundle: dict, contract_index: list[dict]) -> list[str]:
    candidates = build_shared_test_candidates(bundle, contract_index)
    if not candidates:
        return []

    strong_cross_version = [
        candidate
        for candidate in candidates
        if candidate["family"] in {"characterization", "compat"}
        and (candidate["exact_entity"] or candidate["exact_primary"] or candidate["token_overlap"] >= 2)
    ]
    strong_shared_units = [
        candidate
        for candidate in candidates
        if candidate["family"] == "general_unit"
        and (
            candidate["exact_entity"]
            or candidate["exact_primary"]
            or (candidate["token_overlap"] >= 3 and candidate["score"] >= 50)
        )
    ]
    strong_shared_units.sort(
        key=lambda item: (
            0 if item["exact_entity"] or item["exact_primary"] else 1,
            -item["token_overlap"],
            -item["score"],
            item["file_path"],
        )
    )
    if strong_cross_version:
        selected = strong_cross_version + strong_shared_units
    elif str(bundle.get("bundle_type") or "") == "wrapper_entrypoint":
        strong_contract = [
            candidate
            for candidate in candidates
            if candidate["family"] == "module_contract"
            and (candidate["exact_entity"] or candidate["exact_primary"] or candidate["token_overlap"] >= 2)
        ]
        selected = strong_contract + strong_shared_units if strong_contract else candidates
    else:
        selected = candidates

    selected.sort(
        key=lambda item: (
            0 if item["family"] in {"characterization", "compat"} else 1,
            0 if item["exact_entity"] or item["exact_primary"] else 1,
            -item["token_overlap"],
            -item["score"],
            item["file_path"],
        )
    )
    selected_files = dedupe_paths(
        [item["file_path"] for item in strong_shared_units] +
        [item["file_path"] for item in selected]
    )
    return selected_files[:8]


def count_lines_if_exists(repo_root: pathlib.Path, file_path: str) -> int:
    candidate = repo_root / file_path
    if not candidate.exists():
        return 0
    try:
        return len(candidate.read_text(encoding="utf-8").splitlines())
    except OSError:
        return 0


def bundle_note_for_disposition(disposition: str | None) -> str:
    return {
        "equivalent_target_resolved_without_baseline_ancestor": "Grouped lineage family where manifest source ancestry no longer resolves directly, but extracted target helpers resolve cleanly.",
        "equivalent_target_resolution_asymmetry": "Grouped exact-fidelity bundle that required fallback target resolution.",
        "equivalent_whitespace_only": "Grouped normalized-only fidelity bundle with whitespace-only source drift.",
        "equivalent_comment_only": "Grouped AST-only fidelity bundle with comment or roxygen drift.",
        "covered_by_manifest_semantic_review": "Grouped redundant surface drift under the accepted manifest semantic bundle.",
        "helper_equivalent_under_tests": "Grouped helper surface accepted as equivalent under direct tests.",
        "method_equivalent_under_tests": "Grouped S4 method surface accepted as equivalent under direct tests.",
        "entrypoint_shell_equivalent_under_tests": "Grouped wrapper entrypoint accepted as equivalent under contract tests.",
        "canonical_duplicate_retired": "Grouped canonical duplicate retirement bundle.",
        "manual_merge_canonicalized": "Grouped manual merge bundle that retains one canonical implementation.",
    }.get(str(disposition or ""), "Grouped parity evidence bundle.")


def build_bundle_groups(
    *,
    repo_root: pathlib.Path,
    baseline_ref: str,
    target_ref: str,
    exception_records: list[dict],
    manifest_selector_index: dict[tuple[str, str, str], dict],
    inventory_entity_index: dict[tuple[str, str, str], dict],
    contract_index: list[dict],
) -> list[dict]:
    range_inventory_entries: list[dict] = []
    seen_range_inventory_entries: set[tuple] = set()
    for candidate in inventory_entity_index.values():
        candidate_key = (
            candidate.get("side"),
            candidate.get("entity_key"),
            candidate.get("file_path"),
            candidate.get("line_start"),
            candidate.get("line_end"),
        )
        if candidate_key in seen_range_inventory_entries:
            continue
        seen_range_inventory_entries.add(candidate_key)
        if candidate.get("entity_kind") in {"setMethod", "setClass", "setGeneric", "symbol"}:
            range_inventory_entries.append(candidate)

    def infer_active_inventory_entry_for_range(
        side: str,
        file_path: str,
        line_start,
        line_end,
    ) -> dict | None:
        line_start = coerce_int(line_start)
        line_end = coerce_int(line_end)
        if not file_path or line_start is None or line_end is None:
            return None
        matches = []
        for candidate in range_inventory_entries:
            candidate_start = coerce_int(candidate.get("line_start"))
            candidate_end = coerce_int(candidate.get("line_end"))
            if candidate_start is None or candidate_end is None:
                continue
            if str(candidate.get("side") or "") != side:
                continue
            if str(candidate.get("file_path") or "") != file_path:
                continue
            if candidate_start >= line_start and candidate_end <= line_end:
                matches.append(candidate)
        entity_keys = {str(candidate.get("entity_key") or "") for candidate in matches}
        entity_keys.discard("")
        if len(entity_keys) != 1:
            return None
        selected = max(
            matches,
            key=lambda candidate: (
                candidate.get("collate_index") is not None,
                candidate.get("collate_index") or -1,
                candidate.get("line_start") or -1,
            ),
        )
        full_entity_key = str(selected.get("entity_key") or "")
        return inventory_entity_index.get((side, full_entity_key, "<active>")) or selected

    def resolve_inventory_entry(
        side: str,
        entity_key: str,
        file_path: str,
        line_start=None,
        line_end=None,
    ) -> dict | None:
        if entity_key.startswith("manifest::"):
            range_entry = infer_active_inventory_entry_for_range(
                side,
                file_path,
                line_start,
                line_end,
            )
            if range_entry is not None:
                return range_entry
        if entity_key.startswith("selector_lookup::"):
            selector_entry = inventory_entity_index.get((side, entity_key, file_path))
            if selector_entry is None:
                selector_entry = inventory_entity_index.get((side, entity_key, ""))
            if selector_entry is None:
                return None
            full_entity_key = str(selector_entry.get("entity_key") or "")
            active_entry = inventory_entity_index.get((side, full_entity_key, "<active>"))
            return active_entry or selector_entry

        active_entry = inventory_entity_index.get((side, entity_key, "<active>"))
        if active_entry is not None:
            return active_entry
        inventory_entry = inventory_entity_index.get((side, entity_key, file_path))
        if inventory_entry is not None:
            return inventory_entry
        return inventory_entity_index.get((side, entity_key, ""))

    groups: dict[str, dict] = {}
    for record in exception_records:
        disposition = str(record.get("disposition") or "")
        source_paths, target_paths = normalize_exception_file_paths(record)
        selector_surface_key = None
        selector_lookup_key = None
        selector_meta = None
        manifest_key = parse_manifest_exception_entity_key(str(record["entity_key"]))
        if manifest_key is not None:
            selector_meta = manifest_selector_index.get((str(record["run_id"]), manifest_key[0], manifest_key[1]))
            if selector_meta is not None:
                selector_surface_key = selector_meta.get("surface_entity_key")
                selector_payload = selector_meta.get("selector_payload") or {}
                selector_value = selector_payload.get("value")
                selector_kind = selector_meta.get("selector_kind")
                if not selector_surface_key and selector_kind and selector_value:
                    selector_lookup_key = f"selector_lookup::{selector_kind}::{selector_value}"
                if not source_paths and selector_meta.get("source_path"):
                    source_paths = [str(selector_meta["source_path"])]
                if not target_paths and selector_meta.get("target_path"):
                    target_paths = [str(selector_meta["target_path"])]

        canonical_entity_key = selector_surface_key or str(record["entity_key"])
        bundle_type = bundle_type_from_disposition(disposition, selector_surface_key or canonical_entity_key)
        if disposition == "equivalent_target_resolved_without_baseline_ancestor" and source_paths:
            lineage_cluster = infer_lineage_cluster_name(
                manifest_key[1] if manifest_key is not None else None,
                selector_surface_key,
                record.get("exception_key"),
                record.get("entity_key"),
            )
            group_key = f"lineage::{source_paths[0]}::{lineage_cluster}"
            primary_entity_key = f"lineage::{source_paths[0]}::{lineage_cluster}"
        else:
            group_key = f"surface::{canonical_entity_key}"
            primary_entity_key = canonical_entity_key

        bundle = groups.get(group_key)
        if bundle is None:
            bundle = {
                "group_key": group_key,
                "bundle_id": re.sub(r"[^A-Za-z0-9_.:-]+", "_", group_key),
                "bundle_type": bundle_type,
                "primary_entity_key": primary_entity_key,
                "entity_name": extract_surface_entity_name(selector_surface_key or canonical_entity_key),
                "baseline_ref": baseline_ref,
                "target_ref": target_ref,
                "baseline_paths": [],
                "target_paths": [],
                "linked_exception_keys": [],
                "linked_exception_ids": [],
                "linked_exceptions": [],
                "entity_tokens": [],
                "module_family": None,
                "baseline_members": [],
                "target_members": [],
                "scenario_class": scenario_class_for_bundle_type(bundle_type),
                "differential_replay_status": differential_replay_status_for_disposition(disposition),
                "coverage_gate_status": "pending_coverage_capture",
                "review_note": bundle_note_for_disposition(disposition),
            }
            groups[group_key] = bundle

        bundle["baseline_paths"] = sorted(set(bundle["baseline_paths"]) | set(source_paths))
        bundle["target_paths"] = sorted(set(bundle["target_paths"]) | set(target_paths))
        bundle["linked_exception_keys"].append(str(record.get("exception_key") or ""))
        bundle["linked_exception_ids"].append(str(record["exception_id"]))
        bundle["linked_exceptions"].append(record)
        if record.get("exception_key"):
            bundle["entity_tokens"].extend(split_identifier_tokens(str(record["exception_key"])))
        if primary_entity_key:
            bundle["entity_tokens"].extend(split_identifier_tokens(primary_entity_key))
        if not bundle.get("module_family"):
            for candidate in [primary_entity_key, *source_paths, *target_paths]:
                family = infer_module_family(str(candidate))
                if family:
                    bundle["module_family"] = family
                    break
        if selector_meta is not None:
            member_entity_key = selector_surface_key or selector_lookup_key or canonical_entity_key
            append_unique_member(
                bundle["baseline_members"],
                build_range_member(
                    entity_key=member_entity_key,
                    file_path=selector_meta.get("source_path"),
                    line_start=selector_meta.get("source_start_line"),
                    line_end=selector_meta.get("source_end_line"),
                    side="baseline",
                    selector_kind=selector_meta.get("selector_kind"),
                    comparison_status=selector_meta.get("comparison_status"),
                ),
            )
            append_unique_member(
                bundle["target_members"],
                build_range_member(
                    entity_key=member_entity_key,
                    file_path=selector_meta.get("target_path"),
                    line_start=selector_meta.get("target_start_line"),
                    line_end=selector_meta.get("target_end_line"),
                    side="target",
                    selector_kind=selector_meta.get("selector_kind"),
                    comparison_status=selector_meta.get("comparison_status"),
                ),
            )

    bundles: list[dict] = []
    for bundle in groups.values():
        if bundle["bundle_type"] != "lineage_family":
            for side, path_key, member_key in (
                ("baseline", "baseline_paths", "baseline_members"),
                ("target", "target_paths", "target_members"),
            ):
                enriched_members: list[dict] = []
                original_member_paths = {
                    str(member.get("file_path") or "")
                    for member in ensure_list(bundle.get(member_key))
                    if member.get("file_path")
                }
                for existing_member in ensure_list(bundle.get(member_key)):
                    entity_key = str(existing_member.get("entity_key") or bundle["primary_entity_key"] or "")
                    file_path = str(existing_member.get("file_path") or "")
                    lookup_key = (
                        side,
                        entity_key,
                        file_path,
                    )
                    inventory_entry = resolve_inventory_entry(
                        side,
                        entity_key,
                        file_path,
                        existing_member.get("line_start"),
                        existing_member.get("line_end"),
                    )
                    append_unique_member(
                        enriched_members,
                        {
                            **existing_member,
                            "file_path": inventory_entry.get("file_path") if inventory_entry else file_path,
                            "line_start": inventory_entry.get("line_start") if inventory_entry else existing_member.get("line_start"),
                            "line_end": inventory_entry.get("line_end") if inventory_entry else existing_member.get("line_end"),
                        },
                    )
                if not enriched_members:
                    for path_value in bundle.get(path_key, []):
                        entity_key = str(bundle["primary_entity_key"] or "")
                        lookup_key = (
                            side,
                            entity_key,
                            str(path_value),
                        )
                        inventory_entry = resolve_inventory_entry(side, entity_key, str(path_value))
                        append_unique_member(
                            enriched_members,
                            {
                                "entity_key": entity_key,
                                "file_path": inventory_entry.get("file_path") if inventory_entry else path_value,
                                "line_start": inventory_entry.get("line_start") if inventory_entry else None,
                                "line_end": inventory_entry.get("line_end") if inventory_entry else None,
                                "metadata": {"side": side},
                            },
                        )
                represented_paths = {
                    str(member.get("file_path") or "")
                    for member in enriched_members
                    if member.get("file_path")
                }
                for path_value in bundle.get(path_key, []):
                    path_value = str(path_value)
                    if path_value in represented_paths or path_value in original_member_paths:
                        continue
                    entity_key = str(bundle["primary_entity_key"] or "")
                    lookup_key = (
                        side,
                        entity_key,
                        path_value,
                    )
                    inventory_entry = resolve_inventory_entry(side, entity_key, path_value)
                    candidate_path = inventory_entry.get("file_path") if inventory_entry else path_value
                    if candidate_path in represented_paths:
                        continue
                    append_unique_member(
                        enriched_members,
                        {
                            "entity_key": entity_key,
                            "file_path": candidate_path,
                            "line_start": inventory_entry.get("line_start") if inventory_entry else None,
                            "line_end": inventory_entry.get("line_end") if inventory_entry else None,
                            "metadata": {"side": side},
                        },
                    )
                bundle[member_key] = enriched_members
                represented_paths = sorted(
                    {
                        str(member.get("file_path") or "")
                        for member in enriched_members
                        if member.get("file_path")
                    }
                )
                if represented_paths:
                    bundle[path_key] = represented_paths

        primary_path = (
            (bundle["target_paths"][0] if bundle["target_paths"] else None)
            or (bundle["baseline_paths"][0] if bundle["baseline_paths"] else None)
        )
        line_counts = [
            count_lines_if_exists(repo_root, path_value)
            for path_value in sorted(set(bundle["baseline_paths"] + bundle["target_paths"]))
        ]
        bundle["source_file_path"] = primary_path
        bundle["linked_exception_count"] = len(bundle["linked_exceptions"])
        dispositions = Counter(str(item.get("disposition") or "unset") for item in bundle["linked_exceptions"])
        bundle["disposition_counts"] = dict(sorted(dispositions.items()))
        dominant_disposition = dispositions.most_common(1)[0][0] if dispositions else "unset"
        evidence_depth = DISPOSITION_EVIDENCE_DEPTH.get(dominant_disposition, "T2_structural")
        bundle["evidence_depth"] = evidence_depth
        bundle["evidence_depth_rank"] = EVIDENCE_DEPTH_PRIORITY.get(evidence_depth, 1)
        bundle["max_file_lines"] = max(line_counts) if line_counts else 0
        bundle["total_file_lines"] = sum(line_counts)
        bundle["shared_test_files"] = select_shared_tests_for_bundle(bundle, contract_index)
        bundle["shared_test_file_count"] = len(bundle["shared_test_files"])
        if not bundle["baseline_members"] and bundle["bundle_type"] != "lineage_family":
            bundle["baseline_members"] = [
                {
                    "entity_key": bundle["primary_entity_key"],
                    "file_path": path_value,
                    "metadata": {"side": "baseline"},
                }
                for path_value in bundle["baseline_paths"]
            ]
        if not bundle["target_members"]:
            bundle["target_members"] = [
                {
                    "entity_key": bundle["primary_entity_key"],
                    "file_path": path_value,
                    "metadata": {"side": "target"},
                }
                for path_value in bundle["target_paths"]
            ]
        bundle["members"] = (
            list(bundle["target_members"])
            or list(bundle["baseline_members"])
        )
        bundles.append(bundle)

    bundles.sort(
        key=lambda item: (
            -item["linked_exception_count"],
            -BUNDLE_TYPE_PRIORITY.get(item["bundle_type"], 0),
            -item["max_file_lines"],
            item["evidence_depth_rank"],
            item["bundle_id"],
        )
    )
    for index, bundle in enumerate(bundles, start=1):
        bundle["priority_rank"] = index
        bundle["priority_score"] = float(
            item_score(bundle)
        )
    return bundles


def item_score(bundle: dict) -> int:
    return (
        int(bundle.get("linked_exception_count") or 0) * 100000
        + BUNDLE_TYPE_PRIORITY.get(str(bundle.get("bundle_type") or ""), 0) * 10000
        + int(bundle.get("max_file_lines") or 0) * 10
        + (10 - int(bundle.get("evidence_depth_rank") or 0))
    )


def build_priority_families(bundles: list[dict]) -> list[dict]:
    families: dict[str, dict] = {}
    for bundle in bundles:
        family_path = (
            (bundle["baseline_paths"][0] if bundle["baseline_paths"] else None)
            or (bundle["target_paths"][0] if bundle["target_paths"] else None)
            or str(bundle["bundle_id"])
        )
        record = families.setdefault(
            family_path,
            {
                "family_path": family_path,
                "bundle_count": 0,
                "linked_exception_count": 0,
                "max_file_lines": 0,
                "bundle_ids": [],
            },
        )
        record["bundle_count"] += 1
        record["linked_exception_count"] += int(bundle.get("linked_exception_count") or 0)
        record["max_file_lines"] = max(record["max_file_lines"], int(bundle.get("max_file_lines") or 0))
        record["bundle_ids"].append(bundle["bundle_id"])
    ranked = sorted(
        families.values(),
        key=lambda item: (
            -item["linked_exception_count"],
            -item["bundle_count"],
            -item["max_file_lines"],
            item["family_path"],
        ),
    )
    for index, item in enumerate(ranked, start=1):
        item["priority_rank"] = index
    return ranked


def insert_bundle_map_results(
    conn: sqlite3.Connection,
    *,
    run_id: str,
    bundles: list[dict],
) -> None:
    for bundle_index, bundle in enumerate(bundles, start=1):
        coverage_bundle_id = f"{run_id}:bundle:{bundle_index}"
        conn.execute(
            """
            INSERT INTO coverage_bundles (
                coverage_bundle_id, run_id, bundle_id, bundle_type, source_file_path,
                primary_entity_key, baseline_ref, target_ref, linked_exception_count,
                scenario_class, differential_replay_status, coverage_gate_status,
                review_note, priority_rank, priority_score, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                coverage_bundle_id,
                run_id,
                bundle["bundle_id"],
                bundle["bundle_type"],
                bundle.get("source_file_path"),
                bundle.get("primary_entity_key"),
                bundle.get("baseline_ref"),
                bundle.get("target_ref"),
                bundle.get("linked_exception_count"),
                bundle.get("scenario_class"),
                bundle.get("differential_replay_status"),
                bundle.get("coverage_gate_status"),
                bundle.get("review_note"),
                bundle.get("priority_rank"),
                bundle.get("priority_score"),
                stable_json(
                    {
                        "baseline_paths": bundle.get("baseline_paths", []),
                        "target_paths": bundle.get("target_paths", []),
                        "linked_exception_keys": bundle.get("linked_exception_keys", []),
                        "shared_test_files": bundle.get("shared_test_files", []),
                        "evidence_depth": bundle.get("evidence_depth"),
                        "disposition_counts": bundle.get("disposition_counts", {}),
                        "max_file_lines": bundle.get("max_file_lines"),
                        "total_file_lines": bundle.get("total_file_lines"),
                    }
                ),
            ),
        )
        for member_index, member in enumerate(bundle.get("members", []), start=1):
            conn.execute(
                """
                INSERT INTO coverage_bundle_members (
                    coverage_bundle_member_id, run_id, coverage_bundle_id, entity_key,
                    file_path, line_start, line_end, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:bundle_member:{bundle_index}:{member_index}",
                    run_id,
                    coverage_bundle_id,
                    member.get("entity_key"),
                    member["file_path"],
                    member.get("line_start"),
                    member.get("line_end"),
                    stable_json(member.get("metadata", {})),
                ),
            )
        for link_index, test_file_path in enumerate(bundle.get("shared_test_files", []), start=1):
            conn.execute(
                """
                INSERT INTO coverage_test_links (
                    coverage_test_link_id, run_id, coverage_bundle_id, test_file_path, source, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:bundle_test_link:{bundle_index}:{link_index}",
                    run_id,
                    coverage_bundle_id,
                    test_file_path,
                    "bundle_map_heuristic",
                    stable_json({}),
                ),
            )
        for link_index, exception_record in enumerate(bundle.get("linked_exceptions", []), start=1):
            conn.execute(
                """
                INSERT INTO coverage_exception_links (
                    coverage_exception_link_id, run_id, coverage_bundle_id, exception_id,
                    exception_key, audit_layer, entity_key, curation_status, disposition, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    f"{run_id}:bundle_exception_link:{bundle_index}:{link_index}",
                    run_id,
                    coverage_bundle_id,
                    exception_record["exception_id"],
                    exception_record.get("exception_key"),
                    exception_record["audit_layer"],
                    exception_record["entity_key"],
                    exception_record["curation_status"],
                    exception_record.get("disposition"),
                    stable_json(
                        {
                            "severity": exception_record.get("severity"),
                            "exception_type": exception_record.get("exception_type"),
                            "reason": exception_record.get("reason"),
                        }
                    ),
                ),
            )


def build_bundle_manifest(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    closeout_summary: dict,
    bundles: list[dict],
    priority_families: list[dict],
) -> dict:
    return {
        "run_id": run_id,
        "mode": "bundle_map",
        "repo_root": str(repo_root),
        "status": "completed",
        "closeout_run_id": closeout_summary["run_id"],
        "baseline": closeout_summary["baseline"],
        "target": closeout_summary["target"],
        "bundle_count": len(bundles),
        "linked_exception_count": sum(bundle["linked_exception_count"] for bundle in bundles),
        "priority_family_count": len(priority_families),
        "bundles": [
            {
                "bundle_id": bundle["bundle_id"],
                "bundle_type": bundle["bundle_type"],
                "primary_entity_key": bundle["primary_entity_key"],
                "baseline_ref": bundle["baseline_ref"],
                "target_ref": bundle["target_ref"],
                "baseline_paths": bundle["baseline_paths"],
                "target_paths": bundle["target_paths"],
                "linked_exception_keys": bundle["linked_exception_keys"],
                "shared_test_files": bundle["shared_test_files"],
                "shared_test_file_count": bundle["shared_test_file_count"],
                "scenario_class": bundle["scenario_class"],
                "differential_replay_status": bundle["differential_replay_status"],
                "coverage_gate_status": bundle["coverage_gate_status"],
                "review_note": bundle["review_note"],
                "priority_rank": bundle["priority_rank"],
                "priority_score": bundle["priority_score"],
                "linked_exception_count": bundle["linked_exception_count"],
                "evidence_depth": bundle["evidence_depth"],
                "source_file": bundle.get("source_file_path"),
                "members": bundle["members"],
                "baseline_members": bundle.get("baseline_members", []),
                "target_members": bundle.get("target_members", []),
            }
            for bundle in bundles
        ],
        "priority_families": priority_families,
    }


def build_bundle_map_summary(
    *,
    run_id: str,
    repo_root: pathlib.Path,
    closeout_summary: dict,
    bundle_manifest_path: str,
    bundles: list[dict],
    priority_families: list[dict],
) -> dict:
    disposition_counts: Counter[str] = Counter()
    for bundle in bundles:
        for disposition, count in (bundle.get("disposition_counts") or {}).items():
            disposition_counts[str(disposition)] += int(count)
    return {
        "run_id": run_id,
        "mode": "bundle_map",
        "repo_root": str(repo_root),
        "closeout_run_id": closeout_summary["run_id"],
        "baseline": closeout_summary["baseline"],
        "target": closeout_summary["target"],
        "status": "completed",
        "bundle_count": len(bundles),
        "linked_exception_count": sum(bundle["linked_exception_count"] for bundle in bundles),
        "priority_family_count": len(priority_families),
        "bundle_manifest_path": bundle_manifest_path,
        "bundle_type_counts": dict(sorted(Counter(bundle["bundle_type"] for bundle in bundles).items())),
        "disposition_counts": dict(sorted(disposition_counts.items())),
        "scenario_class_counts": dict(sorted(Counter(bundle["scenario_class"] for bundle in bundles).items())),
        "evidence_depth_counts": dict(sorted(Counter(bundle["evidence_depth"] for bundle in bundles).items())),
        "priority_bundles": [
            {
                "bundle_id": bundle["bundle_id"],
                "bundle_type": bundle["bundle_type"],
                "priority_rank": bundle["priority_rank"],
                "linked_exception_count": bundle["linked_exception_count"],
                "shared_test_file_count": bundle["shared_test_file_count"],
                "source_file_path": bundle.get("source_file_path"),
            }
            for bundle in bundles[:25]
        ],
        "priority_families": priority_families[:25],
    }


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


def load_curation_catalog(repo_root: pathlib.Path) -> list[dict]:
    entries: list[dict] = []
    for catalog_path in list_curation_paths(repo_root):
        payload = json.loads(catalog_path.read_text(encoding="utf-8"))
        raw_curations = ensure_list(payload.get("curations"))
        for index, record in enumerate(raw_curations, start=1):
            exception_key = record.get("exception_key")
            if not exception_key:
                raise RuntimeError(f"Curation entry {index} in {catalog_path} is missing `exception_key`")
            entries.append(
                {
                    "catalog_path": relative_label(catalog_path, repo_root),
                    "catalog_id": record.get("id") or f"curation_{index}",
                    "exception_key": str(exception_key),
                    "curation_status": str(record.get("curation_status") or "accepted"),
                    "disposition": record.get("disposition"),
                    "owner": record.get("owner"),
                    "review_note": record.get("review_note"),
                }
            )
    return entries


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


def selector_payload_to_surface_entity_key(selector_kind: str, payload: dict) -> str | None:
    value = payload.get("value")
    if selector_kind == "symbol" and value:
        return f"symbol::{value}"
    if selector_kind == "setMethod" and value and payload.get("signature"):
        return f"setMethod::{value}::{payload['signature']}"
    if selector_kind == "setClass" and value:
        return f"setClass::{value}"
    if selector_kind == "setGeneric" and value:
        return f"setGeneric::{value}"
    return None


def manifest_semantic_overlap_targets(
    conn: sqlite3.Connection,
    *,
    run_ids: list[str],
    include_resolved_manifest: bool = False,
) -> set[str]:
    if not run_ids:
        return set()

    records = fetch_exception_records(
        conn,
        run_ids=run_ids,
        open_only=not include_resolved_manifest,
    )
    manifest_candidates = [
        record
        for record in records
        if record["exception_type"] == "semantic_mismatch"
        and (
            not include_resolved_manifest
            or str(record.get("curation_status") or "") != "rejected"
        )
    ]
    if not manifest_candidates:
        return set()

    entry_keys_by_run: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for record in manifest_candidates:
        parsed = parse_manifest_exception_entity_key(str(record["entity_key"]))
        if parsed is None:
            continue
        entry_keys_by_run[str(record["run_id"])].add(parsed)

    if not entry_keys_by_run:
        return set()

    overlap_targets: set[str] = set()
    for run_id, entry_keys in entry_keys_by_run.items():
        rows = conn.execute(
            """
            SELECT manifest_path, entry_id, selector_kind, selector_payload_json
            FROM manifest_entries
            WHERE run_id = ?
            """,
            (run_id,),
        ).fetchall()
        for row in rows:
            key = (str(row["manifest_path"]), str(row["entry_id"]))
            if key not in entry_keys:
                continue
            payload = parse_json_field(row["selector_payload_json"]) or {}
            surface_key = selector_payload_to_surface_entity_key(str(row["selector_kind"]), payload)
            if surface_key:
                overlap_targets.add(surface_key)
    return overlap_targets


def auto_curate_redundant_surface_manifest_overlaps(
    conn: sqlite3.Connection,
    *,
    run_ids: list[str],
    include_resolved_manifest: bool = False,
) -> list[dict]:
    overlap_targets = manifest_semantic_overlap_targets(
        conn,
        run_ids=run_ids,
        include_resolved_manifest=include_resolved_manifest,
    )
    if not overlap_targets:
        return []

    surface_records = [
        record
        for record in fetch_exception_records(
            conn,
            run_ids=run_ids,
            curation_statuses=["candidate"],
        )
        if record["exception_type"] == "surface_definition_drift"
        and record["entity_key"] in overlap_targets
    ]

    updates: list[dict] = []
    for record in surface_records:
        update_exception_record(
            conn,
            exception_id=str(record["exception_id"]),
            curation_status=REDUNDANT_SURFACE_MANIFEST_RULE["curation_status"],
            disposition=REDUNDANT_SURFACE_MANIFEST_RULE["disposition"],
            review_note=REDUNDANT_SURFACE_MANIFEST_RULE["note"],
            curation_rule=REDUNDANT_SURFACE_MANIFEST_RULE["rule_name"],
        )
        updates.append(
            {
                "exception_id": record["exception_id"],
                "exception_key": record.get("exception_key"),
                "rule_name": REDUNDANT_SURFACE_MANIFEST_RULE["rule_name"],
                "curation_status": REDUNDANT_SURFACE_MANIFEST_RULE["curation_status"],
                "disposition": REDUNDANT_SURFACE_MANIFEST_RULE["disposition"],
            }
        )
    return updates


def apply_curation_catalog(
    conn: sqlite3.Connection,
    *,
    repo_root: pathlib.Path,
    run_ids: list[str],
) -> list[dict]:
    catalog_entries = load_curation_catalog(repo_root)
    if not catalog_entries:
        return []

    open_records = fetch_exception_records(conn, run_ids=run_ids, open_only=True)
    records_by_key: dict[str, list[dict]] = defaultdict(list)
    for record in open_records:
        if record.get("exception_key"):
            records_by_key[str(record["exception_key"])].append(record)

    updates: list[dict] = []
    for entry in catalog_entries:
        matching_records = records_by_key.get(entry["exception_key"], [])
        if not matching_records:
            continue
        curation_rule = f"catalog:{entry['catalog_path']}:{entry['catalog_id']}"
        for record in matching_records:
            update_exception_record(
                conn,
                exception_id=str(record["exception_id"]),
                curation_status=str(entry["curation_status"]),
                disposition=entry.get("disposition"),
                owner=entry.get("owner"),
                review_note=entry.get("review_note"),
                curation_rule=curation_rule,
            )
            updates.append(
                {
                    "exception_id": record["exception_id"],
                    "exception_key": record.get("exception_key"),
                    "curation_rule": curation_rule,
                    "curation_status": entry["curation_status"],
                    "disposition": entry.get("disposition"),
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
        library_paths = dedupe_paths(
            discover_project_library_paths(repo_root) +
            discover_project_library_paths(target_source["path"])
        )
        contract_results: list[dict] = []
        if args.execute:
            for file_record in contract_files:
                execution = run_contract_test_file(
                    target_source["path"],
                    file_record["absolute_path"],
                    load_package=bool(args.load_package),
                    library_paths=library_paths,
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


def coverage_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"coverage-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    side = args.side
    default_ref = "main" if side == "baseline" else "WORKTREE"
    source = resolve_side_source(
        repo_root,
        side=side,
        ref=args.source_ref,
        path=args.source_path,
        default_ref=default_ref,
        allow_worktree=True,
        archive_paths=None,
    )

    bundle_manifest_path, bundles = load_coverage_bundle_manifest(repo_root, args.bundle_manifest)
    bundles = select_coverage_bundles(
        bundles,
        bundle_ids=ensure_list(getattr(args, "bundle_id", None)),
        limit=getattr(args, "limit", None),
    )
    explicit_test_files = dedupe_paths(ensure_list(args.test_file))
    bundle_test_files = [] if getattr(args, "explicit_tests_only", False) else collect_bundle_test_files(bundles)
    if (
        not getattr(args, "explicit_tests_only", False)
        and not ensure_list(getattr(args, "bundle_id", None))
        and getattr(args, "limit", None) is None
    ):
        bundle_test_files = dedupe_paths(
            bundle_test_files + discover_shared_compare_test_files(repo_root)
        )
    requested_test_files = dedupe_paths(explicit_test_files + bundle_test_files)
    if bundles and not requested_test_files:
        raise RuntimeError("No shared coverage tests selected for the requested bundles")

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "coverage",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "side": side,
            "source": source["label"],
            "bundle_manifest_path": bundle_manifest_path,
            "requested_test_file_count": len(requested_test_files),
            "bundle_count": len(bundles),
        },
    )

    workspace = None
    try:
        test_source_root = repo_root if any(bundle.get("shared_test_files") for bundle in bundles) else source["path"]
        active_root = source["path"]
        if test_source_root != source["path"]:
            workspace = materialize_coverage_workspace(
                source["path"],
                test_source_root=test_source_root,
            )
            active_root = pathlib.Path(workspace.name)

        test_paths = list_test_file_paths(active_root, requested_test_files or None)
        if not test_paths:
            raise RuntimeError("No testthat files found for coverage capture")
        selected_test_files = [relative_label(test_path, active_root) for test_path in test_paths]
        explicit_test_paths = (
            list_test_file_paths(active_root, explicit_test_files)
            if explicit_test_files
            else []
        )
        command_test_files = [relative_label(test_path, active_root) for test_path in explicit_test_paths]
        library_paths = dedupe_paths(
            discover_project_library_paths(repo_root) +
            discover_project_library_paths(source["path"])
        )
        coverage_payload = run_coverage_capture(
            active_root,
            test_files=test_paths,
            library_paths=library_paths,
        )
        bundle_results = build_coverage_bundle_results(
            prepare_coverage_bundles_for_side(bundles, side),
            coverage_payload["files"],
        )
    finally:
        if workspace is not None:
            workspace.cleanup()
        source["cleanup"]()

    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=source["label"] if side == "baseline" else None,
            target_ref=source["label"] if side == "target" else None,
            mode="coverage",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "coverage",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        insert_coverage_results(
            conn,
            run_id=run_id,
            side=side,
            scope=side,
            source_label=source["label"],
            bundle_manifest_path=bundle_manifest_path,
            coverage_payload=coverage_payload,
            bundles=prepare_coverage_bundles_for_side(bundles, side),
            bundle_results=bundle_results,
            selected_test_files=selected_test_files,
            command_test_files=command_test_files,
        )
        summary = build_coverage_summary(
            run_id=run_id,
            repo_root=repo_root,
            side=side,
            source_label=source["label"],
            bundle_manifest_path=bundle_manifest_path,
            selected_test_files=selected_test_files,
            coverage_payload=coverage_payload,
            bundle_results=bundle_results,
        )
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
        materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_coverage_summary_markdown(summary))
    write_json(paths["coverage_json"], summary)
    write_text(paths["coverage_md"], render_coverage_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "coverage",
            "status": summary["status"],
            "run_id": run_id,
            "side": side,
            "source": summary["source"],
            "tool_status": summary["tool_status"],
            "test_file_count": summary["test_file_count"],
            "bundle_count": summary["bundle_count"],
            "file_count": summary["file_count"],
            "line_percent": summary["line_percent"],
        },
    )
    print(json.dumps(summary))
    return 0


def coverage_compare_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"coverage-compare-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    baseline_source = resolve_side_source(
        repo_root,
        side="baseline",
        ref=args.baseline_ref,
        path=args.baseline_path,
        default_ref="main",
        allow_worktree=True,
        archive_paths=None,
    )
    target_source = resolve_side_source(
        repo_root,
        side="target",
        ref=args.target_ref,
        path=args.target_path,
        default_ref="WORKTREE",
        allow_worktree=True,
        archive_paths=None,
    )

    bundle_manifest_path, bundles = load_coverage_bundle_manifest(repo_root, args.bundle_manifest)
    bundles = select_coverage_bundles(
        bundles,
        bundle_ids=ensure_list(args.bundle_id),
        limit=args.limit,
    )
    explicit_test_files = dedupe_paths(ensure_list(args.test_file))
    bundle_test_files = [] if getattr(args, "explicit_tests_only", False) else collect_bundle_test_files(bundles)
    if (
        not getattr(args, "explicit_tests_only", False)
        and not ensure_list(args.bundle_id)
        and args.limit is None
    ):
        bundle_test_files = dedupe_paths(
            bundle_test_files + discover_shared_compare_test_files(target_source["path"])
        )
    requested_test_files = dedupe_paths(explicit_test_files + bundle_test_files)
    if bundles and not requested_test_files:
        raise RuntimeError("No shared coverage tests selected for the requested bundles")

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "coverage_compare",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "baseline": baseline_source["label"],
            "target": target_source["label"],
            "bundle_manifest_path": bundle_manifest_path,
            "requested_test_file_count": len(requested_test_files),
            "bundle_count": len(bundles),
        },
    )

    baseline_workspace = None
    target_workspace = None
    try:
        target_workspace = materialize_coverage_workspace(
            target_source["path"],
            test_source_root=target_source["path"],
        )
        baseline_workspace = materialize_coverage_workspace(
            baseline_source["path"],
            test_source_root=target_source["path"],
        )
        target_workspace_root = pathlib.Path(target_workspace.name)
        baseline_workspace_root = pathlib.Path(baseline_workspace.name)

        target_test_paths = list_test_file_paths(target_workspace_root, requested_test_files or None)
        if not target_test_paths:
            raise RuntimeError("No testthat files found for comparative coverage capture")
        baseline_test_paths = list_test_file_paths(baseline_workspace_root, requested_test_files or None)
        if not baseline_test_paths:
            raise RuntimeError("No baseline testthat files found for comparative coverage capture")

        selected_test_files = [relative_label(test_path, target_workspace_root) for test_path in target_test_paths]
        explicit_target_test_paths = (
            list_test_file_paths(target_workspace_root, explicit_test_files)
            if explicit_test_files
            else []
        )
        command_test_files = [
            relative_label(test_path, target_workspace_root)
            for test_path in explicit_target_test_paths
        ]
        library_paths = dedupe_paths(
            discover_project_library_paths(repo_root) +
            discover_project_library_paths(baseline_source["path"]) +
            discover_project_library_paths(target_source["path"])
        )
        baseline_env_overrides = {}
        target_env_overrides = {}
        baseline_fixture = os.getenv("FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_BASELINE", "")
        target_fixture = os.getenv("FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_TARGET", "")
        if baseline_fixture:
            baseline_env_overrides["FIDELITY_COVERAGE_FIXTURE_JSON"] = baseline_fixture
        if target_fixture:
            target_env_overrides["FIDELITY_COVERAGE_FIXTURE_JSON"] = target_fixture

        baseline_payload = run_coverage_capture(
            baseline_workspace_root,
            test_files=baseline_test_paths,
            library_paths=library_paths,
            env_overrides=baseline_env_overrides or None,
        )
        target_payload = run_coverage_capture(
            target_workspace_root,
            test_files=target_test_paths,
            library_paths=library_paths,
            env_overrides=target_env_overrides or None,
        )

        baseline_bundles = prepare_coverage_bundles_for_side(bundles, "baseline")
        target_bundles = prepare_coverage_bundles_for_side(bundles, "target")
        baseline_bundle_results = build_coverage_bundle_results(baseline_bundles, baseline_payload["files"])
        target_bundle_results = build_coverage_bundle_results(target_bundles, target_payload["files"])
        compare_bundle_results = build_coverage_compare_bundle_results(
            bundles,
            baseline_bundle_results=baseline_bundle_results,
            target_bundle_results=target_bundle_results,
            baseline_capture_status=baseline_payload["status"],
            target_capture_status=target_payload["status"],
        )
    finally:
        if baseline_workspace is not None:
            baseline_workspace.cleanup()
        if target_workspace is not None:
            target_workspace.cleanup()
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
            mode="coverage_compare",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "coverage_compare",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        insert_coverage_results(
            conn,
            run_id=run_id,
            side="baseline",
            scope="baseline",
            source_label=baseline_source["label"],
            bundle_manifest_path=bundle_manifest_path,
            coverage_payload=baseline_payload,
            bundles=baseline_bundles,
            bundle_results=baseline_bundle_results,
            selected_test_files=selected_test_files,
            command_test_files=command_test_files,
        )
        insert_coverage_results(
            conn,
            run_id=run_id,
            side="target",
            scope="target",
            source_label=target_source["label"],
            bundle_manifest_path=bundle_manifest_path,
            coverage_payload=target_payload,
            bundles=target_bundles,
            bundle_results=target_bundle_results,
            selected_test_files=selected_test_files,
            command_test_files=command_test_files,
        )
        insert_coverage_compare_results(
            conn,
            run_id=run_id,
            bundles=bundles,
            compare_bundle_results=compare_bundle_results,
        )
        summary = build_coverage_compare_summary(
            run_id=run_id,
            repo_root=repo_root,
            baseline_source_label=baseline_source["label"],
            target_source_label=target_source["label"],
            test_source_label=target_source["label"],
            bundle_manifest_path=bundle_manifest_path,
            selected_test_files=selected_test_files,
            baseline_payload=baseline_payload,
            target_payload=target_payload,
            compare_bundle_results=compare_bundle_results,
        )
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
        materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_coverage_compare_summary_markdown(summary))
    write_json(paths["coverage_compare_json"], summary)
    write_text(paths["coverage_compare_md"], render_coverage_compare_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "coverage_compare",
            "status": summary["status"],
            "run_id": run_id,
            "baseline": baseline_source["label"],
            "target": target_source["label"],
            "test_file_count": summary["test_file_count"],
            "bundle_count": summary["bundle_count"],
            "baseline_line_percent": summary["baseline"]["line_percent"],
            "target_line_percent": summary["target"]["line_percent"],
        },
    )
    print(json.dumps(summary))
    return 0


def coverage_evidence_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"coverage-evidence-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"

    bundle_manifest_arg = args.bundle_manifest or relative_label(paths["bundles_json"], repo_root)
    bundle_manifest_path, bundles = load_coverage_bundle_manifest(repo_root, bundle_manifest_arg)
    manifest_bundle_count = len(bundles)
    bundles = select_coverage_bundles(
        bundles,
        bundle_ids=ensure_list(args.bundle_id),
        limit=args.limit,
    )
    package_gate_applies = len(bundles) == manifest_bundle_count
    coverage_compare_summary = load_coverage_compare_summary(
        repo_root,
        coverage_compare_run_id=args.coverage_compare_run_id,
        audit_dir=args.audit_dir,
        compare_path=args.coverage_compare_path,
    )
    acceptance_policy = load_coverage_acceptance_policy(repo_root)
    evaluated_bundles = build_coverage_evidence_bundles(
        bundles,
        coverage_compare_summary=coverage_compare_summary,
        acceptance_policy=acceptance_policy,
    )

    compare_target_status = str((coverage_compare_summary.get("target") or {}).get("status") or "missing")
    exceptions: list[dict] = []
    if compare_target_status == "completed":
        for bundle in evaluated_bundles:
            record = build_coverage_evidence_exception(run_id, bundle)
            if record is not None:
                exceptions.append(record)
        package_record = build_package_coverage_exception(
            run_id,
            baseline_ref=(coverage_compare_summary.get("baseline") or {}).get("source"),
            target_ref=(coverage_compare_summary.get("target") or {}).get("source"),
            target_line_percent=coerce_float((coverage_compare_summary.get("target") or {}).get("line_percent")),
            threshold_pct=float(acceptance_policy["package_line_coverage_target_pct"]),
            enforce_gate=package_gate_applies,
        )
        if package_record is not None:
            exceptions.append(package_record)

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "coverage_evidence",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "bundle_manifest_path": bundle_manifest_path,
            "coverage_compare_run_id": coverage_compare_summary.get("run_id"),
            "bundle_count": len(evaluated_bundles),
        },
    )

    conn = connect_db(paths["db_path"])
    try:
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=(coverage_compare_summary.get("baseline") or {}).get("source"),
            target_ref=(coverage_compare_summary.get("target") or {}).get("source"),
            mode="coverage_evidence",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "coverage_evidence",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        insert_coverage_evidence_results(
            conn,
            run_id=run_id,
            bundles=evaluated_bundles,
            exceptions=exceptions,
        )
        summary = build_coverage_evidence_summary(
            run_id=run_id,
            repo_root=repo_root,
            bundle_manifest_path=bundle_manifest_path,
            coverage_compare_summary=coverage_compare_summary,
            acceptance_policy=acceptance_policy,
            evaluated_bundles=evaluated_bundles,
            exceptions=exceptions,
            package_gate_applies=package_gate_applies,
        )
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
        materialize_exception_reports(conn, repo_root=repo_root, paths=paths)
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_coverage_evidence_markdown(summary))
    write_json(paths["coverage_evidence_json"], summary)
    write_text(paths["coverage_evidence_md"], render_coverage_evidence_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "coverage_evidence",
            "status": summary["status"],
            "run_id": run_id,
            "bundle_count": summary["bundle_count"],
            "low_coverage_bundle_count": summary["low_coverage_bundle_count"],
            "regression_bundle_count": summary["regression_bundle_count"],
            "exception_count": summary["exception_count"],
        },
    )
    print(json.dumps(summary))
    return 0


def bundle_map_command(args: argparse.Namespace) -> int:
    repo_root = pathlib.Path(args.repo_root).resolve()
    paths = bootstrap_store(repo_root, args.audit_dir, emit_event=False)
    run_id = args.run_id or f"bundle-map-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    closeout_summary = load_closeout_summary(repo_root, args.closeout_run_id, audit_dir=args.audit_dir)
    component_run_ids = sorted(set(closeout_summary.get("component_runs", {}).values()))
    if not component_run_ids:
        raise RuntimeError("Closeout summary does not contain component runs")

    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "bundle_map",
            "status": "started",
            "run_id": run_id,
            "repo_root": str(repo_root),
            "closeout_run_id": closeout_summary["run_id"],
            "component_run_count": len(component_run_ids),
        },
    )

    conn = connect_db(paths["db_path"])
    try:
        manifest_selector_index = build_manifest_selector_index(conn, component_run_ids)
        inventory_entity_index = build_inventory_entity_index(conn, component_run_ids)
        exception_records = fetch_exception_records(
            conn,
            run_ids=component_run_ids,
            curation_statuses=["auto_curated", "accepted"],
        )
        contract_index = build_contract_search_index(repo_root)
        bundles = build_bundle_groups(
            repo_root=repo_root,
            baseline_ref=str(closeout_summary["baseline"]),
            target_ref=str(closeout_summary["target"]),
            exception_records=exception_records,
            manifest_selector_index=manifest_selector_index,
            inventory_entity_index=inventory_entity_index,
            contract_index=contract_index,
        )
        priority_families = build_priority_families(bundles)
        bundle_manifest_path = relative_label(paths["bundles_json"], repo_root)
        bundle_manifest = build_bundle_manifest(
            run_id=run_id,
            repo_root=repo_root,
            closeout_summary=closeout_summary,
            bundles=bundles,
            priority_families=priority_families,
        )
        summary = build_bundle_map_summary(
            run_id=run_id,
            repo_root=repo_root,
            closeout_summary=closeout_summary,
            bundle_manifest_path=bundle_manifest_path,
            bundles=bundles,
            priority_families=priority_families,
        )
        insert_audit_run(
            conn,
            run_id=run_id,
            repo_root=repo_root,
            baseline_ref=str(closeout_summary["baseline"]),
            target_ref=str(closeout_summary["target"]),
            mode="bundle_map",
            status="running",
            summary={
                "run_id": run_id,
                "mode": "bundle_map",
                "repo_root": str(repo_root),
                "status": "running",
            },
        )
        insert_bundle_map_results(conn, run_id=run_id, bundles=bundles)
        update_audit_run(conn, run_id=run_id, status=summary["status"], summary=summary)
        conn.commit()
    finally:
        conn.close()

    write_json(paths["summary_json"], summary)
    write_text(paths["summary_md"], render_bundle_map_summary_markdown(summary))
    write_json(paths["bundles_json"], bundle_manifest)
    write_text(paths["bundles_md"], render_bundle_map_summary_markdown(summary))
    append_jsonl(
        paths["events_path"],
        {
            "timestamp": now_iso(),
            "phase": "bundle_map",
            "status": summary["status"],
            "run_id": run_id,
            "closeout_run_id": closeout_summary["run_id"],
            "bundle_count": summary["bundle_count"],
            "linked_exception_count": summary["linked_exception_count"],
            "priority_family_count": summary["priority_family_count"],
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
        auto_curated_updates.extend(
            auto_curate_redundant_surface_manifest_overlaps(conn, run_ids=component_run_ids)
        )
        catalog_curation_updates = apply_curation_catalog(
            conn,
            repo_root=repo_root,
            run_ids=component_run_ids,
        )
        auto_curated_updates.extend(
            auto_curate_redundant_surface_manifest_overlaps(
                conn,
                run_ids=component_run_ids,
                include_resolved_manifest=True,
            )
        )
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
        closeout_summary["catalog_curated_count"] = len(catalog_curation_updates)
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

    coverage_parser = subparsers.add_parser("coverage", help="Capture line coverage for a source tree and optional bundle manifest")
    coverage_parser.add_argument("--repo-root", default=".")
    coverage_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    coverage_parser.add_argument("--side", default="target", choices=["baseline", "target"])
    coverage_parser.add_argument("--source-ref")
    coverage_parser.add_argument("--source-path")
    coverage_parser.add_argument("--test-file", action="append")
    coverage_parser.add_argument("--explicit-tests-only", action="store_true")
    coverage_parser.add_argument("--bundle-manifest")
    coverage_parser.add_argument("--bundle-id", action="append")
    coverage_parser.add_argument("--limit", type=int)
    coverage_parser.add_argument("--run-id")
    coverage_parser.set_defaults(func=coverage_command)

    coverage_compare_parser = subparsers.add_parser("coverage-compare", help="Capture comparative bundle coverage against baseline and target using the same test subset")
    coverage_compare_parser.add_argument("--repo-root", default=".")
    coverage_compare_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    coverage_compare_parser.add_argument("--baseline-ref")
    coverage_compare_parser.add_argument("--baseline-path")
    coverage_compare_parser.add_argument("--target-ref")
    coverage_compare_parser.add_argument("--target-path")
    coverage_compare_parser.add_argument("--test-file", action="append")
    coverage_compare_parser.add_argument("--explicit-tests-only", action="store_true")
    coverage_compare_parser.add_argument("--bundle-manifest")
    coverage_compare_parser.add_argument("--bundle-id", action="append")
    coverage_compare_parser.add_argument("--limit", type=int)
    coverage_compare_parser.add_argument("--run-id")
    coverage_compare_parser.set_defaults(func=coverage_compare_command)

    coverage_evidence_parser = subparsers.add_parser("coverage-evidence", help="Materialize bundle coverage evidence, deltas, and low-coverage exceptions")
    coverage_evidence_parser.add_argument("--repo-root", default=".")
    coverage_evidence_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    coverage_evidence_parser.add_argument("--bundle-manifest")
    coverage_evidence_parser.add_argument("--coverage-compare-run-id")
    coverage_evidence_parser.add_argument("--coverage-compare-path")
    coverage_evidence_parser.add_argument("--bundle-id", action="append")
    coverage_evidence_parser.add_argument("--limit", type=int)
    coverage_evidence_parser.add_argument("--run-id")
    coverage_evidence_parser.set_defaults(func=coverage_evidence_command)

    bundle_map_parser = subparsers.add_parser("bundle-map", help="Map curated parity evidence onto stable coverage bundles")
    bundle_map_parser.add_argument("--repo-root", default=".")
    bundle_map_parser.add_argument("--audit-dir", default=DEFAULT_AUDIT_DIR)
    bundle_map_parser.add_argument("--closeout-run-id")
    bundle_map_parser.add_argument("--run-id")
    bundle_map_parser.set_defaults(func=bundle_map_command)

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
