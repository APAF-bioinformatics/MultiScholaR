"""Bind PERF-108's inherited obligations to completed evidence, without promotion."""
import hashlib
import json
from pathlib import Path

directory = Path("docs/qa/workflow-execution/perf-108/closeout")
result = json.loads((directory / "result.json").read_text())
assert len(result["checks"]) == 2 and all(x["complete"] for x in result["checks"])
assert result["decision"]["science_outputs_and_measurements_valid"] is True
assert result["decision"]["status"] in ("held", "runtime_non_regression_passed")
ledger_path = Path("docs/qa/workflow-execution/perf-repaired-transfer-20261002.json")
applications = [x for x in json.loads(ledger_path.read_text())["requirement_applications"]
                if x["owner"] == "PERF-108"]
assert len(applications) == 17
check_root = Path(result["checks"][0]["run"])


def evidence(path):
    path = Path(path)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


files = {
    "ledger": ledger_path,
    "protocol": directory / "protocol.md",
    "bindings": directory / "bindings.rds",
    "result": directory / "result.json",
    "members": directory / "members.csv",
    "source_reuse": directory / "source-reuse.json",
    "installed_tests": check_root / "reader-tests.rds",
    "generic_tests": check_root / "generic-tests.rds",
    "compatibility": check_root / "reader-operation.json",
    "rollback": check_root / "wfx-runtime-readback.json",
    "reused_rollup": result["reused_rollup"]["path"],
    "predecessor_acceptance": "docs/qa/workflow-execution/perf-106/repair-v2/acceptance.json",
    "cleanup": directory / "cleanup.json",
}
notes = {
    51: ("protocol bindings reused_rollup", "Exact retained and accepted artifacts, native oracles and populations frozen before timing; rollup reuse checked against source bindings and installed package hashes."),
    52: ("installed_tests generic_tests members", "Installed native parity covers empty, missing, zero, invalid values, duplicate identities, evidence, factors, ordering and annotations; full primary/wide outputs independently verified."),
    53: ("predecessor_acceptance source_reuse installed_tests generic_tests cleanup", "Reuse accepted process-owner lifecycle evidence for unchanged ownership/state components. Recheck changed runtime cancellation, errors, publication recovery and resource cleanup. No claim that old binary measurements qualify new code."),
    54: ("members reused_rollup rollback", "Bounded independent table/summary/ledger/design fingerprints, rollup figures and predecessor consumption. GUI activation and continuing-session rollout remain PERF-107."),
    58: ("installed_tests members reused_rollup", "Exact native scientific comparisons; no mocked scientific computation or relaxed tolerance."),
    59: ("source_reuse installed_tests", "Existing scoped reader, staged writer, checkpoint/history and execution owner retained; public C query ownership replaces the private relation bridge for the retained route."),
    60: ("members reused_rollup installed_tests rollback", "Declared scientific outputs and consumer handoffs verified. Physical and logical selection behavior checked; no whole-GUI claim."),
    61: ("result members", "Science, containment, measured non-regression and default eligibility remain separate. A failed timing pair does not become a pass because science or memory passed."),
    62: ("protocol result", "Pinned Linux x86_64/DuckDB 1.5.5 only; experimental routes remain off. Historical unsigned-extension and other platform holds are not globally cleared."),
    66: ("installed_tests generic_tests compatibility members source_reuse", "Changed runtime/codec/filter paths tested installed; public APIs/S4/history compared to accepted. Unchanged broader scientific and GUI components reuse predecessor acceptance without extending claims."),
    67: ("bindings result source_reuse members", "Exact source, installed binaries, inputs and outputs bound. Earlier failed campaigns remain intact; no full retry or threshold change."),
    68: ("installed_tests generic_tests rollback cleanup", "Input hashes, released scoped resources, managed-process teardown and accepted saved-output compatibility checked."),
    72: ("compatibility installed_tests members reused_rollup", "Public functions, formals, S4 classes/dispatch, axes, metadata, required outputs and native scientific values preserved."),
    73: ("installed_tests rollback", "Selection-aware physical exports and explicit unsupported-format refusals checked; predecessor reads new output and recognizes pending publication."),
    74: ("protocol result", "No GUI/default activation. PERF-107 and PERF-128 retain their respective gates."),
    78: ("bindings source_reuse rollback", "Accepted runtime remains installed and bound; original sources/history preserved. Rollback changes route/runtime selection, never resets the branch or deletes data."),
    79: ("rollback installed_tests", "Actual predecessor bounded readback preserves values/history and recognizes pending output; physical compatibility export remains available."),
}
for item in applications:
    keys, note = notes[item["line"]]
    item["status"] = "satisfied_for_bounded_runtime_scope"
    item["scope"] = note
    item["evidence"] = keys.split()

output = {
    "owner": "PERF-108",
    "case": "CASE-PERF-REL-001",
    "status": "held" if result["decision"]["status"] == "held" else "checks_passed_release_pending",
    "runtime_non_regression": result["decision"],
    "release_authorized_by_this_receipt": False,
    "default_activation": False,
    "requirements": applications,
    "evidence": {key: evidence(path) for key, path in files.items()},
    "limits": ["No GUI/real-study/cross-omics qualification", "No experimental P14 route promotion",
               "No automatic retry, threshold waiver or unsigned-extension approval"],
}
(directory / "acceptance.json").write_text(json.dumps(output, indent=2) + "\n")
