#!/usr/bin/env python3
"""Check authored WFX decision fixtures; this is not the production evaluator."""
import json
import math
from pathlib import Path


def finite(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value > 0)


def expected(case, required_pairs):
    # Missing evidence cannot establish capacity by vacuous all([]) truth.
    pairs = case.get("pairs")
    if not isinstance(pairs, list) or len(pairs) != required_pairs:
        return "held"
    for pair in pairs:
        if not isinstance(pair, dict):
            return "held"
        for arm in ("control", "candidate"):
            if not isinstance(pair.get(arm), dict) or "status" not in pair[arm]:
                return "held"
    if not case["domain_supported"]:
        return "unsupported"
    if any(not case[key] for key in ("science", "lifecycle", "integrity", "outputs", "absolute_limits")):
        return "rejected"
    if any(pair["candidate"]["status"] != "complete" for pair in pairs):
        return "rejected"
    if any(pair["control"]["status"] in ("science", "publication") for pair in pairs):
        return "rejected"
    goal = case["goal"]
    if (goal not in ("benefit", "integration_non_regression", "capacity_only")
            or case["open_holds"] or any(not case[key] for key in
            ("measurement_valid", "bindings_match", "runtime_eligible"))):
        return "held"
    if goal == "integration_non_regression" and not case["predecessor_benefit_bound"]:
        return "held"
    metrics = ("peak_bytes", "elapsed_seconds")
    if goal == "capacity_only":
        valid = all(pair["control"]["status"] in ("memory_guard", "time_guard")
                    and pair["control"].get("input_valid") is True
                    and all(finite(pair["candidate"].get(key)) for key in metrics)
                    for pair in pairs)
        return "capacity_only" if valid else "held"
    if any(pair[arm]["status"] != "complete"
           or not all(finite(pair[arm].get(key)) for key in metrics)
           for pair in pairs for arm in ("control", "candidate")):
        return "held"
    limit = 0.90 if goal == "benefit" else 1.05
    passed = [pair["candidate"]["peak_bytes"] / pair["control"]["peak_bytes"] <= limit
              and pair["candidate"]["elapsed_seconds"] / pair["control"]["elapsed_seconds"] <= 1.20
              for pair in pairs]
    return goal if all(passed) else ("inconclusive" if any(passed) else "rejected")


def main():
    repo = Path(__file__).resolve().parents[2]
    manifest = json.loads((repo / "tests/testdata/workflow-execution/qualification.json").read_text())
    fixtures = manifest["decision_fixtures"]
    assert fixtures["required_pairs"] == 2
    for case in fixtures["cases"]:
        assert expected(case, fixtures["required_pairs"]) == case["expected"], case["id"]
    allocation = manifest["allocations"]["WFX-006"]
    assert allocation["crossover_order_per_size"] == [
        "explicit native", "explicit qualified relational", "explicit qualified relational", "explicit native"]
    assert allocation["crossover_rows"] == [9600, 96000, 960000]
    assert allocation["crossover_members"] == 12
    assert allocation["crossover_seconds"] == 300
    assert allocation["final_members"] == 4
    assert allocation["final_seconds"] == 1800
    print(json.dumps({"status": "passed", "manifest_version": manifest["version"],
                      "decision_cases": len(fixtures["cases"]), "crossover_explicit": True,
                      "scope": "fixture arithmetic and allocation only; no runtime qualification"}))


if __name__ == "__main__":
    main()
