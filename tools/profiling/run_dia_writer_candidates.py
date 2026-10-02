"""Resume the fixed DIA-025 writer screen; reuse completed guarded members."""
import hashlib
import json
import statistics
import subprocess
from pathlib import Path

repo = Path(__file__).resolve().parents[2]
proof = repo.parent / "MultiScholaR-dia-proof-20260911-VhC1n5"
control = proof / "browser-smoke.XsaPvT/library"
candidate = proof / "browser-smoke.U0ZRfy/library"
record_path = candidate.parent / "p14-writer-screen.json"
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
record = json.loads(record_path.read_text())
assert record["status"] == "running"
assert all(sha(Path(path)) == value for path, value in record["sources"].items())
archive = record_path.with_name("p14-writer-screen.interrupted.json")
assert not archive.exists()
archive.write_bytes(record_path.read_bytes())
record["interruption_recovery"] = dict(previous_sha256=sha(archive), script_sha256=sha(Path(__file__)))


def flush():
    record_path.write_text(json.dumps(record, indent=2))


expected = json.loads((proof / "browser-smoke.CHnC1c/writer-comparison.json").read_text())["result"]
try:
    for arm, existing in [("P16", "a6M1rH"), ("P15", None), ("P14", None)]:
        library = control if arm == "P14" else candidate
        assert all(sha(library / path) == value for path, value in record["libraries"][str(library)].items())
        selector = f"gui-bulk-rollup:writer-r:screen:contrast:{arm}:2"
        if existing:
            root = proof / f"browser-smoke.{existing}"
            print("RECOVER", root.name, flush=True)
        else:
            before = set(proof.glob("browser-smoke.*"))
            print("START", selector, flush=True)
            result = subprocess.run(["bash", str(repo / "tools/profiling/run_workflow_checkpoint_proof.sh"),
                "--service", "reader-check", str(library), str(proof), str(repo), selector],
                cwd=repo, capture_output=True, text=True)
            created = set(proof.glob("browser-smoke.*")) - before
            assert len(created) == 1, result.stderr
            root = created.pop()
            if result.returncode:
                record["failed_member"] = root.name
                raise RuntimeError(result.stdout[-3000:] + result.stderr[-1000:])
        payload = json.loads((root / "writer-comparison.json").read_text())
        assert payload["mode"] == selector and payload["status"] == "complete" and payload["exact"]
        result = payload["result"]
        for field in ["output_hashes", "summary", "filtering"]:
            assert result[field] == expected[field], (root.name, field)
        props = dict(line.split("=", 1) for line in (root / "terminal-properties.txt").read_text().splitlines() if "=" in line)
        assert props["Result"] == "success" and props["MemorySwapPeak"] == "0"
        assert len(result["figures"]) == 12 and result["execution"]["intermediate_bytes_written"] == 0
        assert not list(root.rglob("*.tmp-*"))
        record["members"].append(dict(case="contrast", arm=arm, repetition=2, member=root.name,
            exit_code=0, seconds=result["complete_operation_seconds"], peak_bytes=int(props["MemoryPeak"]),
            exact=True, blocks=len(result["output_hashes"]), execution=result["execution"]))
        flush()
        print("PASS", root.name, result["complete_operation_seconds"], flush=True)
    record["medians"] = {arm: {case: {
        key: statistics.median(x[key] for x in record["members"] if x["arm"] == arm and x["case"] == case)
        for key in ["seconds", "peak_bytes"]} for case in ["primary", "contrast"]}
        for arm in ["P14", "P15", "P16"]}
    record["status"] = "complete"
except BaseException as error:
    record["status"] = "failed"
    record["error"] = str(error)
    raise
finally:
    flush()
