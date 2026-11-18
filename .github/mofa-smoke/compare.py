"""MOFA2 smoke lane: compare the result.json files of the runs pairwise (run-mofa-smoke).

Usage: python3 -I compare.py <directory holding mofa-smoke-*/result.json> <output directory>
Writes compare.md (also for the step summary) and compare.json. Exact: every bit of the factors,
weights, variance explained and last ELBO equal. Otherwise the agreement in significant digits:
floor(-log10(max |a - b| / max |a|)) over each block.
"""
import itertools
import json
import math
import pathlib
import sys

root, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
runs = {}
for path in sorted(root.glob("*/result.json")):
    runs[path.parent.name] = json.loads(path.read_text())


def values(result, block):
    model = result["model"]
    if block == "factors":
        return [float.fromhex(v) for v in model["factors"]["values_hex"]]
    if block == "weights":
        return [float.fromhex(v) for view in sorted(model["weights"]) for v in model["weights"][view]["values_hex"]]
    if block == "r2":
        return [float.fromhex(v) for v in model["r2_percent"]["values_hex"]]
    if block == "elbo":
        last = model.get("elbo_last_hex")
        return [float.fromhex(last)] if isinstance(last, str) else []
    raise ValueError(block)


def agreement(a, b):
    if len(a) != len(b):
        return {"exact": False, "note": f"lengths differ ({len(a)} and {len(b)})"}
    if all(x == y for x, y in zip(a, b)):
        return {"exact": True, "max_abs_diff": 0.0, "digits": "all"}
    diff = max(abs(x - y) for x, y in zip(a, b))
    scale = max(max(abs(x) for x in a), 1e-300)
    rel = diff / scale
    return {"exact": False, "max_abs_diff": diff, "max_rel_diff": rel,
            "digits": max(0, math.floor(-math.log10(rel))) if rel > 0 else "all"}


lines = ["### MOFA2 smoke: runs", "",
         "| run | trained | MOFA2 | mofapy2 | Python | numpy | h5py (HDF5) | factors in model / active | iterations | train s | total s | peak MiB | Z1 | Z2 | Z3 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for name, r in runs.items():
    if not r.get("trained"):
        lines.append(f"| {name} | no: {r.get('error')} |" + " |" * 13)
        continue
    pv, m, best = r["python_versions"], r["model"], r["planted"]["best"]
    lines.append(
        f"| {name} | yes | {r['r_versions']['MOFA2']} | {pv['mofapy2']} | {pv['python']} | {pv['numpy']} | "
        f"{pv['h5py']} ({pv['hdf5_h5py']}) | {m['factors_in_model']} / {m['factors_active_1pct']} | {m['iterations']} | "
        f"{r['seconds_train']} | {r['seconds_total']} | {r.get('peak_memory_mib')} | "
        + " | ".join(f"{best[z]['factor']} r={best[z]['r']}" for z in ("Z1", "Z2", "Z3")) + " |")

pairs = {}
lines += ["", "### Pairwise agreement", "", "| pair | factors | weights | variance explained | last ELBO |", "|---|---|---|---|---|"]
for a, b in itertools.combinations([n for n, r in runs.items() if r.get("trained")], 2):
    row = {block: agreement(values(runs[a], block), values(runs[b], block)) for block in ("factors", "weights", "r2", "elbo")}
    pairs[f"{a} vs {b}"] = row
    cell = lambda x: "exact" if x.get("exact") else (x.get("note") or f"{x['digits']} digits (max abs {x['max_abs_diff']:.3g})")
    lines.append(f"| {a} vs {b} | " + " | ".join(cell(row[k]) for k in ("factors", "weights", "r2", "elbo")) + " |")

(out / "compare.md").write_text("\n".join(lines) + "\n")
(out / "compare.json").write_text(json.dumps({"runs": sorted(runs), "pairs": pairs}, indent=2))
print("\n".join(lines))
