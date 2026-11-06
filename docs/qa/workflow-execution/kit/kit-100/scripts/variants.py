# variants.py <member dir>...: every recorded operation and variant from each member's analysis record (KIT-089's variants.txt form).
import json, sys, glob, os
def walk(x, out):
    if isinstance(x, dict):
        if "variant" in x and isinstance(x.get("variant"), str):
            out.append((x.get("state") or x.get("state_name") or x.get("name") or "?", x.get("operation") or x.get("declaration") or "?", x["variant"]))
        for v in x.values(): walk(v, out)
    elif isinstance(x, list):
        for v in x: walk(v, out)
for d in sys.argv[1:]:
    rec = glob.glob(os.path.join(d, "*/project/*/scripts/*/analysis_record.json"))
    out = []
    for r in rec: walk(json.load(open(r)), out)
    seen = set()
    for s, o, v in out:
        if (s, o, v) in seen: continue
        seen.add((s, o, v)); print(os.path.basename(d.rstrip("/")), s, "operation", o, "variant", v)
