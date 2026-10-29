#!/usr/bin/env python3
"""compare-trees.py <before_root> <after_root>: normalised relative path -> content hash; text files hashed with timestamps and temp paths masked."""
import sys, os, re, hashlib, collections
b, a = sys.argv[1], sys.argv[2]
TEXT = {".txt",".tab",".tsv",".csv",".json",".log",".ini",".md",".rmd",".html",".R",".yml",".yaml",".tex"}
norm_path = lambda p: re.sub(r"\d{8}[_-]?\d{6}", "TS", p)
MASK = [re.compile(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}(\.\d+)?( ?[A-Z]{2,4}| ?[+-]\d{4})?"), re.compile(r"\d{4}-\d{2}-\d{2}"), re.compile(r"/home/\S+|/tmp/\S+"), re.compile(r"\d{8}_\d{6}"), re.compile(r"\d{2}:\d{2}:\d{2}")]
def h(path):
    data = open(path, "rb").read()
    if os.path.splitext(path)[1].lower() in TEXT:
        try:
            t = data.decode("utf-8")
            for m in MASK: t = m.sub("X", t)
            data = t.encode()
        except UnicodeDecodeError: pass
    return hashlib.sha256(data).hexdigest()
def walk(root):
    out = collections.defaultdict(list)
    for d, _, fs in os.walk(root):
        for f in fs:
            p = os.path.join(d, f); out[norm_path(os.path.relpath(p, root))].append(h(p))
    return out
B, A = walk(b), walk(a)
keys = set(B) | set(A)
diff = sorted(k for k in keys if k in B and k in A and sorted(B[k]) != sorted(A[k]))
onlyb = sorted(k for k in keys if k not in A); onlya = sorted(k for k in keys if k not in B)
print(f"paths before {len(B)} after {len(A)} same-content {len(keys)-len(diff)-len(onlyb)-len(onlya)} differing {len(diff)} only-before {len(onlyb)} only-after {len(onlya)}")
print("differing by extension:", dict(collections.Counter(os.path.splitext(k)[1] for k in diff)))
for k in diff: print("  DIFF", k)
for k in onlyb: print("  ONLY-BEFORE", k)
for k in onlya: print("  ONLY-AFTER", k)
