# KIT-275 E4 (the scorer's self-test): constructed row files and the predictions of the scorer, derived here in Python (standard library
# only; run as `python3 -I make-selftest.py <out dir>`). No simulation, no emulation of R's random-number generator: every constructed
# data set of an arm carries the same value, so the percentile bootstrap of any statistic collapses onto its estimate and a prediction
# follows by arithmetic; the exact Clopper-Pearson limits are solved here by bisection on the binomial distribution. The rules the
# predictions apply are the ticket's (sections 7.1 to 7.3), written out again below, not read from the scorer.
# Writes <out>/cases/<case>/ (rows, g1.tsv), <out>/index.tsv (case, mode, environment), <out>/expected.tsv (case, key, kind, expected, tolerance).
import csv, gzip, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(os.path.dirname(HERE), "sim-k275.R")
OUT = sys.argv[1]

# ---- the row columns, read from sim-k275.R (a column list, not a prediction) -----------------------------------------------------
TAGS = ["05", "07", "09", "095"]
def row_columns():
    text = open(SIM).read()
    body = text[text.index("ROW_COLUMNS <- c(") + len("ROW_COLUMNS <- c("):text.index("newRow <- function")]
    cols = []
    for m in re.finditer(r'paste0\("([a-z_]+)", TAGS\)|"([A-Za-z0-9_]+)"', body):
        if m.group(1):
            cols += [m.group(1) + t for t in TAGS]
        else:
            cols.append(m.group(2))
    return cols
COLS = row_columns()
assert len(COLS) == len(set(COLS)) and "fa_post_05" in COLS and "smoke" in COLS, COLS

# ---- exact Clopper-Pearson limits (bisection) --------------------------------------------------------------------------------------
def cdf(x, n, p):
    return sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(x + 1))
def cp_upper(x, n, alpha):
    if n <= 0: return None
    if x >= n: return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(x, n, mid) > alpha: lo = mid
        else: hi = mid
    return (lo + hi) / 2
def cp_lower(x, n, alpha):
    if n <= 0: return None
    if x <= 0: return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 1 - cdf(x - 1, n, mid) < alpha: lo = mid
        else: hi = mid
    return (lo + hi) / 2

# ---- the rules of section 7.1 (margins and levels) ------------------------------------------------------------------------------
M_ABS, M_ACT, M_RHO, M_MATCH, M_ANY, M_ALPHA = 0.05, 0.50, 0.95, 0.03, 0.10, 0.0125
TOL = 1e-12
def label_A(lw, up, n, R):
    if n < 0.9 * R: return "not determined (coverage)"
    if up <= M_ABS + TOL: return "absorbed"
    if lw > M_ACT + TOL: return "acts"
    return "mixed"
DESIGNS = [("D3", 3), ("D4", 4), ("D6", 6)]
SIG = "250_31375_5239625"
ARMS = ["LM2", "GK2", "LA", "GA"]
V1G = {"LM2": dict(any=1, p05=0.30, se2=0.81, dot=2.0), "GK2": dict(any=0, p05=0.05, se2=1.0, dot=0.0),
       "LA": dict(any=1, p05=0.30, se2=0.81, dot=2.0), "GA": dict(any=0, p05=0.05, se2=1.0, dot=0.0)}
VD_ARM = {"LM2": "LM1", "GK2": "GK1", "LA": "LA", "GA": "GA"}
R_D = 200
GENUINE_BELOW = {"D3": 150, "D4": 14, "D6": 0}   # K275-r2: R2.4 (the fixed genuine sample below 0.5 in this many of the 200 V1a NO rows of a design)
PC1_FLOOR = 0.80   # K275-r2: 7.1

def blank(**kw):
    row = {c: "NA" for c in COLS}
    row.update({k: str(v) for k, v in kw.items()})
    return row

def d_rows(spec):
    rows = []
    for design, n in DESIGNS:
        for r in range(1, R_D + 1):
            common = dict(tier="D", design=design, cell=design, omic="lipid", r=r, smoke=spec.get("smoke", "FALSE"), status="ok", tested_sig_y=SIG, den_g=2.0,
                          e5_y="TRUE", e5_x="TRUE", seconds=1)
            nonfail = spec.get("pc1_nonfail", {}).get(design, 0)
            gb = GENUINE_BELOW[design]
            rows.append(blank(version="V1a", arm="NO", fa_pre_05=0 if r <= nonfail else 1, fa_post_05=1, any_call_y=0, null_p05_y=0.05, se2_mod_y=1.0, dot_g=0.0,
                              cy_genuine=0.40 if r <= gb else 0.80, **common))
            rows.append(blank(version="V1w", arm="NO", fa_pre_05=1 if r <= 50 else 0, cy_pre_a=0.44, cy_genuine=0.80, any_call_y=0, null_p05_y=0.05, se2_mod_y=1.0, dot_g=0.0, **common))   # K275-r2: R2.4 (the 17th row)
            for arm in ARMS:
                v = V1G[arm]
                fails = spec.get("A_fail", {}).get((design, arm), 0)
                valid = spec.get("valid", {}).get((design, arm), R_D)
                status = "ok" if r <= valid else "crashed"
                c2 = dict(common, status=status)
                rows.append(blank(version="V1a", arm=arm, fa_post_05=1 if r <= fails else 0, any_call_y=v["any"], null_p05_y=v["p05"], se2_mod_y=v["se2"], dot_g=v["dot"], **c2))
                nc1 = spec.get("nc1_fail", {}).get(design, 0) if arm == "GK2" else 0
                rows.append(blank(version="V0", arm=arm, fa_post_05=1 if r <= nc1 else 0, any_call_y=0, null_p05_y=0.05, se2_mod_y=1.0, dot_g=0.0, **common))
                rows.append(blank(version="VD", arm=VD_ARM[arm], fa_post_05=0, any_call_y=0, null_p05_y=0.05, se2_mod_y=1.0, dot_g=0.0, **common))
            derived_sq = (2 * n - 3) / (2 * n - 1)
            se2_lx = spec.get("lx_se2", {}).get(design, derived_sq)
            rows.append(blank(version="V1a", arm="LX", any_call_y=0, null_p05_y=0.05, se2_mod_y=1.0, se2_unmod_y=se2_lx, dot_g=0.1, **common))
            rows.append(blank(version="VD", arm="LX", any_call_y=0, null_p05_y=0.05, se2_mod_y=1.0, se2_unmod_y=1.0, dot_g=0.0, **common))
    return rows

def g1_rows(spec):
    rows = []
    for design, n in DESIGNS:
        for kind, target in (("PC2", math.sqrt((2 * n - 3) / (2 * n - 1))), ("NC2", 1.0)):
            value = spec.get("g1", {}).get((design, kind), target)
            for r in range(1, 21):
                rows.append(dict(design=design, kind=kind, r=r, rho=repr(value), smoke="FALSE"))
    return rows

R_BASE = 2000   # K275-r2: 5.3 (the counted r2 data sets are r = 2001 to 2200; the scorer indexes them 1 to R)
EXPECT_ROWS = {"D": 17, "K": 17, "F": 5, "S": 4, "V": 12, "M": 7, "T": 6, "E": 6, "ID": 5}
def complete(rows):
    """Pad every (tier, design, r) group to the rows a complete data set has (filler rows of arms FILLk), as the scorer's continuation rule counts them."""
    groups = {}
    for r_ in rows:
        if r_["arm"] != "unit":
            groups.setdefault((r_["tier"], r_["design"], r_["r"]), []).append(r_)
    out = list(rows)
    for (t, d, r), g in groups.items():
        for k in range(len(g), EXPECT_ROWS[t]):
            out.append(blank(tier=t, design=d, cell=g[0]["cell"], omic=g[0]["omic"], r=r, version="V0", arm="FILL%d" % k, smoke=g[0]["smoke"], status="ok"))
    return out
def shifted(row):   # K275-r2: 5.3
    out = dict(row)
    if out["r"] not in ("NA", ""): out["r"] = str(int(float(out["r"])) + R_BASE)
    return out
def write_rows(path, rows, gz=True, pad=True):
    if pad: rows = complete(rows)
    rows = [shifted(r_) for r_ in rows]
    opener = gzip.open if gz else open
    with opener(path, "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        for row in rows:
            w.writerow([row[c] for c in COLS])
def write_g1(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["design", "kind", "r", "rho", "smoke"])
        for r in rows:
            w.writerow([r["design"], r["kind"], r["r"], r["rho"], r["smoke"]])

# ---- the predictions of a D-tier case ------------------------------------------------------------------------------------------------
def predict(case, spec, expected):
    E = lambda key, kind, value, tol="1e-8": expected.append((case, key, kind, str(value), tol))
    labels = {}
    pass_ctrl = {}
    for design, n in DESIGNS:
        nonfail = spec.get("pc1_nonfail", {}).get(design, 0)
        x = R_D - nonfail
        lo, hi = cp_lower(x, R_D, 0.01 / 3), cp_upper(x, R_D, 0.01 / 3)
        ok = hi >= PC1_FLOOR - TOL   # K275-r2: 7.1 (F1: the one-sided upper limit at 0.01/3 below 0.80, i.e. fewer than 144 of 200)
        pass_ctrl[("PC1", design)] = ok
        E(f"cell|PC1|D|{design}|V1a|NO|label", "str", "passes" if ok else "FAILS")
        E(f"cell|PC1|D|{design}|V1a|NO|lo", "num", lo)
        E(f"cell|PC1|D|{design}|V1a|NO|hi", "num", hi)
        E(f"cell|PC1|D|{design}|V1a|NO|note", "str", "genuine sample (cy_genuine) below 0.5 in %d of 200 data sets" % GENUINE_BELOW[design])   # K275-r2: R2.4
        E(f"cell|P_pre|D|{design}|V1w|NO|estimate", "num", 0.25)   # K275-r2: R2.4 (V1w is reported, no label)
        E(f"cell|P_pre|D|{design}|V1w|NO|n_valid", "num", 200)
        E(f"cell|q_cy_pre_a|D|{design}|V1w|NO|estimate", "num", 0.44)
        nc1 = spec.get("nc1_fail", {}).get(design, 0)
        up = cp_upper(nc1, R_D, 0.01 / 3)
        ok = up <= M_ABS + TOL
        pass_ctrl[("NC1", design)] = ok
        E(f"cell|NC1|D|{design}|V0|GK2|label", "str", "passes" if ok else "FAILS")
        E(f"cell|NC1|D|{design}|V0|GK2|hi", "num", up)
        for arm in ARMS:
            fails = spec.get("A_fail", {}).get((design, arm), 0)
            valid = spec.get("valid", {}).get((design, arm), R_D)
            # failures are the first `fails` data sets; a data set beyond `valid` is a crashed journey (never a denominator)
            f_eff = min(fails, valid)
            lw, up = cp_lower(f_eff, valid, 0.025 / 12), cp_upper(f_eff, valid, 0.025 / 12)
            lab = label_A(lw, up, valid, R_D)
            labels[(design, arm)] = lab
            E(f"cell|A|D|{design}|V1a|{arm}|label", "str", lab)
            E(f"cell|A|D|{design}|V1a|{arm}|lo", "num", lw)
            E(f"cell|A|D|{design}|V1a|{arm}|hi", "num", up)
            E(f"cell|A|D|{design}|V1a|{arm}|n_valid", "num", valid)
    # the pairs and ratios (every data set of an arm carries the same values, so every interval collapses onto its estimate)
    for design, n in DESIGNS:
        for arm in ("LM2", "GK2"):
            v = V1G[arm]
            rho = math.sqrt(v["se2"] / 1.0)
            if spec.get("valid", {}).get((design, arm), R_D) >= 0.9 * R_D:
                E(f"cell|B|D|{design}|V1a|{arm}|estimate", "num", rho)
                E(f"cell|B|D|{design}|V1a|{arm}|label", "str", "understated" if rho < M_RHO - TOL else "not understated")
                E(f"cell|C|D|{design}|V1a|{arm}|estimate", "num", v["any"] - 0)
                E(f"cell|C|D|{design}|V1a|{arm}|label", "str", "raised (material)" if v["any"] > M_ANY + TOL else "not raised")
                dalpha = v["p05"] - 0.05
                E(f"cell|C'|D|{design}|V1a|{arm}|estimate", "num", dalpha)
                E(f"cell|C'|D|{design}|V1a|{arm}|label", "str", "raised (material)" if dalpha > M_ALPHA + TOL else "not raised")
            pi = (v["dot"] - 0.0) / 2.0
            E(f"cell|pi|D|{design}|V1a|{arm}|estimate", "num", pi)
            vv = spec.get("valid", {}).get((design, arm), R_D)
            E(f"cell|pi|D|{design}|V1a|{arm}|label", "str", "not determined (coverage)" if vv < 0.9 * R_D else "imprint retained" if pi > 0.5 + TOL else "imprint removed" if pi < 0.5 - TOL else "not determined")
        derived = math.sqrt((2 * n - 3) / (2 * n - 1))
        rho_lx = math.sqrt(spec.get("lx_se2", {}).get(design, derived ** 2) / 1.0)
        E(f"cell|B'|D|{design}|V1a|LX|estimate", "num", rho_lx)
        E(f"cell|B'|D|{design}|V1a|LX|label", "str", "as derived" if abs(rho_lx - derived) <= M_MATCH + TOL else f"differs from derived ({round(rho_lx, 4):g})")
        E(f"cell|pi|D|{design}|V1a|LX|label", "str", "imprint removed")
        E(f"cell|B|D|{design}|V1a|LM2|note", "str", "%d pairs dropped for a differing tested set" % spec.get("sig_dropped", {}).get(design, 0)) if spec.get("sig_dropped") else None
    # the outcome row (7.3): the first row that holds
    f1 = not all(pass_ctrl[("PC1", d)] for d, _ in DESIGNS)
    f2 = not all(pass_ctrl[("NC1", d)] for d, _ in DESIGNS)
    g1 = spec.get("g1", {})
    f3 = any(abs(g1.get((d, k), t) - t) > M_MATCH + TOL for d, n in DESIGNS for k, t in (("PC2", math.sqrt((2 * n - 3) / (2 * n - 1))), ("NC2", 1.0)))
    f4 = bool(spec.get("e5_false"))
    lm2 = [labels[(d, "LM2")] for d, _ in DESIGNS]
    gk2 = [labels[(d, "GK2")] for d, _ in DESIGNS]
    if f1: outcome = "Stop (F1)"
    elif f2: outcome = "Stop (F2)"
    elif f3: outcome = "Stop (F3)"
    elif f4: outcome = "Stop (F4)"
    elif all(x == "absorbed" for x in lm2 + gk2): outcome = "R1"
    elif all(x == "absorbed" for x in gk2): outcome = "R2"
    elif all(x == "absorbed" for x in lm2): outcome = "R3"
    elif all(x == "acts" for x in lm2 + gk2): outcome = "R4"
    else: outcome = "R5"
    E("decision|outcome", "str", outcome)
    for d, _ in DESIGNS:
        E(f"cell|PC2|G1|{d}|NA|NA|label", "str", "passes" if abs(g1.get((d, "PC2"), math.sqrt((2 * _ - 3) / (2 * _ - 1))) - math.sqrt((2 * _ - 3) / (2 * _ - 1))) <= M_MATCH + TOL else "FAILS")
    return outcome

CASES = {}
def case(name, **spec): CASES[name] = spec
fail0 = {}
case("r1")
case("r2", A_fail={("D3", "LM2"): 130, ("D4", "LM2"): 100})
case("r3", A_fail={("D3", "GK2"): 2})
case("r4", A_fail={(d, a): 200 for d, _ in DESIGNS for a in ("LM2", "GK2")})
case("r5", A_fail={("D3", "LM2"): 0, ("D4", "LM2"): 100, ("D6", "LM2"): 200, ("D3", "GK2"): 200, ("D6", "GK2"): 100})
case("a_bounds", A_fail={("D3", "LM2"): 1, ("D4", "LM2"): 2, ("D6", "LM2"): 130, ("D3", "GK2"): 120, ("D4", "GK2"): 121})
case("nc1_2", nc1_fail={"D3": 2})
case("nc1_3", nc1_fail={"D3": 3})
case("pc1_144", pc1_nonfail={"D3": 56})   # K275-r2: 7.1 (144 of 200 failed: upper limit 0.8019, passes)
case("pc1_143", pc1_nonfail={"D3": 57})   # K275-r2: 7.1 (143 of 200: upper limit 0.7975, F1)
case("g1_bad", g1={("D3", "PC2"): 1.0})
case("e5_false", e5_false=True)
case("coverage", valid={("D3", "GK2"): 170})
case("bp_off", lx_se2={"D6": 0.64})

expected = []
index = []
os.makedirs(OUT, exist_ok=True)
for name, spec in CASES.items():
    folder = os.path.join(OUT, "cases", name)
    os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
    rows = d_rows(spec)
    if spec.get("e5_false"):
        rows[5]["e5_y"] = "FALSE"
    write_rows(os.path.join(folder, "rows", f"part-{name}.tsv.gz".replace("part-", "").replace("tsv.gz", "tsv.gz")), rows)
    write_g1(os.path.join(folder, "g1.tsv"), g1_rows(spec))
    predict(name, spec, expected)
    index.append((name, "score", "K275_SCORE_R=200,100,20;K275_BOOT_B=200"))

# ---- the tiers named T and F, a plain file beside the gzipped one, and a missing e5 (NA is not FALSE) --------------------------------------------------
name = "tf"
folder = os.path.join(OUT, "cases", name)
os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
write_rows(os.path.join(folder, "rows", "d.tsv.gz"), d_rows({}))
extra = []
for tier, version in (("T", "V1t"), ("F", "V1p")):
    for r in range(1, 11):
        extra.append(blank(tier=tier, design="D3", cell="D3", omic="lipid", r=r, version=version, arm="NO", smoke="FALSE", status="ok", fa_pre_05=0, fa_post_05=1 if r <= 2 else 0,
                           cy_post_a=0.5, e5_y="NA", e5_x="NA", seconds=1))
write_rows(os.path.join(folder, "rows", "part-tf.tsv"), extra, gz=False)
write_g1(os.path.join(folder, "g1.tsv"), g1_rows({}))
predict("tf", {}, expected)
for tier in ("T", "F"):
    E = lambda key, kind, value, tol="1e-8": expected.append(("tf", key, kind, str(value), tol))
    E(f"cell|P_post|{tier}|D3|{'V1t' if tier == 'T' else 'V1p'}|NO|n_valid", "num", 10)
    E(f"cell|P_post|{tier}|D3|{'V1t' if tier == 'T' else 'V1p'}|NO|estimate", "num", 0.2)
    E(f"cell|P_post|{tier}|D3|{'V1t' if tier == 'T' else 'V1p'}|NO|lo", "num", cp_lower(2, 10, 0.025))
    E(f"cell|P_post|{tier}|D3|{'V1t' if tier == 'T' else 'V1p'}|NO|hi", "num", cp_upper(2, 10, 0.025))
    E(f"cell|q_cy_post_a|{tier}|D3|{'V1t' if tier == 'T' else 'V1p'}|NO|estimate", "num", 0.5)
index.append(("tf", "score", "K275_SCORE_R=200,100,20;K275_BOOT_B=200"))

# ---- refusals: the smoke flag, a missing G1 table ------------------------------------------------------------------------------------------------------
name = "smoke"
folder = os.path.join(OUT, "cases", name)
os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
write_rows(os.path.join(folder, "rows", "d.tsv.gz"), d_rows({"smoke": "TRUE"}))
write_g1(os.path.join(folder, "g1.tsv"), g1_rows({}))
expected.append(("smoke", "exit|smoke", "exit", "nonzero", ""))
expected.append(("smoke", "stderr|smoke", "substr", "smoke flag", ""))
index.append(("smoke", "score", "K275_SCORE_R=200,100,20;K275_BOOT_B=200"))
name = "g1_missing"
folder = os.path.join(OUT, "cases", name)
os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
write_rows(os.path.join(folder, "rows", "d.tsv.gz"), d_rows({}))
expected.append((name, "exit|g1_missing", "exit", "nonzero", ""))
expected.append((name, "stderr|g1_missing", "substr", "K275_G1", ""))
index.append((name, "score-nog1", "K275_SCORE_R=200,100,20;K275_BOOT_B=200"))

# ---- blind mode: parse, coverage and E8 on smoke-flagged rows (accepted), D4 against ID: 5 rows compared, 4 equal ------------------------------------------------
name = "blind"
folder = os.path.join(OUT, "cases", name)
os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
brows = []
for arm in ("NO", "LM2", "GK2", "LA", "GA"):
    for tier, cell, omic in (("D", "D4", "lipid"), ("ID", "D4M", "metab")):
        dev = 1e-8 if (tier == "ID" and arm == "GA") else 0
        brows.append(blank(tier=tier, design="D4", cell=cell, omic=omic, r=1, version="V1a", arm=arm, smoke="TRUE", status="ok", seconds=1.5 + (0 if tier == "D" else 0.3),
                           cy_post_a=0.75 + dev, any_call_y=1, fa_post_05=0))
write_rows(os.path.join(folder, "rows", "d.tsv.gz"), brows)
expected.append((name, "scorer|parse|pass", "str", "TRUE", ""))
expected.append((name, "scorer|coverage NO D4|value", "str", "2 journeys, 2 completed", ""))
expected.append((name, "scorer|coverage GA D4|pass", "str", "TRUE", ""))
expected.append((name, "scorer|E8|value", "str", "5 compared, 4 equal, max deviation 1e-08", ""))
expected.append((name, "scorer|E8|pass", "str", "FALSE", ""))
index.append((name, "blind", "K275_BOOT_B=200"))

# ---- the label rules alone, at their boundaries (a limit "equal to" a margin) ------------------------------------------------------------------------------------
name = "labels"
folder = os.path.join(OUT, "cases", name)
os.makedirs(folder, exist_ok=True)
L = []   # id, kind, lo, hi, n, R, derived, est, margin, expected label, shown
def lab_case(i, kind, lo, hi, n, R, label, derived="NA", est="NA", margin="NA", shown=None):
    L.append((i, kind, lo, hi, n, R, derived, est, margin, label, shown))
# family A: upper limit equal to 0.05 is absorbed; lower limit equal to 0.5 is not acts; coverage at exactly 0.90 R is enough
lab_case("A1", "A", 0.0, 0.05, 200, 200, "absorbed")
lab_case("A2", "A", 0.0, 0.0507, 200, 200, "mixed")
lab_case("A3", "A", 0.5, 0.9, 200, 200, "mixed")
lab_case("A4", "A", 0.5000001, 0.9, 200, 200, "acts")
lab_case("A5", "A", 0.0, 0.0, 179, 200, "not determined (coverage)")
lab_case("A6", "A", 0.0, 0.0, 180, 200, "absorbed")
# family B: upper limit equal to 0.95 is not determined; below it understated; lower limit equal to 0.95 not understated
lab_case("B1", "B", 0.9, 0.95, 200, 200, "not determined")
lab_case("B2", "B", 0.9, 0.9499, 200, 200, "understated")
lab_case("B3", "B", 0.95, 1.0, 200, 200, "not understated")
lab_case("B4", "B", 0.9, 1.0, 200, 200, "not determined")
lab_case("B5", "B", 0.9, 0.94, 150, 200, "not determined (coverage)")
# B': the interval within derived +/- 0.03 (both ends equal to the bounds still matches)
d6 = math.sqrt(9 / 11)
lab_case("P1", "Bp", d6 - 0.03, d6 + 0.03, 100, 100, "as derived", derived=d6, est=d6)
lab_case("P2", "Bp", d6 - 0.03, d6 + 0.0301, 100, 100, f"differs from derived ({round(d6 + 0.01, 4):g})", derived=d6, est=d6 + 0.01)
# C (margin 0.10): lower limit equal to the margin is not raised(material); upper limit equal to the margin is not raised; any rise shown needs lower > 0
lab_case("C1", "C", 0.10, 0.2, 200, 200, "not determined", margin=0.10, shown="TRUE")
lab_case("C2", "C", 0.1000001, 0.2, 200, 200, "raised (material)", margin=0.10, shown="TRUE")
lab_case("C3", "C", 0.0, 0.10, 200, 200, "not raised", margin=0.10, shown="FALSE")
lab_case("C4", "C", 0.000001, 0.10, 200, 200, "not raised", margin=0.10, shown="TRUE")
lab_case("C5", "C", 0.0, 0.2, 200, 200, "not determined", margin=0.10, shown="FALSE")
lab_case("C6", "C", 0.0125, 0.02, 200, 200, "not determined", margin=0.0125, shown="TRUE")
lab_case("C7", "C", -0.01, 0.0125, 200, 200, "not raised", margin=0.0125, shown="FALSE")
# pi at 0.5
lab_case("I1", "pi", 0.5, 0.9, 200, 200, "not determined")
lab_case("I2", "pi", 0.5001, 0.9, 200, 200, "imprint retained")
lab_case("I3", "pi", 0.1, 0.5, 200, 200, "not determined")
lab_case("I4", "pi", 0.1, 0.4999, 200, 200, "imprint removed")
with open(os.path.join(folder, "table.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["id", "kind", "lo", "hi", "n", "R", "derived", "est", "margin"])
    for row in L:
        w.writerow([row[0], row[1], repr(row[2]), repr(row[3]), row[4], row[5], row[6] if isinstance(row[6], str) else repr(row[6]),
                    row[7] if isinstance(row[7], str) else repr(row[7]), row[8] if isinstance(row[8], str) else repr(row[8])])
for row in L:
    expected.append((name, f"labels|{row[0]}|label", "str", row[9], ""))
    if row[10] is not None:
        expected.append((name, f"labels|{row[0]}|shown", "str", row[10], ""))
index.append((name, "labels", ""))


# ---- S16: the cases the audit asked for (a) to (l), (g) non-degenerate ------------------------------------------------------------------------------------
def new_case(name, rows, g1=None, files=None, mode="score"):
    folder = os.path.join(OUT, "cases", name)
    os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
    write_rows(os.path.join(folder, "rows", "d.tsv.gz"), rows)
    write_g1(os.path.join(folder, "g1.tsv"), g1_rows({}))
    for fname, text in (files or {}).items():
        open(os.path.join(folder, "rows", fname), "w").write(text)
    index.append((name, mode, "K275_SCORE_R=200,100,20;K275_BOOT_B=200"))
    return lambda key, kind, value, tol="1e-8": expected.append((name, key, kind, str(value), tol))
def sel(rows, **kw):
    return [r for r in rows if all(r[k] == str(v) for k, v in kw.items())]

# (a) unlabelled cells; R-2: the emulation against VD at the SAME k (tier K's VD LM2, not tier D's VD LM1)
rows = d_rows({})
for r in range(1, 11):
    rows.append(blank(tier="E", design="D3", cell="D3", omic="lipid", r=r, version="V1a", arm="PEL2", smoke="FALSE", status=["uncorrected", "dropped", "refused_after_exclusion"][r % 3] if r > 6 else "uncorrected",
                      pe_a_excluded="TRUE" if r <= 5 else "FALSE", pe_n_excluded=1, pe_excluded_x=0, tested_sig_y=SIG, se2_mod_y=0.64, fa_post_05=0, any_call_y=0, null_p05_y=0.05, dot_g=0.0, den_g=2.0))
    rows.append(blank(tier="K", design="D3", cell="D3", omic="lipid", r=r, version="VD", arm="LM2", smoke="FALSE", status="ok", tested_sig_y=SIG, se2_mod_y=0.5, any_call_y=0, null_p05_y=0.05, dot_g=0.0, den_g=2.0))
E = new_case("unlabelled", rows)
predict("unlabelled", {}, expected)
E("cell|rho_mod|D|D3|V1a|LM2|estimate", "num", 0.9)
E("cell|pi|D|D3|V1a|LM2|estimate", "num", 1.0)
E("cell|pi|D|D3|V1a|LM2|label", "str", "imprint retained")
E("cell|delta_any_call_y_V0|D|D3|V1a|LM2|estimate", "num", 1.0)
E("cell|delta_any_call_y_VD|D|D3|V1a|LM2|estimate", "num", 1.0)
E("cell|rho_mod|E|D3|V1a|PEL2|estimate", "num", math.sqrt(0.64 / 0.5))   # VD LM2 (same k), not VD LM1 (which would give 0.8)
E("cell|rho_mod|E|D3|V1a|PEL2|estimand", "str", "rho_mod against VD LM2 (K)")
# (e) H2 and H3: statuses r=7..10 -> r%3: 7:dropped? computed below by the same rule
st = [(["uncorrected", "dropped", "refused_after_exclusion"][r % 3] if r > 6 else "uncorrected") for r in range(1, 11)]
E("cell|H3_uncorrected|E|D3|V1a|PEL2|estimate", "num", st.count("uncorrected") / 10)
E("cell|H3_dropped|E|D3|V1a|PEL2|estimate", "num", st.count("dropped") / 10)
E("cell|H3_refused_after_exclusion|E|D3|V1a|PEL2|estimate", "num", st.count("refused_after_exclusion") / 10)
E("cell|H3_refused_after_exclusion|E|D3|V1a|PEL2|hi", "num", cp_upper(st.count("refused_after_exclusion"), 10, 0.025))
E("cell|H2|E|D3|V1a|PEL2|estimate", "num", 0.5)

# (b) pairs with differing tested sets: 20 data sets of V1a LM2 at D3 carry another signature (dropped from the pair, counted)
rows = d_rows({})
n_drop = 0
for r_ in rows:
    if r_["tier"] == "D" and r_["design"] == "D3" and r_["version"] == "V1a" and r_["arm"] == "LM2" and int(r_["r"]) <= 20:
        r_["tested_sig_y"] = "249_1_1"; n_drop += 1
E = new_case("sigdiff", rows)
predict("sigdiff", {"sig_dropped": {"D3": 20}}, expected)

# (c) tier K: LG2's pi and its reference VD LG1 across tiers
rows = d_rows({})
for r in range(1, 101):
    rows.append(blank(tier="K", design="D3", cell="D3", omic="lipid", r=r, version="V1a", arm="LG2", smoke="FALSE", status="ok", tested_sig_y=SIG, dot_g=2.0, den_g=2.0, se2_mod_y=1.0))
    rows.append(blank(tier="K", design="D3", cell="D3", omic="lipid", r=r, version="VD", arm="LG1", smoke="FALSE", status="ok", tested_sig_y=SIG, dot_g=0.0, den_g=2.0, se2_mod_y=1.0))
E = new_case("tierk", rows)
predict("tierk", {}, expected)
E("cell|pi|K|D3|V1a|LG2|estimate", "num", 1.0)
E("cell|pi|K|D3|V1a|LG2|label", "str", "imprint retained")

# (d) E8 in scoring mode: 3 data sets of D4 against the ID rows, one value differing by 0.5 (a finding, not a stop)
rows = d_rows({})
ids = []
for r_ in rows:
    if r_["tier"] == "D" and r_["design"] == "D4" and r_["version"] == "V1a" and r_["arm"] in ("NO", "LM2", "GK2", "LA", "GA") and int(r_["r"]) <= 3:
        c = dict(r_); c["tier"] = "ID"; c["cell"] = "D4M"; c["omic"] = "metab"
        ids.append(c)
ids[0]["se2_mod_y"] = "0.5" if ids[0]["se2_mod_y"] != "0.5" else "0.25"
dev = abs(float(ids[0]["se2_mod_y"]) - float(sel(rows, tier="D", design="D4", version="V1a", arm=ids[0]["arm"], r=ids[0]["r"])[0]["se2_mod_y"]))
def write_sub(name, parts, g1=None):
    folder = os.path.join(OUT, "cases", name)
    for sub_, rows_ in parts.items():
        os.makedirs(os.path.join(folder, "rows", sub_), exist_ok=True)
        write_rows(os.path.join(folder, "rows", sub_, "d.tsv.gz"), rows_, pad=False)
    write_g1(os.path.join(folder, "g1.tsv"), g1_rows({}))
    index.append((name, "score", "K275_SCORE_R=200,100,20;K275_BOOT_B=200"))
    return lambda key, kind, value, tol="1e-8": expected.append((name, key, kind, str(value), tol))
E = write_sub("e8score", {"j1": complete(rows), "j2": complete(ids)})   # H01: the two sides come from different sources
predict("e8score", {}, expected)
E("decision|E8", "str", "15 compared, 14 equal, max deviation %s" % ("%g" % dev))

# (R-3) tier E's pass-1 NO against tier D's NO: identical -> reported; one differing -> Stop (F4)
for name, change in (("crossno_ok", False), ("crossno_bad", True)):
    rows = d_rows({})
    es = []
    for r_ in rows:
        if r_["tier"] == "D" and r_["design"] == "D3" and r_["version"] == "V1a" and r_["arm"] == "NO" and int(r_["r"]) <= 5:
            c = dict(r_); c["tier"] = "E"; es.append(c)
    if change: es[0]["cy_post_a"] = "0.123"
    E = write_sub(name, {"j1": complete(rows), "j2": complete(es)})   # H01: tier D and tier E from different sources
    predict(name, {}, expected)
    if change:
        expected[:] = [e for e in expected if not (e[0] == name and e[1] == "decision|outcome")]
        E("decision|outcome", "str", "Stop (F4)")
    E("decision|crossNO", "str", "5 compared, %d differ" % (1 if change else 0))

# (f) no_da, dropped, uncorrected stay in the denominators of family A
rows = d_rows({})
k = 0
for r_ in rows:
    if r_["tier"] == "D" and r_["design"] == "D3" and r_["version"] == "V1a" and r_["arm"] == "LM2" and int(r_["r"]) <= 30:
        r_["status"] = ["no_da", "dropped", "uncorrected"][(int(r_["r"]) - 1) // 10]
E = new_case("statuses", rows)
predict("statuses", {}, expected)

# (h) rows outside the declared range; a duplicate journey row
rows = d_rows({})
bad = dict(rows[0]); bad["r"] = "201"
E = new_case("range_stop", rows + [bad]); expected += [("range_stop", "exit|range_stop", "exit", "nonzero", ""), ("range_stop", "stderr|range_stop", "substr", "RANGE-1", "")]
E = new_case("dup_stop", rows + [dict(rows[0])]); expected += [("dup_stop", "exit|dup_stop", "exit", "nonzero", ""), ("dup_stop", "stderr|dup_stop", "substr", "appears twice", "")]

# (i) an instrument status, and a partial-match guard class
rows = d_rows({}); rows[3]["status"] = "instrument"
E = new_case("instrument", rows); expected += [("instrument", "exit|instrument", "exit", "nonzero", ""), ("instrument", "stderr|instrument", "substr", "instrument status", "")]
rows = d_rows({}); rows[3]["class"] = "partial_match_guard: Partial match of 'abc' to 'abcdef' in data frame"
E = new_case("guard_class", rows); expected += [("guard_class", "exit|guard_class", "exit", "nonzero", ""), ("guard_class", "stderr|guard_class", "substr", "partial-match guard", "")]

# (j) pi without coverage, and pi with a missing reference
rows = [r_ for r_ in d_rows({}) if not (r_["tier"] == "D" and r_["design"] == "D4" and r_["version"] == "VD" and r_["arm"] == "LM1")]
rows = [r_ for r_ in rows if not (r_["tier"] == "D" and r_["design"] == "D3" and r_["version"] == "V1a" and r_["arm"] == "LM2" and int(r_["r"]) > 100)]
E = new_case("picov", rows)
E("cell|pi|D|D3|V1a|LM2|label", "str", "not determined (coverage)")
E("cell|pi|D|D4|V1a|LM2|label", "str", "not determined (coverage)")
E("cell|pi|D|D4|V1a|LM2|n_valid", "num", 0)

# (k) F4 through a failing E1 row in a checks file; (l) an E12 failure is reported, not a stop
hdr = "check\tcell\trep\tcondition\tvalue\tthreshold\tpass\tnote\n"
E = new_case("e1_fail", d_rows({}), files={"checks-x.tsv": hdr + "E1\tD3\t1\tdeterminism\tNA\tbyte-identical\tFALSE\t\n"})
E("decision|outcome", "str", "Stop (F4)")
E = new_case("e12_fail", d_rows({}), files={"checks-x.tsv": hdr + "E12\tD3\t1\tV1a LG2\tNA\tpositional\tFALSE\t\n"})
E("decision|outcome", "str", "R1")
E("decision|E12 failing checks", "str", "1")

# (g) non-degenerate data: se2_mod of V1a LM2 at D3 alternates 0.5 and 1.5 (VD 1.0): the estimate and the jackknife limits are predicted exactly,
# the bootstrap limits are bounded by the data range
rows = d_rows({})
for r_ in rows:
    if r_["tier"] == "D" and r_["design"] == "D3" and r_["version"] == "V1a" and r_["arm"] == "LM2":
        r_["se2_mod_y"] = "0.5" if int(r_["r"]) % 2 else "1.5"
E = new_case("nondegen", rows)
vals = [0.5 if r % 2 else 1.5 for r in range(1, 201)]
tot = sum(vals); n = 200
est = tot / n
jack = [(tot - v) / (n - 1) for v in vals]
mj = sum(jack) / n
se = math.sqrt((n - 1) / n * sum((j - mj) ** 2 for j in jack))
from statistics import NormalDist
z = NormalDist().inv_cdf(1 - 0.025 / 6)
E("cell|B|D|D3|V1a|LM2|estimate", "num", math.sqrt(est))
E("cell|B|D|D3|V1a|LM2|jk_lo", "num", math.sqrt(est - z * se))
E("cell|B|D|D3|V1a|LM2|jk_hi", "num", math.sqrt(est + z * se))
E("cell|B|D|D3|V1a|LM2|lo", "between", "%r:%r" % (math.sqrt(0.5), math.sqrt(est)))
E("cell|B|D|D3|V1a|LM2|hi", "between", "%r:%r" % (math.sqrt(est), math.sqrt(1.5)))


# ---- the unit paths of the guard (coordinator's ruling after smoke 2): an instrument unit row, a crashed unit row, a unit that did not complete ----------------------
rows = d_rows({}); u = blank(tier="D", design="D3", r=9101, arm="unit", smoke="FALSE", status="instrument", cls="")
u["class"] = "partial_match_guard: partial-match guard: Partial match of 'abc' to 'abcdef' in data frame"
E = new_case("unit_instrument", rows + [u]); expected += [("unit_instrument", "exit|unit_instrument", "exit", "nonzero", ""), ("unit_instrument", "stderr|unit_instrument", "substr", "instrument status", "")]
u2 = blank(tier="D", design="D3", r=1, arm="unit", smoke="FALSE", status="crashed"); u2["class"] = "child process exit 1, no result"
E = new_case("unit_crashed", rows + [u2]); expected += [("unit_crashed", "exit|unit_crashed", "exit", "nonzero", ""), ("unit_crashed", "stderr|unit_crashed", "substr", "left no journey rows", "")]
E = new_case("unit_completed_false", rows, files={"checks-x.tsv": hdr + "unit_completed\tD\t1\tD3\tNA\tran to the end\tFALSE\tchild process exit 1, no result\n"})
expected += [("unit_completed_false", "exit|unit_completed_false", "exit", "nonzero", ""), ("unit_completed_false", "stderr|unit_completed_false", "substr", "did not complete", "")]
# ---- the continuation rule (RG4) and tier V's share (RG6) ---------------------------------------------------------------------------------------------------
# (1) an interrupted data set: D3 r = 7 lost three rows (V1a LM2, V0 LM2, VD LM1); it is dropped whole and listed; coverage counts it missing
rows = complete(d_rows({}))
drop = [r_ for r_ in rows if r_["tier"] == "D" and r_["design"] == "D3" and r_["r"] == "7" and ((r_["version"], r_["arm"]) in (("V1a", "LM2"), ("V0", "LM2"), ("VD", "LM1")))]
rows = [r_ for r_ in rows if not any(r_ is x for x in drop)]
E = write_sub("interrupted", {"": rows})
predict("interrupted", {"valid": {}}, expected)
expected[:] = [e for e in expected if not (e[0] == "interrupted" and e[1].startswith("cell|") and e[1].split("|")[3] == "D3")]
E("decision|interrupted data sets", "str", "D|D3|7|")
E("cell|A|D|D3|V1a|LM2|n_valid", "num", 199)
E("cell|A|D|D3|V1a|LM2|label", "str", "absorbed")
# (2) a data set complete in two sources stops scoring
rows = complete(d_rows({}))
write_sub("double_source", {"j1": rows, "j2": [r_ for r_ in rows if r_["design"] == "D3" and r_["r"] == "5"]})
expected += [("double_source", "exit|double_source", "exit", "nonzero", ""), ("double_source", "stderr|double_source", "substr", "two sources", "")]
# (3) a continuation that completes a job: j1 holds D3 r = 1..151 (r = 151 interrupted), j2 (the continuation) holds D3 r = 151..200
rows = complete(d_rows({}))
j1 = [r_ for r_ in rows if not (r_["design"] == "D3" and int(r_["r"]) >= 152)]
cut = [r_ for r_ in j1 if r_["design"] == "D3" and r_["r"] == "151" and r_["tier"] == "D"][:3]
j1 = [r_ for r_ in j1 if not any(r_ is x for x in cut)]
j2 = [r_ for r_ in rows if r_["design"] == "D3" and int(r_["r"]) >= 151]
E = write_sub("continuation", {"j1": j1, "j2": j2})
predict("continuation", {}, expected)
E("decision|interrupted data sets", "str", "D|D3|151|j1")
E("cell|A|D|D3|V1a|LM2|n_valid", "num", 200)
# (4) tier V (H04): the share of journeys with RUV not applied in Y, by KIT-242's reason (legacy: ruv_not_applied_reason; declared: factor status none)
rows = d_rows({})
vs = []
for r in range(1, 11):
    kind = "corrected" if r <= 5 else ("nuc" if r <= 8 else "nfs")
    vs.append(blank(tier="V", design="D3", cell="D3", omic="lipid", r=r, version="V1a", arm="LM2", smoke="FALSE", seconds=1,
                    status="uncorrected" if kind == "nuc" else "ok", ruv_status_y="corrected" if kind == "corrected" else "uncorrected",
                    ruv_not_applied_reason_y={"corrected": "NA", "nuc": "no_usable_controls", "nfs": "no_factor_supported"}[kind]))
    vs.append(blank(tier="V", design="D3", cell="D3", omic="lipid", r=r, version="V1a", arm="GK2", smoke="FALSE", seconds=1, status="ok", ruv_status_y="none" if r <= 4 else "ok"))
E = new_case("tier_v", rows + vs)
predict("tier_v", {}, expected)
E("cell|V_ruv_not_applied|V|D3|V1a|LM2|estimate", "num", 0.5)
E("cell|V_ruv_not_applied|V|D3|V1a|LM2|hi", "num", cp_upper(5, 10, 0.025))
E("cell|V_ruv_not_applied_no_usable_controls|V|D3|V1a|LM2|estimate", "num", 0.3)
E("cell|V_ruv_not_applied_no_factor_supported|V|D3|V1a|LM2|estimate", "num", 0.2)
E("cell|V_ruv_not_applied_no_factor_supported|V|D3|V1a|LM2|lo", "num", cp_lower(2, 10, 0.025))
E("cell|V_ruv_not_applied|V|D3|V1a|GK2|estimate", "num", 0.4)   # the declared arm: factor status none counts as not applied
E("cell|V_ruv_not_applied|V|D3|V1a|GK2|lo", "num", cp_lower(4, 10, 0.025))

# (H05) the controls-only mode: PC1 and NC1 only; no outcome row, no cell file, no 'outcome' in the console
for name, spec in (("controls_ok", {}), ("controls_f1", {"pc1_nonfail": {"D3": 57}})):   # K275-r2: 7.1 (143 of 200 failed: F1)
    folder = os.path.join(OUT, "cases", name)
    os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
    write_rows(os.path.join(folder, "rows", "d.tsv.gz"), d_rows(spec))
    index.append((name, "controls", "K275_CONTROLS_ONLY=1;K275_SCORE_R=200,100,20;K275_BOOT_B=200"))
    expected += [(name, f"exit|{name}", "exit", "0" if name == "controls_ok" else "1", ""), (name, f"stdout|{name}", "nosubstr", "outcome", ""),
                 (name, f"stdout|{name}", "substr", "PC1 D3 " + ("passes" if name == "controls_ok" else "FAILS"), ""),
                 (name, "absent|cells.tsv", "absent", "", ""), (name, "absent|decision.tsv", "absent", "", ""), (name, "present|controls-judged.tsv", "present", "", "")]

# ---- K275-r2 (K2-07, 8.1): a V1w NO row and NO V1a row gives "not determined (coverage)" for PC1 and never a pass or fail; V1w is read only by its own reported cell ---------
rows = [r_ for r_ in d_rows({}) if r_["version"] != "V1a"]
E = new_case("v1w_only", rows)
E("decision|outcome", "str", "Stop (PC1 not determined (coverage))")
for d_, _n in DESIGNS:
    E(f"cell|PC1|D|{d_}|V1a|NO|label", "str", "not determined (coverage)")
    E(f"cell|PC1|D|{d_}|V1a|NO|n_valid", "num", 0)
    E(f"cell|P_pre|D|{d_}|V1w|NO|estimate", "num", 0.25)
    E(f"cell|P_pre|D|{d_}|V1w|NO|n_valid", "num", 200)
    E(f"cell|A|D|{d_}|V1a|LM2|label", "str", "not determined (coverage)")
name = "controls_v1w_only"
folder = os.path.join(OUT, "cases", name)
os.makedirs(os.path.join(folder, "rows"), exist_ok=True)
write_rows(os.path.join(folder, "rows", "d.tsv.gz"), rows)
index.append((name, "controls", "K275_CONTROLS_ONLY=1;K275_SCORE_R=200,100,20;K275_BOOT_B=200"))
expected += [(name, f"exit|{name}", "exit", "1", ""), (name, f"stdout|{name}", "substr", "PC1 D3 not determined (coverage)", ""), (name, f"stdout|{name}", "nosubstr", "FAILS", ""),
             (name, "present|controls-judged.tsv", "present", "", "")]
# the genuine sample's counts beside PC1 (R2.4): printed with the labels of the controls-only read
expected += [("controls_ok", "stdout|controls_ok", "substr", "PC1 D3 passes (genuine below 0.5: %d of 200)" % GENUINE_BELOW["D3"], "")]

# ---- C02: the --g1 mode judges the G1 controls alone (exit 1 on F3); C05 (c): pi with a reference present and no matching pair -> coverage -----------------------
for name, spec in (("g1_pass", {}), ("g1_fail", {"g1": {("D3", "PC2"): 1.0}})):
    folder = os.path.join(OUT, "cases", name)
    os.makedirs(folder, exist_ok=True)
    write_g1(os.path.join(folder, "g1.tsv"), g1_rows(spec))
    index.append((name, "g1", "K275_BOOT_B=200"))
    for d, n_ in DESIGNS:
        for k_, t_ in (("PC2", math.sqrt((2 * n_ - 3) / (2 * n_ - 1))), ("NC2", 1.0)):
            val = spec.get("g1", {}).get((d, k_), t_)
            expected.append((name, f"g1|{d}|{k_}|label", "str", "passes" if abs(val - t_) <= M_MATCH + TOL else "FAILS", ""))
            expected.append((name, f"g1|{d}|{k_}|estimate", "num", val, "1e-8"))
    expected.append((name, f"exit|{name}", "exit", "0" if name == "g1_pass" else "1", ""))
rows = d_rows({})
for r_ in rows:
    if r_["tier"] == "D" and r_["design"] == "D3" and r_["version"] == "V1a" and r_["arm"] == "LM2":
        r_["tested_sig_y"] = "249_1_1"
E = new_case("pi_nopair", rows)
E("cell|pi|D|D3|V1a|LM2|label", "str", "not determined (coverage)")
E("cell|pi|D|D3|V1a|LM2|n_valid", "num", 0)
with open(os.path.join(OUT, "index.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["case", "mode", "env"])
    for row in index: w.writerow(row)
with open(os.path.join(OUT, "expected.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["case", "key", "kind", "expected", "tolerance"])
    for row in expected: w.writerow(row)
print(f"make-selftest: {len(index)} cases, {len(expected)} predictions, {len(COLS)} row columns")
