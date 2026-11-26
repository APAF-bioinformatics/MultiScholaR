# E16 (7.4, r1.23): the constructed row set for score-k188.R and the values predicted for it, written before the scorer runs.
# Standard library only. The rows are constructed (no simulated values): every safety value is 0 or a planted constant, so each
# prediction below follows from the ticket's 38.2 and 38.3 and from independent arithmetic (a t quantile and a beta quantile written
# here), never from the scorer. Usage: python3 -I make-e16.py <out dir>. Writes <out>/cases/<case>/rows.tsv.gz, cases.tsv, expected.tsv.
import csv, gzip, math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = open(os.path.join(HERE, "..", "sim-k188.R")).read()
BLIND = re.findall(r'"([^"]+)"', SIM[SIM.index("BLIND <- c("):SIM.index(")\n\nwarning_count")])
ASSAY = re.findall(r'"([^"]+)"', SIM[SIM.index("ASSAY_COLUMNS <- c("):SIM.index('ROW_COLUMNS <- c("part"')])
COLUMNS = ["part", "cell", "job", "scenario", "mask", "rep", "arm", "arm_scale", "scoreset", "assay", "contrast", "status", "class", "step",
           "seconds", "warnings", "samples_removed", "k", "pct", "n_controls", "complete_controls", "k_estimate", "e1_pass", "n_T", "n_tested",
           "lost"] + BLIND + ["q_method", "miss_share", "miss_after", "tie_rule_match"] + ASSAY + ["smoke"]
SCEN = ["N0", "N1", "E0", "E1", "N2", "E2"]
MASKS = ["none", "vd5", "vd20"]
ARMS = ["LA", "GA", "NO", "OR"]
FAMILY = 720.0

# ---- independent special functions ------------------------------------------------------------------------------------------
def betacf(a, b, x):
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; d = d if abs(d) > tiny else tiny
        c = 1 + aa / c; c = c if abs(c) > tiny else tiny
        d = 1 / d; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; d = d if abs(d) > tiny else tiny
        c = 1 + aa / c; c = c if abs(c) > tiny else tiny
        d = 1 / d; delta = d * c; h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return h

def betainc(a, b, x):
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
    if x < (a + 1) / (a + b + 2):
        return bt * betacf(a, b, x) / a
    return 1 - bt * betacf(b, a, 1 - x) / b

def bisect(f, lo, hi):
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(mid) < 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2

def qbeta(p, a, b):
    return bisect(lambda x: betainc(a, b, x) - p, 0.0, 1.0)

def qt(p, df):  # p above 0.5
    def cdf(t):
        tail = 0.5 * betainc(df / 2, 0.5, df / (df + t * t))
        return 1 - tail
    return bisect(lambda t: cdf(t) - p, 0.0, 200.0)

def clopper(x, n, level):
    return (0.0 if x == 0 else qbeta(1 - level, x, n - x + 1), 1.0 if x == n else qbeta(level, x + 1, n - x))

def mean_sd(v):
    n = len(v); m = sum(v) / n
    return m, (math.sqrt(sum((a - m) ** 2 for a in v) / (n - 1)) if n > 1 else float("nan"))

def classify(v, est, edge):
    """5.7: meets when the one-sided 97.5 % upper limit is below the edge, fails when the Bonferroni lower limit is at or above it."""
    n = len(v)
    if est == "S5":
        up = clopper(sum(v), n, 0.975)[1]; lo = clopper(sum(v), n, 1 - 0.025 / FAMILY)[0]
    else:
        m, sd = mean_sd(v); se = sd / math.sqrt(n)
        up = m + qt(0.975, n - 1) * se; lo = m - qt(1 - 0.025 / FAMILY, n - 1) * se
    return "meets" if up < edge else ("fails" if lo >= edge else "not determined")

def excess_class(diff):
    if len(diff) < 2: return "not determined (excess over OR)"
    m, sd = mean_sd(diff); se = sd / math.sqrt(len(diff))
    up = m + qt(0.975, len(diff) - 1) * se; lo = m - qt(1 - 0.025 / FAMILY, len(diff) - 1) * se
    return "meets (excess over OR)" if up < 0.025 else ("fails (excess over OR)" if lo >= 0.025 else "not determined (excess over OR)")

# ---- rows ------------------------------------------------------------------------------------------------------------------
def row(cell, job, scen, mask, rep, arm, **over):
    r = {c: "NA" for c in COLUMNS}
    nul = scen.startswith("N")
    r.update(part="L", cell=cell, job=job, scenario=scen, mask=mask, rep=rep, arm=arm, arm_scale="1", scoreset="T", assay="LCMS_Neg",
             contrast="B", status="ok", n_T=125, n_tested=125, lost=0, n_null=125 if nul else 100, null_p05=0, null_p01=0,
             any_call=0 if nul else 1, calls=0 if nul else 10, false_calls=0, true_calls=0 if nul else 10, n_planted=0 if nul else 10,
             bh_calls=0, bh_false=0, bh_any=0, bh_fdp=0, gross1=0, n_common=125, s5=0, or_s5=0, gross1_common=0, or_gross1_common=0,
             err_max=0.1, logfc_max=0.2, samples_removed=0, itsd_ok="TRUE" if arm != "OR" else "NA", norm_ok="TRUE" if arm != "OR" else "NA",
             da_status="ok", corr_removed_one_assay="")
    if arm == "LA":
        r.update(ruv_status="corrected", ruv_status_product="corrected", ruv_changed="TRUE", ruv_k_found=1, ruv_k_sweep=1, ruv_k_used=1,
                 ruv_k_cap=2, k=1, n_controls=20, complete_controls=18, ruv_n_controls=20, ruv_record_match="TRUE")
    if arm == "GA":
        r.update(ruv_status="ok", k=1, k_estimate=1)
    r.update({k: v for k, v in over.items()})
    return r

def crashed(cell, job, scen, mask, rep, arm, **over):
    r = row(cell, job, scen, mask, rep, arm)
    for c in BLIND + ["n_T", "n_tested", "lost", "n_null"]:
        r[c] = "NA"
    r.update(status="crashed", **{"class": "planted"}, da_status="failed")
    r.update(over)
    return r

def conds(cell, job="deciding"):
    return [(s, m) for s in (SCEN if cell == "L12" else SCEN[:4]) for m in MASKS]

def build(cell, n, plants=None, unit_rep=None, tweak=None):
    """All data rows of one (cell, deciding) set at R = n; `plants` maps (scenario, mask, arm) to a function (rep -> row or None)."""
    rows = []
    for s, m in conds(cell):
        for arm in ARMS:
            for rep in range(1, n + 1):
                if rep == unit_rep: continue
                make = (plants or {}).get((s, m, arm))
                r = make(rep) if make else None
                if r is None:
                    r = row(cell, "deciding", s, m, rep, arm)
                if tweak: tweak(r)
                rows.append(r)
    return rows

def unit_row(cell, rep):
    r = {c: "NA" for c in COLUMNS}
    r.update(part="L", cell=cell, job="deciding", scenario="all", mask="all", rep=rep, arm="unit", arm_scale="1", status="crashed",
             **{"class": "child process exit 3, no result"})
    return r

def write(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([r[c] for c in COLUMNS])

out = sys.argv[1]
expected = []   # (case, file, key, column, value)
def expect(case, file, key, column, value):
    expected.append((case, file, key, column, str(value)))
def pk(*a): return "|".join(a)
def cond_key(cell, job, scen, mask, arm, est): return "|".join([cell, job, scen, mask, arm, est])

# ---- case main: L12 (R = 400) and L24 (R = 100), unit row at data set 5 -----------------------------------------------------
R = 400
plants = {}
def mk(s, m, arm, f):
    plants[(s, m, arm)] = lambda rep: f(rep, s, m)
# P1: 19 crashed + the unit row = 20 of 400 (5 %, not above it): not "no result"
mk("N0", "none", "LA", lambda rep, s, m: crashed("L12", "deciding", s, m, rep, "LA", ruv_status="not_reached") if 11 <= rep <= 29 else None)
# P2: 20 crashed + the unit row = 21 of 400 (above 5 %): no result
mk("N0", "vd5", "LA", lambda rep, s, m: crashed("L12", "deciding", s, m, rep, "LA", ruv_status="not_reached") if 11 <= rep <= 30 else None)
# P3: a phi event (uncorrected, no k) whose DA then fails: gate first, not an event; 10 + unit
mk("E0", "none", "LA", lambda rep, s, m: crashed("L12", "deciding", s, m, rep, "LA", ruv_status="uncorrected", ruv_k_found="NA") if 41 <= rep <= 50 else None)
# P4: 30 phi events (every assay uncorrected, no k, journey uncorrected): not in the gate
mk("N1", "none", "LA", lambda rep, s, m: row("L12", "deciding", s, m, rep, "LA", status="uncorrected", ruv_status="uncorrected", ruv_status_product="uncorrected",
    ruv_k_found="NA", ruv_changed="FALSE", ruv_n_controls=0) if 61 <= rep <= 90 else None)
# P5: 5 refusals at the RUV step (not events, in the gate)
mk("E1", "none", "LA", lambda rep, s, m: crashed("L12", "deciding", s, m, rep, "LA", status="refused", ruv_status="refused", **{"class": "multischolar_refusal"}) if 101 <= rep <= 105 else None)
# P5b: 21 dropped journeys (a DA table exists; in the gate; a phi event by 38.2 (6) but gate first) + unit = 22
mk("E1", "vd5", "LA", lambda rep, s, m: row("L12", "deciding", s, m, rep, "LA", status="dropped", ruv_status="dropped", ruv_status_product="dropped",
    ruv_k_found="NA") if 111 <= rep <= 131 else None)
# P5c: 19 dropped journeys + the unit row = 20 of 400 (5 %, not above it): the stand-alone "dropped" condition just below the gate
mk("E1", "vd20", "LA", lambda rep, s, m: row("L12", "deciding", s, m, rep, "LA", status="dropped", ruv_status="dropped", ruv_status_product="dropped",
    ruv_k_found="NA") if 141 <= rep <= 159 else None)
# P6: the oracle fails S1 in N1 vd20 (0.2), GA returns a table at data set 1 only (not_run elsewhere): the excess path with one pair
mk("N1", "vd20", "OR", lambda rep, s, m: row("L12", "deciding", s, m, rep, "OR", null_p05=0.2))
mk("N1", "vd20", "GA", lambda rep, s, m: row("L12", "deciding", s, m, rep, "GA", null_p05=0.2) if rep == 1 else
    row("L12", "deciding", s, m, rep, "GA", status="not_run", **{"class": "planted not_run"}))
rows = build("L12", R, plants, unit_rep=5) + [unit_row("L12", 5)]
rows += build("L24", 100)
def block(job, arm, n, over=None, or_over=None):
    """The lp or las job of L12 (base scenarios, all masks) at data sets 1..n: the arm and the oracle. Its oracle rows are the deciding job's."""
    out_rows = []
    for s_, m_ in conds("L12"):
        if s_ not in SCEN[:4]: continue
        for rep in range(1, n + 1):
            f = (over or {}).get((s_, m_))
            out_rows.append(f(rep) if f else row("L12", job, s_, m_, rep, arm))
            g = plants.get((s_, m_, "OR")); r = g(rep) if g else row("L12", "deciding", s_, m_, rep, "OR")
            r["job"] = job
            if or_over and (s_, m_) in or_over: r.update(or_over[(s_, m_)])
            out_rows.append(r)
    return out_rows
rows += block("lp", "LP", 100, {("N0", "none"): lambda rep: row("L12", "lp", "N0", "none", rep, "LP", null_p05=0.04)})
rows += block("las", "LAs", 100, {("N0", "vd20"): lambda rep: row("L12", "las", "N0", "vd20", rep, "LAs", null_p05=0.02)})
write(os.path.join(out, "cases", "main", "rows.tsv.gz"), rows)
# M-repair (LA against LP) and the LAs-LA table: paired per data set on the data sets where both returned a DA table
# LA (N0, none): crashed at 19 data sets and the unit row at 5 -> 100 - 20 = 80 paired; LP's S1 is 0.04 and LA's 0: mean difference -0.04
expect("main", "partm-repair.tsv", pk("L12", "N0", "none", "S1", "LA", "LP"), "n", 80)
expect("main", "partm-repair.tsv", pk("L12", "N0", "none", "S1", "LA", "LP"), "mean_diff", -0.04)
# LAs (N0, vd20): LA unplanted there, the unit row removes data set 5 -> 99 paired; LAs' S1 is 0.02: mean difference +0.02
expect("main", "partm-levels.tsv", pk("L12", "N0", "vd20", "S1", "LAs", "LA"), "n", 99)
expect("main", "partm-levels.tsv", pk("L12", "N0", "vd20", "S1", "LAs", "LA"), "mean_diff", 0.02)

def exp_gate(case, cell, scen, mask, arm, n_gate, r_out, estims, nres, r=400):
    for e in estims:
        k = cond_key(cell, "deciding", scen, mask, arm, e)
        expect(case, "conditions.tsv", k, "n_gate", n_gate); expect(case, "conditions.tsv", k, "r_out", r_out)
        expect(case, "conditions.tsv", k, "r_used", r - r_out); expect(case, "conditions.tsv", k, "no_result", "TRUE" if nres else "FALSE")
        if nres:
            expect(case, "conditions.tsv", k, "classification", "no result")

# P1 (S1, S2, S5 of N0): 380 values at 0; best adds 20 zeros, worst 20 ones
base = [0.0] * 380
for e, edge in (("S1", .075), ("S2", .075), ("S5", .01)):
    k = cond_key("L12", "deciding", "N0", "none", "LA", e)
    expect("main", "conditions.tsv", k, "n_gate", 20); expect("main", "conditions.tsv", k, "r_used", 380); expect("main", "conditions.tsv", k, "r_out", 20)
    expect("main", "conditions.tsv", k, "no_result", "FALSE")
    expect("main", "conditions.tsv", k, "classification", classify(base, e, edge))
    best = classify(base + [0.0] * 20, e, edge); worst = classify(base + [1.0] * 20, e, edge)
    expect("main", "conditions.tsv", k, "class_best", best); expect("main", "conditions.tsv", k, "class_worst", worst)
    expect("main", "conditions.tsv", k, "sensitive", "TRUE" if best != worst else "FALSE")
exp_gate("main", "L12", "N0", "vd5", "LA", 21, 21, ["S1", "S2", "S5"], True)
for e in ("S1", "S2", "S5"):
    k = cond_key("L12", "deciding", "N0", "vd5", "LA", e)
    expect("main", "conditions.tsv", k, "class_best", "no result"); expect("main", "conditions.tsv", k, "class_worst", "no result")
    expect("main", "conditions.tsv", k, "sensitive", "FALSE")
for e in ("S3", "S5"):
    k = cond_key("L12", "deciding", "E1", "vd20", "LA", e)
    expect("main", "conditions.tsv", k, "n_gate", 20); expect("main", "conditions.tsv", k, "r_out", 1)
    expect("main", "conditions.tsv", k, "no_result", "FALSE")
expect("main", "phi.tsv", pk("L12", "deciding", "E1", "vd20", "LA"), "events", 0)
expect("main", "phi.tsv", pk("L12", "deciding", "E1", "vd20", "LA"), "events_also_gate", 19)
# phi: P3 gate first; P4 events; P5 refusals; P5b dropped
expect("main", "phi.tsv", pk("L12", "deciding", "E0", "none", "LA"), "events", 0)
expect("main", "phi.tsv", pk("L12", "deciding", "E0", "none", "LA"), "events_also_gate", 10)
expect("main", "phi.tsv", pk("L12", "deciding", "E0", "none", "LA"), "r_reached", 399)
expect("main", "phi.tsv", pk("L12", "deciding", "N1", "none", "LA"), "events", 30)
expect("main", "phi.tsv", pk("L12", "deciding", "N1", "none", "LA"), "events_also_gate", 0)
expect("main", "phi.tsv", pk("L12", "deciding", "N1", "none", "LA"), "r_declared", 400)
expect("main", "phi.tsv", pk("L12", "deciding", "N1", "none", "LA"), "r_reached", 399)
expect("main", "phi.tsv", pk("L12", "deciding", "N1", "none", "LA"), "phi", 30 / 399)
expect("main", "phi.tsv", pk("L12", "deciding", "N1", "none", "LA"), "phi_all", 30 / 400)
expect("main", "phi.tsv", pk("L12", "deciding", "N1", "none", "LA"), "lan0_mark", "TRUE")
expect("main", "phi.tsv", pk("L12", "deciding", "E1", "none", "LA"), "events", 0)
expect("main", "phi.tsv", pk("L12", "deciding", "E1", "vd5", "LA"), "events", 0)
expect("main", "phi.tsv", pk("L12", "deciding", "E1", "vd5", "LA"), "events_also_gate", 21)
exp_gate("main", "L12", "E1", "vd5", "LA", 22, 1, ["S3", "S5"], True)   # r_out: dropped journeys have a table; only the unit row has none
exp_gate("main", "L12", "E1", "none", "LA", 6, 6, ["S3", "S5"], False)
# P6: the excess path with one pair; GA's other 398 data sets with an oracle table put at 0 (best) and 1 (worst)
k = cond_key("L12", "deciding", "N1", "vd20", "GA", "S1")
expect("main", "conditions.tsv", k, "oracle_check", "TRUE"); expect("main", "conditions.tsv", k, "classification", excess_class([0.0]))
expect("main", "conditions.tsv", k, "n_paired_dropped", 399)
best = excess_class([0.0] + [0 - 0.2] * 398); worst = excess_class([0.0] + [1 - 0.2] * 398)
expect("main", "conditions.tsv", k, "class_best", best); expect("main", "conditions.tsv", k, "class_worst", worst)
expect("main", "conditions.tsv", k, "sensitive", "TRUE" if best != worst else "FALSE")
# L24: reported, R = 100, S5 cannot meet
k = cond_key("L24", "deciding", "N0", "none", "LA", "S5")
up0 = clopper(0, 100, 0.975)[1]
expect("main", "conditions.tsv", k, "role", "reported"); expect("main", "conditions.tsv", k, "R", 100)
expect("main", "conditions.tsv", k, "label", "cannot meet S5 at this R (zero events give an upper limit of %.4f)" % up0)
expect("main", "conditions.tsv", k, "classification", "not determined")
# statuses (journeys, by arm, status and class)
expect("main", "statuses.tsv", pk("L", "L12", "deciding", "LA", "1", "crashed", "planted"), "n", 19 + 20 + 10)
expect("main", "statuses.tsv", pk("L", "L12", "deciding", "LA", "1", "crashed", "unit: child process exit 3, no result"), "n", 18)
expect("main", "statuses.tsv", pk("L", "L12", "deciding", "LA", "1", "refused", "multischolar_refusal"), "n", 5)
# GA's not_run journeys in (N1, vd20): data sets 2 to 400 are planted not_run (399); the unit row at data set 5 replaces every L12 deciding
# journey of data set 5 (all arms and conditions) by a crashed one, so data set 5 is not among them: 399 - 1 = 398. (First written as 399, which
# forgot the unit row; corrected after send 2's observation, the scorer's count following from the planted rows.)
expect("main", "statuses.tsv", pk("L", "L12", "deciding", "GA", "1", "not_run", "planted not_run"), "n", 398)
# L24 is in no verdict
expect("main", "conditions.tsv", "absent:L24", "deciding_rows", 0)

# ---- outcome cases: L12 only, R = 400, no unit row, no plants beyond the case's own -------------------------------------------
def outcome_case(case, plants, line_start, also=()):
    write(os.path.join(out, "cases", case, "rows.tsv.gz"), build("L12", R, plants))
    expect(case, "outcome.txt", "line", "startswith", line_start)
    for part in also: expect(case, "outcome.txt", "line", "contains", part)
def uncorrected_all(arm="LA"):
    p = {}
    for s, m in conds("L12"):
        p[(s, m, arm)] = (lambda s, m: lambda rep: row("L12", "deciding", s, m, rep, arm, status="uncorrected", ruv_status="uncorrected",
            ruv_status_product="uncorrected", ruv_k_found="NA", ruv_changed="FALSE", ruv_n_controls=0) if 61 <= rep <= 90 else None)(s, m)
    return p
def bad_s1(arm):
    return {("N0", "none", arm): (lambda rep: row("L12", "deciding", "N0", "none", rep, arm, null_p05=0.5))}
outcome_case("k1", {}, "lipid_metab: K1 ")
p = uncorrected_all(); outcome_case("k2phi", p, "lipid_metab: K2 ", ["K1's non-correction clause failed", "GA's largest upper limit"])
p = uncorrected_all(); p.update(bad_s1("GA")); outcome_case("val248", p, "lipid_metab: no outcome by rule: K1 blocked by non-correction")
p = bad_s1("LA"); p.update(bad_s1("GA")); outcome_case("k3", p, "lipid_metab: K3 ")
p = bad_s1("LA"); p.update(bad_s1("GA")); p.update(bad_s1("NO")); outcome_case("k4", p, "lipid_metab: K4 ")
for case, safe in (("k1", ("safe", "safe", "safe")), ("k2phi", ("safe", "safe", "safe")), ("val248", ("safe", "unsafe", "safe")),
                   ("k3", ("unsafe", "unsafe", "safe")), ("k4", ("unsafe", "unsafe", "unsafe"))):
    for arm, v in zip(("LA", "GA", "NO"), safe):
        expect(case, "verdicts.tsv", "lipid_metab|" + arm, "verdict", v)
# phi of the k2phi case: 30 events of 400 reached in every LA condition; K1's clause fails by the upper limit
up = clopper(30, 400, 0.975)[1]
expect("k2phi", "phi.tsv", pk("L12", "deciding", "N0", "none", "LA"), "events", 30)
expect("k2phi", "phi.tsv", pk("L12", "deciding", "N0", "none", "LA"), "phi_upper", up)
expect("k1", "phi.tsv", pk("L12", "deciding", "N0", "none", "LA"), "events", 0)

# ---- a column of the letters T and F must stay text (smoke check, defect 1): LA's class is "T", every step "F", scoreset is "T" throughout
def coerce(r):
    if r["arm"] == "LA":
        r["class"] = "T"
    r["step"] = "F"
write(os.path.join(out, "cases", "tf", "rows.tsv.gz"), build("L12", R, {}, tweak=coerce))
expect("tf", "statuses.tsv", pk("L", "L12", "deciding", "LA", "1", "ok", "T"), "n", 18 * R)
expect("tf", "outcome.txt", "line", "startswith", "lipid_metab: K1 ")
write(os.path.join(out, "cases", "smokeflag", "rows.tsv.gz"), [row("L12", "deciding", "N0", "none", 1, "LA", smoke="1")])
# ---- E16-1 (script audit 2): the cases the first set did not exercise --------------------------------------------------------
def over(arm, scen, mask, **kw):
    return {(scen, mask, arm): (lambda rep: row("L12", "deciding", scen, mask, rep, arm, **kw))}
# k2la: the canonical K2 row, LA not safe (S1 0.5), GA and NO safe
outcome_case("k2la", bad_s1("LA"), "lipid_metab: K2 ")
for arm, v in (("LA", "unsafe"), ("GA", "safe"), ("NO", "safe")): expect("k2la", "verdicts.tsv", "lipid_metab|" + arm, "verdict", v)
# k1eff: LA and GA safe, LA much lower than GA in one E condition (power 0.7 against 1: F1 0.3 >= 0.20) -> K2 with the efficiency clause failed;
# the F1 and F2 words: F1 0.3 "much lower", 0.1 "lower", 0 "comparable"; F2 0.01 "less precise", 0.06 "much less precise", 0 "comparable"
p = {}
p.update(over("LA", "E0", "none", true_calls=7)); p.update(over("LA", "E1", "none", true_calls=9))
p.update(over("LA", "N0", "vd5", gross1_common=0.01)); p.update(over("LA", "N1", "vd5", gross1_common=0.06))
outcome_case("k1eff", p, "lipid_metab: K2 ", ["K1's efficiency clause failed"])
for (sc, mk, e), word in ((("E0", "none", "F1"), "much lower"), (("E1", "none", "F1"), "lower"), (("E0", "vd5", "F1"), "comparable"),
                          (("N0", "vd5", "F2"), "less precise"), (("N1", "vd5", "F2"), "much less precise"), (("N0", "none", "F2"), "comparable")):
    expect("k1eff", "conditions.tsv", cond_key("L12", "deciding", sc, mk, "LA", e), "classification", word)
expect("k1eff", "verdicts.tsv", "lipid_metab|LA", "verdict", "safe"); expect("k1eff", "verdicts.tsv", "lipid_metab|GA", "verdict", "safe")
# nss: one LA condition "not determined" (140 of 400 data sets at 0.2: mean 0.07, upper limit above 0.075, Bonferroni lower limit below it)
vals = [0.2] * 140 + [0.0] * 260
assert classify(vals, "S1", .075) == "not determined"
p = {("N0", "none", "LA"): (lambda rep: row("L12", "deciding", "N0", "none", rep, "LA", null_p05=0.2 if rep <= 140 else 0))}
outcome_case("nss", p, "lipid_metab: K2 ")
expect("nss", "conditions.tsv", cond_key("L12", "deciding", "N0", "none", "LA", "S1"), "classification", "not determined")
expect("nss", "verdicts.tsv", "lipid_metab|LA", "verdict", "not shown safe")
# s5fail: 15 of 400 data sets with an attributable extreme error (N0, vd20): S5 fails
ev = [1.0] * 15 + [0.0] * 385
assert classify(ev, "S5", .01) == "fails"
p = {("N0", "vd20", "LA"): (lambda rep: row("L12", "deciding", "N0", "vd20", rep, "LA", s5=1 if 201 <= rep <= 215 else 0))}
outcome_case("s5fail", p, "lipid_metab: K2 ")
expect("s5fail", "conditions.tsv", cond_key("L12", "deciding", "N0", "vd20", "LA", "S5"), "classification", "fails")
expect("s5fail", "verdicts.tsv", "lipid_metab|LA", "verdict", "unsafe")
# excess: the oracle fails S1 (0.2) in (N1, vd20); NO and GA at 0 have an excess of -0.2 (meets), LA at 0.5 an excess of +0.3 (fails)
p = {("N1", "vd20", "OR"): (lambda rep: row("L12", "deciding", "N1", "vd20", rep, "OR", null_p05=0.2)),
     ("N1", "vd20", "LA"): (lambda rep: row("L12", "deciding", "N1", "vd20", rep, "LA", null_p05=0.5))}
outcome_case("excess", p, "lipid_metab: K2 ")
for arm, want in (("NO", excess_class([-0.2] * 400)), ("GA", excess_class([-0.2] * 400)), ("LA", excess_class([0.3] * 400))):
    expect("excess", "conditions.tsv", cond_key("L12", "deciding", "N1", "vd20", arm, "S1"), "classification", want)
assert excess_class([-0.2] * 400) == "meets (excess over OR)" and excess_class([0.3] * 400) == "fails (excess over OR)"
# tiny stop cases at the ranges 2,2 (K188_SCORE_R): a norm_ok FALSE row; a partial-match guard row; a data set outside the declared range;
# an oracle row of the lp job that differs from the deciding job's (E5)
write(os.path.join(out, "cases", "normstop", "rows.tsv.gz"), [row("L12", "deciding", "N0", "none", 1, "LA", norm_ok="FALSE")])
write(os.path.join(out, "cases", "guardstop", "rows.tsv.gz"), [row("L12", "deciding", "N0", "none", 1, "LA", status="instrument",
    **{"class": "partial_match_guard: partial-match guard: partial match of 'da' to 'da_reached'"})])
# a unit-level guard stop: the child ended outside any journey; the parent's unit row is `instrument` with the guard's text
write(os.path.join(out, "cases", "guardunit", "rows.tsv.gz"), [{**row("L12", "deciding", "N0", "none", 1, "unit"), "scenario": "N0", "status": "instrument", "scoreset": "NA", "arm": "unit",
    "class": "child process exit 1, no result Error: partial-match guard: partial match of 'partial_probe' to 'partial_probe_value'"}])
write(os.path.join(out, "cases", "rangestop", "rows.tsv.gz"), [row("L12", "deciding", "N0", "none", 3, "LA")])
tiny = build("L12", 2)
for s_, m_ in conds("L12"):
    if s_ in SCEN[:4]:
        for rep in (1, 2):
            tiny.append(row("L12", "lp", s_, m_, rep, "LP"))
            tiny.append(row("L12", "lp", s_, m_, rep, "OR", null_p05=0.5 if (s_, m_, rep) == ("N0", "none", 1) else 0))
write(os.path.join(out, "cases", "e5diff", "rows.tsv.gz"), tiny)
# ---- stop cases (tiny) --------------------------------------------------------------------------------------------------------
write(os.path.join(out, "cases", "itsd", "rows.tsv.gz"), [row("L12", "deciding", "N0", "none", 1, "LA", itsd_ok="FALSE")])
write(os.path.join(out, "cases", "gap", "rows.tsv.gz"), [row("L12", "deciding", "N0", "none", 1, "LA")])
with open(os.path.join(out, "cases.tsv"), "w") as f:
    f.write("case\trowsdir\tcells\tranges\texpect\n")
    f.write("main\tmain\tL12:deciding,L24:deciding,L12:lp,L12:las\t400,100\tok\n")
    for case in ("k1", "k2phi", "val248", "k3", "k4"):
        f.write("%s\t%s\tL12:deciding\t400,100\tok\n" % (case, case))
    for case in ("k2la", "k1eff", "nss", "s5fail", "excess"):
        f.write("%s\t%s\tL12:deciding\t400,100\tok\n" % (case, case))
    f.write("normstop\tnormstop\tL12:deciding\t2,2\tstop:norm_ok\n")
    f.write("guardstop\tguardstop\tL12:deciding\t2,2\tstop:partial-match guard\n")
    f.write("guardunit\tguardunit\tL12:deciding\t2,2\tstop:partial-match guard\n")
    f.write("rangestop\trangestop\tL12:deciding\t2,2\tstop:outside the declared range\n")
    f.write("e5diff\te5diff\tL12:deciding,L12:lp\t2,2\tstop:E5 differs\n")
    f.write("tf\ttf\tL12:deciding\t400,100\tok\n")
    f.write("smokeflag\tsmokeflag\tL12:deciding\t400,100\tstop:smoke flag\n")
    f.write("itsd\titsd\tL12:deciding\t400,100\tstop:itsd_ok\n")
    f.write("gap\tgap\tL12:deciding\t400,100\tstop:absent without a unit row\n")
    f.write("cellgap\tk1\tL12:deciding,L24:deciding\t400,100\tstop:L24 deciding\n")   # the k1 rows, the declared set naming L24 (SA-P14)
with open(os.path.join(out, "expected.tsv"), "w") as f:
    f.write("case\tfile\tkey\tcolumn\tvalue\n")
    for e in expected:
        f.write("\t".join(e) + "\n")
print("cases written; %d predictions" % len(expected))
