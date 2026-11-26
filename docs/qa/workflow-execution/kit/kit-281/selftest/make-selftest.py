# KIT-281 scorer self-test (the manner of KIT-188's E16): constructed row files for score-k281.R and the values PREDICTED for them,
# written before the scorer ever runs on them. Python standard library only (run as `python3 -I make-selftest.py <selftest dir>`).
# No simulation and no emulation of R's random-number generator: every constructed table is DEGENERATE (every replicate of a setting
# alike, so every bootstrap resample and every jackknife replicate equals the estimate and every limit equals the estimate), except the
# one case (`partial`) whose only unpredictable quantity (the count of non-finite resamples) is asserted as "> 0".
# Every prediction is derived below from the planted counts by arithmetic with the SAME operations as the ticket's definitions
# (a ratio of sums, ratios of those ratios), never from the scorer. IEEE doubles are correctly rounded, so a ratio of totals equals the
# ratio of the per-replicate-mask counts: the boundary cases are decided exactly as written (no tolerance).
# Writes: <dir>/cases/<case>/<job>/{A,B,C,D}.tsv.gz, metrics.tsv, range.tsv, provenance-sim.json (+ determinism.tsv, instrument-stop.txt),
#         <dir>/cases.tsv, <dir>/expected.tsv, <dir>/g3-cases.tsv, <dir>/requests/*.json and <dir>/requests.tsv (the G3 request-check cases).
import csv, gzip, io, json, math, os, shutil, sys
from fractions import Fraction

HERE = sys.argv[1]
CELLS = ["A", "B", "C", "D"]
MASKS = ["vd5", "vd10", "vd25"]
REPS = [401, 402, 403, 404]
ARMS = ["S", "B", "B1", "N"]
NAN = float("nan")
RHO_MARGIN, FDP_LIMIT, GUARD_MARGIN, POWER_MARGIN, PC_MARGIN = 1.10, 0.05, 0.0025, 0.005, 1.20
B_BOOT = 1000
ROW_COLUMNS = ["cell", "mask", "r", "arm", "subset", "class", "stratum", "ell", "band", "ell0", "n", "n_tested", "n_p05", "n_call", "n_bh",
               "sum_err", "sum_s_err", "sum_sgn_err"]
METRIC_COLUMNS = ["cell", "mask", "r", "arm", "kept", "added", "calls", "t_r", "mono_ok", "qmethod", "threshold", "viol", "e10_dev", "e3_pass",
                  "n_planted", "n_down_all", "n1_band_all", "n_onoff_all", "calls_both_s", "calls_new_s", "gained", "lost", "seconds"]


# ---- the planted counts of one replicate-mask (all in "per replicate-mask" units) --------------------------------------------------
def default(cell, mask, r):
    return dict(
        Bn=1000, Bp=50,                       # B's whole list, null proteins: tested, p < 0.05 (alpha_B(H_B) = 0.05)
        B_bhV=10, B_bhT=200,                  # B: BH calls among null (false) and among planted (true; split half up, half down)
        B1n=1000, B1p=50, B1_bhV=10, B1_bhT=200,
        flag_n=20, flag_p=1, flag_err=0,      # B1's whole list: null proteins flagged one-sided (one_treated, ell 1): tested, p < 0.05, sum of signed error
        gA0n=1900, gA0p=95, gA1n=100, gA1p=5,  # B1's g_added (H_B1 without H_B): null, ell 0 and ell 1: tested, p < 0.05
        add={a: (100, 5) for a in ARMS},      # `added` (against O), null, ell 1: (tested, p05), per arm (Phi's numerator)
        kb={a: (1000, 50) for a in ARMS},     # `kept_both`, null: (tested, p05), per arm (Phi's denominator)
        planted=600, g_added_rows=True)


def build_spec(over):
    spec = {}
    for cell in CELLS:
        for mask in MASKS:
            for r in REPS:
                d = default(cell, mask, r)
                over(d, cell, mask, r)
                spec[(cell, mask, r)] = d
    return spec


# ---- writing the tables --------------------------------------------------------------------------------------------------------
def row(cell, mask, r, arm, subset, cls, stratum, ell, band, ell0, n, tested, p05, call, bh, err=0, serr=0, sgn=0):
    return dict(cell=cell, mask=mask, r=r, arm=arm, subset=subset, **{"class": cls}, stratum=stratum, ell=ell, band=band, ell0=ell0, n=n,
                n_tested=tested, n_p05=p05, n_call=call, n_bh=bh, sum_err=err, sum_s_err=serr, sum_sgn_err=sgn)


def unit_rows(cell, mask, r, d):
    out = []
    for arm in ARMS:
        n_whole, p_whole = (d["B1n"], d["B1p"]) if arm == "B1" else (d["Bn"], d["Bp"])
        bhV, bhT = (d["B1_bhV"], d["B1_bhT"]) if arm == "B1" else (d["B_bhV"], d["B_bhT"])
        up, down = bhT - bhT // 2, bhT // 2
        out.append(row(cell, mask, r, arm, "whole", "null", "none", 0, 0, -1, n_whole, n_whole, p_whole, 0, bhV))
        if arm == "B1":
            out.append(row(cell, mask, r, arm, "whole", "null", "one_treated", 1, 0, -1, d["flag_n"], d["flag_n"], d["flag_p"], 0, 0,
                           serr=d["flag_err"]))
        half = d["planted"] // 2
        out.append(row(cell, mask, r, arm, "whole", "up", "none", 0, 0, -1, half, half, half, half // 2, up))
        out.append(row(cell, mask, r, arm, "whole", "down", "none", 0, 0, -1, half, half, half, half // 2, down))
        n1, p1 = d["add"][arm]
        out.append(row(cell, mask, r, arm, "added", "null", "none", 1, -1, -1, n1, n1, p1, 0, 0))
        nk, pk = d["kb"][arm]
        out.append(row(cell, mask, r, arm, "kept_both", "null", "none", 0, -1, -1, nk, nk, pk, 0, 0))
        if arm == "B1" and d["g_added_rows"]:
            if d["gA0n"] or d["gA0p"]:
                out.append(row(cell, mask, r, arm, "g_added", "null", "none", 0, -1, -1, d["gA0n"], d["gA0n"], d["gA0p"], 0, 0))
            if d["gA1n"] or d["gA1p"]:
                out.append(row(cell, mask, r, arm, "g_added", "null", "one_treated", 1, -1, -1, d["gA1n"], d["gA1n"], d["gA1p"], 0, 0))
    return out


def unit_metrics(cell, mask, r, arm, d, over_metric):
    m = dict(cell=cell, mask=mask, r=r, arm=arm, kept=3000, added=0, calls=0, t_r=0, mono_ok="TRUE", qmethod="qvalue", threshold=12.5, viol=0,
             e10_dev=0, e3_pass="TRUE", n_planted=d["planted"], n_down_all=d["planted"] // 2, n1_band_all=0, n_onoff_all=0, calls_both_s=0,
             calls_new_s=0, gained=0, lost=0, seconds="NA" if arm != "S" else "1.5")
    over_metric(m)
    return m


def write_tsv(path, rows, columns, gz=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    f = io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0), newline="") if gz else open(path, "w", newline="")   # mtime 0: reproducible bytes
    with f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow([r[c] for c in columns])


PROVENANCE_OK = {"ticket": "KIT-281", "status": "complete", "smoke": False, "counted": False}


def write_case(case, spec, jobs=None, reps=REPS, over_metric=lambda m: None, row_hook=None, provenance=None, extra_files=None, drop=None,
               range_complete=None, determinism=False):
    """One case: `jobs` maps a job name to the replicates it holds (default one job `job1` with all `reps`)."""
    jobs = jobs or {"job1": reps}
    for job, jreps in jobs.items():
        base = os.path.join(HERE, "cases", case, job)
        for cell in CELLS:
            rows, metrics = [], []
            for mask in MASKS:
                for r in jreps:
                    if drop and (cell, mask, r) in drop:
                        continue
                    d = spec[(cell, mask, r)]
                    rows += unit_rows(cell, mask, r, d)
                    metrics += [unit_metrics(cell, mask, r, arm, d, lambda m: over_metric(m)) for arm in ARMS]
            if row_hook:
                row_hook(cell, rows, metrics)
            write_tsv(os.path.join(base, cell + ".tsv.gz"), rows, ROW_COLUMNS, gz=True)
            write_tsv(os.path.join(base, "metrics.tsv") if cell == "A" else os.path.join(base, "metrics-" + cell + ".part"), metrics, METRIC_COLUMNS)
        # one metrics.tsv per job: concatenate the per-cell parts
        parts = [os.path.join(base, "metrics.tsv")] + [os.path.join(base, "metrics-%s.part" % c) for c in CELLS[1:]]
        lines = open(parts[0]).read().splitlines()
        for p in parts[1:]:
            lines += open(p).read().splitlines()[1:]
            os.remove(p)
        open(parts[0], "w").write("\n".join(lines) + "\n")
        rng = []
        for cell in CELLS:
            have = {r for r in jreps if not (drop and any((cell, m, r) in drop for m in MASKS))}
            rng.append(dict(cell=cell, first=min(jreps), last=max(jreps), expected=len(jreps),
                            complete=len(have) if range_complete is None else range_complete))
        write_tsv(os.path.join(base, "range.tsv"), rng, ["cell", "first", "last", "expected", "complete"])
        json.dump(provenance if provenance is not None else PROVENANCE_OK, open(os.path.join(base, "provenance-sim.json"), "w"), indent=2)
        if determinism:
            write_tsv(os.path.join(base, "determinism.tsv"), [dict(cell="A", r=r, rows_identical="TRUE", metrics_identical="TRUE") for r in jreps[:5]],
                      ["cell", "r", "rows_identical", "metrics_identical"])
        for name, text in (extra_files or {}).items():
            open(os.path.join(base, name), "w").write(text)


# ---- the independent predictions -----------------------------------------------------------------------------------------------------
def div(a, b):
    return a / b if b != 0 else NAN


def agg(spec, keys):
    """Sums of the planted counts over (cell, mask) in `keys` and all replicates (one replicate-mask per entry)."""
    s = dict(units=0, Bn=0, Bp=0, B_bhV=0, B_bhT=0, B1_bhV=0, B1_bhT=0, B1n=0, B1p=0, flag_n=0, flag_err=0, gAn=0, gAp=0, gA1n=0, gA1p=0, gA0n=0,
             gA0p=0, planted=0, add={a: [0, 0] for a in ARMS}, kb={a: [0, 0] for a in ARMS}, gadd_units=0)
    for (cell, mask, r), d in spec.items():
        if (cell, mask) not in keys:
            continue
        s["units"] += 1
        for k in ("Bn", "Bp", "B_bhV", "B_bhT", "B1_bhV", "B1_bhT", "B1n", "B1p", "flag_n", "flag_err", "planted"):
            s[k] += d[k]
        if d["g_added_rows"]:
            s["gA0n"] += d["gA0n"]; s["gA0p"] += d["gA0p"]; s["gA1n"] += d["gA1n"]; s["gA1p"] += d["gA1p"]
        for a in ARMS:
            s["add"][a][0] += d["add"][a][0]; s["add"][a][1] += d["add"][a][1]
            s["kb"][a][0] += d["kb"][a][0]; s["kb"][a][1] += d["kb"][a][1]
    s["gAn"] = s["gA0n"] + s["gA1n"]
    s["gAp"] = s["gA0p"] + s["gA1p"]
    return s


def stats(s):
    """The deciding quantity, the guards, F5', Phi from summed counts: ratios of sums, with the ticket's operations in the ticket's order."""
    alpha_whole = div(s["Bp"], s["Bn"])
    out = {}
    out["rhog.B1"] = div(div(s["gAp"], s["gAn"]), alpha_whole)
    out["rhog_l1.B1"] = div(div(s["gA1p"], s["gA1n"]), alpha_whole)
    out["rhog_l0.B1"] = div(div(s["gA0p"], s["gA0n"]), alpha_whole)
    out["B1.fdr_bh"] = div(s["B1_bhV"], s["B1_bhV"] + s["B1_bhT"])
    out["B.fdr_bh"] = div(s["B_bhV"], s["B_bhV"] + s["B_bhT"])
    out["dfdrbh.B1"] = out["B1.fdr_bh"] - out["B.fdr_bh"]
    out["B.powbh"] = div(s["B_bhT"], s["planted"])
    out["B1.powbh"] = div(s["B1_bhT"], s["planted"])
    out["dpowbh.B1"] = out["B.powbh"] - out["B1.powbh"]
    out["B1.m6_l_H"] = div(s["flag_err"], s["flag_n"])
    out["B1.gadd_per_r"] = div(s["gAn"], s["units"])
    out["B.alpha_whole"] = alpha_whole
    rho = {a: div(div(s["add"][a][1], s["add"][a][0]), div(s["kb"][a][1], s["kb"][a][0])) for a in ARMS}
    for a in ("S", "B", "B1"):
        out["Phi." + a] = div(rho[a], rho["N"])
    return out


def fmt(x):
    if isinstance(x, bool):
        return "TRUE" if x else "FALSE"
    if isinstance(x, float) and math.isnan(x):
        return "NA"
    return repr(x) if isinstance(x, float) else str(x)


def state(upper, limit):
    return "pass" if (not math.isnan(upper) and upper <= limit) else "fail"


def label(lower, upper):
    if not math.isnan(upper) and upper <= RHO_MARGIN:
        return "within margin"
    if not math.isnan(lower) and lower > RHO_MARGIN:
        return "exceeds margin"
    return "not determined"


EXPECTED = []
OUTCOMES = {}


def expect(case, file, key, column, value):
    EXPECTED.append((case, file, key, column, fmt(value) if not isinstance(value, str) else value))


def predict_decision(case, spec, outcome=None, jack_finite=True, boot_nonfinite_pool=0):
    """Predicted decision.tsv, falsifier.tsv, positive-control.tsv and outcome.tsv rows of a degenerate (or declared) case."""
    keys = {(c, m) for c in CELLS for m in MASKS}
    s = agg(spec, keys)
    st = stats(s)
    items = [("rho_g", "rhog.B1", RHO_MARGIN), ("G1a", "B1.fdr_bh", FDP_LIMIT), ("G1b", "dfdrbh.B1", GUARD_MARGIN), ("G2", "dpowbh.B1", POWER_MARGIN)]
    states = {}
    for id_, stat, limit in items:
        est = st[stat]
        up = est  # degenerate: the bootstrap 0.975 quantile of identical resamples is the estimate
        states[id_] = state(up, limit)
        expect(case, "decision.tsv", id_, "estimate", est)
        expect(case, "decision.tsv", id_, "upper", up)
        expect(case, "decision.tsv", id_, "limit", limit)
        expect(case, "decision.tsv", id_, "state", states[id_])
        expect(case, "decision.tsv", id_, "upper_jack", up if jack_finite else NAN)
        expect(case, "decision.tsv", id_, "state_jack", state(up if jack_finite else NAN, limit))
        expect(case, "decision.tsv", id_, "boot_nonfinite", boot_nonfinite_pool if id_ == "rho_g" else 0)
    m6 = st["B1.m6_l_H"]
    fires = (not math.isnan(m6)) and (m6 > 0 or m6 < 0)
    expect(case, "falsifier.tsv", "F5'", "estimate", m6)
    expect(case, "falsifier.tsv", "F5'", "fires", fires)
    pc = (not math.isnan(st["Phi.S"])) and st["Phi.S"] > PC_MARGIN
    expect(case, "positive-control.tsv", "Phi.S", "estimate", st["Phi.S"])
    expect(case, "positive-control.tsv", "Phi.S", "passes", pc)
    if outcome is None:
        outcome = "Withheld" if fires else ("Calibrated" if all(v == "pass" for v in states.values()) else "Not calibrated")
    expect(case, "outcome.tsv", "", "outcome", outcome)
    expect(case, "outcome.tsv", "", "rho_g_state", states["rho_g"])
    expect(case, "outcome.tsv", "", "g1a_state", states["G1a"])
    expect(case, "outcome.tsv", "", "g1b_state", states["G1b"])
    expect(case, "outcome.tsv", "", "g2_state", states["G2"])
    expect(case, "outcome.tsv", "", "f5_fires", fires)
    expect(case, "outcome.tsv", "", "pc_pass", pc)
    expect(case, "outcome.tsv", "", "rho_g", st["rhog.B1"])
    expect(case, "outcome.tsv", "", "kind", "outcome")
    OUTCOMES[case] = outcome
    return st, outcome


CASES = []  # (case, rowsdir, env, expect)


def add_case(case, env="K281_RANGE=401,404;K211_B_BOOT=%d" % B_BOOT, expect_text="ok", rowsdir=None):
    CASES.append((case, rowsdir or case, env, expect_text))


def noop(d, cell, mask, r):
    pass


shutil.rmtree(os.path.join(HERE, "cases"), ignore_errors=True)

# ---- outcome cases: the deciding limit just below, at and just above 1.10 -----------------------------------------------------------------
def rho_case(case, n_add, p_add, split=None, extra=None, expect_outcome=None):
    def over(d, cell, mask, r):
        if split:
            d["gA1n"], d["gA1p"], d["gA0n"], d["gA0p"] = split
        else:
            d["gA1n"], d["gA1p"] = 0, 0
            d["gA0n"], d["gA0p"] = n_add, p_add
        if extra:
            extra(d)
    spec = build_spec(over)
    write_case(case, spec)
    add_case(case)
    st, outcome = predict_decision(case, spec, expect_outcome)
    return spec, st


rho_case("cal105", 2000, 105)                               # rho_g = 1.05 (0.0525 / 0.05), guards pass: Calibrated
rho_case("bnd_below", 100000, 5499)                         # 1.0997999999999999: pass
rho_case("bnd_at", 2000, 110)                               # 550 / 10,000 over 500 / 10,000 = 1.0999999999999999: pass (as written)
rho_case("bnd_above", 100000, 5501)                         # 1.1002: fail
rho_case("rho112", 10000, 560)                              # 1.12: Not calibrated
spec_split, st_split = rho_case("split", 0, 0, split=(1000, 90, 99000, 4950))  # flagged 90/1000 (1.8), unflagged 4950/99000 (1.0)
# the split as the scorer reports it (pooled; the per-design and per-setting split is in `fields`)
expect("split", "pooled.tsv", "pool|twelve|vd|rhog_l1.B1", "estimate", st_split["rhog_l1.B1"])
expect("split", "pooled.tsv", "pool|twelve|vd|rhog_l0.B1", "estimate", st_split["rhog_l0.B1"])
expect("split", "pooled.tsv", "pool|twelve|vd|rhog.B1", "estimate", st_split["rhog.B1"])
assert abs(st_split["rhog_l1.B1"] - 1.8) < 1e-12 and abs(st_split["rhog_l0.B1"] - 1.0) < 1e-12

# ---- guards -----------------------------------------------------------------------------------------------------------------------------------
def fdp_over(b_v, b_t, b1_v, b1_t, planted=600):
    def f(d):
        d["B_bhV"], d["B_bhT"], d["B1_bhV"], d["B1_bhT"], d["planted"] = b_v, b_t, b1_v, b1_t, planted
    return f


rho_case("g1a_fail", 2000, 100, extra=fdp_over(26, 474, 26, 474, 3000))      # FDP_BH(B1) = 26/500 = 0.052 > 0.05; paired difference 0
rho_case("g1a_pass_edge", 2000, 100, extra=fdp_over(25, 475, 25, 475, 3000))  # 25/500 = 0.05 exactly: pass (<= as written)
rho_case("g1b_fail", 2000, 100, extra=fdp_over(23, 477, 99, 1901, 3000))      # B 23/500 = 0.046, B1 99/2000 = 0.0495: difference 0.0035 > 0.0025
rho_case("g1b_pass", 2000, 100, extra=fdp_over(23, 477, 24, 476, 3000))       # difference 24/500 - 23/500 = 0.002: pass
rho_case("g2_fail", 2000, 100, extra=fdp_over(10, 300, 10, 294, 1000))        # power cost (300 - 294) / 1000 = 0.006 > 0.005
rho_case("g2_below", 2000, 100, extra=fdp_over(10, 5000, 10, 4951, 10000))    # 49 / 10,000 = 0.0049: pass
rho_case("g2_at", 2000, 100, extra=fdp_over(10, 200, 10, 195, 1000))          # 0.2 - 0.195 = 0.0050000000000000044 in doubles: fail as written
rho_case("g2_above", 2000, 100, extra=fdp_over(10, 5000, 10, 4949, 10000))    # 51 / 10,000: fail

# ---- positive control and F5' --------------------------------------------------------------------------------------------------------------
rho_case("pc_pass", 2000, 105, extra=lambda d: d["add"].__setitem__("S", (100, 15)))   # Phi_S = 3: the positive control passes; the outcome is the same as cal105
rho_case("pc_fail_ok", 2000, 105)                                                      # Phi_S = 1: the positive control fails; the outcome is unchanged (Calibrated)
rho_case("f5_fires_pos", 2000, 105, extra=lambda d: d.__setitem__("flag_err", 4))      # bias +0.2: F5' fires: Withheld (although every other row says Calibrated)
rho_case("f5_fires_neg", 2000, 105, extra=lambda d: d.__setitem__("flag_err", -4))     # bias -0.2: fires (upper < 0)
rho_case("f5_zero", 2000, 105, extra=lambda d: d.__setitem__("flag_err", 0))           # 0: does not fire: Calibrated
rho_case("f5_and_fail", 10000, 560, extra=lambda d: d.__setitem__("flag_err", 4))      # F5' fires and rho_g fails: Withheld (the first row that holds)

# ---- a pool with no added protein (B1 adds nothing): the deciding quantity is not finite; unlabelled settings ----------------------------------------
def none_over(d, cell, mask, r):
    d["g_added_rows"] = False
spec_none = build_spec(none_over)
write_case("noadded", spec_none)
add_case("noadded")
predict_decision("noadded", spec_none, "Not calibrated", boot_nonfinite_pool=B_BOOT)
for cell in CELLS:
    for mask in MASKS:
        expect("noadded", "labels-settings.tsv", "%s|%s" % (cell, mask), "label", "not determined")
        expect("noadded", "labels-settings.tsv", "%s|%s" % (cell, mask), "estimate", NAN)
        expect("noadded", "settings.tsv", "setting|%s|%s|rhog.B1" % (cell, mask), "boot_nonfinite", B_BOOT)
        expect("noadded", "settings.tsv", "setting|%s|%s|rhog.B1" % (cell, mask), "estimate", NAN)

# ---- non-finite resamples dropped and counted: only replicate 401 of A at vd5 carries added proteins ----------------------------------------------------
def partial_over(d, cell, mask, r):
    if (cell, mask, r) == ("A", "vd5", 401):
        d["gA0n"], d["gA0p"], d["gA1n"], d["gA1p"] = 100, 10, 0, 0     # rate 0.1 against 0.05: 2.0
    else:
        d["g_added_rows"] = False
spec_partial = build_spec(partial_over)
write_case("partial", spec_partial)
add_case("partial")
# pool: the finite resamples all hold the same ratio (n = 100 k, p = 10 k), so every limit is 2.0; the count of non-finite resamples is > 0
est_partial = div(div(10, 100), div(50 * 48, 1000 * 48))
expect("partial", "decision.tsv", "rho_g", "estimate", est_partial)
expect("partial", "decision.tsv", "rho_g", "upper", est_partial)
expect("partial", "decision.tsv", "rho_g", "state", "fail")
expect("partial", "decision.tsv", "rho_g", "boot_nonfinite", ">0")
expect("partial", "decision.tsv", "rho_g", "upper_jack", NAN)      # deleting replicate 401 of A leaves no added protein: the jackknife SE is not finite
expect("partial", "outcome.tsv", "", "outcome", "Not calibrated")
expect("partial", "settings.tsv", "setting|A|vd5|rhog.B1", "estimate", est_partial)
expect("partial", "settings.tsv", "setting|A|vd5|rhog.B1", "lo95", est_partial)
expect("partial", "settings.tsv", "setting|A|vd5|rhog.B1", "hi95", est_partial)
expect("partial", "settings.tsv", "setting|A|vd5|rhog.B1", "boot_nonfinite", ">0")
expect("partial", "labels-settings.tsv", "A|vd5", "label", "exceeds margin")
expect("partial", "labels-settings.tsv", "A|vd5", "label_jack", "not determined")
expect("partial", "labels-settings.tsv", "B|vd5", "label", "not determined")

# ---- every field in its place: distinct values in every setting (the defect of positional arguments landing in the wrong field) ----------------------------
def fields_over(d, cell, mask, r):
    i, j = CELLS.index(cell) + 1, MASKS.index(mask) + 1
    d["Bn"], d["Bp"] = 1000 + 37 * i + 5 * j, 50 + 3 * i + j
    d["gA0n"], d["gA0p"] = 1000 * i + 100 * j, 50 * i + 6 * j        # rates differ per setting: some within the margin, some above
    d["gA1n"], d["gA1p"] = 10 * i + j, 1 + (i + j) % 5
    d["B1_bhV"], d["B1_bhT"] = 7 * i + j, 100 * i + 10 * j + 3
    d["B_bhV"], d["B_bhT"] = 5 * i + 2 * j, 90 * i + 9 * j + 1
    d["planted"] = 600 + 100 * i + 10 * j
    d["flag_n"], d["flag_err"] = 11 * i + j, 0
spec_fields = build_spec(fields_over)
write_case("fields", spec_fields)
add_case("fields")
predict_decision("fields", spec_fields)
nrep = len(REPS)
for cell in CELLS:
    for mask in MASKS:
        s = agg({k: v for k, v in spec_fields.items()}, {(cell, mask)})
        st = stats(s)
        for stat in ("rhog.B1", "rhog_l1.B1", "rhog_l0.B1", "B1.fdr_bh", "dfdrbh.B1", "dpowbh.B1", "B1.gadd_per_r", "B.alpha_whole"):
            key = "setting|%s|%s|%s" % (cell, mask, stat)
            if stat in ("rhog.B1", "B1.gadd_per_r"):
                expect("fields", "settings.tsv", key, "R", nrep)
                expect("fields", "settings.tsv", key, "estimate", st[stat])
                expect("fields", "settings.tsv", key, "lo95", st[stat])
                expect("fields", "settings.tsv", key, "hi95", st[stat])
                expect("fields", "settings.tsv", key, "jack_se", 0.0)
                expect("fields", "settings.tsv", key, "boot_nonfinite", 0)
            else:
                expect("fields", "settings.tsv", key, "estimate", st[stat])
        expect("fields", "labels-settings.tsv", "%s|%s" % (cell, mask), "R", nrep)
        expect("fields", "labels-settings.tsv", "%s|%s" % (cell, mask), "estimate", st["rhog.B1"])
        expect("fields", "labels-settings.tsv", "%s|%s" % (cell, mask), "lower", st["rhog.B1"])
        expect("fields", "labels-settings.tsv", "%s|%s" % (cell, mask), "upper", st["rhog.B1"])
        expect("fields", "labels-settings.tsv", "%s|%s" % (cell, mask), "label", label(st["rhog.B1"], st["rhog.B1"]))
        expect("fields", "labels-settings.tsv", "%s|%s" % (cell, mask), "label_jack", label(st["rhog.B1"], st["rhog.B1"]))
    s = agg(spec_fields, {(cell, m) for m in MASKS})
    st = stats(s)
    for stat in ("rhog.B1", "rhog_l1.B1", "rhog_l0.B1", "B1.gadd_per_r", "B1.fdr_bh"):
        expect("fields", "designs.tsv", "design|%s|vd|%s" % (cell, stat), "estimate", st[stat])
        expect("fields", "designs.tsv", "design|%s|vd|%s" % (cell, stat), "R", nrep)
s = agg(spec_fields, {(c, m) for c in CELLS for m in MASKS})
st = stats(s)
for stat in ("rhog.B1", "rhog_l1.B1", "rhog_l0.B1", "B1.gadd_per_r", "B1.fdr_bh", "dfdrbh.B1", "dpowbh.B1"):
    expect("fields", "pooled.tsv", "pool|twelve|vd|%s" % stat, "estimate", st[stat])
expect("fields", "pooled.tsv", "pool|twelve|vd|rhog.B1", "R", nrep * len(CELLS))

# ---- instrument and falsifier Stops ------------------------------------------------------------------------------------------------------------------------
def base_spec():
    return build_spec(noop)


def stop_case(case, phrase, kind="instrument", env=None, **kw):
    write_case(case, base_spec(), **kw)
    add_case(case, env=env or "K281_RANGE=401,404;K211_B_BOOT=%d" % B_BOOT, expect_text="stop:%s" % phrase)
    expect(case, "outcome.tsv", "", "outcome", "Stop")
    expect(case, "outcome.tsv", "", "kind", kind)


stop_case("miss_rep", "a replicate is absent", drop={("C", m, 404) for m in MASKS})          # one replicate absent for one cell
stop_case("smoke", "smoke flag", provenance=dict(PROVENANCE_OK, smoke=True))
stop_case("smoke_absent", "smoke flag", provenance={"ticket": "KIT-281", "status": "complete"})
stop_case("range_out", "outside the declared range", env="K281_RANGE=401,403;K211_B_BOOT=%d" % B_BOOT)    # replicate 404 lies outside 401 to 403
stop_case("qmethod_bh", "q-method other than qvalue", over_metric=lambda m: m.__setitem__("qmethod", "BH") if (m["cell"], m["mask"], m["r"], m["arm"]) == ("B", "vd25", 403, "B1") else None)
stop_case("guardstop", "instrument-stop.txt", extra_files={"instrument-stop.txt": "partial-match guard: partial match of 'ab' to 'abcd'\n"})
stop_case("probe", "partial-match guard", env="K281_RANGE=401,404;K211_B_BOOT=%d;K281_GUARD_PROBE=1" % B_BOOT)   # a valid input; the scorer's own deliberate `$` prefix match
stop_case("f1", "F1", kind="falsifier", over_metric=lambda m: m.__setitem__("viol", 1) if (m["cell"], m["mask"], m["r"], m["arm"]) == ("D", "vd5", 402, "B1") else None)
stop_case("f2_e10", "F2", kind="falsifier", over_metric=lambda m: m.__setitem__("e10_dev", "1e-09") if (m["cell"], m["mask"], m["r"], m["arm"]) == ("A", "vd10", 401, "B") else None)
stop_case("f2_e3", "F2", kind="falsifier", over_metric=lambda m: m.__setitem__("e3_pass", "FALSE") if (m["cell"], m["mask"], m["r"], m["arm"]) == ("A", "vd10", 401, "B1") else None)
stop_case("incomplete_range", "range is incomplete", range_complete=3)
# a numeric column holding text, NA; labels T, F, NA in the cell, arm and mask columns
def row_cell_hook(label_):
    def hook(cell, rows, metrics):
        if cell == "B":
            rows[0]["cell"] = label_
    return hook
for lab in ("T", "F", "NA"):
    stop_case("cell_%s" % lab, "unexpected label", row_hook=row_cell_hook(lab))
def arm_hook(cell, rows, metrics):
    if cell == "C":
        rows[3]["arm"] = "NA"
stop_case("arm_NA", "unexpected label", row_hook=arm_hook)
def mask_hook(cell, rows, metrics):
    if cell == "D":
        rows[2]["mask"] = "mcar5"
stop_case("mask_mcar", "unexpected label", row_hook=mask_hook)
def text_hook(cell, rows, metrics):
    if cell == "A":
        rows[1]["n_tested"] = "abc"
stop_case("text_num", "not a number", row_hook=text_hook)
def na_hook(cell, rows, metrics):
    if cell == "A":
        rows[1]["n_p05"] = "NA"
stop_case("na_num", "not a number", row_hook=na_hook)
def logical_na_hook(cell, rows, metrics):
    if cell == "A":
        metrics[0]["mono_ok"] = "NA"
stop_case("logical_na", "not TRUE, FALSE, T or F", row_hook=logical_na_hook)
# T and F are accepted spellings of a logical: the scorer's result equals the base's (an UNLABELLED-column check as well: no value is lost to a logical guess)
def tf_hook(cell, rows, metrics):
    for m in metrics:
        m["mono_ok"] = "T" if m["mono_ok"] == "TRUE" else "F"
        m["e3_pass"] = "T"
spec_tf = rho_case("tf", 2000, 105)[0]
write_case("tf", spec_tf, row_hook=tf_hook)   # rewritten with T / F spellings; the prediction is cal105's

# ---- the join of two ranges and a duplicated replicate ------------------------------------------------------------------------------------------------------
spec_join = base_spec()
write_case("join", spec_join, jobs={"job1": [401, 402], "job2": [403, 404]})
add_case("join")
predict_decision("join", spec_join)
write_case("dup", spec_join, jobs={"job1": [401, 402, 403], "job2": [403, 404]})
add_case("dup", expect_text="stop:duplicated")
expect("dup", "outcome.tsv", "", "outcome", "Stop")
expect("dup", "outcome.tsv", "", "kind", "instrument")

# ---- the F2 determinism failure and the counted mode ------------------------------------------------------------------------------------------------------------
spec_det = base_spec()
write_case("f2_det", spec_det, extra_files={"determinism.tsv": "cell\tr\trows_identical\tmetrics_identical\nA\t401\tTRUE\tFALSE\n"})
add_case("f2_det", expect_text="stop:F2")
expect("f2_det", "outcome.tsv", "", "outcome", "Stop")
expect("f2_det", "outcome.tsv", "", "kind", "falsifier")
G3_OK = "cell\tr_cell\tfirst_r\tlast_r\tn_hat\tbranch\n" + "".join("%s\t4\t401\t404\t276300\tselftest\n" % c for c in CELLS)
G3_FIVE = "cell\tr_cell\tfirst_r\tlast_r\tn_hat\tbranch\n" + "".join("%s\t5\t401\t405\t276300\tselftest\n" % c for c in CELLS)
os.makedirs(os.path.join(HERE, "cases", "g3tables"), exist_ok=True)
open(os.path.join(HERE, "cases", "g3tables", "g3-ok.tsv"), "w").write(G3_OK)
open(os.path.join(HERE, "cases", "g3tables", "g3-five.tsv"), "w").write(G3_FIVE)
# counted mode (the default B = 10,000, the range read from the G3 table): determinism present and absent, a table whose R differs from the rows', no table
spec_c = base_spec()
write_case("counted_ok", spec_c, determinism=True)
add_case("counted_ok", env="K281_G3_TABLE=@HERE@/cases/g3tables/g3-ok.tsv")
predict_decision("counted_ok", spec_c)
expect("counted_ok", "outcome.tsv", "", "counted", "TRUE")
expect("counted_ok", "outcome.tsv", "", "range_source", "g3-table")
expect("counted_ok", "outcome.tsv", "", "bootstrap_B", 10000)
write_case("counted_nodet", spec_c)
add_case("counted_nodet", env="K281_G3_TABLE=@HERE@/cases/g3tables/g3-ok.tsv", expect_text="stop:no determinism table")
write_case("counted_g3five", spec_c, determinism=True)
add_case("counted_g3five", env="K281_G3_TABLE=@HERE@/cases/g3tables/g3-five.tsv", expect_text="stop:a replicate is absent")
write_case("counted_g3none", spec_c, determinism=True)
add_case("counted_g3none", env="K281_G3_TABLE=@HERE@/cases/g3tables/absent.tsv", expect_text="stop:no declared replicate range")
for case in ("counted_nodet", "counted_g3five", "counted_g3none"):
    expect(case, "outcome.tsv", "", "outcome", "Stop")
    expect(case, "outcome.tsv", "", "kind", "instrument")
# checks.tsv of one ordinary case, as written
expect("cal105", "checks.tsv", "complete", "pass", True)
expect("cal105", "checks.tsv", "qmethod", "pass", True)
expect("cal105", "checks.tsv", "range_declared", "pass", False)     # K281_RANGE is a mechanics run
expect("cal105", "checks.tsv", "bootstrap_B", "pass", False)        # B = 1000
expect("cal105", "replicates.tsv", "A", "R", 4)
expect("cal105", "outcome.tsv", "", "counted", "FALSE")
# (the T / F case "tf" carries cal105's predictions: rho_case("tf", ...) predicted them before its files were rewritten with T and F spellings)

# ---- the outcomes the ticket's rules give, written by hand (a coding mistake in the predictions above stops here) --------------------------------------------
HAND = {"cal105": "Calibrated", "bnd_below": "Calibrated", "bnd_at": "Calibrated", "bnd_above": "Not calibrated", "rho112": "Not calibrated",
        "split": "Calibrated", "g1a_fail": "Not calibrated", "g1a_pass_edge": "Calibrated", "g1b_fail": "Not calibrated", "g1b_pass": "Calibrated",
        "g2_fail": "Not calibrated", "g2_below": "Calibrated", "g2_at": "Not calibrated", "g2_above": "Not calibrated", "pc_pass": "Calibrated",
        "pc_fail_ok": "Calibrated", "f5_fires_pos": "Withheld", "f5_fires_neg": "Withheld", "f5_zero": "Calibrated", "f5_and_fail": "Withheld",
        "noadded": "Not calibrated", "tf": "Calibrated", "join": "Calibrated", "counted_ok": "Calibrated"}
for case, want in HAND.items():
    assert OUTCOMES[case] == want, (case, OUTCOMES[case], want)

# ---- the G3 function's table: Fractions, not floats --------------------------------------------------------------------------------------------------------
def g3(n_hat):
    n_hat = Fraction(n_hat)
    if n_hat >= 276220:
        return 200
    if n_hat == 0:
        return 571
    q = Fraction(200 * 276220) / n_hat
    c = -((-q.numerator) // q.denominator)       # ceiling
    return int(min(571, max(200, c)))


g3_values = [0, 1, 96786, 100000, 138110, 199444, 199445, 200000, 276099, 276100, 276200, 276219, 276220, 276221, 276300, 300000, 1000000, 96786.5,
             96785, 138109, 138110.5, 55244000 // 571, 55244000 // 571 + 1, 55244000 // 571 - 1]
with open(os.path.join(HERE, "g3-cases.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["n_hat", "r"])
    for v in g3_values:
        w.writerow([repr(v) if isinstance(v, float) else v, g3(v)])
assert g3(96786) == 571 and g3(200000) == 277 and g3(276220) == 200 and g3(276219) == 201 and g3(0) == 571

# ---- the G3 request checks (check-request-k281.py): a request whose R differs from the rule's is refused --------------------------------------------------
os.makedirs(os.path.join(HERE, "requests"), exist_ok=True)


def cite_of():
    import hashlib
    return hashlib.sha256(open(os.path.join(HERE, "requests", "g3-277.tsv"), "rb").read()).hexdigest()[:16]


def request(jobs, label="selftest", cite=True):
    scripts = []
    for cells, first, last in jobs:
        scripts.append({"path": "docs/qa/workflow-execution/kit/kit-281/sim-k281.R", "args": ["k281-out", str(first), str(last)],
                        "env": {"K211_CELLS": cells}, "out": "k281-out", "out_pattern": "\\.(tsv|gz|json|txt|md)$", "timeout_minutes": 75})
    return {"requested_by": "selftest; G3 table sha256 " + (cite_of() if cite else "none"), "label": label, "package_commit": "0" * 40, "base_commit": "none", "tests": [], "scripts": scripts}


G3_R = "cell\tr_cell\tfirst_r\tlast_r\tn_hat\tbranch\n" + "".join("%s\t277\t401\t677\t200000\tselftest\n" % c for c in CELLS)
open(os.path.join(HERE, "requests", "g3-277.tsv"), "w").write(G3_R)
req_cases = []
def req(name, jobs, expect_code, cite=True):
    json.dump(request(jobs, cite=cite), open(os.path.join(HERE, "requests", name + ".json"), "w"), indent=2)
    req_cases.append((name, expect_code))
# R = 277 for every cell: two jobs per cell (200 + 77): accepted
req("req_ok", [(c, f, l) for c in CELLS for f, l in ((401, 600), (601, 677))], 0)
req("req_r276", [(c, f, l) for c in CELLS for f, l in ((401, 600), (601, 676))], 1)            # R = 276 for every cell: differs from the rule's 277
req("req_cellshort", [(c, f, l) for c in CELLS[:3] for f, l in ((401, 600), (601, 677))] + [("D", 401, 600)], 1)   # one cell stops at 600
req("req_overlap", [(c, f, l) for c in CELLS for f, l in ((401, 600), (600, 677))], 1)         # replicate 600 run twice
req("req_toolong", [(c, f, l) for c in CELLS for f, l in ((401, 677),)], 1)                    # one job of 277 replicates (above 200)
req("req_nocite", [(c, f, l) for c in CELLS for f, l in ((401, 600), (601, 677))], 1, cite=False)   # a request that does not cite the table
req("req_range", [(c, f, l) for c in CELLS for f, l in ((402, 601), (602, 678))], 1)           # not 401 to 677
with open(os.path.join(HERE, "requests.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["request", "g3_table", "exit"])
    for name, code in req_cases:
        w.writerow([name + ".json", "g3-277.tsv", code])

# ---- cases.tsv and expected.tsv -------------------------------------------------------------------------------------------------------------------------------
with open(os.path.join(HERE, "cases.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["case", "rowsdir", "env", "expect"])
    for c in CASES:
        w.writerow(c)
with open(os.path.join(HERE, "expected.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["case", "file", "key", "column", "value"])
    for e in EXPECTED:
        w.writerow(e)
print("make-selftest: %d cases, %d predicted values, %d G3 values, %d request cases" % (len(CASES), len(EXPECTED), len(g3_values), len(req_cases)))
