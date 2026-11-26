"""KIT-262 M3 E0 (local, standard library only): the copied blocks of m3-k262.R, score-m3-k262.R and reestimate-m3-k262.R equal their
sources except the marked lines; the stress generator's floor order is `after`; the row format is the same in the simulation and the
scorer. Run from the repository root:
    python3 -I docs/qa/workflow-execution/kit/kit-262/e0-m3-k262.py
Sources: KIT-190's effect.R, decide.R, bounds.R, reestimate.R; tests/testthat/helper-kit172.R; KIT-262's m2-k262.R (all as on the branch,
which holds the binding's tree). Rule: every line of a copied function in the destination is either a line of the source function (same
text after stripping, in the source's order) or carries a mark `# [X]`; every source line the destination lacks is counted against the
marked lines (each marked line may replace at most one source line). Prints one line per check; exits 1 on a failure."""
import re
import sys
from pathlib import Path

root = Path(".")
here = root / "docs/qa/workflow-execution/kit/kit-262"
k190 = root / "docs/qa/workflow-execution/kit/kit-190"


def lines(path):
    return Path(path).read_text().split("\n")


def span(text, header):
    start = next(i for i, l in enumerate(text) if l.startswith(header))
    end = start
    if text[start].rstrip().endswith("}") and text[start].count("{") == text[start].count("}"):
        return text[start : start + 1]
    while text[end] != "}":
        end += 1
    return text[start : end + 1]


def norm(seq):
    return [l.strip() for l in seq if l.strip()]


failures = 0


def check(name, ok, detail=""):
    global failures
    print(("PASS " if ok else "FAIL ") + name + (" " + detail if detail else ""))
    if not ok:
        failures += 1


def compare(name, src, dst, rename_ok=False):
    """dst lines: unmarked ones must be source lines in order; marked ones are allowed; missing source lines <= marked lines."""
    s, d = norm(src), norm(dst)
    marked = [l for l in d if re.search(r"#\s*\[[A-Z]\]", l)]
    unmarked = [l for l in d if not re.search(r"#\s*\[[A-Z]\]", l)]
    cursor, bad = 0, []
    for l in unmarked:
        try:
            cursor = s.index(l, cursor) + 1
        except ValueError:
            bad.append(l)
    kept = len(unmarked) - len(bad)
    missing = len(s) - kept
    ok = not bad and missing <= len(marked)
    check(name, ok, "%d source lines, %d kept verbatim, %d marked, %d source lines replaced%s" % (len(s), kept, len(marked), missing,
          "; NOT IN SOURCE: " + " | ".join(bad[:3]) if bad else ""))


m3, sc, re_ = lines(here / "m3-k262.R"), lines(here / "score-m3-k262.R"), lines(here / "reestimate-m3-k262.R")
m2 = lines(here / "m2-k262.R")
eff = lines(k190 / "effect.R")
dec, bnd, rst = lines(k190 / "decide.R"), lines(k190 / "bounds.R"), lines(k190 / "reestimate.R")
hk = lines(root / "tests/testthat/helper-kit172.R")

# ---- m3-k262.R against KIT-190's effect.R ----
compare("m3 e2PickCells = effect.R", span(eff, "e2PickCells <- function"), span(m3, "e2PickCells <- function"))
compare("m3 e2Values = effect.R (marked [L])", span(eff, "e2Values <- function"), span(m3, "e2Values <- function"))
compare("m3 e2Data = effect.R (marked [L])", span(eff, "e2Data <- function"), span(m3, "e2Data <- function"))
check("m3 E2_CELLS = effect.R", span(eff, "E2_CELLS <- ")[0].strip() == span(m3, "E2_CELLS <- ")[0].strip())
check("m3 ALIGNED_CELLS = effect.R", span(eff, "ALIGNED_CELLS <- ")[0].strip() == span(m3, "ALIGNED_CELLS <- ")[0].strip())
check("m3 CUTOFF = effect.R", span(eff, "CUTOFF <- ")[0].strip() == span(m3, "CUTOFF <- ")[0].strip())
compare("m3 k190Oracle = effect.R", span(eff, "k190Oracle <- function"), span(m3, "k190Oracle <- function"))
compare("m3 k190Candidate = effect.R", span(eff, "k190Candidate <- function"), span(m3, "k190Candidate <- function"))
# the seed rule of e2Data in the source is the one the ticket names
check("e2Data seeds 190000000 + 10000 x idx + r", any("seed <- 190000000 + 10000 * idx + r" in l for l in span(m3, "e2Data <- function")))

# ---- m3-k262.R against M2's script and the harness ----
compare("m3 m2Design = m2-k262.R", span(m2, "m2Design <- function"), span(m3, "m2Design <- function"))
compare("m3 m2Object = m2-k262.R", span(m2, "m2Object <- function"), span(m3, "m2Object <- function"))
check("m3 M2_RUV_PARAMS = m2-k262.R", norm(span(m2, "M2_RUV_PARAMS <- ")) == norm(span(m3, "M2_RUV_PARAMS <- ")) if False else
      norm(m2[next(i for i, l in enumerate(m2) if l.startswith("M2_RUV_PARAMS <- ")) :][:2]) == norm(m3[next(i for i, l in enumerate(m3) if l.startswith("M2_RUV_PARAMS <- ")) :][:2]))
compare("m3 m2Prepare = m2-k262.R", span(m2, "m2Prepare <- function"), span(m3, "m2Prepare <- function"))
compare("m3 m2Sweep = m2-k262.R", span(m2, "m2Sweep <- function"), span(m3, "m2Sweep <- function"))
compare("m3 m3Object = kit172Object (marked [D])", span(hk, "kit172Object <- function"), span(m3, "m3Object <- function"))
compare("m3 m3StateMatrix = kit172StateMatrix (marked [D])", span(hk, "kit172StateMatrix <- function"), span(m3, "m3StateMatrix <- function"))

# ---- the stress generator: its lines are lines of m2Values or e2Values, or marked; the floor order is `after` ----
def code(l):
    return re.sub(r"\s+#.*$", "", l.strip())


pool = set(code(l) for l in norm(span(m2, "m2Values <- function")) + norm(span(eff, "e2Values <- function")))
gen = span(m3, "m3StressValues <- function")
unmarked_new = [l.strip() for l in gen if l.strip() and code(l) not in pool and not re.search(r"#\s*\[[A-Z]\]", l)]
check("m3StressValues: every line not in m2Values/e2Values is marked", not unmarked_new, "; ".join(unmarked_new[:3]))
text = [l.strip() for l in gen]
i_mask = next(i for i, l in enumerate(text) if l.startswith("values[sample.int(length(values), floor(rate * length(values)))] <- NA_real_"))
i_pick = next(i for i, l in enumerate(text) if l.startswith("ones <- e2PickCells(!is.na(values), groups, 2L, round(0.02 * sum(!is.na(values))))"))
i_one = next(i for i, l in enumerate(text) if l.startswith("values[ones] <- 1"))
check("m3StressValues: the mask comes first, then the floor values on observed cells (KIT-172's order)", i_mask < i_pick < i_one)
check("m3StressValues: floor order named `after` and asserted", any('stopifnot(identical(floor_order, "after"))' in l for l in gen) and
      any(l.startswith('M3_STRESS_FLOOR_ORDER <- "after"') for l in m3))
check("m3StressValues: default floor_order is the named constant", any("floor_order = M3_STRESS_FLOOR_ORDER" in l for l in gen[:2]))
check("stress seeds 262200000 + 10000 x idx + r", any("262200000 + 10000 * idx + r" in l for l in m3))
# the stream of draws up to the mask is m2Values'
seq_m2 = [l for l in norm(span(m2, "m2Values <- function")) if l.startswith(("level <-", "logs <-", "effects <-", "loadings <-", "changed <-", "w <-", "pattern <-"))]
seq_m3 = [l for l in norm(gen) if l.startswith(("level <-", "logs <-", "effects <-", "loadings <-", "changed <-", "w <-", "pattern <-"))]
check("m3StressValues: the draws before the mask are M2's (the effect-share line aside)", [re.sub(r"\s*#.*$", "", l) for l in seq_m3 if not l.startswith("changed <-")] ==
      [re.sub(r"\s*#.*$", "", l) for l in seq_m2 if not l.startswith("changed <-")])

# ---- score-m3-k262.R against decide.R and bounds.R ----
for name in ("tangoZ", "bergerBoos", "skewness", "effectTest", "nullTest"):
    compare("score %s = decide.R" % name, span(dec, name + " <- function"), span(sc, name + " <- function"))
src_b = span(bnd, "tangoZ <- function")
dst_b = span(sc, "tangoZb <- function")
check("score tangoZb = bounds.R tangoZ (renamed)", norm(src_b)[1:] == norm(dst_b)[1:] and norm(dst_b)[0].startswith("tangoZb <- function(b, c, n, m)"))
check("score alpha, margins and the sign-flip seed", any(l.startswith("ALPHA <- 0.05") for l in sc) and any(l.startswith("M_NULL <- 0.10") for l in sc) and
      any(l.startswith("M_EFFECT <- 0.05") for l in sc) and any(l.startswith("SIGNFLIP_SEED <- 262000002L") for l in sc) and any(l.startswith("BOOT_SEED <- 262000001L") for l in sc) and
      any(l.startswith("BOOT_B <- 10000L") for l in sc))
check("decide.R's margins and alpha are the same", any(l.startswith("ALPHA <- 0.05") for l in dec) and any(l.startswith("M_NULL <- 0.10") for l in dec) and any(l.startswith("M_EFFECT <- 0.05") for l in dec))

# ---- reestimate-m3-k262.R against reestimate.R ----
compare("reestimate estimate() = reestimate.R (marked [R])", span(rst, "estimate <- function"), span(re_, "estimate <- function"))
compare("reestimate writeEstimate() = reestimate.R", span(rst, "writeEstimate <- function"), span(re_, "writeEstimate <- function"))
check("reestimate SCHEMA = reestimate.R", span(rst, "SCHEMA <- ")[0].strip() == span(re_, "SCHEMA <- ")[0].strip())
check("reestimate has no ceiling loop (R is never reduced)", not any("t - 0.0005" in l or "while (core" in l for l in re_))

check("m3Tag and m3Seed are the same function in the simulation and the scorer", span(m3, "m3Tag <- function")[0].split("<-", 1)[1].strip() == span(sc, "m3Tag <- function")[0].split("<-", 1)[1].split("#")[0].strip() and
      span(m3, "m3Seed <- function")[0].split("<-", 1)[1].strip() == span(sc, "m3Seed <- function")[0].split("<-", 1)[1].strip())

# ---- the row format is one thing ----
def columns(text):
    start = next(i for i, l in enumerate(text) if l.startswith("M3_COLUMNS <- c("))
    end = start
    while not text[end].rstrip().endswith(")"):
        end += 1
    return re.findall(r'"([^"]+)"', " ".join(text[start : end + 1]))


c3, cs = columns(m3), columns(sc)
check("M3_COLUMNS in m3-k262.R = score-m3-k262.R (%d columns)" % len(c3), c3 == cs and len(c3) == len(set(c3)))

# ---- the self-test's expected values all exist in the Python-worked table; the environment variables read are the ones documented ----
expected_keys = set(l.split("\t")[0] for l in lines(here / "e5hand-m3-expected.tsv")[1:] if l)
used_keys = set(re.findall(r'ev\("([^"]+)"\)', "\n".join(sc)))
check("score self-test: every ev() key is in e5hand-m3-expected.tsv (%d keys used, %d available)" % (len(used_keys), len(expected_keys)), used_keys <= expected_keys,
      "missing: " + ", ".join(sorted(used_keys - expected_keys)))
header = "\n".join(l for l in m3 if l.startswith("#"))
documented = set(re.findall(r"K262_[A-Z0-9_]+", header))
read = set(re.findall(r'Sys\.getenv\("(K262_[A-Z0-9_]+)"', "\n".join(m3)))
check("m3-k262.R: every environment variable it reads is documented in its header", read <= documented, "undocumented: " + ", ".join(sorted(read - documented)))
check("m3-k262.R: every documented environment variable is read", documented <= read,
      "documented, not read: " + ", ".join(sorted(documented - read)))
# no global handler that stops: the partial-match guard is installed only through withCallingHandlers
check("no globalCallingHandlers in the M3 scripts", not any("globalCallingHandlers(" in l and not l.lstrip().startswith("#") for l in m3 + sc + re_))
check("every stop of the guard is reached through a local calling handler (withCallingHandlers present around journeys and units)", sum("withCallingHandlers(" in l for l in m3) >= 3)

sys.exit(1 if failures else 0)
