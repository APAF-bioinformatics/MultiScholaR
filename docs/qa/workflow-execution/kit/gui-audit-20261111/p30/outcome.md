# P-30 route B: the report step on the two large corpora (2026-11-11)

**Outcome: route B does what it was approved for. After the report, the application's level falls by 0.92 GiB (metab) and 0.86 GiB (lipid), to within 0.11 and 0.06 GiB of its level before the report. The report-step peak does not move (−0.07 and −0.11 GiB, which counts as no effect for single runs), as the diagnosis predicted.**

## Runs

| Member | Corpus | Build | Decision | Exact |
|---|---|---|---|---|
| `108v` (KIT-108, route A only) | `metab_large_missing` | `kit_108` (`X6BzNV`, `7b52f535`) | accepted | yes |
| `108w` (KIT-108, route A only) | `lipid_large_missing` | `kit_108` | accepted | yes |
| **`920v`** (route B) | `metab_large_missing` | **`kit_920`** (`kH2HPy`, installed at `afcb3854`) | accepted | yes |
| **`920w`** (route B) | `lipid_large_missing` | `kit_920` | accepted | yes |

- **How they ran:** each once, in sequence, through `kit/scripts-audit-20261111/p30/member.sh`. No other `multischolar*` unit was live. AC power, low-power profile, load about 2, 41 GiB available ([power-state.txt](power-state.txt)). Stage summaries: `audit-kit-journey-920{v,w}-stage.json`.
- **Instrument:** I1 only (KIT-108's `analyse-traces.R`, the harness's controller trace split by the driver's marks). The same script on `108v`/`108w` gives KIT-108's published figures to the third decimal. The tables are [p30-i1.tsv](p30-i1.tsv) and [p30-i1-summary.tsv](p30-i1-summary.tsv). I2 (the per-PID sampler) was not attached.
- **Declared before the run:** one attempt each; compare against `108v`/`108w`; with single runs, a difference under about 0.1 GiB counts as no effect.

## The report step (GiB, controller anon)

| | `108v` → **`920v`** (metab) | `108w` → **`920w`** (lipid) |
|---|---|---|
| Level in the 5 s before the report | 1.527 → 1.310 | 1.435 → 1.391 |
| **Report-step peak** | 2.458 → **2.385** (−0.073) | 2.543 → **2.437** (−0.106) |
| Level after the copy | 2.342 → 1.236 | 2.309 → 1.269 |
| Level after the render | 2.342 → 1.415 | 2.308 → 1.451 |
| **Anon at teardown** | 2.342 → **1.418** (−0.924) | 2.308 → **1.453** (−0.855) |
| Slice `resident_end_bytes` | 2.702 → 1.775 (−0.927) | 2.666 → 1.810 (−0.856) |

After the trim, the render re-grows the heap by about 0.3 GiB, because it now starts from returned pages instead of retained ones. The level is then flat to teardown.

## D8's two criteria, per corpus

| Criterion | `920v` (metab) | `920w` (lipid) |
|---|---|---|
| Teardown ≤ the pre-report level + 0.25 | 1.418 against 1.560: **met** | 1.453 against 1.641: **met** |
| Report-step peak ≤ DA-step peak + 0.25 | 2.385 against 1.718: not met (+0.667) | 2.437 against 1.752: not met (+0.685) |

The user retired the peak criterion with route B (parked-decisions answers, P-30): the report step is not the workflow's peak, which is the slice's (about 4.0 to 4.1 GiB at KIT-108). Removing it would need route C (the copy on the worker) or route D (another reader).

## The comparison's limits

- **The two builds differ by more than route B.** `kit_108` → `kit_920` includes KIT-109 to KIT-113, the remediation, and the interactive volcano page at the DA step for every omic. Route B is the only change that touches the copy's exit. The teardown fall (0.86 to 0.92 GiB) matches the replica's prediction (`malloc_trim(0)` returned about 94% of the retained memory; [m4-replica.md](../../kit-108/memory/m4-replica.md)). The peak change is within single-run noise, and the report started 0.04 to 0.22 GiB lower.
- **Wall time.** Both members took about 38 minutes against about 21 at KIT-108. Every step is about 1.7 times slower, at a similar load: for example the DA step on lipid took 82.9 → 142.2 s and the save took 7.4 → 12.8 s. The copy took 47.1 → 84.8 s (metab) and 48.6 → 81.8 s (lipid), which is the same ratio. So the slowdown does not come from the trim. These are not timing evidence, and the trim's own duration was not isolated, because I2 was not attached.
- **Single runs on each side.** There is no spread.

**Decision:** P-30 is closed with route B. The level after the report meets D8 on both corpora. The peak stays as measured, and the user retired that criterion.
