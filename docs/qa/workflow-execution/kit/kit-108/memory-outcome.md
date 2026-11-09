# KIT-108: the report-step memory after the fix (M4's GUI members, 2026-11-09)

**Outcome: D8's threshold is not met. The fix is kept, and P-30 stands.**

On the two large corpora, route A brings the application's report-step peak down by 0.52 GiB (metab) and 0.27 GiB (lipid). Its level after the report comes down by 0.60 and 0.36 GiB. The report step still rises with data size:
- **The peak** misses D8's limit by 0.68 GiB (metab) and 0.79 GiB (lipid).
- **The teardown level** misses it by 0.57 and 0.62 GiB.

Both misses are as the diagnosis predicted:
- one `openxlsx::read.xlsx` of a 45 MB DA workbook still needs about 1 GB of transient heap;
- the combined workbook still holds about 0.41 GB per sheet until it is saved;
- glibc still keeps the freed pages after the copy.

Route A cannot remove either cause; only P-30's routes B and C can ([report-memory-diagnosis.md](memory/report-memory-diagnosis.md) §6, [m4-replica.md](memory/m4-replica.md)). Both post-fix members are valid, exact and accepted, with their baselines' dispatches. No worker or slice peak is worse by more than 0.1 GiB.

## Runs

| Member | Corpus | Build | Started | Decision | Exact |
|---|---|---|---|---|---|
| `108p` (baseline, M2) | `metab_large_missing` | `kit_107` (`66sPmE`) | 2026-11-09 17:11 | accepted | yes |
| `108q` (baseline, M2) | `lipid_large_missing` | `kit_107` | 17:34 | accepted | yes |
| **`108v`** (post-fix, M4) | `metab_large_missing` | **`kit_108` (`X6BzNV`, `7b52f535`)** | 22:31 | accepted | yes |
| **`108w`** (post-fix, M4) | `lipid_large_missing` | `kit_108` | 22:52 | accepted | yes |

- **How they ran:** each once, in sequence, through `member.sh` (wait 7,200 s), with I2 attached (`sampler.sh`). No other `multischolar*` unit was live. AC power, low-power profile unchanged, 29 GiB available at start ([power-state.txt](power-state.txt)).
- **Validity:** for each member, 0 sample gaps above 250 ms, swap 0, phases complete, verification complete.
- **I2 (the R process):** its largest tick gap inside the report step was 1.13 s (`108v`). All 1,065 of `108v`'s rows are within 10% of the harness's controller anon (max relative error 0.030), and the controller ran as one R PID. The stage summaries are `kit108-kit-journey-108{v,w}-stage.json`.
- **Raw outputs:** in `/home/doktersmol/.cache/k108-mem` (`108{v,w}-i1*`, `108{v,w}-i2*`).
- **Dispatches:** `108v` has the same dispatches as `108p`, and `108w` the same as `108q`.

## The report step (I1, the harness's controller trace split by the driver's marks; GiB)

| | `108p` → **`108v`** (metab) | `108q` → **`108w`** (lipid) |
|---|---|---|
| DA-step peak | 1.540 → 1.527 | 1.398 → 1.499 |
| Level in the 5 s before the report | 1.308 → 1.527 | 1.398 → 1.435 |
| **Report-step peak** | 2.980 → **2.458** (−0.522) | 2.815 → **2.543** (−0.272) |
| **Anon at teardown** | 2.945 → **2.342** (−0.603) | 2.669 → **2.308** (−0.361) |
| **Slice `resident_end_bytes`** | 3.304 → **2.702** (−0.602) | 3.029 → **2.666** (−0.363) |
| Copy time (s) | 49.6 → 47.1 | 43.0 → 48.6 |

**The copy window (I2, the R process's RssAnon, entry → peak):**
- **metab:** 1.115 → 2.879 (+1.76) before, 1.407 → 2.458 (+1.05) after;
- **lipid:** 1.423 → 2.788 (+1.37) before, 1.251 → 2.493 (+1.24) after.

The render adds nothing after the fix (−0.0 to +0.002 GiB); before, it added +0.065 and +0.078. The level is flat from the render's end to teardown.

## D8's threshold, per corpus

| Criterion | `108v` (metab) | `108w` (lipid) |
|---|---|---|
| Report-step peak ≤ DA-step peak + 0.25 | 2.458 against 1.777: **not met** (+0.681) | 2.543 against 1.749: **not met** (+0.794) |
| Teardown ≤ the pre-report level + 0.25 | 2.342 against 1.777: **not met** (+0.565) | 2.308 against 1.685: **not met** (+0.623) |

**Why not, from the attribution:**
- **The peak:** route A's floor is one sheet's read over the earlier sheets' cell store, about 1.54 GB of R heap ([m4-replica.md](memory/m4-replica.md), "the floor of this mechanism"). Removing it means changing the reader (route D, a new dependency) or moving the copy out of the application (route C).
- **The teardown:** the copy's freed heap stays with the allocator. In the replica, `malloc_trim(0)` returned 94% of it (route B).

Both routes are P-30 and need the user's decision.

## No regression (D8's second rule; the slice and the worker steps)

The slice's anon peak, against the baseline:
- **metab:** 4.130 → 4.043 (−0.087 GiB);
- **lipid:** 4.085 → 4.104 (+0.019 GiB).

Per step, against the baseline:
- **The slice's largest per-step rise** is +0.052 GiB (`108w`, DA). Every other step is within ±0.06 GiB, or lower.
- **The controller's DA step on `108w`** is +0.101 GiB (1.398 → 1.499); on `108v` it is −0.013. KIT-108 changes nothing the DA step runs, and the DA step is not one of the rule's peaks (worker and slice). It is reported here as observed, with single runs on each side.

**No regression.**

## The comparison's limits

- **`108p`/`108q` → `108v`/`108w` spans two builds.** `kit_107` → `kit_108` is the whole of KIT-108: the shared summary frame and core (W3, W4), F104-01, and route A. Route A's own effect is isolated by the replica at the W4 and W5 trees ([m4-replica.md](memory/m4-replica.md)): the copy's peak rise fell 0.32 GiB on `108p`'s project. The GUI's fall on metab (−0.52 GiB) is larger than that. The members alone do not separate the frame's share from route A's, and single runs do not give the spread.
- **`108p` ran while the W3 host-rule breach happened** ([findings.md](findings.md)). The figures here are the application's own process, so the comparison holds; its timings had company.
- **`108v` entered the report 0.22 GiB higher than `108p`** (1.527 against 1.308), and the copy started from 1.407 against 1.115. A process that already holds free pages can reuse them, so part of the lower peak may come from that. It was not measured separately.

## Content

- **Both post-fix members are exact** against the same oracles as their baselines, and their dispatches are identical. So their states and products are equal.
- **The replica's content checks** on the three projects are identical at W4 and W5 ([m4-replica.md](memory/m4-replica.md), "D8's content checks"): every `results_summary` file, every workbook sheet, and the report text.
