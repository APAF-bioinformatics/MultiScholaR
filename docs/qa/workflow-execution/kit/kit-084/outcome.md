# KIT-084 outcome: how the stage-3b numbers scale, and the close of stage 3b

**Synthetic corpora only: this supports no claim about real data. No default changed.**

**At 200 samples by 20,000 features per assay, both omics complete exactly on `kit_084`, and the application stays at about 1 GiB through the worker steps. Its largest event is now the report step, at 3.33 GiB on both omics, and it stays near that afterwards. The whole slice peaks at 3.79 to 3.97 GiB in normalisation. There is no old-versus-new result at this size: the `kit_032` controls were dropped at the user's direction.**

Stage 3b closes here (user decision, 2026-11-05): real data moves to a final gate after integration (KIT-063).

## Members

| Member | Build | Corpus | Result | Time | Exactness |
|---|---|---|---|---|---|
| Lipid record | `kit_082` | lipid large | valid | 4,134 s | the `lipid84` oracle |
| `kit-journey-084s` (small test) | `kit_032` | metab screen | accepted | 360 s | exact |
| `kit-journey-084a` (member 2) | `kit_082` | metab large | **invalid** | stopped at a worker deadline | none |
| `kit-journey-084c` (member 2 again) | `kit_084` | metab large | **accepted** | 1,235 s | exact against `metab84`, 199 files |
| `kit-journey-084b` (control) | `kit_032` | metab large | **invalid** | 8,831 s, machine suspended | none |
| `kit-journey-084f` (control again) | `kit_032` | metab large | **stopped by the user** | 8,661 s, still in normalisation | none |
| `kit-journey-084d` | `kit_084` | lipid large | **accepted** | 2,136 s | exact against `lipid84`, 208 files |
| `kit-journey-084e` (lipid control) | `kit_032` | lipid large | **not run** (user) | | |
| `kit-084p`, `kit-084q`, `kit-084r` | `kit_084` | both large | valid, table below | 3 to 8 min per stage | all 12 members exact |

Every worker dispatch on `084c` and `084d` completed. Neither member used swap.

## From screen to large

The screen assays are about 300 rows by 24 samples, and the large ones 20,000 rows by 200 samples, so about 555 times the cells. Screen members ran on `kit_082`, large members on `kit_084`; the two builds differ only in the Pearson method.

**Application (controller), peak anonymous GiB:**

| Member | import | design | intensity | duplicates | QC finalise | normalisation | correlation | export | DA | report | peak |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lipid screen | .51 | .53 | .64 | .66 | .67 | .59 | .59 | .76 | .87 | 1.02 | **1.02** |
| lipid large | .83 | .88 | .92 | 1.00 | 1.01 | .98 | .99 | 1.15 | 1.54 | 3.33 | **3.33** |
| metab screen | .51 | .53 | .65 | .65 | .65 | .58 | .58 | .61 | .83 | .94 | **.94** |
| metab large | .82 | .92 | .97 | .98 | .99 | .94 | .96 | 1.11 | 1.56 | 3.33 | **3.33** |

**Whole slice (application plus worker), peak anonymous GiB:**

| Member | import | design | intensity | duplicates | QC finalise | normalisation | correlation | export | DA | report | peak |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lipid screen | .85 | .87 | 1.24 | 1.28 | 1.29 | 1.77 | 1.38 | 1.57 | 1.83 | 1.59 | **1.83** |
| lipid large | 1.18 | 1.23 | 3.32 | 2.72 | 2.37 | 3.79 | 2.23 | 2.34 | 3.05 | 3.69 | **3.79** |
| metab screen | .85 | .87 | 1.24 | 1.28 | 1.28 | 1.76 | 1.36 | 1.38 | 1.86 | 1.54 | **1.86** |
| metab large | 1.16 | 1.26 | 3.25 | 2.22 | 2.40 | 3.97 | 2.21 | 2.32 | 3.13 | 3.69 | **3.97** |

**Time, seconds (reported, not judged):**

| Member | import | design | intensity | duplicates | QC finalise | normalisation | correlation | export | DA | report | total |
|---|---|---|---|---|---|---|---|---|---|---|---|
| lipid screen | 7 | 16 | 21 | 20 | 20 | 109 | 22 | 15 | 52 | 11 | **323** |
| lipid large | 17 | 87 | 89 | 89 | 28 | 1,343 | 41 | 30 | 137 | 126 | **2,136** |
| metab screen | 6 | 16 | 17 | 15 | 15 | 110 | 27 | 14 | 59 | 15 | **327** |
| metab large | 12 | 75 | 56 | 25 | 19 | 789 | 26 | 17 | 79 | 59 | **1,235** |

The memory figures above are for one run of each member.

- **The application scales far less than the matrix.** With about 555 times the cells, the application's normalisation peak rises about 1.6 times, and the whole slice's peak about 2.1 times.
- **The worker steps stay off the application.** It holds about 1 GiB through normalisation and correlation.
- **Whole-slice peak including file cache:** 5.63 GiB lipid and 5.36 GiB metab, against 1.94 and 1.96 at screen.

## What still costs what

1. **The report step, in the application: 3.33 GiB on both omics, up from about 1 GiB at screen.** It is the application's largest event at this size, and the only one that scales strongly. Resident memory at the end of the journey is 3.64 to 3.65 GiB (1.21 to 1.33 at screen), so it stays. Its cause was not measured.
2. **DA, in the application: about 1.55 GiB,** up from about 0.85.
3. **Normalisation, in the worker: 3.79 to 3.97 GiB and the most time,** 789 s on metab and 1,343 s on lipid. **The RUV worker is 602 s and 1,065 s of that, 76 to 79%** (`worker-times.txt`). It runs a sweep of 20 negative-control percentages, each with an ANOVA over every feature and a cancor plot per assay, then the third-party `RUVIIIC::RUVIII_C_Varying` once per assay on about 18,000 features. Proteomics uses the same third-party function, and its normalisation step is also its largest in time. The split between sweep and `RUVIII_C_Varying` was not timed, and nor was why lipid's RUV took about 1.8 times metab's.
4. **The intensity filter, in the worker: 3.25 to 3.32 GiB on `whole`.** `partitioned` lowers this (below).

## Variants

`whole` against `partitioned`, on `kit_084`, every member exact and inside its 1,800 s ceiling. Ratios are candidate over control. Decisions as recorded by `kitOperationDecision`, per stage:

| Stage | Decision | Operation | Peak | Anon peak | Time | Not worse (≤ 1.02) | Reduction (≤ .95) |
|---|---|---|---|---|---|---|---|
| `kit-084p` | accepted | lipid intensity filter | .813 | .808 | .637 | yes | **yes** |
| | | metab intensity filter | .769 | .759 | .649 | yes | **yes** |
| `kit-084q` | accepted | lipid log2 | .925 | .928 | .606 | yes | **yes** |
| | | metab log2 | .954 | .944 | .851 | yes | no |
| `kit-084r` | **rejected** (`not_worse`) | lipid duplicates | .972 | .973 | .751 | yes | no |
| | | metab duplicates | **1.042** | **1.052** | .867 | **no** | no |

- The metab duplicates result is valid and unfavourable. It stands and was not retried.
- **No default changed; `whole` stays.** KIT-080's rule would switch a step only with `reduction` true and an accepted decision, and it was written for real data. On these numbers it would select the intensity filter on both omics, and lipid log2 if log2's rule is read per omic. That choice is the user's.
- Each pair ran once, control first, so the time ratios may favour the candidate through the file cache.
- The harness recorded no batch counts for these members. The variant each ran is confirmed from its `timed-member.rds`.

## What this does not establish

- **No old-versus-new claim at this size.** The rules `whole_slice_anon` and `app_process` were not applied, because no control completed. `pair-decisions.R` is committed and tested on a self-pair, and was not used for a decision. The stage's old-versus-new result remains the screen one from KIT-082: application peak 1.22 to 1.02 GiB on lipid and 1.16 to 0.94 on metab, with the whole slice not lower.
- **A partial reading only:** on `kit_032` the application reached 6.05 GiB (`084b`) and 5.89 GiB (`084f`) in an unfinished normalisation, a lower bound. `kit_084`'s application stayed under 1 GiB in that step. No ratio is taken from it.
- **No real-data claim of any kind.**
- **Screen journeys and proteomics were not rerun on `kit_084`.** Its only code change since `kit_082` is the lipid and metab Pearson method. Both large journeys exercise it exactly, the 1:1 test covers it at small size, and the full suite covers the rest.
- **The time comparison is rough.** Other work was on the host throughout (`power-state.txt` records each member's power and load). All `kit_084` members ran on mains in `low-power`.

## Deviations from the declaration

1. **The candidate moved from `kit_082` to `kit_084`, at the user's direction.** `084a` failed at the application's 1,800 s worker deadline in `metab_qc_figures`. The cause was the per-pair Pearson loop: about 8.6 minutes per assay, repeated per assay in every figure pass. The user directed the repair ("i think you can just fix it now by refactoring directly, test its 1:1 on a smaller matrix"), an explicit exception to "kernels are never rewritten" for this one function, carried under this ticket with no separate one. `084a` stands as a valid result of `kit_082`.
2. **`kit_084` is judged on `kit_082`'s oracles,** with no new records. That is the 1:1 proof at size.
3. **The gate's third member never completed, and the gate was not satisfied.** `084b` was stopped when KDE suspended the machine at 06:50:09 for 18m30s, a machine event and not a result of `kit_032`. `084f`, the stage's one correction, was stopped by the user at 8,661 s, inside the ceiling and with no error ("i dont want the fstupid fucking slow kit"). `084e` was not run. Everything after the gate ran at the user's direction, and the `metab84c` and `lipid84c` oracles went unused.
4. **The log2 populations and oracles come from the records' stored states** (`norm-states.R`).
5. **Sleep lock added to the run scripts** after the suspend: `systemd-inhibit`, released on battery at 15%. No power setting or profile changed.
6. **Early launches of `084b` were refused twice** by a dirty worktree from my script. No member ran and no attempt was used.

## The Pearson repair

`pearsonCorForSamplePairs` for lipid and metab now calls one `stats::cor` over the sample matrix, as proteomics does (commit 2c4a453e, helper `calculateSamplePairCorrelations`).

- **1:1 test:** `test-workflow-omics-pearson-matrix.R`, `expect_identical` against the per-pair kernels on 24 samples, with NA, wide scale, zero-variance and all-NA samples, duplicate ids and one-feature assays. It failed its mock check before the change.
- **At size:** 1.97 s for all 19,900 pairs on a simulated 200 × 18,243 matrix, against about 8.6 minutes per 9,900 pairs, extrapolated from 81 pairs.
- **Changed behaviour:** non-numeric sample columns now give a warning and no result for that assay, where the old pivot stopped with an error.
- **Recorded pin change:** `test-metab-01j-s4-pearson-correlation-characterization.R`.
- **Left alone:** the exported per-pair functions, and the repeated Pearson call per assay in the figure code, now about 2 s each.

## Other notes

- **Copy warnings:** 2 to 3 files (heatmaps, significance plots, composite QC) fail to copy into the results summary on the in-process records. This is a headless recorder limit, also present on the 50 × 5,000 trial.
- **Not touched:** the proteomics composite still uses the old reader (`R/mod_prot_norm_qc_support_helpers.R:225`).

## Checks

- **Full suite (one process per file) on `kit_084`:** 358 of 360 files reported, 31,319 pass, 0 fail, 2 error, 313 skip (`full-suite.txt`). The two errors are the known diagnostic errors already in KIT-078's suite, in `test-workflow-perf101-diagnostic.R` and `test-workflow-ruv-pearson-diagnostic.R`. They passed at KIT-075, and neither is in the Pearson change. The two files with no result line pass alone (15 and 29).
- **Dependents of the Pearson change:** 56 files, 11,038 pass, 0 fail, 0 error, 54 skip (`matrix-pearson-dependents.txt`).
- **`R CMD check`:** at baseline, 3 WARNINGs and 3 NOTEs (`check-matrix-pearson.log`).
- **Corpus and trial:** generated and hashed in the bindings; the stepped trial is in `size-trial.txt`.

## Hand-off

Real data is at the final gate: KIT-066, KIT-080 and KIT-083, after KIT-063. Integration (KIT-060) needs synthetic corpora on shared samples, which do not exist yet.

**Decision open for the user:** whether any default moves to `partitioned` on this synthetic evidence.
