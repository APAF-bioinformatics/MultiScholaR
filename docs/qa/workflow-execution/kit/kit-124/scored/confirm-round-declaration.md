# KIT-124 confirmation round: declaration (2026-11-19, written before any job of the round; the lifecycle run below is the instrument check only)

Authority: the user's decision M29 of 2026-11-18 (`design-review-20261118.md` section 6 on `performance-v2` b35d8089): "Every option's figure is re-measured on the shipped code in a final confirmation round, since only the recommended design was measured under the exact method". With M28 (k limited by biological units; the count removes them), M30 (200 permutations; the classic test applies its trend and robust settings under the new designs and No RUV), M31 to M33 (typed k above the limit stops lipid/metab; k below the number of controls everywhere; the count recorded beside a typed k). Simulations run on GitHub-hosted runners (M13), with the equivalence read on the values the rules use (M19). Nothing under `R/`, `tests/`, `DESCRIPTION` or `tools/` changes on this branch; the package is installed in each job from the bound commit.

Every choice below that is not carried over from an earlier declaration (`exact-round-declaration.md`, the confirming study, `final-20261118`) is marked **[for the user to approve]**.

## Purpose

To give, for every RUV option the app will offer in each omic, the figures measured on the code that ships (the installed package's new-design kernels and its DA entry points), by the ticket's rules 1 to 7 unchanged, and to show for each option that has one whether the shipped code reproduces the simulator arm its earlier figure came from. The figures shown beside the options (`omicGuiRuvTable*` in `R/utils_omic_gui_config.R`) are replaced by this round's where the user accepts them.

## Package under test

**BOUND COMMIT (filled at launch): `kit-129-integration` @ `________________________________________` [for the user to approve: the commit that carries the wiring and M31 to M33].**

- The lifecycle run used `33cfc99b` (`kit-129-integration`, 2026-11-19: M28 to M30 and the vocabulary pass; M31 to M33 and the wiring not yet in it).
- How it is installed (as module-ci and KIT-134 do, no new mechanism): params `package: "installed"`, `package_commit` = `tools_ref` = the bound commit; the job checks out that commit, `.github/actions/snapshot-library` installs R 4.6.1, the KIT-134 snapshot library, the two internal packages built in the job from the read token `DEPS_READ_TOKEN` (never printed, cached or uploaded), and MultiScholaR from that checkout; the copy is kept and the simulators load it (`K124_PACKAGE=installed`). The job checks that the checkout is the bound commit and that the installed namespace holds the new-design kernels.
- The bound commit must be reachable on `origin`: on this run-only branch it is merged with the `ours` strategy (the branch's tree does not change). The lifecycle commit `33cfc99b` was made reachable that way (9c8e6e76).

## Simulators, cells, designs

| Simulator | File | Data generation, seeds, scoring |
|---|---|---|
| S1 classic proteomics | `s1_classic.R` of `sim-k124-exact` 2ba5cf87, with this branch's additions | unchanged |
| S2 DIA with limpa (`fastlimpadpc` engine) | `s2_limpa.R` of `sim-k124-exact` 2ba5cf87, with additions | unchanged |
| S3 transcriptomics | `s3_trx.R` of `sim-k129-exact` 2d03b103 (technical replicates summed at the design step, as the shipped package does), with additions **[for the user to approve: the request named `sim-k124-exact`; the two files differ in C6 and C8 only, and only the `sim-k129-exact` file matches the shipped design step]** | unchanged |
| S4 lipidomics and metabolomics | `s4_lipid_metab.R` of `sim-k129-exact` 2d03b103 (the extract declared as technical replicate), with additions | unchanged |

The additions: the arms 34 to 40 below, the switch `K124_PACKAGE=installed`, and `K124_FAULT` (a deliberate crash, lifecycle only). The existing arms are not edited; under `K124_PACKAGE=installed` they run on the installed package's functions instead of the work tree's.

Cells (carried over): core C0 to C7 (C6 technical duplicates: 4 units per group, 2 injections or libraries each) and C8 (C6 with no unwanted factor; reported, judges no design). Designs with 3 per group exist only among the extra-hard cells (S1 3 per group; S11 with the factor tied to groups at 0.3; S12 with three factors; S13 with no factor). The extra-hard cells S1 to S14 (S3: no S6) are included because the notes shown beside the options quote them (for example lipid/metab "13% at a correlation of 0.6" from S3, "10% at 3 per group" from S1, "14% with an unbalanced batch" from S8) **[for the user to approve: all extra-hard cells, or only those the notes quote: S1, S3, S8, S11, S12, S13]**.

## Options measured, per omic

Arm numbers 34 to 40 are new (17 to 33 were taken). Each new arm calls the installed package's kernels as `kit-129/final-method/prove-arms.R` calls them; each runs on its own random stream; the count's seed is each simulator's E seed, so the shipped count equals the simulators' E at 200 permutations (proven for run 1 of C0 and C6 by `prove-arms.R count`).

| Option (app id) | Classic proteomics (S1) | DIA with limpa (S2) | Transcriptomics (S3) | Lipid/metab (S4) |
|---|---|---|---|---|
| Recommended, `groups_all_features`, estimated count | 34 (beside 32) | 34 (beside 32; in C6/C8 also the technical-replicate design) | 34 (beside 31) | 34 (beside 32) |
| Technical replicates, estimated count (C6, C8) | 38 (beside 33): the kernel handed the unit as replicate set; the kernel documents the groups as its only structure **[gap: no shipped kernel for this option; see "Before the launch"]** | = 34 (one computation on this route) | unavailable (summed at the design step) | 38 (beside 33): `ruv_grouping_variable = "extract"` |
| Supplied controls, estimated count | 35 (oracle list; arm 8's) | 35 | 35 | 35 (the internal standards) |
| Supplied controls, contaminated list | 36 (arm 9's set; not C4) | 36 (not C4) | 36 | 36 (2 of 10 standards changing) |
| Supplied controls, separation rule | existing arm 22 on the installed package **[the kernel takes this k from the wiring: no shipped computation to call]** | (not in the request; the kernel computes it as for 40) | existing arm 22, as classic | existing arm 22, as classic |
| First-pass nulls, estimated count | **gap: no shipped kernel** (existing arm 23 on the installed package reported) | 39 | 39 (beside 23, which has B = 20) | **gap: no shipped kernel** (existing arm 23 reported) |
| First-pass nulls, separation rule | gap (existing arm 5 reported) | 40 (the kernel's default for this rule: today's sweep k) | existing arm 5 / 24 (the kernel takes the k from the wiring) | gap (existing arm 5 reported) |
| Legacy | existing arm 3 on the installed package (today's sweep, `ruvIII_C_Varying`, `runTestsContrasts`: the package's own legacy functions) **[for the user to approve: legacy measured by the simulator's composition of the legacy functions, not through a new entry point]** | existing arms 4 (covariates) and 3 (correction) | existing arm 3 | existing arm 3 |
| None | 37 (beside 1; the final method, M30, would be arm 27: reported beside it) | 37 (beside 1) | 37 (beside 1) | 37 (beside 1) |
| User's k | not simulated (no rule to simulate); the option's figure is the estimated count's, as KIT-129 section 1a states **[for the user to approve]** | as classic | as classic | as classic |

Each run of every cell also runs arms 1 and 2 (the rules read them: rule 3's limit, rule 4's pairing, the void test) and the simulator arms named "beside".

k rules per option, as the bound code defines them: the estimated count (`num_sv_be`, 200 permutations, the groups removed; on technical-replicate designs the biological units removed, M28) wherever defined, pre-selected (M15); the separation rule only where a separate control set exists (refused with all features as controls); a user's k. A count of 0 is no RUV (status ok, k = 0, the shipped no-RUV table). A k above the limit (M28, M31, M32) is a refusal: a classed refusal of the package (`ruv_k_refused`, `ruv_factor_count_unavailable`, `prot_ruv_factors_unavailable`, `omics_ruv_no_residual_df`, the block-correlation and no-replication classes; in S3 also `trx_error`; in S2 a run the kernel leaves without RUV for another reason than a count of 0; in S4 an assay the kernel refuses) is `refused` and scored as arm 1 (R-07, carried over). Any other error is `crashed` and counts in rule 6.

## Runs and seeds (carried over)

1,000 runs per core cell (runs 1 to 1,000), 500 per extra-hard cell (runs 1 to 500); seed = 100000 + 10000 x cell index + run; the data under `set.seed(seed)`; arm a under `set.seed(seed * 100 + a)` (S2: `100 * seed + a`), arms 34 to 40 on their own number; E's seed as each simulator has it. The rows of an arm are the same whichever other arms run.

## Pass rules (carried over unchanged)

Rules 1 to 7 of `KIT-124-rev1.4-ticket-20261118.md` (W3), computed by `aggregate.R` as on `sim-k124-exact` (sha256 `089d852c...`, C8 a reported no-factor cell), with their limits unchanged: rule 1 null rate at 0.05 and 0.01 (upper bounds 0.075 and 0.015), rule 2 mean false discovery proportion (0.075), rule 3 share of runs above 0.10 against arm 2's share + 0.05, rule 4 the reduction against arm 1 in cells with a factor, rule 5 power loss against arm 1 in cells without one (0.02), rule 6 crash share (0.01), rule 7 global null C4 (0.075). Void cells (arm 2 failing rule 1 or 2), the design verdict over the judged core cells, the primary subset (`full`) and the extra-hard cells listing their failing rules without failing a design: as before. `aggregate.R`'s recommendation table (`k124-omics.tsv`) does not know arms 34 to 40 and is not read; the per-arm verdicts (`k124-arms.tsv`, `k124-cells.tsv`) are. The figure per option is read as before: mean false discovery proportion in the worst core cell with true effects, and power in C0.

## Gate and agreement

1. **Gate (carried over, M19 reading).** Before a cell is scored, arms 1 and 2 of the round equal the earlier hosted-runner rows of the same seeds under `packaging/verify/compare-results.R` (relative 1e-6, absolute floor 1e-9, discrete values exact), read on the values the rules use (status, k, counts, null rates, false discovery proportion, power); bias and rmse are reported. References: `exact-gha` (run 44233215700) for S1, S2 C0 to C8 and S3, S4 C0 to C5, C7; `exact-k128-gha` (44233219582) for S3, S4 C6, C8; the runner rows of `final-20261118` (44201244793) and the laptop rows of `s*-core` for the extra-hard cells where they ran. A cell that does not meet the gate is not scored and is reported.
2. **Agreement of the shipped code with the arm it replaces.** Per pair (shipped arm beside simulator arm, same cell, run and subset), the same rule, read the same way: S1 34~32, 38~33, 37~1; S2 34~32, 37~1; S3 34~31, 37~1; S4 34~32, 38~33, 37~1. Computed in the workflow's compare job (`.github/k124-sim/compare-arms.R` then the bound commit's `compare-results.R`). Declared now, before any row exists: S4 34~32 and 38~33 are **not expected** to agree (S4 forms W on the features complete in study and pooled-QC injections; the package holds the study injections only); S1 37~27 is reported and **not expected** to agree (M30's No RUV is not wired on the bound commit; its DA is `runTestsContrasts`). **[for the user to approve]** What follows from the reading: where a pair agrees, the earlier figure of the simulator arm is confirmed on the shipped code; where it does not, the shipped arm's own figure is the one reported (it is what ships) and the difference is stated with its cause before any figure goes into the app. Options without a counterpart (35, 36, 39, 40) have no agreement check: their figures are the shipped code's.

## Runner budget and attempts

One job per simulator and cell (S1 cells split in run ranges where one job would not finish in 300 minutes, the `chunk` entry; the rows are the same). Per-job limit 300 minutes for the simulator (the job's own limit 345). One attempt per job; one bounded correction: a job that fails for an infrastructure reason may be rerun once as the same command (`gh run rerun <run> --failed`, the same rows). A job that reaches its limit, or a crash of a shipped-kernel arm (the job fails after its upload, `check-arms.R`), is reported, not relaunched with other settings. The size is in "Sizing" below; above about 120 core-hours the round is not launched without the user's word (carried over).

## What may be uploaded

The simulators' result tables (`s<i>-runs.tsv`, `s2-data.tsv`, `s2-errors.tsv`, `s3-events.tsv`, `s3-timing.tsv`, `s4-notes.tsv`, `s4-timing.tsv`, the per-cell restart files), their meta and log files, the runner's library and job tables, and the compare job's tables. No console output (in this mode the simulator's output stays in the job log), no install log, nothing built with the token; every file is scanned for build output, internal source locations and token shapes before the upload (`check-output.R`, carried over). Artifacts kept 30 days; the per-run tables are not committed (sha256 recorded), as before.

## Before the launch (what must be bound or changed)

| # | Item | Why |
|---|---|---|
| 1 | Bind the commit (the line above), merge it into this branch with the `ours` strategy, set `package_commit` = `tools_ref` to it in `params.json` | The package under test; it must be reachable on `origin` |
| 2 | **Technical replicates in classic proteomics and lipid/metab: the kernels refuse the design.** Handed the unit as the replicate set, both count the 8 units as 8 groups and leave k at most 0 (lifecycle: arm 38 refused in every run; classic `prot_ruv_factors_unavailable`, lipid `ruv_k_above_limit`, "limit of 0"). The app lists the option as measured (arm 16). Needed: a kernel that takes the units as RUV-III's replicate sets while counting the groups for the limit, or the option withdrawn there **[the user's decision; wiring]** | Otherwise the round measures a refusal, not the option |
| 3 | **First-pass nulls have no kernel in classic proteomics and lipid/metab** (the app offers them in every omic, M7). Until one exists the round can only report the simulators' arms 5 and 23 run on the installed package's functions **[the user's decision]** | The option's figure would not be the shipped code's |
| 4 | **The separation rule is computed by the wiring, not by the kernels**, in classic proteomics, lipid/metab and transcriptomics (the kernels take that k as given); only DIA's kernel computes it (arm 40). Either the wiring's computation is callable at the bound commit, or the separation-rule figures stay those of arms 22, 5 and 24 **[the user's decision]** | As 3 |
| 5 | **No RUV in classic proteomics is not yet M30's** at `33cfc99b`: arm 37 equals arm 1 (today's `runTestsContrasts`) and differs from arm 27 (M30) in the null rate at 0.01 (2 of 4 rows), as "found and left" (1) of the M29/M30 pass says. The bound commit should carry it | M30 |
| 6 | M31, M32, M33 are not in `33cfc99b`. M31 and M33 concern a typed k, which the round does not simulate; M32 changes the limit in lipid/metab and so its refusals | The bound commit's limit is the one measured |
| 7 | The sizing below is far above the carried-over limit of about 120 core-hours: the user chooses the variant **[the user's decision]** | Carried-over limit |

## Lifecycle run (2026-11-19, the one run of this preparation)

Run 44341242746 (branch `sim-k124-confirm`, commit 240c97d1, label `confirm-lc`): package installed in each job from `kit-129-integration` 33cfc99b; one cell per simulator (S1 C6, S2 C6, S3 C0, S4 C6), runs 1 and 2, 2 cores; S4 arm 37 crashed on purpose in run 1 (`K124_FAULT`).

| Check | Result |
|---|---|
| Install | Every job built the library with the internal packages from the read token, kept the bound commit's MultiScholaR 0.5.0, confirmed the checkout was 33cfc99b and that the kernels exist; the library check against the laptop's versions passed |
| New arms ran to completion | S1 34 to 37, S2 34 to 37, 39, 40, S3 34 to 37, 39, S4 34 to 36: every row `ok`. S1 38 and S4 38 `refused` in every run (item 2). S4's lipid and metab tables agreed in every arm |
| Failure is reported | S4 arm 37's deliberate crash: rows `crashed` with the message in `s4-notes.tsv`, the table uploaded, then the job failed at `check-arms.R`; the run concluded `failure`. The other three simulator jobs and the compare job succeeded |
| Table assembles | `assemble.sh 44341242746 confirm-lc` laid out all four simulators; the compare job's artifact `k124-confirm-lc-compare` holds the pairs |
| Agreement (fixed rule, M19 values) | S1 34~32, 37~1; S2 34~32, 37~1; S3 34~31, 37~1: identical (largest relative difference 0). S1 37~27: outside on the null rate at 0.01 only (item 5, declared). S1 and S4 38~33: refused against computed (item 2). S4 34~32: run 1 identical, run 2 different calls (the declared cause). S4 37~1: run 2 identical in every measure (k 0 against arm 1's NA, a convention); run 1 the deliberate crash |

Seconds per run inside a worker (mean of the two runs; each shipped arm runs its own count with 200 permutations, and S2's arm 40 its own sweep):

| Simulator, cell | Shipped arms | Counterparts |
|---|---|---|
| S1 C6 | 34, 35, 36: 2.5 to 2.7; 37: 1.4; 38: 0.05 (refuses) | 32: 1.6; 33: 1.1; 1: 0.8; 27: 0.8 |
| S2 C6 | 34, 35, 36: 8.6 to 9.0; 37: 2.5; 39: 10.2; 40: 36.0 | 32: 4.4; 1: 1.8 |
| S3 C0 | 34, 35, 36, 39: 13.3 to 14.8; 37: 0.5 | 31: 7.3; 1: 0.4 |
| S4 C6 | 34, 35, 36: 2.6 (both subsets); 37: 0.6; 38: 2.3 | 32, 33: 0.3; 1: 1.2 |

## Sizing (estimates from the lifecycle; not launched)

Core-hours = seconds per run x runs / 3600, at the lifecycle's per-run cost, with the earlier rounds' cost for arms not in the lifecycle (S1's sweep 23 to 27 s for arm 3; S2's sweep about 30 s for arms 3 to 5; S3's grid). One cell was measured per simulator, so a cell's cost may differ by about half either way.

| Variant | S1 | S2 | S3 | S4 | Core-hours, about | Jobs (300-minute limit, 4 cores) | Wall time at 20 concurrent jobs |
|---|---|---|---|---|---|---|---|
| A. Core C0 to C8, 1,000 runs; every option arm, arms 1 and 2, the counterparts, legacy recomputed | 105 | 250 | 210 | 35 | 600 | about 45 | about 8 h |
| B. As A without legacy recomputed (its code is unchanged, M3) and without S2's arm 40 | 40 | 120 | 160 | 35 | 355 | about 30 | about 5 h |
| C. B plus the extra-hard cells at 500 runs | 70 | 215 | 280 | 65 | 630 | about 55 | about 9 h |
| D. A plus the extra-hard cells at 500 runs | 185 | 445 | 365 | 65 | 1,060 | about 80 | about 14 h |

Every variant is above the carried-over limit of about 120 core-hours; none is launched without the user's word. No launch file is written: `params.json` keeps the lifecycle's parameters and a push that leaves it alone starts nothing.
