# KIT-209 outcome: in-process memory of the release candidate per omic, size and step, with b468a885 beside it

A memory discovery (VAL-151): figures with their run-to-run spread, no verdicts. In process on GitHub-hosted ubuntu-22.04
runners (VAL-145; 4 cores, 16 GB, processor by job: AMD EPYC 7763, 9V74 or 9V45 (`runner.txt`); image ubuntu22 20261116.365.1, kernel 6.8). Ticket revision 1.3 and the coordinator's
rulings of 2026-11-24 (declaration addendum 1, 55fbe19c). Declaration f6b21504; dated section d8fbd7f2 (committed alone
before M1). Runs and minutes: `runs.tsv`. G0: `g0.md`.

- **Candidate (RC):** `91e3cdcb` (staging head 2026-11-24, after KIT-207). **Base:** `b468a885`. Integration base:
  build C `324e14b4`.
- **Every figure describes 91e3cdcb;** staging has since changed files under `R/`, `inst/` and `DESCRIPTION` (ticket 9.4's
  re-measure map applies if the release is cut later).
- **Build labels in the tables:** in `units.tsv`, `cells.tsv`, `findings.tsv`, `journeys-counted.tsv` and
  `steps-counted.tsv`, build `C` is the candidate 91e3cdcb; in `integration-m3.tsv`, `m3-*-rows.tsv` and G0i, `build_C` / `C`
  is the integration base 324e14b4 and `RC` the candidate; in `g0.md`, `C` is the staging head c03cb02a.
- **Every base figure:** base on the candidate's internal dependency versions (`fastduckdb` 0.2.x, `fast_limpa_dpc`
  0ad8c527). The two upgrades since b468a885 are outside what the comparison shows (ruling A4).
- **Settings:** bound identically on both builds (`declaration.md` section 3), except where stated: U11 and U12 with RUV
  **skipped** on both builds (ruling A2); U4 and U5 compared through the session export only (ruling A1).
- **Declared method differences** beside the figures they touch: KIT-177's plain log2 in the LFQ/TMT design step; KIT-187's
  count rule (feature counts differ from `protein_intensity` on S1, S2, U4); the LFQ/TMT base figures include four protein
  QC steps the base's GUI bypassed by default (`g0.md`, B5); on the candidate only, and not behind a setting: KIT-128 (b396a2e6; NEWS; ledger M6, DAS-048) sums transcriptomics' declared technical replicates at the design step (U3: 12 summed samples against the base's 24 runs; U8: 91 against 200, see Finding 1), and fits lipidomics' and metabolomics' DA with the biological unit as block (`duplicateCorrelation`; 12 blocks for 24 samples at screen size, 100 for 200 at size); KIT-172: the lane leaves the normalisation tab's log offset at each build's UI default, log2(x + 1) on the base and plain log2 on the candidate (lipidomics and metabolomics; at most 1.44e-4 log2 on these intensity tables). Integration (candidate against build C 324e14b4): KIT-151 INT-14 trains five models per run on the candidate (`training_seeds` 5, provisional) against one on C; KIT-155 pages the alignment.

Labels: OBSERVED = in the returned records named; INFERRED = reasoning from them.

## Tables (all OBSERVED)

- `units.tsv`: per unit, the journey peak, the retained-heap maximum and the journey wall, each build's runs, median,
  range and n, the difference (MiB and %), the listing aid, the order contrast; and for U4 and U5 the rows "through
  session_export" (ruling A1).
- `cells.tsv`: the same per unit and step (step process peak, retained R heap at step entry, step wall) with the
  committed-state check per build.
- `findings.tsv`: the rows the listing aid selects (ticket 7.3), largest first. `nonc.tsv`: every non-completion with its
  class.
- `journeys-counted.tsv`, `steps-counted.tsv`: every journey and step record (M1, J1 to J6, D, X10, P3).
- `integration-m3.tsv` (and `m3-M-rows.tsv`, `m3-L-rows.tsv`): M3 per size and stage.

### Per unit (median [min to max], n runs per build)

Journey peak = the largest step `VmHWM` (MiB); retained heap = the largest R heap after a full collection at a step
boundary (MiB); wall = baseline to the last step (s), including 13 to 27 s per journey of the lane's own reads between steps (11 to 28 % of a screen journey); step seconds alone are in `cells.tsv`. Counters: the step peak is the R process's resident-set high-water mark (`VmHWM`): native allocations and resident file pages included, child processes (pandoc) and page cache excluded. The cgroup's `memory.peak`, which includes both, is 1.06 to 1.98 times the journey peak (`journeys-counted.tsv`, `cgroup_peak_mib`; U9 candidate 8,008.9 MiB). States committed: every mapped step committed its expected state
on both builds in every completed journey (`cells.tsv`).

| Unit | Journey peak, base | candidate | difference of medians | Paired differences (MiB) | Retained heap, base / candidate | Wall s, base / candidate | cgroup peak, median, base / candidate |
|---|---|---|---|---|---|---|---|
| U1 lipid screen | 1,273.7 [1,271.0 to 1,275.5], 3 | 1,274.6 [1,270.5 to 1,278.0], 3 | +0.9 (+0.1 %) | +4.3, +3.6, −5.0 | 535.7 / 537.5 | 138.0 / 140.6 | 1,574.5 / 1,576.3 |
| U2 metab screen | 1,224.9 [1,224.9 to 1,230.9], 3 | 1,249.2 [1,248.9 to 1,249.2], 3 | +24.3 (+2.0 %) | +18.0, +24.3, +24.3 | 530.2 / 532.2 | 131.1 / 130.7 | 1,503.0 / 1,505.8 |
| U3 trx screen (a) | 1,606.4 [1,605.5 to 1,608.6], 3 | 1,641.1 [1,637.2 to 1,650.3], 3 | +34.7 (+2.2 %) | +31.7, +43.9, +32.5 | 938.7 / 929.9 | 78.3 / 75.9 | 1,794.8 / 1,811.7 |
| U4 DIA screen, through the session export | 1,246.5 [1,221.2 to 1,248.0], 3 | 1,212.0 [1,201.7 to 1,233.6], 3 | −34.5 | −12.9, −9.2, −46.3 | 336.2 / 337.0 | 269.8 / 278.4 (step seconds) | (c) |
| U4, whole journey, base only (A1) | 1,346.4 [1,335.0 to 1,380.4], 3 | base only; the candidate refuses this design (one biological unit per group) | | | 659.7 / — | 309.9 / — | 1,941.0 / — |
| U5 DIA_LIMPA screen, through the session export | 1,222.9 [1,200.1 to 1,229.7], 3 | 1,225.9 [1,181.5 to 1,258.2], 3 | +3.0 | −48.2, +25.8, +35.3 | 336.0 / 335.8 | 254.0 / 256.5 (step seconds) | (c) |
| U5, whole journey, base only (A1) | 1,281.1 [1,264.5 to 1,283.9], 3 | base only; as U4 | | | 657.4 / — | 286.6 / — | 1,718.7 / — |
| S1 LFQ fixture (5 proteins) | 1,012.1, 1 | 965.0, 1 | −47.1 | −47.1 | 628.3 / 627.7 | 65.2 / 65.5 | 1,302.0 / 1,250.0 |
| S2 TMT fixture (5 proteins) | 965.6, 1 | 964.6, 1 | −1.0 | −1.0 | 628.5 / 627.8 | 66.3 / 66.2 | 1,255.7 / 1,250.4 |
| U6 lipid large, missing | 3,847.0 [3,820.6 to 3,873.4], 2 | 4,096.8 [4,009.1 to 4,184.4], 2 | +249.8 (+6.5 %) | +135.7, +363.8 | 1,333.0 / 1,326.4 | 1,038.9 / 1,211.9 | 5,426.0 / 5,817.5 |
| U7 metab large, missing | 4,020.4 [3,733.3 to 4,249.8], 3 | 4,243.5 [4,143.3 to 4,420.4], 3 | +223.1 (+5.5 %) | −106.5, +400.0, +510.2 | 1,331.8 / 1,325.7 | 808.9 / 924.4 | 5,503.8 / 5,759.9 |
| U8 trx 60k × 200 (a) | 5,353.0 [5,337.7 to 5,384.8], 3 | 3,790.4 [3,780.0 to 3,822.7], 3 | −1,562.6 (−29.2 %) | −1,515.0, −1,562.6, −1,604.8 | 2,257.6 / 1,567.5 | 409.1 / 242.9 | 5,983.3 / 4,366.3 |
| U9 DIA full (300) (b) | 3,987.8 [3,935.2 to 4,237.6], 3 | 4,100.1 [3,791.5 to 4,223.8], 3 | +112.3 (+2.8 %) | −143.7, −13.8, +112.3 | 712.5 / 712.9 | 2,730.7 / 2,671.9 | 7,889.2 / 8,008.9 |
| U10 DIA_LIMPA full, `fastlimpadpc` | 4,061.4 [3,822.5 to 4,163.8], 3 | 4,080.3 [3,909.0 to 4,205.6], 3 | +18.9 (+0.5 %) | +86.5, −83.5, +144.2 | 694.0 / 694.0 | 2,715.2 / 2,737.8 | 5,717.1 / 5,844.8 |
| U11 LFQ 300, RUV skipped | 2,818.6 [2,817.8 to 2,819.4], 2 | 2,734.2 [2,733.1 to 2,735.3], 2 | −84.4 (−3.0 %) | −82.5, −86.3 | 767.3 / 767.7 | 511.6 / 493.3 | 3,977.9 / 3,895.4 |
| U12 TMT 300, RUV skipped | 2,809.9 [2,805.1 to 2,821.0], 3 | 2,810.0 [2,783.8 to 2,811.5], 3 | +0.1 (0.0 %) | +6.4, −26.1, −11.0 | 767.5 / 767.9 | 917.4 / 933.7 | 3,754.2 / 3,753.7 |

Paired differences: candidate minus base, in pair order B C, C B, B C. (a) U3 and U8: not like for like: 12 / 91 summed
samples against 24 / 200 runs (Finding 1). (b) U9: the difference of the medians (+112.3) is not the paired difference
(median −13.8). (c) The cgroup peak is recorded per whole journey; for U4 and U5 the candidate's journey ends at its
refused DA step (1,523.2 and 1,477.9 MiB), so no through-export cgroup figure is given.

Where each journey's peak sits (OBSERVED, both builds): lipid and metab screen at normalisation; trx screen at the RUVIII
export; DIA and DIA_LIMPA full at the precursor rollup; lipid and metab large at the report or normalisation; trx large
at normalisation; LFQ/TMT 300 at the report.

### Run-to-run spread (OBSERVED)

- Journey peak ranges within a build: screen units 0.3 to 76.7 MiB (U5 candidate through the session export); size units 1.6 to 516.5 MiB (U7 base). The large ranges sit at the rollup, normalisation and report steps.
- Retained heap: 0.0 to 16.4 MiB within a build in every unit.
- Pairs and order: per pair, candidate minus base, B C / C B / B C: U6 +135.7 / +363.8; U7 −106.5 / +400.0 / +510.2; U10 +86.5 / −83.5 / +144.2. The two orders do not separate: within one order the pairs differ by more than between orders, and the retained heap is identical in both orders. Each job ran both builds on one processor; processors differed between jobs (AMD EPYC 7763, 9V74, 9V45), which moves wall times by up to about 40 %, so differences are stated per pair.
- Third pairs (ticket 7.2) were run for U7, U8, U9, U10 and U12; U6 and U11 did not meet the rule.

### Integration (M3; KIT-150's lane, capped; median, n = 3 at M, 2 at L): candidate against build C (324e14b4; not b468a885), both on the candidate's snapshot library (G0i)

| Size, stage | Build C (324e14b4) peak MiB | Candidate peak MiB | Build C (324e14b4) wall s | Candidate wall s | R heap peak, build C (324e14b4) / candidate |
|---|---|---|---|---|---|
| M write | 933.0 | 927.8 | 3.9 | 3.9 | 110.3 / 110.4 |
| M import | 1,155.4 | 1,155.1 | 43.0 | 42.9 | 471.8 / 471.9 |
| M align | 1,185.1 | 1,066.8 | 17.3 | 18.6 | 411.7 / 335.0 |
| M train | 1,210.7 | 1,238.4 | 312.6 | 1,344.1 | 540.6 / 463.9 |
| M results | 1,098.8 | 1,175.5 | 5.5 | 6.7 | 382.9 / 494.0 |
| L write | 1,766.2 | 1,766.8 | 18.3 | 18.2 | 126.3 / 126.4 |
| L import | 1,565.8 | 1,565.7 | 145.0 | 146.2 | 588.9 / 589.2 |
| L align | 1,557.8 | 1,181.8 | 63.8 | 82.0 | 744.9 / 456.7 |
| L train | 1,547.7 | 1,342.2 | 803.1 | 3,323.2 | 847.1 / 558.9 |
| L results | 2,391.9 | 2,413.8 | 18.2 | 19.6 | 1,511.5 / 1,523.3 |

Every stage completed and was valid. Pipeline peak (the largest stage peak): M 1,238.4 (candidate, training) against
1,210.7 MiB (C, training); L 2,413.8 against 2,391.9 MiB (both at the results stage). Peak ranges within a build at most
28.9 MiB.

### The candidate's shipped defaults and the other engine (candidate only, one run each; reported)

| Run | Outcome | Journey peak MiB (step) | Wall s | Beside |
|---|---|---|---|---|
| D1 lipid large | completed | 3,843.2 (report) | 902.7 | U6 candidate 4,096.8 |
| D2 metab large | completed | 3,847.2 (normalisation) | 712.9 | U7 candidate 4,243.5 |
| D3 trx 60k × 200 | completed | 3,811.7 (normalisation) | 241.0 | U8 candidate 3,790.4 |
| D4 DIA full | completed | 4,280.7 (rollup) | 2,811.7 | U9 candidate 4,100.1 |
| D5 DIA_LIMPA full (normalisation none, covariates) | refused at `ruv_covariates`: "No valid percentage found" | 4,195.3 (rollup), 10 of 14 steps | 2,272.1 | U10 |
| D6 LFQ 300 (automatic sweep) | refused at `normalisation_ruv`: "at every control percentage tried (1-20%), only 0 were [observed in every sample] (k = 2)" | 2,001.4 (accession), 6 of 11 | 968.4 | U11 |
| D7 TMT 300 (automatic sweep) | completed | 2,919.8 (report) | 1,098.1 | U12 (RUV skipped) 2,810.0 |
| X10 DIA_LIMPA full, engine `limpa` | completed | 6,103.9 (dpc_quant, 7,925 s) | 10,366.4 | U10 candidate (`fastlimpadpc`, n 3): dpc_quant 2,778.2 MiB [2,740.4 to 2,917.3], 280.8 s [279.9 to 281.4] |

## Findings for the owner

Listed by ticket 7.3's aid (a difference of the medians above both max(0.25 GiB, 10 % of the base) and twice the larger
range), plus non-completions and other differences the owner should see. An aid for reading, not a judgement.

1. **Transcriptomics (U8 `trx_60k_200`; U3 `trx_screen`): the two builds analyse different numbers of samples, so the memory figures are not like for like.** On the candidate the design step sums declared technical replicates (KIT-128, b396a2e6; NEWS; ledger M6, DAS-048); its logs read "the counts of 200 runs summed into 91 biological samples" (U8; the corpus declares 100, cause not determined) and "24 runs summed into 12" (U3); the base analyses every run. The expression filter's code is unchanged (`git log b468a885..91e3cdcb -- R/func_trx_filter.R` is empty); it keeps 22,880 genes on the candidate against 21,247 (U3: 5,095 against 4,565) because its input differs. Observed on U8, candidate against base: journey peak 3,790.4 against 5,353.0 MiB (paired differences −1,515.0, −1,562.6, −1,604.8; ranges 42.7 and 47.1), retained heap 1,567.5 against 2,257.6 MiB, journey wall 242.9 against 409.1 s (step seconds 218.2 against 383.9), store after the session export 1.09 against 2.20 GB. These differences follow the change in the analysed design; they are not an engineering difference at equal input.
   **OPEN (K209-X02):** the corpus declares 100 biological samples, the candidate's log says 91. Unexplained here; another
   agent is investigating it. Until then the candidate's U8 and D3 figures carry "the candidate's design summed the 200
   runs into 91 units against the corpus's 100 (cause not determined)".
2. **RUV at 300 samples refuses on both builds.** With the bound manual RUV (control percentage 5, k 3 from the tab's
   empty default) both builds refuse at `normalisation_ruv`: "RUV-III cannot estimate 3 factors from 0 control proteins
   observed in every sample" (`lfq_300`) and "from 2" (`tmt_300`) (probe run 45383701044 for the base; G0p for the staging head c03cb02a; the candidate 91e3cdcb was not run under manual RUV). The candidate's shipped automatic sweep also refuses on `lfq_300` ("at every control percentage tried
   (1-20%), only 0 were", D6) but completes on `tmt_300` (D7). INFERRED: the corpora's random missingness (8 % and 2 % of
   values) leaves almost no protein observed in all 300 samples (0.92^300 ≈ 1e-11 and 0.98^300 ≈ 0.002 per protein), so
   any rule that takes controls only from complete proteins finds none or very few at this sample count. This bears on how
   RUV behaves at large sample counts with scattered missing values, not on memory. The corpora draw missing values independently per value (`kitProteinRecipe`, `missingness = "random"`); with intensity-dependent missingness, abundant proteins are usually complete. The rule and its remedies are owned by KIT-114 (d), KIT-126 and KIT-210. Separately, the DIA_LIMPA shipped default (D5: normalisation none, RUV covariates) stops at `ruv_covariates` with "No valid percentage found. Please check your data and parameters.", which gives no reason; DPC quantification leaves no missing values, so this is not the complete-control case (INFERRED).
3. **Integration (candidate against build C 324e14b4, both on the candidate's snapshot library): the candidate trains five models per run, C one; each model takes the same time.** The training stage's marks show five rounds on the candidate (`training_seeds` defaults to 5 since 4f22808c, KIT-151 INT-14; provisional; no NEWS entry) and one on C. Per round: S 27.5 to 43.0 s against 32.1 and 42.1 s; M 248.8 to 266.9 s against 256.3 to 260.7 s; L 608.4 to 642.6 s against 619.8 and 653.5 s. Stage walls: M 1,344.1 against 312.6 s (n 3), L 3,323.2 against 803.1 s (n 2). The reported model is the first seed's on both builds (identical factors and iterations). The lane imports numpy before training, so the single-thread maths pin (INT-17, 2dbfef61) applied on neither build; in the application it applies, so a user's training time on the candidate was not measured here and is expected to be at least these figures. Peaks: at M the candidate's training stage is 1,238.4 against 1,210.7 MiB and its results stage 1,175.5 against 1,098.8 MiB; at L its alignment (1,181.8 against 1,557.8) and training (1,342.2 against 1,547.7) are lower; the candidate's Python peak is the fifth round's and is not comparable with C's.
4. **DIA and DIA_LIMPA at screen size (U4, U5): the candidate refuses DA** on the screen design (one biological unit per
   group; a refusal by design, KIT-128, 40ae84bd; ruling A1). Through the session export U4 is paired −12.9, −9.2, −46.3 MiB
   (ranges 26.8 and 31.9) and U5 −48.2, +25.8, +35.3 MiB (ranges 29.6 and 76.7); DA and the report are measured on the base only.
5. **Not listed by the aid.** U6 (lipid large, 2 pairs): journey peak paired differences +135.7 and +363.8 MiB (ranges 52.8 and 175.3); journey wall +180.9 and +165.2 s, of which the DA step +149 s (the candidate's block fit, KIT-128). U7 (metab large, 3 pairs): journey peak paired differences −106.5, +400.0, +510.2 MiB (ranges 516.5 and 277.1), within its spread; journey wall +109.0, +115.5, +177.3 s, of which the DA step +100 s. U11 (2 pairs): −82.5 and −86.3 MiB (ranges 1.6 and 2.2). U2 and U3 at screen size: +18 to +44 MiB in every pair. Step detail in `cells.tsv`.
6. **The other DPC engine (X10, candidate):** `limpa` peaks at 6,103.9 MiB in DPC quantification (5.96 GiB) and takes 10,366 s for
   the whole journey (one run); `fastlimpadpc` (U10) 2,778.2 MiB there (median of 3) and 2,738 s for the journey.
7. **Corpora:** the runner's `screen` and `full` DIA `report.tsv` and the trx `truth.json` differ from the bound laptop
   files' sha256, identically in every job (G0a-1, G0a-2; `full` in J3, J4 and X10): deterministic on the runner, not
   run-to-run; which part of the generator causes it was not determined. Within each job both builds read one generated
   copy (`declaration.md` addendum 1, A5). Bindings not regenerated.

## Runner time

| Rows | Ticket expected / ceiling min | Used min |
|---|---|---|
| G0 (G0a with two instrument runs, G0i, G0p) | 549 / 1,350 | 269.1 |
| Probes after G0 (base manual RUV, U6/U8, G0i repeat) | (G0 re-run and G0p rows) | 93.5 |
| M1 | 115 / 240 | 112.6 |
| M2 (J1 to J6) with D1 to D5 appended | 1,280 / 1,710 | 798.1 |
| D6, D7 | 366 / 595 | 55.7 |
| X10 | 193 / 330 | 187.3 |
| M3 | 220 / 510 | 278.6 |
| P3 (third pairs) | 127 / 1,620 | 317.8 |
| Small jobs (parameters, summaries) | | about 3 |
| **Total** | **2,850 / 13,125** | **about 2,116 min = 141 core-hours** (ticket 190 expected, 875 ceiling) |

## Not checked

- Why the candidate summed U8's 200 runs into 91 samples where the corpus declares 100 (K209-X02; OPEN, under investigation elsewhere).
- The training time a user gets with the maths pin in effect (INT-17), and the time per round with `training_seeds = 1`.
- The per-step sample count (not recorded).
- Whether the base's own dependency pins change any figure (G0a' not run).
- The GUI layer, the worker path, enrichment, normalisation re-runs and exports (VAL-145; seen by nothing here).
- Spread beyond 2 or 3 runs.
- Why the runner's DIA reports differ in bytes from the bound files.

Every figure describes 91e3cdcb; staging has since changed files under `R/`, `inst/` and `DESCRIPTION` (ticket 9.4).

## For the change report

No behaviour changed; measurement only.
