# KIT-209 outcome: in-process memory of the release candidate per omic, size and step, with b468a885 beside it

A memory discovery (VAL-151): figures with their run-to-run spread, no verdicts. In process on GitHub-hosted ubuntu-22.04
runners (VAL-145; 4 cores, 16 GB; image ubuntu22 20261116.365.1, kernel 6.8). Ticket revision 1.3 and the coordinator's
rulings of 2026-11-24 (declaration addendum 1, 55fbe19c). Declaration f6b21504; dated section d8fbd7f2 (committed alone
before M1). Runs and minutes: `runs.tsv`. G0: `g0.md`.

- **Candidate (RC):** `91e3cdcb` (staging head 2026-11-24, after KIT-207). **Base:** `b468a885`. Integration base:
  build C `324e14b4`.
- **Every base figure:** base on the candidate's internal dependency versions (`fastduckdb` 0.2.x, `fast_limpa_dpc`
  0ad8c527). The two upgrades since b468a885 are outside what the comparison shows (ruling A4).
- **Settings:** bound identically on both builds (`declaration.md` section 3), except where stated: U11 and U12 with RUV
  **skipped** on both builds (ruling A2); U4 and U5 compared through the session export only (ruling A1).
- **Declared method differences** beside the figures they touch: KIT-177's plain log2 in the LFQ/TMT design step; KIT-187's
  count rule (feature counts differ from `protein_intensity` on S1, S2, U4); the LFQ/TMT base figures include four protein
  QC steps the base's GUI bypassed by default (`g0.md`, B5); transcriptomics' expression filter keeps a different number of
  genes on the two builds (cause not determined, below).

Labels: OBSERVED = in the returned records named; INFERRED = reasoning from them.

## Tables (all OBSERVED)

- `units.tsv`: per unit, the journey peak, the retained-heap maximum and the journey wall, each build's runs, median,
  range and n, the difference (MiB and %), the listing aid, the order contrast.
- `cells.tsv`: the same per unit and step (step process peak, retained R heap at step entry, step wall) with the
  committed-state check per build.
- `findings.tsv`: the rows the listing aid selects (ticket 7.3), largest first. `nonc.tsv`: every non-completion with its
  class.
- `journeys-counted.tsv`, `steps-counted.tsv`: every journey and step record (M1, J1 to J6, D, X10, P3).
- `integration-m3.tsv` (and `m3-M-rows.tsv`, `m3-L-rows.tsv`): M3 per size and stage.

### Per unit (median [min to max], n runs per build)

Journey peak = the largest step `VmHWM` (MiB); retained heap = the largest R heap after a full collection at a step
boundary (MiB); wall = baseline to the last step (s). States committed: every mapped step committed its expected state
on both builds in every completed journey (`cells.tsv`).

| Unit | Journey peak, base | candidate | difference | Retained heap, base / candidate | Wall s, base / candidate |
|---|---|---|---|---|---|
| U1 lipid screen | 1,273.7 [1,271.0 to 1,275.5], 3 | 1,274.6 [1,270.5 to 1,278.0], 3 | +0.9 (+0.1 %) | 535.7 / 537.5 | 138.0 / 140.6 |
| U2 metab screen | 1,224.9 [1,224.9 to 1,230.9], 3 | 1,249.2 [1,248.9 to 1,249.2], 3 | +24.3 (+2.0 %) | 530.2 / 532.2 | 131.1 / 130.7 |
| U3 trx screen | 1,606.4 [1,605.5 to 1,608.6], 3 | 1,641.1 [1,637.2 to 1,650.3], 3 | +34.7 (+2.2 %) | 938.7 / 929.9 | 78.3 / 75.9 |
| U4 DIA screen, through the session export | 1,246.5 [1,221.2 to 1,248.0], 3 | 1,212.0 [1,201.7 to 1,233.6], 3 | −34.5 | 336.2 / 337.0 | 269.8 / 278.4 |
| U4, whole journey, base only (A1) | 1,346.4 [1,335.0 to 1,380.4], 3 | base only; the candidate refuses this design (one biological unit per group) | | 659.7 / — | 309.9 / — |
| U5 DIA_LIMPA screen, through the session export | 1,222.9 [1,200.1 to 1,229.7], 3 | 1,225.9 [1,181.5 to 1,258.2], 3 | +3.0 | 336.0 / 335.8 | 254.0 / 256.5 |
| U5, whole journey, base only (A1) | 1,281.1 [1,264.5 to 1,283.9], 3 | base only; as U4 | | 657.4 / — | 286.6 / — |
| S1 LFQ fixture (5 proteins) | 1,012.1, 1 | 965.0, 1 | −47.1 | 628.3 / 627.7 | 65.2 / 65.5 |
| S2 TMT fixture (5 proteins) | 965.6, 1 | 964.6, 1 | −1.0 | 628.5 / 627.8 | 66.3 / 66.2 |
| U6 lipid large, missing | 3,847.0 [3,820.6 to 3,873.4], 2 | 4,096.8 [4,009.1 to 4,184.4], 2 | +249.8 (+6.5 %) | 1,333.0 / 1,326.4 | 1,038.9 / 1,211.9 |
| U7 metab large, missing | 4,020.4 [3,733.3 to 4,249.8], 3 | 4,243.5 [4,143.3 to 4,420.4], 3 | +223.1 (+5.5 %) | 1,331.8 / 1,325.7 | 808.9 / 924.4 |
| U8 trx 60k × 200 | 5,353.0 [5,337.7 to 5,384.8], 3 | 3,790.4 [3,780.0 to 3,822.7], 3 | −1,562.6 (−29.2 %) | 2,257.6 / 1,567.5 | 409.1 / 242.9 |
| U9 DIA full (300) | 3,987.8 [3,935.2 to 4,237.6], 3 | 4,100.1 [3,791.5 to 4,223.8], 3 | +112.3 (+2.8 %) | 712.5 / 712.9 | 2,730.7 / 2,671.9 |
| U10 DIA_LIMPA full, `fastlimpadpc` | 4,061.4 [3,822.5 to 4,163.8], 3 | 4,080.3 [3,909.0 to 4,205.6], 3 | +18.9 (+0.5 %) | 694.0 / 694.0 | 2,715.2 / 2,737.8 |
| U11 LFQ 300, RUV skipped | 2,818.6 [2,817.8 to 2,819.4], 2 | 2,734.2 [2,733.1 to 2,735.3], 2 | −84.4 (−3.0 %) | 767.3 / 767.7 | 511.6 / 493.3 |
| U12 TMT 300, RUV skipped | 2,809.9 [2,805.1 to 2,821.0], 3 | 2,810.0 [2,783.8 to 2,811.5], 3 | +0.1 (0.0 %) | 767.5 / 767.9 | 917.4 / 933.7 |

Where each journey's peak sits (OBSERVED, both builds): lipid and metab screen at normalisation; trx screen at the RUVIII
export; DIA and DIA_LIMPA full at the precursor rollup; lipid and metab large at the report or normalisation; trx large
at normalisation; LFQ/TMT 300 at the report.

### Run-to-run spread (OBSERVED)

- Journey peak ranges within a build: screen units 0.3 to 45.4 MiB; size units 1.6 to 432.3 MiB (largest on U9, U10 and
  U7, whose peaks sit at the rollup or normalisation and vary between runs by 170 to 430 MiB).
- Retained heap: 0.0 to 16.4 MiB within a build in every unit.
- Order contrast (the candidate-minus-base journey peak in pairs run B C against C B; `units.tsv`): U7 +201.8 (mean of 2)
  against +400.0 (1); U6 +135.7 against +363.8; U10 +115.4 against −83.5; elsewhere within the spread.
- Third pairs (ticket 7.2) were run for U7, U8, U9, U10 and U12; U6 and U11 did not meet the rule.

### Integration (M3; KIT-150's lane, capped; median, n = 3 at M, 2 at L)

| Size, stage | Build C peak MiB | Candidate peak MiB | Build C wall s | Candidate wall s | R heap peak C / candidate |
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
| X10 DIA_LIMPA full, engine `limpa` | completed | 6,103.9 (dpc_quant, 7,925 s) | 10,366.4 | U10 (`fastlimpadpc`): dpc_quant about 2,750 MiB, 272 s |

## Findings for the owner

Listed by ticket 7.3's aid (a difference of the medians above both max(0.25 GiB, 10 % of the base) and twice the larger
range), plus non-completions and other differences the owner should see. An aid for reading, not a judgement.

1. **Transcriptomics at size (U8, `trx_60k_200`): the candidate uses less memory.** Journey peak 3,790.4 against
   5,353.0 MiB (−1,562.6, −29.2 %; ranges 42.7 and 47.1); retained heap 1,567.5 against 2,257.6 MiB (−690.1; ranges 0).
   The difference begins at the expression filter (step peak 1,400.1 against 1,916.1) and holds through normalisation
   (3,790.4 against 5,338.2), DA, the RUVIII export, the report and the session export (`findings.tsv`, 12 rows). Beside
   it: the candidate's expression filter keeps 22,880 genes against the base's 21,247 (the trx screen corpus: 5,095
   against 4,565), and the store after the session export holds 1.09 GB against 2.20 GB. Journey wall 242.9 against
   409.1 s. Which merge made the change was not determined (attribution is the coordinator's follow-up).
2. **RUV at 300 samples refuses on both builds.** With the bound manual RUV (control percentage 5, k 3 from the tab's
   empty default) both builds refuse at `normalisation_ruv`: "RUV-III cannot estimate 3 factors from 0 control proteins
   observed in every sample" (`lfq_300`) and "from 2" (`tmt_300`) (probe run 45383701044 for the base; G0p for the
   candidate). The candidate's shipped automatic sweep also refuses on `lfq_300` ("at every control percentage tried
   (1-20%), only 0 were", D6) but completes on `tmt_300` (D7). INFERRED: the corpora's random missingness (8 % and 2 % of
   values) leaves almost no protein observed in all 300 samples (0.92^300 ≈ 1e-11 and 0.98^300 ≈ 0.002 per protein), so
   any rule that takes controls only from complete proteins finds none or very few at this sample count. This bears on how
   RUV behaves at large sample counts with scattered missing values, not on memory. The DIA_LIMPA shipped default (D5:
   normalisation none, RUV covariates) also stops, at `ruv_covariates`, "No valid percentage found".
3. **Integration: the candidate's MOFA training takes about 4 times as long at equal or lower memory.** S: 197.6 and
   149.1 s against 48.2 and 35.8 s (G0i and its repeat; identical factors, 2,056 iterations both); M: 1,344.1 against
   312.6 s (n 3 each, ranges under 26 s); L: 3,323.2 against 803.1 s (n 2). In the same jobs both builds ran on the same
   system library and R installation, the same bundle Python, thread variables unset (records `runner.txt`,
   `mofa-library-versions.tsv`, each stage's `threads`). Peaks: the candidate's alignment and training are lower at L
   (1,181.8 against 1,557.8 MiB; 1,342.2 against 1,547.7), the pipeline maximum (the results stage) about equal.
4. **DIA and DIA_LIMPA at screen size (U4, U5): the candidate refuses DA** on the screen design (one biological unit per
   group; KIT-128, correct product behaviour by ruling A1). Through the session export the two builds are within their
   spread; DA and the report are measured on the base only.
5. **Not listed by the aid, but larger than the spread:** U6 (lipid large) journey peak +249.8 MiB (+6.5 %) with journey
   wall +173 s (+16.7 %); U7 (metab large) +223.1 MiB (+5.5 %), wall +115.5 s; U11 −84.4 MiB. Step detail in `cells.tsv`
   (U6: normalisation +289.0, intensity −167.8, duplicates +166.4 MiB).
6. **The other DPC engine (X10, candidate):** `limpa` peaks at 6,103.9 MiB in DPC quantification and takes 10,366 s for
   the whole journey; `fastlimpadpc` (U10) about 2,750 MiB there and 2,738 s for the journey.
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

- Which merges cause the U8 difference, the integration training time and the trx expression-filter gene counts.
- The base on its own dependency pins (G0a' not run).
- The GUI layer, the worker path, enrichment, normalisation re-runs and exports (VAL-145; seen by nothing here).
- Spread beyond 2 or 3 runs; whether U7's and U6's order contrasts persist.
- Why the runner's DIA reports differ in bytes from the bound files.

## For the change report

No behaviour changed; measurement only.
