# KIT-275: lifecycle record (scripts written; nothing run except the one smoke send)

Branch `kit-275`, from `origin/v0.9-staging` at `ec27f69e0` (the bound commit of every tier but V: it holds KIT-255 with its I01 repair, merge aa119308f, VAL-358).
Ticket r1.2 is the frozen protocol; nothing in it was changed. No R was run on this machine; no simulation, no emulation of R's random-number generator; Python
(`python3 -I`) was used for arithmetic, the self-test's predictions, the E0 diff and bracket and quote balance checks of the R scripts only.
**This file is the implementer's record for the script audit. It is not a result: no estimand has been produced.**

**Bounded work item (stated before editing).** Outcome: the scripts of KIT-275's measurement (simulation, scorer, G1 fixture, scorer self-test, E0 check), their
exactness records, the job plan and the running-ceiling ledger, one smoke send. Files: `docs/qa/workflow-execution/kit/kit-275/` only. Exit: this record, handed to the
script audit. Not done: L1, G1, G2, any counted job, any result, `declaration.md`, the ledger row (`VAL-<merge>`, written at the end of the ticket), the KIT-242-dependent tier V run.

## 1. Step 0: what the ticket cites, against the bound commit and against 0866982f

- **Package lines (the ticket's r1.2 cites 75acd421)**, compared by text at `ec27f69e0` (`step0-lines.txt`): unmoved: `R/utils_workflow_omics_declarations.R:692-723, 965, 988, 1036-1044, 1047`;
  every cited line of `R/func_omics_ruv_factors.R`, `R/func_prot_norm_optimization_helpers.R`, `R/func_lipid_s4_ruv_methods.R`, `R/func_lipid_da_model.R`, `R/func_general_ruv_factor_count.R`,
  `R/func_lipid_s4_ruv_control_methods.R`, and `R/func_lipid_s4_correlation_methods.R:12` (the pairs kernel). **Moved:** `omicsRouteSteps` `:1447-1453` is now **`:1710-1716`** (the RUV operation at `:1714`, was `:1451`;
  KIT-255 added 266 lines above it); `R/func_lipid_s4_correlation_methods.R:341-420` (the selection) is **changed** by KIT-255: the method now starts at `:344`, returns
  `omicsCorrelationFilterWithinAssay` (`R/utils_workflow_omics_declarations.R:1080-1215`) unless the replay marker holds its earlier token, and the earlier body, whose keep rule was `:388`, is at `:397`;
  `R/utils_workflow_analysis_record.R:693` (the replay-marker doc) is at `:726`. The keep rule of the new function is the same (`correlation >= threshold`, at least one pair: `:1115-1118`).
  The code is as the ticket derives it: nothing of RUV-III-C, rule (d), the cap, the DA kernels or the pairs kernel moved.
- **`sim-k188.R` at 0866982f**: every marked line was applied by exact text (the script that builds `sim-k275.R` asserted each old string); `generate` is `:164-309`, `CELLS` `:66-80`, `cellDef` `:62`,
  `ARM_MARKERS` `:137-139`, `maskOf` `:315-319`, `journeyL` `:613-`, the child and parent `:2046-2283` (not carried). `generator-diff.md` lists every hunk line by line.
- **G0** (`g0.txt`): the only changes on the ticket's G0 list since 75acd421 are KIT-255's filter (aa119308f) and the figure-axis hunks (7f17547f9) and two heatmap hunks; no change to RUV-III-C, rule (d), the cap, the DA moderation or the pairs
  and selection kernels beyond KIT-255's: G0 does not stop the ticket. KIT-255 with I01 is in the bound commit. **KIT-210 is not on staging** (`origin/kit-210-ruviiic-remedy` is not an ancestor): recorded; the deciding tiers and the other reported tiers
  bind this commit. **KIT-242 is not on staging**: tier V waits (section 5).

## 2. Requirement table: each ticket requirement and the script lines that meet it

`S` = `sim-k275.R`, `C` = `score-k275.R`, `F` = `fixture-k275.R`, `T` = `selftest/`, line numbers at the committed tip.

| Ticket requirement | Where met |
|---|---|
| 5.1 the generator copy with exactly the ticket's marks (seed, sig, cells, C, noIS, plant, ret, mask, W, markers) | S `:158` (sig), `:160` (seed), `:177` (noIS), `:225` (C), `:280-301` (plant, after `lbl_ids`), `:324` (ret), `:331-335` (mask), `:698` (W), `:130-133` (markers), `:51-78` (cellDef and the cells D3, D4, D6, D4M, X3, X4, X6); `generator-diff.md`, `generator-diff.patch`, `e0-diff.py` (E0) |
| 5.3 the planting block, own stream, session RNG state restored, V1g/p/s/t/V2 share draws, sparse form = the 25 smallest `u_feat` of Y | S `:280-301` (as the ticket's text; the loop variable is `s`) |
| 5.3 seeds 275M, planting 276M, bootstrap 277M | S `:160`, `:282`; C `:155-158` (`277000000 + 100 c + t`, t = 1..8 by tier in table order; ID 9; **G1 uses t = 9 too: the ticket gives none, see section 7**) |
| 4 versions V0, V1g, V1p, V1s, V1t, V2, VD (V0 with a removed from design and both assays), VE (built as VD) | S `makeVersion` `:965`, `dropSamples` `:954`, `runEmulation` `:1276` |
| 4 arms NO, LMk, LM2@20, LG1/LG2 through `ARM_MARKERS`, LA, GKk, GA, GSC, LX, PEL2, PEG2 | S `armSpec` `:909`, `requestK275` `:922`, `ARM_MARKERS` `:128-133` (the `earlier` token asserted by `armMarkers` in every journey), LX `lxRow` `:1189` |
| 4 tiers D, K, F, S, V, M, T, E, ID with their versions, arms, designs | S `TIER_PLAN` `:892`, `tierDesigns`, `runEmulation` (E), cell `D4M` (ID) |
| 6.1 A: C_s from the filter's own pairs table before and after RUV, P_pre/P_post at 0.5, 0.7, 0.9, 0.95, C_a, q_a, false fails pre and post (Y, and X at 0.5), kappa_c, PR, pi | S `journeyL` reads (K275-pairs `:676-678`, `:702-704`), `bestCorrelation` `:995`, `journeyRow` `:1107-1188`, `kappaControls` `:1028`, `pseudoReplicateRatio` `:1013`; pi's dot products `daNumbers` (`dot_g/p/t/two`), `plantedShares` `:981`, denominators `den_*`; C `:393-` (the cells) |
| 6.1 B: moderated SE = abs(logFC/t); unmoderated by `lm.fit` on the matrix and model the DA read (legacy: corrected matrix; declared: uncorrected with W; supplied controls left out); residual df | S `daNumbers`, `unmoderatedReading` `:1039`, `journeyRow` (the `mat`, `w` lines), W kept by K275-W |
| 6.1 C: any call, per-feature null rate (Y and X), paired against V0 and VD, V0's and VD's own rates printed beside | S `journeyRow` (`any_call_y/x`, `null_p05_y/x`); C `:366-386` (labelled), the reported cells in the descriptive loop |
| 6.1 H: pre-RUV emulation, H1 to H3, `refused_after_exclusion` counted as a result | S `runEmulation` `:1276`, `runArm` `:1260` (status rewrite when samples were excluded), `pe_*` columns; C `VALID` includes it |
| 7.1 margins, levels, Clopper-Pearson, paired percentile bootstrap (type 7, B = 10,000, shared index matrices), jackknife-normal limit beside, label rules, coverage rule, no tie-break | C `:44-83` (constants, label functions), `:153-180` (`cpUpper`, `cpLower`, `bootIndex`, `ratioStat`), `:312-409` (the labelled cells) |
| 7.2 falsifiers F1 to F4, stop rules | C `:319-332` (PC1, NC1), `:394-425` (G1 PC2/NC2, E5 and E1/E2 from the rows and the counted jobs' checks) |
| 7.3 outcome table, first row that holds | C `:436-` |
| 8.3 E0 (copy diff; requests build; defaults; markers' `earlier`; hashes) | `e0-diff.py`; S `armMarkers` per journey; S `checkE3` `:449` and the parent's E3; provenance hashes S `:1517`, `:1612` |
| 8.3 E1 (determinism, planting line, V1g/V1t/V1s/VD exact, X untouched) | S `checkVersionsE1` `:1307` (run in every tier D unit, also counted: F4) |
| 8.3 E2 (DA of X identical across V0 and V1g under NO) | S `checkVersionsE1` |
| 8.3 E3 (route defaults) | S `checkE3` `:449`; the parent runs it before any unit (a counted job stops on a mismatch) |
| 8.3 E4 (the scorer self-test) | `selftest/` (section 4) |
| 8.3 E5 (filter identity: fail list from the pairs table = the package filter's removed set) | S `journeyRow` (`e5_y`, `e5_x`, from `corr_removed` per assay by K275-pairs `:714-721`), `runArm` (check rows), C F4 |
| 8.3 E6 (SE reading on the G1 fixtures) | F (E6 rows: moderated SE, unmoderated against the kernel's sigma and stdev.unscaled, LX on the hand-made matrix) |
| 8.3 E7 (a recorded crash and a refused journey: GK3 at D3 on VD; neither a denominator) | S `runUnit` `:1400`, `:1405`; `runArm` |
| 8.3 E8 (metabolomics identity) | C `e8Compare` `:109` (blind and counted: the scorer compares tier ID with tier D at D4 to 1e-10; a difference is a finding) |
| 8.3 E9 (a data set run twice gives identical rows) | S `runParentMain` `:1465` (lifecycle: the first data set of each (tier, design), rows without their seconds) |
| 8.3 E10 (kappa_c, PR, pi on a constructed matrix) | F (E10 rows) |
| 8.3 E11 (emulation: pass 1 fail list = the NO journey's filter set; VE exact; refusal handling) | S `runEmulation` |
| 8.3 E12 (LG's positional controls) | S `checkE12` `:1238`, called by `runArm` for LG arms in every lifecycle journey and in counted data sets 1 to 20 |
| 8.4 G0 (reading) | `g0.txt`, section 1 |
| 8.4 G1 (X3, X4, X6, N0, 2,000 features, 20 data sets, NC2, PC2) | F |
| 8.4 L1 blind: only checks, timing by arm class, hashes, scorer parse and coverage, console lines; row files under `tempdir()`; failure prints class and message only | S `runParentMain` `:1432-1620` (lifecycle branch), C `--blind` |
| 8.4 G2 and the ceiling stop | section 3 and 4 of this record (the job plan is derived from L1's `timing-<job>.tsv`) |
| 10 running-ceiling rule | `runs.tsv` (columns used, request ceiling, still planned, sum, limit), section 4 |
| Smoke flag written to provenance and rows, refused by the scorer for an outcome line, refused in a counted request | S `:1438` (`SMOKE_counted_refusal`), S rows column `smoke`, provenance `smoke`; C `:143` (outcome mode refuses), `F` (the G1 table carries it; C refuses it); T case `smoke` |
| Lessons 1 to 3 (read tables as text; no `$` prefix match, also inside the quieting wrapper; progress files and provenance before any unit) | C `readText` `:88` (`colClasses = "character"`), T cases `tf`; S header `:18-28` and `quiet()` (K275-guard `:99`), journey fields renamed (K275-prefix), `[[` with exact names in every new read; S `appendTable`/`flushChecks` `:1222-1235`, the start provenance and `checks-<job>-start.tsv` before the first unit |
| Tier V written and marked as waiting | S `:1447-1449` (tier V stops unless `K275_TIER_V_BOUND` is a 40-digit commit; `all` excludes V without it) |

## 3. Job plan, every timeout re-derived from section 10 on the measured basis (Automatic 22.8 s, Manual 17.8 s)

Derived with Python from the tier plans in `TIER_PLAN` (journey counts per data set: manual = NO, LMk, LGk, LM2@20, PEL2 and the emulation's pass-1 NO; automatic = LA, GA, GKk, GSC, PEG2). Seconds per data set at 12 samples:
D 294.2, K 327.6, F 104.0, S 91.2, V 248.6, M 139.6, T 121.8, E 116.8, ID 104.0 (all equal the ticket's). Plan = seconds x data sets x 1.25 / 4 / 60; `timeout_minutes` = ceil(1.5 x plan); ceiling = timeout + 10.
Every figure equals the ticket's section 10 table.

| Job | Request | Plan (min) | `timeout_minutes` | Ceiling (job-min) | Expected (min) |
|---|---|---|---|---|---|
| D | 12 jobs: D3, D4, D6 x seeds 1-50, 51-100, 101-150, 151-200 | 76.61 | 115 | 125 | 68.2 (D3, D4), 74.3 (D6) |
| K | 6 jobs: 3 designs x seeds 1-50, 51-100 | 85.31 | 128 | 138 | 74.4, 81.3 |
| F+S | 6 jobs (tiers `F,S` in one call): 3 designs x seeds 1-50, 51-100 | 50.83 | 77 | 87 | 49.6, 53.7 |
| V (waits) | 6 jobs: 3 designs x seeds 1-50, 51-100 | 64.74 | 98 | 108 | 59.6, 64.8 |
| M | 1 job: D6 seeds 1-100 | 72.71 | 110 | 120 | 71.2 |
| ID | 1 job: D4 seeds 1-20 | 10.83 | 17 | 27 | 20.8 |
| T | 3 jobs: 3 designs x seeds 1-100 | 63.44 | 96 | 106 | 58.7, 63.8 |
| E | 3 jobs: 3 designs x seeds 1-100 | 60.83 | 92 | 102 | 56.8, 61.7 |
| L1 (V waits: tiers `all` without V; 3 data sets + E9's second run per design, scripts: one call over `D3,D4,D6`) | 13,641.6 s | 71.05 | **107** | 117 | about 60 |
| L1 (if V is bound with the others: all nine tiers) | 16,624.8 s | 86.59 | 130 | 140 | 77.7 |
| L1-V (tier V alone at its own commit) | 2,983.2 s | 15.54 | 24 | 34 | 24.6 |
| G1 | the fixtures, 20 data sets x 3 designs | | 20 | 30 | 15 |
| SC | scoring, B = 10,000 | | 30 | 40 | 20 |

Counted jobs: 38 (with V), 24,200 journeys; every plan is at most 85.31 minutes (G2's 150-minute split does not act). Counted expected 2,424.8 job-minutes, ceiling 4,269 (= 284.60 core-hours).
G2 re-plans every job from L1's `timing-<job>.tsv` (mean seconds per arm class and design) before the declaration.

## 4. Running-ceiling ledger (ticket 10 (2): used + this request's ceiling + everything still planned at most 350.73 core-hours before any request)

Core-hours = job-minutes x 4 / 60. Planned, everything but the smoke send (worst case: V waits, so L1 without V and L1-V both run): counted 4,269 + L1 117 + L1-V 34 + G1 30 + SC 40 + reserve 150 +
correction (largest recomputed job K 138 + rescoring 40 = 178) = 4,818 job-minutes = **321.20 core-hours** (the ticket's own total is 320.47 with L1 at 130 and no L1-V; 322.74 with it). The smoke declaration printed 321.21 (a 0.01 rounding slip); `runs.tsv` has 321.20.

| Request | Ceiling (core-h) | Used before | Still planned (all later requests) | Sum | Limit | Sent? |
|---|---|---|---|---|---|---|
| smoke-1 | 3.87 | 0.00 | 321.20 | 325.07 | 350.73 | yes (section 6) |
| L1 | 7.80 | smoke's actual | 321.20 minus L1's own 7.80 (and the smoke has gone) | written at the time | 350.73 | no |
| later requests | | | | written in `runs.tsv` before each is sent | 350.73 | no |

The rule binds at every later request (counted job, continuation, instrument rerun, correction, L1-V): the implementer of the next step writes the sum in `runs.tsv` and sends only if it is at most 350.73.

## 5. What waits for KIT-242 (tier V)

Tier V (missing values, `vd5`, reported) needs a commit holding KIT-242 Parts A and B (and binds after KIT-210 if that is on staging first; ticket 13). Written and marked: the tier plan, the mask construction (`maskedLog2`
with the limit of detection solved on V0's values; VD's mask is V0's minus a's columns), the job's refusal without `K275_TIER_V_BOUND` (a 40-digit commit). Nothing assumes KIT-242's code: the arms are the same requests;
what KIT-242 changes (what a legacy RUV arm returns when missing values leave too few complete controls) shows in the rows' `status` and `ruv_status_y`. At that commit: its own G0, L1-V (tier V's arms, 3 data sets per design and E9's second run; the same scripts with
`K275_TIER_V_BOUND=<commit>`), then its counted jobs. Its numbers never decide the outcome row.

## 6. The smoke send

Declared in `smoke/declaration.md` (what it must show, ten items). One request, `kit275-smoke-1e381cab` on `run-tests`, package commit `1e381cab29e7d13e68f059c02ed6c26040cd19e1`, ceiling 58 job-minutes = 3.87 core-hours, charged to L1's line.
Result: see the end of this section (filled in when the run has returned; if it shows a defect it is recorded here with its file and line and **not repaired or re-sent**).

**Result, run 45742592872 (package 1e381cab, scripts at 892 s of 48 min declared; about 25 job-minutes with set-up = 1.67 core-hours against the 3.87 ceiling). Raw: `smoke/result/`.** Nothing is a result; no repair was made and none was re-sent.

| Declared item | Observed |
|---|---|
| 1 all scripts parse | met (`parse: all scripts parse`) |
| 2 E0 | met (29 hunks, 0 unmarked) |
| 3 E4 | met: 19 cases, 1,859 predictions (not 1,860 as declared: the count is the generator's), 0 differ; scorer exit codes 1 for `smoke` and `g1_missing` as predicted |
| 4 G1 fixture | met: 4 readings written (NA), 6 checks (E6, E10), 0 failed |
| 5 job A (D,K,S,E,T at D3) | **not met**: 2 failed checks. (a) `L1_status` 49 of 50: journey `VD NO` ended `crashed` at step `da` ("Lipidomics DA formula cannot produce a model matrix: contrasts can be applied only to factors with 2 or more levels"): the pre-RUV filter at 0.5 removed enough samples of VD (5 samples) that a group had no values left; this is the package's behaviour on a genuine data set (ticket P04, KIT-255 I01), not a script fault, but L1 as written would flag it: the audit should decide whether such a journey is a result (`no_da`) in L1's status check. (b) **Defect: unit `E-D3-9101` crashed** with `invalid 'type' (list) of argument` inside `journeyAssayColumns` (the kept KIT-188 function, called from `journeyL`) on a pass-2 journey of tier E (a `PEL2`/`PEG2` journey on VE, `correlation = FALSE`); trace `smoke/result/3/progress/unit-trace-E-D3-9101.txt`. Only V0's PEL2 and PEG2 rows exist; V1g's pass 2 did not run. Checks E1, E2, E3, E5 (92 of 92), E7_crash, E7_refusal, E9, E11, E12 all passed; blind scorer exit 0 |
| 6 job B (D,ID at D4) | met: 6 checks, 0 failed; the blind scorer's E8 line compared journeys (see `smoke/result/4`) |
| 7 SMOKE_counted_refusal | met (exit 1 after 4 s with that message) |
| 8 no partial-match stop | met (no occurrence in any console, progress file or table) |
| 9 provenance and progress | met (`provenance-smoke-*.json`, `checks-*-start.tsv`, `progress/` present) |
| 10 timing | recorded: mean seconds per journey at D3: NO 18.7, LM 22.8, LG 21.6, GK 24.7, GA 24.2, LA 40.1, GSC 31.9, PEL2 22.9, PEG2 25.8 (4 workers busy; for information only; LA is above the 22.8 s basis, L1 decides) |

**Desk-review findings made after the send, before the result was read, NOT repaired (for the script audit):** (1) `score-k275.R` lines 279-280, 284 and 290 pass `st$jk_lo, st$jk_hi` positionally to `emit`, so they land in `label` and `jk_lo` (and line 290 also names `label`, an error "formal argument matched by multiple actual arguments" when a pi cell is reached): the E4 self-test did not reach this because the predictions cover the labelled cells only; the fix is to name them. (2) `sim-k275.R` K275-pairs block (`corr_removed`, about line 718): if the legacy RUV output drops an assay (`m_out` without it), `before` is NULL and the journey crashes at step `correlation`; use `m_in` for that assay. (3) The E4 predictions say 1,860 in the smoke declaration; the generator writes 1,859.

## 7. Drafts of the declaration files (not sent; `declarations/`)

`declarations/l1-draft.md` (the L1 request and its checks), `declarations/counted-draft.md` (arms, tiers, designs, R, seeds, endpoints, margins, levels, outcome table by reference to the accepted
ticket sections, the bound commit and the script hashes as placeholders, the job list). They are drafts: the bound commit, the hashes and G2's re-planned timeouts are filled by the step that sends them.

## 8. What the script audit should know (interpretations and additions; none changes what is measured)

1. **"Copy"**: the generator and the journeys are copied; five line ranges of KIT-188's file that this ticket does not use are removed and the child and parent mechanics are written new (generator-diff.md, tables A to C). Marks beyond the ticket's list:
   K275-guard, K275-prefix, K275-pairs, K275-REFUSALS, K275-E3 (table A of generator-diff.md); the first three are the lead's lessons or the ticket's 8.2; K275-REFUSALS and K275-E3 are mine and are listed for acceptance.
2. **E3**: reduced to the route defaults of 5.2 and the marker check; the ITSD import-pattern clause of KIT-188 would fail on staging (KIT-260 not merged) and is not part of this route (ITSD off).
3. **Boundaries** ("equal to the margin", "lower limit exactly 0") are compared at 1e-12, since 0.95, 0.05, 0.10 and 0.0125 are not binary numbers; E4 tests the rule on constructed intervals (`labels` case).
4. **G1's bootstrap stream**: PC2 and NC2 are two-sided at 0.01/3 on 20 data sets; the ticket names the percentile bootstrap for differences and ratios and gives stream indices t = 1 to 8 only. G1 uses t = 9 (`277000000 + 100 c + 9`), the one ID also has (ID has no bootstrap). If the audit prefers a t interval, it is a change in `C` `:402-411` only.
5. **The NO arm** is the journey with `ruv = NULL` (RUV skipped, k = 0); its "post-RUV" pairs table is its pre-RUV one (the filter reads the normalised state).
6. **E5** reads the package filter's removed set from the matrices before and after the step (samples with a value before and none after, per assay) rather than from `correlation_excluded_per_assay`; the two are the same set by KIT-255's code (`R/utils_workflow_omics_declarations.R:1080-1135`).
7. **E8** is evaluated by the scorer on the rows (tier ID against tier D at D4, same seed index; the cell `D4M` is D4 with the metabolomics class, so its draws are the same); L1 therefore needs tiers D and ID at D4 in one job (`all`).
8. **LX** uses the DA kernel pair without `technical_replicate_id` (an unblocked `lmFit`, as the DA does on a design without technical replicates).
9. **Pairing across tiers**: per-data-set scalars (`se2_mod_y`, `dot_g/p/t/two`, `den_*`, `tested_sig_y`) make the paired estimands (B, pi, C) computable between tier D's references and the other tiers' rows without per-feature files; a pair whose tested sets differ is dropped and counted in the cell's note.
10. **ID rows and E8's L1 coverage**: the blind scorer prints, for E8, only the counts of rows compared and equal and the maximum deviation.
11. **Not verified by any run**: every R script (parse included). The smoke send is the first parse.

## 9. Rule slips

None known. (The smoke declaration's 321.21 is a rounding slip, corrected in `runs.tsv`; recorded in section 4.)


## 10. After the script audit (AUDIT-KIT-275-scripts-20261126): repairs S01 to S16 and the coordinator's rulings

Applied exactly as the audit's section 9 writes them (marks K275-assaycols and K275-nogroup added to E0 and `generator-diff.md`; K275-REFUSALS moved to the end of its statement): S01 (list-safe control count), S02 (`refused_after_exclusion` only
when too few samples are left and the class says so), S03 (a lost group is `no_da`, class `filter_removed_group:<groups>`), S04 (L1 status accepts `ok`, a lost group, and the emulation's results after an exclusion), S05, S06 (named `emit` arguments at :279-280 and :284 as corrected; :290 was right),
S07 (guard stop is `instrument`, the scorer stops on instrument rows, rows outside the declared range and duplicate rows; G2 timing over completed journeys only with a `G2_coverage` check), S08 (case-insensitive guard), S09, S10 (E3 asserts the marker names and
`correlation_exclusion` at `within_assay`; `K275_REQUIRED_MARKERS` for tier V), S11 (E8 and tier E's pass-1 NO against tier D's NO in scoring mode), S12 (H2 and H3 cells), S13 (pi coverage and missing reference; the dropped-pair count in every paired cell's note),
S14 (tier V's mask clauses of E1), S15, S16 (34 self-test cases, 2,745 predictions: the twelve of the audit's list plus the existing ones; non-degenerate data with exact jackknife limits and bounded bootstrap limits). The declaration's count is 1,859 (send 1); `runs.tsv` has the smoke's 25 job-minutes = 1.67 core-hours and run id 45742592872.

**Rulings.**
- **R-1.** A journey in which the filter removes a whole group is a RESULT: `no_da`, class `filter_removed_group:<groups>`, accepted by L1's status check and counted in the denominators (7.1). The tier E statuses (`uncorrected`/`dropped` for PEL2; `refused_after_exclusion` for PEG2 only when too few samples are left and the class says so) are results; any other refusal after an exclusion is an instrument defect (E11). The product finding (the DA stops with a technical message when the filter removes a whole group; with three or more groups the other comparisons are lost too; pre-existing) is filed by the coordinator; the scripts do not work round it.
- **R-2 (coordinator's ruling, told to the owner).** The emulation arms are compared with VD at the SAME k: PEL2 with VD LM2, PEG2 with VD GK2, both from tier K at r up to 100 (`REF_ARM` in `score-k275.R`), not at k - 1, because in VE the anomalous sample is excluded before RUV, so no factor is spent on it (the k - 1 rule of section 4 exists because one factor is spent on a). This changes H3's rho and the emulation's pi. **Sentence for the ticket (sections 4 and 6.1), not edited here:** "The reference arm for the emulation arms PEL2 and PEG2 is VD at the same k (PEL2 against VD LM2, PEG2 against VD GK2, both in tier K), not k - 1: in VE the anomalous sample is excluded before RUV, so no factor is spent on it; H3's rho_mod and the emulation's pi are read against that reference." (Section 4's list of k - 1 pairs ends with this exception; 6.1 H3 gains "against VD at the same k".)
- **R-3.** A cross-job identity mismatch (tier E's pass-1 NO rows against tier D's NO rows, r up to 100) STOPS scoring as F4 does (`crossNO`); E8 stays a finding.
- **R-4.** G1's bootstrap stream is t = 10; the scorer emits no bootstrap cell for tier ID (E8 only).
- **Routed to KIT-242, not ours:** on the route where the legacy control method skips an assay (a group with fewer than two samples), the package's record says RUV succeeded with k = 2 for an assay that got no controls.
- Interpretation 9 (pairing by per-data-set scalars and a tested-set signature) is accepted with the dropped-pair count printed in every paired cell; in tier V pi's denominator runs over all 250 features while its numerator runs over the tested ones: a stated limit.

## 11. Smoke send 2 (run 45748396833, package dfd14ad2; scripts 88 + 333 + 355 s, about 23 job-minutes = 1.53 core-hours against 3.33; raw in `smoke/result2/`)

Nothing is a result; nothing was repaired; no third send.

| Declared item (smoke/declaration-2.md) | Observed |
|---|---|
| 1 parse, E0 with the new marks | met: 36 hunks, 0 unmarked (K275-assaycols, K275-nogroup, the moved K275-REFUSALS listed) |
| 2 self-test | met: 34 cases, 2,745 predictions, 0 differ |
| 3 unit E-D3-9101 completes | met: V1g PEL2 `uncorrected` (class `filter_removed_group:B`: the VE design lost group B), V1g PEG2 `refused_after_exclusion` (`ruv_k_above_limit`); E11 4 of 4 |
| 4 VD NO `no_da` accepted | met: `no_da`, class `filter_removed_group:B`, `L1_status` 34 of 34 |
| 5 planted guard journey ends `instrument`, the blind scorer stops on it | **MISSED**: the unit holding the planted guard (`D-D3-9101`, K275_FAULT_GUARD) died after 13.9 s with exit 1 and no result (`child process exit 1, no result`), before any row, check or trace was written, so no `GUARD_planted` check exists and the blind scorer saw 0 instrument rows (`instrument_stop` item passes, exit 0). The cause is not visible (no trace: the failure is outside the child's `tryCatch`, or the process was ended); the other units ran without a guard stop. Not repaired |
| 6 E3 lists `correlation_exclusion` at `within_assay` | met |
| 7 completed-journey means at D3 and D6, `G2_coverage` passes | met (8 of 8 and 5 of 5; tier D at D6 has NO, LM, LA, GA, GK) |
| 8 no other failed check | one: `unit_completed` for the guard unit (item 5) |

Measured seconds per completed journey (4 workers): D3: NO 15.7, LM 19.4, LG 17.4, LA 31.7, GA 20.1, GK 21.0, PEL2 16.3, PEG2 15.9; D6: NO 16.6, LM 21.7, LA 29.8, GA 22.5, GK 22.3.

## 12. G2's job plan re-derived on the times of smoke 2 (`budget-replan.py`; "measured" = D3 and D6 above; "scaled" = D4 as the midpoint of D3 and D6; LG, PEL2, PEG2 at D4 and D6 and GSC (send 1's D3 31.9 s) scaled by LM's time at that design)

Plan = seconds x data sets x 1.25 / 4 / 60; timeout = ceil(1.5 x plan); ceiling = timeout + 10; expected = seconds x data sets / 4 / 60 + 13.

| Job (tier, design: jobs) | s per data set | Plan | `timeout_minutes` | Ceiling | Expected |
|---|---|---|---|---|---|
| D D3, D4, D6: 4 each | 308.0, 315.0, 322.1 | 80.2, 82.0, 83.9 | 121, 124, 126 | 131, 134, 136 | 77.2, 78.6, 80.1 |
| K D3, D4, D6: 2 each | 328.1, 344.1, 360.1 | 85.4, 89.6, 93.8 | 129, 135, 141 | 139, 145, 151 | 81.4, 84.7, 88.0 |
| F+S: 2 each | 224.6, 233.4, 242.2 | 58.5, 60.8, 63.1 | 88, 92, 95 | 98, 102, 105 | 59.8, 61.6, 63.5 |
| V (waits): 2 each | 256.2, 263.0, 269.8 | 66.7, 68.5, 70.3 | 101, 103, 106 | 111, 113, 116 | 66.4, 67.8, 69.2 |
| M D6 (R 100) | 148.6 | 77.4 | 117 | 127 | 74.9 |
| ID D4 (R 20) | 110.4 | 11.5 | 18 | 28 | 22.2 |
| T D3, D4, D6 (R 100) | 125.3, 128.8, 132.4 | 65.3, 67.1, 68.9 | 98, 101, 104 | 108, 111, 114 | 65.2, 66.7, 68.2 |
| E D3, D4, D6 (R 100) | 95.8, 100.5, 105.2 | 49.9, 52.4, 54.8 | 75, 79, 83 | 85, 89, 93 | 52.9, 54.9, 56.8 |
| L1 (V waits; 3 data sets + E9 repeat per design) | 14,499 s | 75.5 | 114 | 124 | about 70 |
| L1-V | 3,156 s | 16.4 | 25 | 35 | about 28 |
| G1, SC, reserve, correction | | | 20, 30, 150 pool | 30, 40, 150, 191 (largest job 151 + 40) | 15, 20 |

Counted jobs: expected 2,690.1 job-min, ceiling 4,519 (301.27 core-hours). **Total (with L1, L1-V, G1, SC, reserve, correction): ceiling 5,089 job-minutes = 339.27 core-hours; expected 2,813 = 187.5 core-hours.** The limit is 350.73: **339.27 is under it** (11.46 core-hours of headroom); nothing is trimmed and R is unchanged. Every D-job and K-job timeout is at or above 120 minutes, so those requests set `timeout_minutes` explicitly. Every plan is at most 93.8 (G2's 150 does not act).
Measured: D3 and D6 for NO, LM, LA, GA, GK; D3 for LG, PEL2, PEG2. Scaled: D4 everything; LG, PEL2, PEG2 at D4 and D6; GSC (also from send 1). With send 1's times the auditor's figure was 386 to 424: LA, the largest class, was 40 s there and is 30 to 32 s here (the second send ran fewer simultaneous heavy units); L1 measures it for G2 and the stop at 350.73 stays in force. Not marked "awaiting the owner" because the total does not exceed 350.73 on these times.

## 13. Smoke send 3 (run 45751802889, package d3eb29d2; scripts 89 + 347 s, about 17.3 job-minutes = 1.15 core-hours against 1.87; raw in `smoke/result3/`)

Cause found by reading (not observed): the guard's `stop()` was raised from a GLOBAL calling handler, which R runs with every local handler disabled, so it reached no `tryCatch` and ended the child at top level (exit 1, nothing written). Repair: `guardWarning` as a local calling handler around each journey and each unit, the global handler kept as a net with an `options(error=)` trace writer; the planted fault now also exists at unit level. Details in `smoke/declaration-3.md`.

| Declared item | Observed |
|---|---|
| 1 parse; E0 | met (38 hunks, 0 unmarked) |
| 2 self-test | met: 37 cases, 2,751 predictions, 0 differ |
| 3 planted journey-level guard | met: `GUARD_planted` 1 of 1 passes (journey status `instrument`, class `partial_match_guard: ...`), trace `progress/unit-trace-D-D3-9101-journey.txt`, the unit went on to complete |
| 4 planted unit-level guard | met: unit `T-D3-9101` ended after 9 s with error `partial-match guard: partial match of 'abc' to 'abcdef'`, trace `progress/unit-trace-T-D3-9101.txt`, a failed `unit_completed`, a unit row of status instrument (the blind scorer counts 2 instrument rows, `instrument_stop` false) and exit 3 of the blind scorer; the job continued |
| 5 control unit `S-D3-9101` | completes; none of its checks fails |
| 6 no other failed check | **one**: `E9` fails for 2 of 3 comparisons (units D and T): the planted fault is in the first run only (journey-level in D, unit-level in T), so the repeat run differs from it by construction; the control pair (S) is equal. An artefact of the plants, not a defect; recorded as a miss of the declared wording |

Running ceiling (core-hours; `runs.tsv`): smoke 1 used 1.67; smoke 2 used 1.53; smoke 3 used 1.15 (ceiling 1.87): used 4.35; planned (accepted re-plan) 339.27; sum 343.62 against 350.73 (headroom 7.11). L1 re-measures every arm class; G2 re-plans from L1's own times and stops at 350.73.

## 14. Closure check (C01 to C05) applied

C01: every plant (the fault, the E7 refusal, both guard plants) sits on a data set E9 does not repeat (`units$r != first_rep`: r = 9102 in L1). E7's fault and refusal sit on the first unit that E9 does not repeat. C02: `score-k275.R --g1 <g1.tsv> <out>` judges PC2 and NC2 alone (`g1-judged.tsv`, exit 1 on F3); full scoring calls the same `judgeG1`; self-test cases `g1_pass`, `g1_fail`. C03: the header's guard lines carry `# K275-guard`; E0 exempts only hunks made of comments and blank lines. C04: the guard stop is re-raised in the analysis-record and E12 handlers. C05: the L1 draft is in line with the re-plan (114, 15, 139 job-minutes = 9.27 core-hours); `budget-replan.py` has no 0.8 factor (planned ceiling 340.27); `labelPi` applies the coverage rule when n = 0 (case `pi_nopair`).
**Continuation rule as the scorer applies it:** a counted job that did not complete is replaced by its continuation in the scoring input; the scorer stops on unit rows, on a false `unit_completed` and on duplicate rows; the failed job's directory stays in the records but is not scored.
Self-test now: 40 cases, 2,779 predictions.

## 15. G1 (run 45755206721, package 6075959d; not blind; run once; fixture 39 s, judging 1 s, about 11 job-minutes = 0.73 core-hours against 2.33; `g1/`)

F3 did not fire: all six controls pass (two-sided 0.01/3 bootstrap interval of the mean over 20 data sets, within 0.03 of the value).

| Design | Kind | Estimate | Interval | Value | Label |
|---|---|---|---|---|---|
| D3 | PC2 | 0.7753 | 0.7702 to 0.7809 | 0.7746 | passes |
| D3 | NC2 | 1.0017 | 0.9953 to 1.0080 | 1 | passes |
| D4 | PC2 | 0.8485 | 0.8429 to 0.8538 | 0.8452 | passes |
| D4 | NC2 | 1.0037 | 0.9987 to 1.0091 | 1 | passes |
| D6 | PC2 | 0.9055 | 0.9021 to 0.9090 | 0.9045 | passes |
| D6 | NC2 | 1.0003 | 0.9963 to 1.0043 | 1 | passes |

E6 (moderated SE reading, unmoderated reading against the kernel, LX on the hand-made matrix: deviations at most 3.2e-16, 0) and E10 (kappa_c 0.1785714286, PR 3.16227766, dot product 2): 12 of 12 pass.
