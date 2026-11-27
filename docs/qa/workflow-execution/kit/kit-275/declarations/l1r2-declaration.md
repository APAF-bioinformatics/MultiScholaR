# KIT-275 L1-r2 (the blind lifecycle, repeated once) and the conditional G1 re-run: declaration, committed alone, before the request (ticket 8.4, 8.5; VAL-401 (3), VAL-402)

## 1. Nature and purpose

L1-r2 is **blind** (VAL-11): it writes only checks, timings by arm class and design, hashes, the scorer's parse success and cell coverage, and the self-test; no failure rate, correlation, standard-error ratio, any call or rate. It proves the instrument of revision 2.1 through completion and failure at the **single bound commit**, on **its own seeds r = 9701, 9702, 9703** in each design (disjoint from every other range: `pilot/seeds-disjoint.py`, 0 shared), with the second run of the first data set (9701) of every (tier, design) for E9. It repeats r1's L1 scope (tiers D, K, F, S, T, E, M at D6, ID at D4) **and covers tier V's arms** (the vd5 mask, E1's mask clauses, the tier V acceptance rule, `K275_TIER_V_BOUND`, `K275_REQUIRED_MARKERS`), so **L1-V is dropped** (audit RS1). It exercises every changed path: the plant's "add" branch in `generate()`, the version builders, tier D's 17 rows per data set (V1w's NO journey), `journeyRow`'s version list and share key (C-2), `cy_genuine`, E1's exact identities and the reported mean and SD readings, the provenance's tree status (C-3), the scorer's labels, `EXPECT_ROWS`, the r offset and the K2-07 case, and the cross-implementation identity check. **It runs once.** The pilot (run 45923367873) accepted V1a (50 of 50 in every design); V1b is not a counted arm.

## 2. Binding

- **Bound commit `bfbcb34dfc888637d879e89029a9ee86505b476a`** (origin/v0.9-staging; merge `ddbdaaff134aa56f6c251f68255975ed9fc96a00` into kit-275; not re-merged). The request's `package_commit` is this declaration's own commit. **Tree check, re-done at this commit:** `git diff --stat bfbcb34dfc88 <this commit> -- . ':!docs/qa/workflow-execution/kit/kit-275'` is **empty**; the provenance records the same diff on the runner as `package_tree_vs_bound` with **git's exit status beside it** (`git_status` 0 and an empty `diff_stat` = equal; an empty string with another status is a failed git, C-3). G0 was re-read at this commit (`g0-r2.txt`, with the two inert hunks added).
- **Scripts (SHA-256, first 16 hex, at this commit; the provenance records the full hashes):** `sim-k275.R` d972252dee0bb81a, `score-k275.R` ee99cbb40e5c621d, `fixture-k275.R` 407949a9a6c0f90f (unchanged against r1), `crosscheck-k275.R` 1a0a515cc4c31f6a, `e0-r2.py` a2cef82db3b6a691, `selftest/make-selftest.py` a8653847bc7391ba, `selftest/check-selftest.py` 40a2a88c6c58e300 (unchanged), `selftest/run-selftest.R` 253daedb83c17777, `pilot/seeds-disjoint.py` 19757b91f0562d76. `e0-r2.txt` lists every hunk of the changed scripts against r1's audited versions, each tagged with the ticket section or finding it implements, **0 untagged**.
- **What r2.1 changed in the counted scripts (all tagged in `e0-r2.txt`):** 5.3's "add" branch inside `generate()` exactly as written (V1a: kind "add", s = 4, scope "all"; V1s scope "spared"; V1t scope "tenth"; V2 two samples; V1w and V1p kind "replace", c = 0 and 0.5); r1's V1g is **V1w**, a reported arm with its NO journey only (tier D: 17 rows per data set); every deciding reading uses **V1a** and V1w is never read as deciding (8.1's label check: 8.1 diff of every version label; the self-test case `v1w_only`); `kappa_c`'s version list and the share-key switch include V1a (C-2); `cy_genuine` in every row (the fixed genuine sample's best correlation, R2.4); E1 keeps ONLY exact identities (V1a's deviation = 4 z to 1e-12 on every planted feature and labelled-control row, z redrawn from the planting stream; V1w's = r1's formula with c = 0) and the mean (4 SE) and SD (4 SE = 4/sqrt(2 (n - 1))) readings are the reported checks `R1_dev_mean` and `R1_dev_sd`, outside E1 and F4 (C-1); the scorer's PC1 is the single floor (the one-sided upper limit at 0.01/3 below 0.80, i.e. fewer than 144 of 200, is F1) with the genuine sample's count beside it, no V1a row is "not determined (coverage)" and never a pass or fail, `EXPECT_ROWS` D = 17, the data-set index offset 2000 (r = 2001 to 2200); the provenance field of C-3. The diagnostic and pilot scripts assert r1's `sim-k275.R` hash and are history: they are not re-run.
- **Cross-implementation identity (audit Q1 (b)).** `crosscheck-k275.R`, run by the runner self-test, builds at r = 9601 in each design the counted `makeVersion(design, 9601, "V1a", "none", v0)` and the pilot's `pilotPlant()` V1a and requires `identical()` on y, on the values passed to the pipeline and on `planted_delta`: **the pilot's control is the counted control**. A mismatch stops (exit 1).

## 3. G1: do its identity conditions hold? No: the re-run is part of this request

G1 (run 45755206721) stands only if `fixture-k275.R` and the scorer's G1 path (`judgeG1`) are byte-identical to r1's **and** the bound package tree equals L1's (f3c1358b, equivalently 577a6f7b, over `R/`, `src/`, `DESCRIPTION`, `NAMESPACE`, `inst/`, `packaging/`). `fixture-k275.R` is byte-identical (`e0-r2.txt`); `judgeG1` is byte-identical (no hunk of `score-k275.R` falls in its lines; the scorer's only other change on the `--g1` path is the added `R_BASE` constant, unused there). **The tree condition fails:** 62 files under `R/` and `inst/` differ between 577a6f7b and bfbcb34d (`git diff --stat 577a6f7b bfbcb34dfc88 -- R/ src/ DESCRIPTION NAMESPACE inst/ packaging/`; `g0-r2.txt`), and it is not shown to hold. **So G1 is re-run once in this request** (35 job-minutes = 2.33 core-hours at r1's G1 ceiling): `fixture-k275.R D3,D4,D6 1 20 k275-g1` (limit 20) and `score-k275.R --g1 k275-g1/g1.tsv k275-g1` (limit 5); not blind; its readings are a result (reported); F3 (any of the six controls outside its tolerance) stops the ticket, the owner is told; the scripts of the request run on regardless (the request is one), and a G1 stop makes this L1-r2's checks moot only for the counted step that follows.

## 4. The request (one), the budget and the running line

| # | Script | Arguments and environment | `timeout_minutes` |
|---|---|---|---|
| 1 | `selftest/run-selftest.R` | out `k275-l1r2-selftest`, `K275_SELFTEST_OUT=k275-l1r2-selftest` | 15 |
| 2 | `fixture-k275.R` | `D3,D4,D6 1 20 k275-g1` | 20 |
| 3 | `score-k275.R` | `--g1 k275-g1/g1.tsv k275-g1` | 5 |
| 4 | `sim-k275.R` | `all D3,D4,D6 9701 9703 k275-l1r2`; `K275_LIFECYCLE=1 K275_FAULT=1 K275_JOB=l1r2 K275_DECLARATION=docs/qa/workflow-execution/kit/kit-275/declarations/l1r2-declaration.md K275_TIER_V_BOUND=bfbcb34dfc888637d879e89029a9ee86505b476a K275_REQUIRED_MARKERS=normalisation_fallback,itsd_outcome,ruv_fallback OMP_NUM_THREADS=1`; no smoke flag, no planted guard | **168** |

Request fields as the counted requests: `blas` `reference`, `os` `ubuntu-22.04`, `extra_packages` `["sva"]`, `tests` `[]`, `check` `false`, `tools_ref` `""`, `base_commit` `none`, `fetch_refs` `[]`; upload patterns `\.(tsv|txt|json)$` (the self-test, G1, L1-r2). Every `timeout_minutes` is explicit and at most 280.
- **L1-r2's timeout by the ticket's rule (K3-01).** Timeout = ceil(1.5 x compute), compute from L1's measured seconds per data set (`budget-replan-l1.py`'s arithmetic, 3 data sets plus E9's repeat in every tier and design, tier V included, V1w's NO journey included): compute 111.72 minutes, **timeout 168**, ceiling 168 + self-test 15 + set-up 10 = **193 job-minutes = 12.87 core-hours** (the audit's figure, RS1).
- **G1's row:** 35 job-minutes = 2.33 core-hours. **The request's own ceiling:** 15 + 20 + 5 + 168 + one set-up of 10 = 218 job-minutes = 14.53 core-hours, inside the two rows' sum 228 = 15.20.
- **Running line (core-hours; the ceiling of VAL-402 is 527.50).** Used: 112.7 before the pilot (coordinator's figure) + the pilot's about 2.3 = about 115.0. This request: 15.20 (the rows' sum, conservative). Still planned after it: the remaining plan of R2.7 (411.60 with the pilot at its ceiling 5.00; the pilot is now spent, so 406.60) with L1-r2 repriced to 193 (12.87 against 9.27: +3.60), L1-V dropped (-2.73), G1 added (+2.33) = 409.80 less this request's two rows (12.87 + 2.33 = 15.20) = 394.60. **115.0 + 15.20 + 394.60 = 524.80, at most 527.50** (527.50 exactly if the pilot is counted at its ceiling of 5.00 instead of its actual): met. The coordinator reports to the owner when the spend passes 406.90 (VAL-401 (4)).
- **P-3:** the pilot's and the diagnostic's timings are not used for G2, any timeout or any ceiling; G2 uses L1-r2's timings only.

## 5. Predictions for every mechanical check (nothing about any outcome is predicted)

Counts are derived by counting (units: 23 (tier, design) cells x 3 data sets + 23 repeats = **92 units**; journeys and rows as `TIER_PLAN`: D 17, K 17, F 5, S 4, V 12, M 7, T 6, E 6, ID 5 rows per data set), as in r1's L1 (80 units, 697 rows, 695 `L1_status`). A count below the prediction is a finding, not repaired.

| Check | Predicted |
|---|---|
| Self-test step | the bound hashes as declared; the five scripts parse; `crosscheck`: the counted V1a equals the pilot's at D3, D4 and D6 (exit 0); E0: **41 hunks, 0 unmarked**; `e0-r2`: **0 untagged or unexpected hunks**; `seeds-disjoint`: **0 shared** (pilot, L1-r2, counted); scorer self-test **48 cases, 3,519 predictions, 0 differ** (77 row columns) |
| G1 (run again) | the six controls pass (PC2 0.775, 0.845, 0.905 and NC2 1 within 0.03); E6 and E10 12 of 12 (r1's readings; the fixture, the DA kernels and `judgeG1` are unchanged, the 62 changed files are not on the kernels' path, G0) |
| Units | 92 units, none with an error; no `unit_completed` false row; no instrument row; no guard stop |
| Rows | **853 rows** (852 + E7's refusal row) |
| E1 | **120 of 120** (tier D units, 12 x 8 = 96; tier V units, 12 x 2 mask clauses = 24) |
| E2 | 12 of 12 |
| E3 | 2 of 2 (lipid, metab), with `correlation_exclusion` at `within_assay` and the three tier V markers present |
| E5 | **1,558 of 1,558** (779 journeys that reach the filter x 2 assays: 852 rows less 24 LX rows and 48 emulation pass-2 rows, less the one E7-fault journey) |
| E7 | crash 1 of 1; refusal 1 of 1 |
| E9 | **23 of 23** (every repeat equal) |
| E11 | 48 of 48 |
| E12 | 48 of 48 |
| `L1_status` | **851 of 851** (853 rows less the fault journey and the refusal row, which have none; LX rows included) |
| `G2_coverage` | 27 of 27 |
| Reported readings (not E1, not F4) | `R1_dev_mean` 12 of 12 and `R1_dev_sd` 12 of 12 TRUE expected (chance of a failure about 6e-5 each); a FALSE is reported, not a stop |
| Blind scorer | parse; **853 rows**; every arm class covered at every design; **E8: 15 compared, 15 equal, max deviation 0**; 0 instrument rows; exit 0 |
| Provenance (each job file) | `lifecycle` true, `smoke` false, `final` true, `package_tree_vs_bound`: `git_status` 0, `diff_stat` empty, `bound_commit` as above, the hashes of this section |
| Written | `checks-l1r2.tsv`, `timing-l1r2.tsv`, `timing-units-l1r2.tsv`, `hashes-l1r2.tsv`, `scorer-l1r2.tsv`, `provenance-l1r2.json`, `console` lines of failed checks only |

## 6. Decision and hand-back

G2 re-plans every counted job from L1-r2's own completed-journey means per arm class and design (8.4; the 1.5 x rule; a job planned above 150 minutes is split at the same R; section 10's total is recomputed; **if it exceeds the owner's ceiling the matter returns to the owner**). A failed check is recorded and handed back without repair; a script or mechanics failure with no result is an instrument defect (VAL-54: four conditions, at most three further sends, `runs.tsv`, the ceiling unchanged); a result keeps one correction. F3 on the re-run G1 stops the ticket. The counted declaration follows (committed alone: the pilot's choice V1a and its expectation as R2.3 records them, the single F1 floor, the bound commit, the running ceiling). The coordinator: this declaration's commit and the request file (committed locally on `run-tests`, not pushed); the script audit reads everything before the send.
