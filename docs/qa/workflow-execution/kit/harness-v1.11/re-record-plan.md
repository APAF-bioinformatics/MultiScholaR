# The single re-record after harness v1.11: plan for the laptop (DRAFT, 2026-11-19)

Nothing here has been run. The plan follows KIT-129 section 4: one build per cause, store against store with the build
before, unit references re-frozen with their mechanism, and journey and entry oracles at the end. It adds the items that
other tickets deferred to the same revision (`inventory.md`).

R stays at 4.6.1 (IgnorePkg). The proof library is `browser-smoke.y5ordS/library` or a library before it on `R_LIBS`,
with Bioconductor 3.23 `sva` 3.60.0, `genefilter` 1.94.0 and `annotate` 1.90.0.

## 1. Preconditions (each a user go)

1. The user approves `protocol-v1.11-draft.md` item 6 (items 1 to 5 were approved as U11).
2. The user gives the go for the laptop's machine time (section 5). The machine is quiet, on the balanced power profile
   (KIT-123 lost four measurements to the power-saver profile).
3. `kit-harness-v1.11` is reviewed. The branch is not merged. Its runner evidence is `runner-evidence.md` (runs 44461036735 and 44467185692).
4. The open decisions D1 to D4 in section 7 are answered. D2 decides whether P4 joins.

## 2. The builds and what each one adds

Each build is one `R CMD INSTALL` of a clean tree (`run_workflow_kit.sh install <previous library> <core dir>`), with
the source commit, build identity and exact-sum identity recorded in the bindings before any member uses it (protocol
1, Builds).

B1 to B3 do not exist as branches on top of `performance-v2` 791cc5f3. They are built on the laptop day by merging,
onto B0, the commits KIT-129's B1 rehearsal merged cleanly. Each merge is tree-checked with `git merge-tree` before it
is built.

| Build | Tree | Adds | Gates |
|---|---|---|---|
| B0 | `performance-v2` 791cc5f3 | baseline (KIT-052, 134, 137, 145 code, 146 part 1) | both off |
| B1 | B0 + `kit-128` ace1f4de | technical replicates as blocks; trx sums technical replicates | off |
| B2 | B1 + `kit-126-audit` 95f9d1d3 | classic: factors in the model (default unchanged) | off |
| B3a | B2 + the KIT-127 commit that restores the feature levels on the legacy path | lipid/metab: levels restored | off |
| B3b | B3a + the rest of `kit-127-audit` b3bd74c3 | lipid/metab: factors in the model (default unchanged) | off |
| B4 | `kit-harness-v1.11` tip (= `kit-129-integration` 039ce71d + this change set) | KIT-129 presets and wiring, the harness revision, the folder names | off |
| B4p | B4 + P1 (no `D66_debug.log`) + P3 (trx template paragraph) + P4 (`kit-135-integration-inputs`, if D2 = before) | report text, integration-input files | off |
| B4g | B4p + G2 (`cyclic_loess_minimum` on) + P2 (KIT-145 fixtures and manifest) | the refusal under 10 features; the four e2e feature-table fixtures | G2 on |
| B5 | B4g + G1 (`ruv_design_panel` on) | the default switch: Recommended pre-selected; technical replicates on designs that declare them | both on |

The B4 build carries two changes against B3b that are not scientific: the KIT-129 tab UI and the folder names. B3b and
B4 are compared store against store through the v1.11 readers, which take both layouts. An earlier build's store is
therefore compared by content, not by path.

## 3. Journeys per build (members: one at a time, nothing else on the machine)

Journey keys follow KIT-129 section 4:
- a, b: lipid screens;
- s, v: metab screens;
- c: trx screen;
- m: DIA_LIMPA covariates;
- r: DIA_LIMPA correction;
- l: LFQ;
- t: TMT;
- f: LFQ with RUV;
- q: DIA full.

They are bound under new single-journey ticket keys in `bindings-generator.R` (for example `kit-journey-129[a-z]`).
Each binding carries its build, its corpus, its oracles (7.1) and, from B4 on, `ruv_design`.

| Step | Build | Members | Judged against | Expected |
|---|---|---|---|---|
| 0 | B0 | a b s v c m r l t f q | (baseline stores; kit_940's are reused where `R/` is unchanged for the route) | |
| 1 | B1 | a b s v c; m l t as controls | B0, store against store | lipid/metab: DA on; trx: design on; proteomics identical |
| 2 | B2 | q f; l t as controls | B1 | classic: DA on; states before it identical |
| 3 | B3a, B3b | a b s v each | B2, then B3a | B3a: corrected state by a constant per feature, DA following; B3b: DA only |
| 4a | B4 | lifecycle pair first (`kit:lifecycle:...:complete`, `:fail`, through the new `ruv_design` input), then all eleven bound to Legacy (no `ruv_design`, or `legacy`/`legacy_correction`) | B1 (classic, trx, DIA_LIMPA), B3a (lipid, metab) | no scientific difference; UI and record fields only; store paths now `.n-`/`.we`/`.wx` (content compared) |
| 4b | B4p | a b s v c (report and integration); m l as controls | B4 | `D66_debug.log` gone; `integration/<omic>/` present (P4); trx report paragraph; nothing else |
| 4c | B4g | a b s v | B4p | identical (the kit corpora have 14+ features per assay; KIT-145 found no stored reference outside the four e2e lanes) |
| 5 | B5 | all eleven on `ruv_design = "default"`; then a b c m r f again on `legacy` | B4g | differences from the RUV step on, each route; the Legacy members equal to B4g |

l and t: at B5 they bind `legacy` with `ruv_mode = skip`. That reproduces today's route, which "verified unchanged"
requires (D1).

**Attribution.** Each step's store-against-store table assigns every differing component, product or file to that
step's mechanism. Tools: `kit-118/compare-content.R`, `compare-da.R` and `compare-records.R`, run as in KIT-118 and
KIT-123. A difference no mechanism explains stops the sequence before the re-record (section 6).

## 4. The re-record, on B5 (one commit per reference family, each with its attribution comment)

| Order | Family | Recorder | Notes |
|---|---|---|---|
| 1 | In-process oracles (7.1) for m, r, a, b, s, v, c, f | the route's `oracle-record*.R` (KIT-053, KIT-098 `oracle-record.R`, `kit-118/oracle-record-limpa.R`) with `R_LIBS=<B5 library>` | then the B5 journeys verified exact against them |
| 2 | Bindings regenerated once | `Rscript docs/qa/workflow-execution/kit/bindings-generator.R` | protocol hash (after `protocol.md` section 13), harness hashes, step-map hashes, `core_pin_lifted`, B-builds, the 129 keys |
| 3 | GUI-entry oracles, ten | `MULTISCHOLAR_GUI_ENTRY_FREEZE=1` on `test-gui-entry-{proteomics,lipidomics,metabolomics,transcriptomics}.R`, then a verify run without it | causes: H6 keys, KIT-128, KIT-129 UI, G1, P1, P2+G2, P4 |
| 4 | W1 captures and frozen dumps | `MULTISCHOLAR_KIT1xx_W1_DUMP` in `test-omic-{norm,summary,da,enrich}-module.R`, `test-prot-norm-core.R`, `test-trx-norm-module.R`, `test-workflow-saved-da.R` | the 335/159/49 folder keys (H6), KIT-128's 14 + 12, `omicSummaryW1Frozen` (P4), the stale S4 digest |
| 5 | Unit digests and fixtures | per file, in each test's own freeze procedure | `test-omics-da-core.R`, `test-trx-kernels.R`, `test-trx-publication.R`, `test-omics-summary-core.R`, `kit073-omics-declarations-frozen.rds`, `kit094-characterisation.rds`, MCI-027.2, MCI-028.2, the integration-input tests, `test-workflow-figure-paths.R:747-751` (H6: 64/125 become about 77/138) |
| 6 | Module CI | `tests/testdata/mci/refs`, `helper-module-ci-oracles.R` | lipid, metab, proteomics, trx module-CI files |
| 7 | e2e lanes and report-template doubles | the seven lanes (`test-e2e-*`) | expected counts unchanged by G2 (KIT-145, 6.); P3 moves the trx report double |
| 8 | Journey records | the B5 members re-run against the new oracles | exact |

## 5. Machine time and load (from recorded times; KIT-129 section 4, KIT-118)

| Part | Wall | Peak memory on the laptop |
|---|---|---|
| Nine builds and installs | about 2 h (seven took 1.5 h) | install only |
| Steps 0 to 4c journeys, serial (screen journeys 140 to 320 s; DIA full 2,130 s) | about 5.5 h | screen journeys 1 to 2.5 GiB; DIA full q about 5.7 GiB (the DIA route's peak, KIT-053); every member under the 16 GiB slice cap, 12 GiB `MemAvailable` at start |
| Step 5 journeys and oracle records | about 2.5 h | as above |
| Re-record (entry oracles record and verify, digests, module CI, e2e) | about 4.5 h (`test-gui-entry-proteomics.R` 1,465 s; `test-prot-norm-core.R` 1,790 s; e2e lanes 1,490 s) | the suite six files at a time |
| Exit: full suite in both modes, `R CMD check`, seven e2e lanes | about 3 h | six at a time |
| **Total** | **about 17.5 h of machine time; the re-record and exit from B5 on about 10 h** | one heavy process at a time during journeys |

R0 (KIT-124 C6 arms) is not counted: it ran on GitHub runners under the user's platform decision.

## 6. Stop rules

- One correction pass for the whole sequence (KIT-129). A failed comparison is recorded and reported, never retried
  under another name.
- Before step 4a, the lifecycle pair must show completion and forced failure through the new `ruv_design` input. If it
  does not, stop.
- If any step has a difference its mechanism does not explain, stop before section 4 and report it with the table.
- A member over its ceiling, an OOM or a guard stop is recorded as it stands (protocol 1). Only a demonstrated
  instrumentation fault allows one rerun.
- If a step's reference does not exist, stop and ask; do not substitute (B3a: see D8).
- Nothing is merged into `performance-v2` until the exit is clean and the user decides.

## 7. Open decisions (each answerable in one line)

| | Question | Recommendation |
|---|---|---|
| D1 | At B5, do LFQ/TMT screens (l, t) bind `legacy` + skip (today's route, identical) or `none` (the new option set's No RUV, which applies trend/robust: M30)? | `legacy` + skip, so "verified unchanged" holds |
| D2 | KIT-135/136 (integration inputs for lipid, metab, trx): merge before B5 so their oracle rows join this re-record, or after it? | before (B4p); otherwise the five entry oracles are re-frozen twice |
| D3 | U6: keep the KIT-040 notebook as the Legacy check on a new corpus variant without technical replicates? Then the direct `RUVSeq::RUVs` + edgeR check of item 5 gets written before step 5. | keep (as recommended in KIT-129) |
| D4 | `.dia-native-import-`, `.dia-native-legacy-` and `.dia-gui-adoption-` keep their long names (KIT-052 named only four). Shorten them as well in this revision? | not now; record the limit in the Windows notes |
| D5 | The GUI-entry recorder's mask (`helper-gui-entry.R:60`) hides runs of 10 or more digits (KIT-134). Fix it in this re-record (every entry oracle moves anyway)? | yes, as its own commit before step 3 of section 4, if wanted; otherwise leave it |
| D6 | Record the platform-sensitive unit references on a runner (section 8)? | B: runner for unit references, laptop for journeys |
| D8 | B3a needs a tree with KIT-127's level repair and without its factors. Read here, `kit-127-audit` has one kernel commit (abf30703: levels and factors together, on f0a09a30's tests). Split it into two builds by a scratch commit, or judge B3a and B3b as one step? | one step (B3) with both attributions in one table, unless the level repair is found as its own commit |
| D7 | Push CI: with `MULTISCHOLAR_TRX_DIGEST_RULE=values` workflow-wide, the impacted module-CI jobs note the trx digests without judging values (only the full gate runs `compare-results.R`). Also add the compare step to the impacted job? | yes, a two-line copy of the step |

## 8. What could be recorded on GitHub runners without breaking the protocol

The protocol binds memory and counted members to the laptop: a managed slice, systemd, 50 ms sampling. Exactness
references are not memory claims.

- **Stay on the laptop:**
  - every journey member and lifecycle pair (sections 3, 4.1, 4.8);
  - the in-process oracles they are verified against;
  - the bindings.
- **Could move to runners (`run-tests`, reference BLAS, R 4.6.1, same snapshot library):**
  - the full suite in `load_all` mode, `R CMD check`;
  - the seven e2e lanes (KIT-134 ran them on runners);
  - the GUI-entry lanes and every unit-level freeze (section 4, rows 3 to 7).
  - This saves about 7 of the 10 hours from B5 on.

**Platform-sensitive references.** Laptop-recorded hashes differ from runner results in the last digits. Causes: R
built with `-march=native` (fused multiply-add, vectorised sums) against generic binaries, and OpenBLAS kernel choice by
CPU where it is used. Which references are sensitive:

| Reference | Form | Sensitive? |
|---|---|---|
| GUI-entry oracles | values rounded to 10 significant digits, then sha256; long hex/digit runs masked | mostly stable (KIT-134: prot, lipid and trx lanes identical on runners); metab cyclic loess on 3 features was not (fixed by P2) |
| trx digests (nine), `test-trx-kernels.R`, `test-trx-publication.R` | exact sha256 of serialised doubles | yes (KIT-134 B-01; item 4 handles them off the laptop) |
| Frozen digests in the omic module tests, `test-prot-norm-core.R`, `test-omics-da-core.R`, `test-omics-summary-core.R`, KIT-073 dump, `kit094-characterisation.rds`, MCI refs | exact sha256 of doubles | yes, wherever doubles come from a fit (limma, RUV, loess) |
| Folder-name keys, file lists, UI digests | text | no |
| Journey and in-process oracles | exact fingerprints | yes; laptop-only by protocol |

**Proposal for the user (D6):**
- (A) Keep everything on the laptop, as before. The runners then go on comparing under the value rule, with digests as
  notes.
- (B) Record the unit-level references (section 4, rows 3 to 7) on a runner image with reference BLAS, so that CI
  verifies them exactly thereafter. The laptop's exact gate for those files then becomes the value rule, the mirror of
  item 4. That needs one more protocol sentence before step 3 of section 4. Journeys and their oracles stay on the
  laptop.

Recommendation: (B). The runners agree among themselves under reference BLAS (KIT-134), and CI is where these references
are checked most often.
