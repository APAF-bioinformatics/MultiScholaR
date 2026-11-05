# KIT-093 audit (audit_only), 2026-11-06

**Mode:** `audit_only` (skill `audit-evidence-contracts`). No ticket, code or doc was edited and no test, install or journey was run. The only write is this file.

**Target:** `tickets/kit/KIT-093.md` revision 1.0, sha256 `3db64c38…c91d999` (the tickets directory is gitignored, `.gitignore:81`, so the ticket has no commit identity). Repository HEAD `bbe1fae9` (branch `performance-v2`).

**Authorities read:**
- `AGENTS.md` and `conventions/R.md`, `R-shiny.md`, `R-golem.md`, `R-bioinformatics.md`.
- `tickets/kit/00-KIT-PLAN.md` (plan 3.8, rows 231-238 and decision record line 305) and the approved plan `~/.claude/plans/write-a-series-of-atomic-dragon.md`.
- The neighbour tickets KIT-040 (sha256 `6ac5ea0e…99f53`) and KIT-094 to 099.
- The superseded TRX-004/005 and the TRX series audit, `TRX-SERIES-AUDIT-20261007.json`.
- The notebook `dev/distributable/plans/transcriptomics_workflow.rmd` (793 lines by `awk`, no trailing newline; sha256 `44d5f643…0da652a`), the runbook and its audit JSON.
- `docs/qa/workflow-execution/kit/protocol.md`.

**Mechanical validity:** `next-openai-ticket.sh KIT-093` parses the ticket (frontmatter explicit, `needs_planning: false`, three acceptance items). The index entry matches (dependency KIT-040, 8h).

**Environment facts checked (read-only):**
- R 4.6.1, Bioconductor 3.23.
- Installed: edgeR 4.10.1, limma 3.68.4, ruv 0.9.7.1, Biobase 2.72.0, Glimma 2.15.2, statmod 1.5.2.
- **Not installed:** RUVSeq, EDASeq, and EDASeq's dependencies ShortRead, Rsamtools, aroma.light and biomaRt.
- BLAS/LAPACK: reference `/usr/lib/libblas.so.3.12.0` / `liblapack.so.3.12.0`.

## Finding register

Severity: **blocking** = the ticket cannot be executed or closed as written. **should-fix** = it can run, but the claim or closure would be wrong, vacuous or ambiguous. **minor** = clarity or record.

### KIT093-F01: blocking, dimension 1 (completeness) and 4 (references): RUVSeq and EDASeq are not installed, and the ticket has no dependency-acquisition step or decision

- **Evidence:**
  - `KIT-093.md:35` adds the packages to `Suggests` only.
  - `packageVersion()`: RUVSeq and EDASeq are absent, and so are EDASeq's dependencies ShortRead, Rsamtools, aroma.light and biomaRt.
  - The notebook needs them: `rmd:53` `p_load(RUVSeq)`, `:406` `newSeqExpressionSet` (EDASeq), `:410` `makeGroups`, `:413` `RUVs`.
  - R is pinned at 4.6.1 because a library change broke every stored oracle (memory `r-version-pinned`). A Bioconductor install that touches shared dependencies risks the lipid/metab/DIA oracles.
  - KIT-040's reference render needs the same packages and does not say how they arrive (`p_load` would try to install them at render time).
- **Impact:**
  - Neither the KIT-040 reference nor the KIT-093 RUVs and edgeR-with-W kernels can run.
  - Tests guarded by `requireOptionalPackage()` or `skip_if_not_installed()` would skip, so the parity acceptance would close vacuously (see F11).
- **Repair:**
  - A user decision recorded in the ticket covering: how the packages are installed (for example `BiocManager::install(c("RUVSeq","EDASeq"), update = FALSE)`); into which library; the resulting versions; and a check that no already-installed package version changed.
  - Add it as a precondition to KIT-040 and KIT-093.
- **Closure evidence:**
  - `installed.packages()` diff before and after: new packages only, no version changes.
  - Recorded versions of RUVSeq and EDASeq.
  - One existing oracle replay (for example the lipid frozen dump) unchanged after the install.
- **Blocked state:** decision, then execution.

### KIT093-F02: should-fix, dimension 4 (references) and 3 (practice): `ruv` is already a hard Import; moving it to Suggests is wrong, and Biobase/BiocGenerics are not named

- **Evidence:**
  - `KIT-093.md:35` says "edgeR, RUVSeq, EDASeq and ruv enter `Suggests`".
  - But `DESCRIPTION:522` lists `ruv` under Imports, and `NAMESPACE:816` has `importFrom(ruv,ruv_cancorplot)`, used unqualified in `R/func_prot_s4_norm_methods.R:179` and `R/func_pept_s4_norm_methods.R:758,831`. limma is also in Imports.
  - The notebook also calls `pData` (Biobase, already in Suggests, `DESCRIPTION:409`) and `counts()` (BiocGenerics/EDASeq).
  - edgeR is absent from DESCRIPTION, so adding it is correct.
- **Impact:** a literal implementation either breaks proteomics RUV and R CMD check, or contradicts "No existing class, generic or public function changed". `pData` would be called without its optional-package guard.
- **Repair:** state that ruv and limma stay in Imports, called `ruv::` and `limma::` directly with no optional guard; edgeR, RUVSeq and EDASeq go into Suggests; Biobase (already in Suggests) and the `counts` generic are guarded and qualified too.
- **Closure evidence:** DESCRIPTION diff adds only edgeR, RUVSeq and EDASeq to Suggests; `NAMESPACE:816` unchanged; R CMD check at baseline.
- **Blocked state:** implementation.

### KIT093-F03: should-fix, dimension 2 (science): the filter pre-pass is misdescribed, and the zero-variance filter and the two TMM passes are missing

- **Evidence:**
  - `KIT-093.md:28` says "`filterByExpr` with library sizes as a declared global pre-pass".
  - In edgeR 4.10.1, `filterByExpr.DGEList` uses `lib.size <- y$samples$lib.size * y$samples$norm.factors`. The notebook calls it on `dge_raw` after `calcNormFactors(method="TMM")` on the *unfiltered* matrix (`rmd:166-170`).
  - It then applies a zero-variance filter (`rmd:174-175`) and recomputes TMM on the filtered matrix (`rmd:178-180`).
  - The runbook says the same: "Preserve initial TMM, filtering arguments/group levels, and the separate post-filter TMM calculation".
- **Impact:** a kernel fed raw column sums gives a different keep mask from the notebook, and the error is caught only late by parity. KIT-095 inherits the wrong pre-pass definition for its partitioned filter (`KIT-095.md`, `trx_expression_filter` row).
- **Repair:** reword the stage as follows:
  - the pre-pass is the effective library size `colSums × TMM norm.factors` on the unfiltered, design-ordered counts;
  - then `filterByExpr(group = design group)`;
  - then the zero-variance filter on the filtered counts;
  - then post-filter TMM and logCPM.

  Name both factor vectors as products.
- **Closure evidence:** separate tests of the pre-pass vector, keep mask, zero-variance mask and post-filter factors, each `identical` to a named reference product.
- **Blocked state:** implementation.

### KIT093-F04: should-fix, dimension 2 (science) and 1 (completeness): notebook stages missing from the kernel list

- **Evidence:** `KIT-093.md:27-34` against the notebook:
  - **Preliminary non-robust QL ranking** (`rmd:232-237`: `estimateDisp`, `glmQLFit`, `glmQLFTest(coef=2)`, `topTags(sort.by="none")`, `order(PValue, decreasing=TRUE)`). It was in the approved plan ("preliminary QL ranking, control sets, `ruv_cancorplot`") and was dropped from the ticket.
  - **Control-set construction** `max(10, floor(n × pct/100))` (`rmd:248-250`).
  - **"Raw" diagnostic logCPM** without TMM, on filtered counts (`rmd:490`; the runbook says the "raw" diagnostic is logCPM of already filtered counts).
  - **Contrast parsing:** the first non-header line of `contrast_strings.tab`, then `makeContrasts` (`rmd:540-548`).
  - **Annotation join:** GeneName fallback (`rmd:144-148`), `DisplaySymbol` and the column order (`rmd:558-562`), plus the top-10 up/down tables (`rmd:569-583`).
  - **Glimma** (`rmd:673-693`) and the provenance summary (`rmd:754-779`). No ticket in KIT-093 to 099 names them; grep finds "Glimma" in none.
- **Impact:** these products have no owning kernel. Either they are reimplemented ad hoc in KIT-095/097 adaptors without parity, or they are silently lost. The runbook's preservation rules for them (unrobust versus robust fits, tie order, one annotation match per GeneID, FDR populations) have no executor.
- **Repair:**
  - Add the preliminary ranking and control sets to the RUV-grid kernel's stated scope, raw logCPM to the QC-inputs kernel, and the contrast parse plus annotation join to the DA kernel.
  - Name explicitly which ticket owns Glimma and the provenance summary (or record them as out of scope with the user's decision).
- **Closure evidence:** a kernel-to-notebook-line map in `kit-093/` covering every computational chunk at `rmd:115-779`.
- **Blocked state:** implementation.

### KIT093-F05: should-fix, dimension 1 (completeness) and producer-consumer: no map from kernel to reference product; KIT-040 freezes too little for the acceptance criterion

- **Evidence:**
  - `KIT-093.md:52`: "Every kernel output identical to its KIT-040 reference product".
  - KIT-040 freezes only: filter mask, "TMM factors" (it does not say which pass), logCPM, grid scores, best k and percentage, W, the edgeR table, RUVs-corrected logCPM, RUVIII `newY`, "PCA and RLE inputs" and file hashes.
  - It does not freeze: the pre-filter effective library sizes, the zero-variance mask, the preliminary ranking, per-percentage control sets, the cancor curve data (the input to the CCA-diagnostics kernel), the design/contrast matrices, raw logCPM (unless it counts as a "PCA input"), the annotated table, or the RLE statistics.
  - The RLE statistics source was already a TRX finding (TRX-A08: use `ggplot2::layer_data`).
  - Plot objects and PNG/PDF hashes cannot be `identical` across processes (environments, timestamps).
- **Impact:** for several kernels, the acceptance criterion cannot be closed, or is closed against nothing.
- **Repair:**
  - Add a table to the ticket: kernel → output → KIT-040 product name → comparison (`identical`, or the data layer for figures).
  - Name the missing products as a KIT-040 dependency requirement (an addendum to KIT-040 before it runs).
  - State that plot kernels compare plot data or layer data, never objects or file hashes.
- **Closure evidence:** every row of the table has a reference product, and the test file has one expectation per row.
- **Blocked state:** implementation (and KIT-040 execution).

### KIT093-F06: should-fix, dimension 1 (completeness) and 2 (science): the preflight refusals and grid-selection rules required by TRX-005 and the runbook are missing

- **Evidence:**
  - `KIT-093.md:25-26,53` cover only design faults and invalid count fixtures.
  - TRX-005's acceptance criteria and runbook line 115 also required classed errors for:
    - too few features (`max(10, …)` above the available genes; in R, `ctl_idx[NA] <- TRUE` with a scalar is silently ignored, so the notebook runs on fewer controls);
    - no successful grid member (the notebook's `purrr::compact` then `which.max` of an empty vector makes `grid_results[[integer(0)]]` an uninformative error);
    - insufficient rank or W confounded with group;
    - an unparseable or non-estimable contrast;
    - zero libraries, and an empty filter result (runbook line 117).
  - The runbook (line 42) and TRX-005 require failed grid members to stay visible. The notebook drops them (`rmd:281`).
  - The ticket also does not state that the selection uses the *unrounded* composite score and the first maximum (`rmd:295`), whereas `df_summary` is rounded (`rmd:288-290`).
- **Impact:** degenerate inputs either crash without a class or silently select from fewer controls. If an implementer keeps failed members, the product differs from the notebook without a declared divergence.
- **Repair:**
  - List each preflight refusal with its condition class.
  - Decide and record whether failed grid members are retained. If they are, record it as a declared divergence, with an extra product alongside the notebook-identical winner.
  - State the unrounded first-maximum rule.
- **Closure evidence:** one fixture and one `expect_error(class=)` per refusal; a test with a tied composite score selects the first.
- **Blocked state:** implementation.

### KIT093-F07: should-fix, dimension 2 (science): the count policy is under-specified

- **Evidence:**
  - `KIT-093.md:26,40` say only "NA, negative or fractional".
  - The test TRX-005 proposed, `x == round(x)`, accepts `Inf` (`Inf == round(Inf)` is `TRUE`, checked here).
  - The ticket does not say:
    - where the check runs (S4 validity, constructor, or the entry of each kernel) or that it runs before any statistics;
    - that integer-valued doubles pass;
    - whether the notebook's `round()` before `newSeqExpressionSet` (`rmd:407`) is kept;
    - the storage type of the count matrix. `identical` distinguishes integer from double; `vroom` reads doubles; `round(1L)` returns a double.
  - TRX-A07 required asserting that the final-model counts equal the filtered counts.
- **Impact:** non-finite counts could reach edgeR, and parity could fail, or pass for the wrong reason, on type alone.
- **Repair:**
  - Refuse NA, NaN, ±Inf, negative and non-whole values.
  - Accept integer-valued doubles.
  - Run the check once at construction and validation, before any statistics.
  - Keep the `round()` call verbatim (a no-op under the policy).
  - Declare double storage.
  - Assert `counts(set_ruvs)` is identical to the filtered counts.
- **Closure evidence:** fixtures for each refused class including Inf; one for an integer-valued double; the identity assertion.
- **Blocked state:** implementation.

### KIT093-F08: should-fix, dimension 2 (science) and 1 (completeness): design validation omits ordering, the source of group levels and the coupling to the group-column name

- **Evidence:**
  - `KIT-093.md:25` lists the checks only.
  - TRX-004 (after TRX-A01) required the design in design-file order plus a count-column index (`counts[, index]`, as at `rmd:152`), and refusal against *declared* levels.
  - The notebook hard-codes `levels = c("RPMI","Sera")` (`rmd:139`). Level order fixes `coef = 2` in the preliminary fit and the column order of `~ 0 + group + W` (`rmd:529-533`).
  - The contrast string is written against the model column names (`groupSera - groupRPMI`), which depend on the design column being named `group`. The class has a configurable `group_id`.
  - The existing corpus design writer uses groups `WT`/`KO` (`tools/profiling/kit_omics_corpus.R:26`). If KIT-040 reuses it, the unchanged notebook turns every group into `NA`.
- **Impact:**
  - Sample order of every product (W rows, PCA, `newY` rows) can diverge.
  - Changing the level order changes floating-point results, so the `identical` check fails.
  - A non-`group` `group_id` breaks the contrast.
  - The reference may be degenerate.
- **Repair:**
  - Specify: input levels declared (default from the contrast or design, in the notebook's order); output in design-file order with an index; the model frame column is named `group` regardless of `group_id`.
  - Add a shuffled-column test.
  - Flag to KIT-040 that the trx corpus must use the levels the notebook hard-codes (or the driver must supply them without editing the `.rmd`).
- **Closure evidence:** a shuffled-column test, a test that a non-`group` `group_id` gives the same products, and KIT-040's manifest showing the levels.
- **Blocked state:** implementation.

### KIT093-F09: should-fix, dimension 2 (science): the technical-replicate RUVIII acceptance is not falsifiable and its source column is undefined

- **Evidence:**
  - `KIT-093.md:33,54`: "runs and is recorded as a distinct product". The notebook has no reference for this path.
  - The ticket does not say which design column feeds `technical_replicate_id`. In the corpus generator, `replicates` is the *biological* replicate number within a group, and `tech_rep_group` is the biological sample that technical replicates share (`kit_omics_corpus.R:15-16,46,142`).
  - Building `M` from `replicates` would pair samples across treatment groups and remove biological signal.
  - The runbook (line 114) says: "Do not silently substitute the `Replicate` column or relabel biological replicates as technical ones."
  - KIT-040's corpus does not say it contains technical replicates.
- **Impact:** a wrong `M` passes "runs", and the product can remove the treatment effect undetected.
- **Repair:**
  - Define the source (`tech_rep_group`-type column) and refuse with a classed error when it is absent or NA, when a replicate set spans groups, or when no set has ≥ 2 members.
  - Make the oracle a direct `ruv::RUVIII(Y = t(logCPM), M = ruv::replicate.matrix(tech_rep_id), ctl, k)` on the reference's logCPM and controls.
  - Require that the corpus used carries technical replicates.
  - Label the default product "M from treatment group (notebook)".
  - Record the pseudo-replication caveat when technical replicates also enter edgeR as independent samples.
- **Closure evidence:** `identical` to the direct call; the three refusals tested; the product label in its metadata.
- **Blocked state:** implementation.

### KIT093-F10: should-fix, dimension 4 (references) and 3 (practice): the class-declaration entry changes an asserted set; the test is not in the allowlist, and the file is also claimed by KIT-095

- **Evidence:**
  - `tests/testthat/test-workflow-omics-owner.R:8-9` asserts `expect_setequal(workflowDeclaredClassNames("assay"), c(four classes))`. A fifth assay class fails it.
  - The owner's metadata reader takes that list (`R/utils_workflow_dia_state.R:229,256`).
  - KIT-093's `files_to_modify` does not include the test.
  - KIT-095 also lists `R/utils_workflow_class_declarations.R`.
  - `validateWorkflowClassDeclaration` (`:139-161`) is satisfiable: `count_data` in `lists`, `design_matrix` in `tables`, five roles.
- **Impact:** either the targeted run fails, or the test is edited outside scope. The ownership of the entry is ambiguous.
- **Repair:**
  - Add the owner test to `files_to_modify`, with the expected-set change named as the only intended observable.
  - State that KIT-093 owns the entry and KIT-095 only reads it.
  - Name the new R files' addition to the DESCRIPTION `Collate` field (`DESCRIPTION:35`).
- **Closure evidence:** the owner test diff changes that one vector only; `validateWorkflowClassDeclaration("TranscriptomicsAssayData")` passes.
- **Blocked state:** implementation.

### KIT093-F11: should-fix, dimension 1 (completeness): the parity tests can pass by skipping

- **Evidence:**
  - The reference lives under `docs/qa/workflow-execution/kit/kit-040/` (`KIT-040.md`), outside `tests/testthat/fixtures`.
  - The package's frozen-dump test pattern skips when its file is absent (`test-workflow-omics-declarations.R:9`, `skip_if_not(file.exists(...))`).
  - Optional packages that are missing also skip (F01).
  - TRX-005 required an installed-package test. KIT-093 has only "install" (`KIT-093.md:58`).
- **Impact:** "identical" can be reported from a run that compared nothing.
- **Repair:**
  - State where the tests read the reference (a path resolved from the repo, or a fixture copy with its sha256).
  - Require zero skips in `test-trx-kernels.R` as closure evidence.
  - Run the kernels once from the installed build.
- **Closure evidence:** a test log showing `[ FAIL 0 | WARN x | SKIP 0 | PASS n ]` for the file; the reference sha256 matches KIT-040's.
- **Blocked state:** implementation.

### KIT093-F12: should-fix, dimension 3 (practice): the divergences from the conventions are not recorded

- **Evidence:**
  - `conventions/R-bioinformatics.md:28-29` says omics containers **MUST** inherit `SummarizedExperiment`; `:31-33` requires `setValidity` plus a validating constructor and accessor-only access.
  - `AGENTS.md` requires recording a divergence when a compatibility contract wins.
  - TRX-005 decision 1 recorded exactly this (the kit's class declarations are the container contract). KIT-093 (`:24,36`) does not.
  - The approved plan's "≲ 50 lines per function" is also dropped (`conventions/R.md:87`).
- **Impact:** a reviewer cannot tell a deliberate divergence from an oversight.
- **Repair:**
  - One line recording the SummarizedExperiment divergence and its reason (parity with the lipid and metab classes and the kit declarations).
  - Require validity plus a constructor.
  - Restate ≲ 50 lines and ≤ 5–7 parameters.
- **Closure evidence:** the divergence recorded in the implementation summary; a validity test.
- **Blocked state:** implementation.

### KIT093-F13: minor, dimension 2 (science): determinism and numeric environment are not declared

- **Evidence:**
  - The ticket says nothing on RNG.
  - The path is deterministic: `ruv::ruv_cancorplot` uses `svd` and `cancor` (body read); RUVs and RUVIII with a fixed `k`, and edgeR, are deterministic.
  - Plan rule (`00-KIT-PLAN.md:41`): a kernel needing randomness takes an explicit `seed`.
  - The worker pins `OMP`/`OPENBLAS=1` (`R/utils_workflow_execution.R:1163`). The reference render does not say what it pins.
- **Impact:** no stated basis for `identical`; a later BLAS or thread change would show up as a mysterious parity break.
- **Repair:**
  - State "no RNG; no seed parameter; determinism shown by two runs giving identical results".
  - Record the BLAS/LAPACK identity, thread environment and package versions alongside KIT-040's.
- **Closure evidence:** both runs identical; versions and BLAS recorded in `kit-093/`.
- **Blocked state:** implementation.

### KIT093-F14: minor, dimension 2 (science): the tolerance mechanism and the PCA routine are not pinned

- **Evidence:**
  - `KIT-093.md:40,52`: "tolerances … declared first", with no location, time or approver.
  - The approved plan offers "`pcaWithMissing`, or `prcomp` as the notebook does". NIPALS is not `identical` to `prcomp`.
- **Impact:** a tolerance could be introduced after a failure.
- **Repair:**
  - Tolerances live in a committed `kit-093/parity-criteria.md` before the first comparison.
  - A tolerance is admissible only for a declared divergence (for example retained failed grid members). Any non-identical result from an unchanged call on the same inputs is a defect, not a tolerance.
  - Pin `stats::prcomp(t(mat))`.
- **Closure evidence:** the criteria file's commit precedes the first parity commit.
- **Blocked state:** decision.

### KIT093-F15: minor, dimension 2 (science): the inferred-controls label has no check

- **Evidence:** `KIT-093.md:40` states the label; no acceptance item tests it.
- **Impact:** the provenance claim is unverified.
- **Repair:** the grid and RUVIII products carry `control_source = "inferred_least_significant"` (or similar), with a test.
- **Closure evidence:** the test.
- **Blocked state:** implementation.

### KIT093-F16: minor, dimension 3 (practice) and 1 (completeness): one large item without a planning gate; the kernel interfaces are undefined

- **Evidence:**
  - The ticket has `needs_planning: false`, 7 new R files, a class, validation, a policy and about 8 kernels, in 8h.
  - The kernel signatures and product names that KIT-095 declares as steps are not given.
  - TRX-005 handled this with a planning step.
  - The budget (`:44`) gives no runtime ceiling. AGENTS.md asks for one, though the in-process cost is small.
- **Impact:** KIT-095's step table and record outcomes are designed against an unknown interface; correcting it later costs a ticket.
- **Repair:** add a kernel interface table (input, parameters, output products) to the ticket, or set `needs_planning: true` with that as the planning output; add a per-file ceiling.
- **Closure evidence:** the table is present before code.
- **Blocked state:** implementation.

## References verified accurate

- `findBestKElbow` exists (`R/func_prot_norm_optimization_helpers.R:64`) and differs from the notebook's `find_best_k`.
- `requireOptionalPackage(package, purpose)` exists (`R/utils_optional_packages.R:6`, class `multischolar_package_missing`).
- `workflowClassDeclarations()` exists, and the slot shape the ticket proposes mirrors `LipidomicsAssayData` (`R/func_lipid_s4_assay_data_class.R:47-59`).
- No `R/func_trx_*` or `test-trx-*` exist yet, so all are new files.
- The `R CMD check` baseline of 3 WARNINGs and 3 NOTEs matches `kit-089/check.log:150`.
- The notebook has 793 lines.
- edgeR is genuinely absent from DESCRIPTION; Glimma and Biobase are already in Suggests.

The inaccurate references are F02 (`ruv`) and F10 (the test and file ownership).

## Verdict per dimension

1. **Completeness: not executable as written.**
   - Blocking: the dependency-acquisition gap (F01).
   - Missing relative to TRX-004/005 and the runbook: the preflight refusals (F06), design order and levels (F08), the count-policy details (F07), the kernel-to-reference map with the products KIT-040 must freeze (F05), and protection against skip-vacuous closure (F11).
2. **Scientific rigor: the contract is sound in intent but not yet falsifiable in four places.**
   - The filter pre-pass is misdescribed (F03).
   - Several notebook stages have no kernel (F04).
   - The technical-replicate RUVIII acceptance is "runs" (F09).
   - The tolerance and PCA routine are unpinned (F14).
   - RNG: the path is deterministic, which should be declared (F13).
   - RUVIII `M` defaulting to group is correctly carried as the user's decision but should be labelled.
   - edgeR QL with W and the k caps (3 below 15 samples, else 4) match the notebook.
3. **Implementation practice: mostly aligned** (reuse of the S4 pattern, `requireOptionalPackage`, no new framework).
   - The `ruv` move to Suggests is wrong (F02).
   - The SummarizedExperiment and validity divergences are unrecorded (F12).
   - The work is one large item without an interface table (F16).
4. **Reference accuracy: mostly accurate.**
   - `ruv` is already a hard Import (F02).
   - The class entry breaks an unlisted test and its file is co-owned with KIT-095 (F10).
   - RUVSeq and EDASeq are assumed available but are not installed (F01).

**Overall:** KIT-093 should not start until F01 is decided and F02–F11 are repaired in a ticket revision (KIT-040 needs a matching addendum for F05 and F08).
