# KIT-040 audit (audit_only), 2026-11-06

- **Target:** `tickets/kit/KIT-040.md` revision 2.0 (sha256 `6ac5ea0e249c89e155e9020380939b5148bd237669668686e99a05b996799f53`), status `pending`, depends on KIT-092.
- **Mode:** `audit_only`. Nothing was repaired. No ticket, code or plan was edited, and no tests, installs, renders or journeys were run. The only R used was quick `packageVersion`/`sessionInfo`/`DGEList` probes.
- **Authorities read:** `AGENTS.md`; `conventions/R.md`, `R-shiny.md`, `R-golem.md`, `R-bioinformatics.md`; `tickets/kit/00-KIT-PLAN.md` (plan 3.8 rows, lines 231-242, 266, 301, 305); the approved plan `~/.claude/plans/write-a-series-of-atomic-dragon.md`; `docs/qa/workflow-execution/kit/protocol.md` (v1.6); the notebook `dev/distributable/plans/transcriptomics_workflow.rmd` (793 lines, sha256 `44d5f643946212c18cee10b5ebf77b1e7ebbe8d22ca37b40bafd08b5a0da652a`); the runbook and its audit JSON; `tickets/transcriptomics/TRX-001..005` with the TRX series audit; the neighbours KIT-092 to KIT-099, KIT-088 and KIT-050; `tools/profiling/kit_omics_corpus.R`; `tests/testthat/test-workflow-kit-harness.R`; `docs/qa/workflow-execution/kit/bindings-generator.R`.
- **Validator:** `next-openai-ticket.sh KIT-040` reports the schema valid and `status_consistent: true`. Mechanical validity only.
- **Environment facts (probed read-only):** R 4.6.1; Bioconductor 3.23; edgeR 4.10.1; limma 3.68.4; ruv 0.9.7.1; Glimma 2.15.2; pacman 0.5.1; here 1.0.2; pandoc 3.11. **RUVSeq and EDASeq are not installed in any library path.** BLAS is netlib `blas 3.12.1` (single-threaded).

## Finding register

| ID | Severity | Dimension | Blocked state |
|---|---|---|---|
| KIT040-F01 | blocking | completeness / best practice | execution |
| KIT040-F02 | should-fix | completeness / references | execution |
| KIT040-F03 | should-fix | best practice / scientific rigor | execution |
| KIT040-F04 | should-fix | scientific rigor | decision (the downstream identity claims) |
| KIT040-F05 | should-fix | scientific rigor | decision |
| KIT040-F06 | should-fix | scientific rigor | decision |
| KIT040-F07 | should-fix | scientific rigor / consumer trace | implementation (KIT-093) |
| KIT040-F08 | should-fix | scientific rigor | decision |
| KIT040-F09 | should-fix | completeness / scientific rigor | implementation |
| KIT040-F10 | should-fix | completeness / best practice | execution |
| KIT040-F11 | should-fix | completeness / best practice | execution |
| KIT040-F12 | should-fix | completeness / references | implementation (closure) |
| KIT040-F13 | should-fix | completeness | decision |
| KIT040-F14 | minor | scientific rigor | none |
| KIT040-F15 | minor | consumer trace | implementation (KIT-093) |
| KIT040-F16 | minor | graph and documentation integrity | none |
| KIT040-F17 | minor | references | none |

### KIT040-F01: the notebook cannot run "unchanged" because RUVSeq and EDASeq are missing, and its setup chunk would install packages at render time

- **Evidence:** the notebook's setup chunk calls `p_load(... RUVSeq, ruv, ...)` (rmd:45-60). If Glimma is absent it also calls `BiocManager::install("Glimma", ...)` (rmd:62-66). `packageVersion("RUVSeq")` and `packageVersion("EDASeq")` fail in both library paths (`~/R/x86_64-pc-linux-gnu-library/4.6`, `/usr/lib/R/library`). KIT-040:30 requires rendering the notebook unchanged. KIT-040:15 lists no library or install step, and KIT-040:65 says "No package code changes". The superseded TRX-001 §1 had an isolated library with a version receipt, and TRX-002 §1 forbade runtime installation. The runbook's setup row says "No runtime package installation ... freeze versions outside timed work". The approved plan puts RUVSeq/EDASeq into `Suggests` only at KIT-093.
- **Impact:** the first render either fails at `p_load` or makes pacman install RUVSeq, EDASeq and their dependencies from the network into the user's default library mid-ticket. Versions would then be chosen at render time, and the dependency set the series pins (R 4.6.1, memory note) could change silently. The reference would not be bound to a known environment. AGENTS.md requires explaining the necessity of a dependency and getting approval before adding one.
- **Smallest repair:** add a pre-run step to the ticket. Install RUVSeq and EDASeq for Bioconductor 3.23 once, before any render, into a declared library (isolated, or the default library with the user's approval). Use `update = FALSE` so existing edgeR, limma and ruv stay unchanged. Record a version receipt (R, Bioconductor, edgeR, limma, RUVSeq, EDASeq, ruv, Glimma, BLAS/LAPACK) in `kit-040/`. Make the driver refuse to render if any `p_load` target is missing, so the notebook never installs at render time.
- **Closure evidence:** the receipt in `kit-040/`; the driver's missing-package refusal test; edgeR, limma and ruv versions identical before and after the install.

### KIT040-F02: the notebook's fixed input contract is not written down, and the generator as described would break it

- **Evidence:**
  - The notebook reads `data/transcriptomics/GSE152840_gene_counts.txt` (rmd:117).
  - It reads `scripts/transcriptomics/design_matrix.tab` and `contrast_strings.tab` (rmd:124, 540), with the base directory taken from `here::here()` (rmd:78-81).
  - It forces `group = factor(group, levels = c("RPMI", "Sera"))` (rmd:139).
  - The contrast parser takes the first non-header line, `name=expression`, and the `make.names` design columns are `groupRPMI` and `groupSera` (rmd:540-550).
  - KIT-040:25 says to reuse "the design writer". That writer (`kitOmicsDesign`, `kit_omics_corpus.R:43-52`) takes groups from the recipe, which are `c("WT", "KO")` (`kit_omics_corpus.R:26`), and writes `sample`, not `Run`.
  - Probe: `factor(c("WT","KO"), levels = c("RPMI","Sera"))` gives all `NA`.
  - KIT-040:27 lists the files but not their names, locations, group labels or contrast string.
- **Impact:** if the generator follows the ticket's instruction to reuse the design writer, every group becomes NA and the render fails. If the corpus is renamed or remapped ad hoc, the mapping goes unrecorded. Either way the "once each" reference budget (F10) can be spent on input plumbing.
- **Smallest repair:** state the notebook input contract in the ticket:
  - the exact file names and layout under the prepared base directory;
  - the group labels `RPMI` and `Sera` set in the `trx` recipe;
  - `design_matrix.tab` columns `Run, group, replicates`;
  - the contrast line, for example `Sera_vs_RPMI=groupSera-groupRPMI` under a `contrasts` header;
  - the driver refuses to render when the design file is absent, so the positional fallback at rmd:126-133 can never run (see F06).
- **Closure evidence:** a harness test that checks the generated files against this contract; a driver refusal test for a missing design.

### KIT040-F03: "renders in a prepared working directory" leaves path resolution, product capture and repository writes unspecified

- **Evidence:**
  - `here::here()` (rmd:78-81) finds its root by searching upward from the working directory for `.here`, `DESCRIPTION`, `.git` and similar markers.
  - `rmarkdown::render` defaults `knit_root_dir` to the `.rmd`'s own directory, which is `dev/distributable/plans/`, inside the repository.
  - Rendered in place, the notebook would read from and write to `<repo>/data/transcriptomics`, `<repo>/scripts/transcriptomics` and `<repo>/results/...`. None of these is gitignored (`.gitignore` covers only `data-raw/*`). The `.nb.html` would be written into `dev/`.
  - The proof root has no root markers, so a copy placed there falls back to its working directory.
  - Of the products KIT-040:31-34 lists, the notebook writes to disk only `best_res` (rmd:771), the DE table and two text TSVs. The filter mask, both TMM factor sets, logCPM, W, `design_final` and the PCA/RLE inputs exist only in the knit environment. The TSVs (`vroom_write`, rmd:428, 721) are decimal text and cannot be compared bit for bit.
- **Impact:** an implementer can dirty the repository, read the wrong inputs, or freeze products from text round-trips, which loses the bits that every downstream `identical` comparison (KIT-093:46, KIT-098:51) depends on.
- **Smallest repair:** specify the mechanism. The driver runs a fresh `Rscript` that:
  - copies the `.rmd` byte for byte into a prepared directory outside the repository (under the proof root) containing a `.here` sentinel, and checks the copy's sha256 equals `44d5f643…`;
  - sets `knit_root_dir`, `output_dir` and `intermediates_dir` to that directory;
  - calls `render(envir = e)` with a new environment;
  - reads the named objects from `e` and freezes them with `saveRDS`, never from TSV.

  Assert that `git status` in the repository is unchanged after a render.
- **Closure evidence:** the driver header; a test or lifecycle record showing the copy hash and an unchanged repository tree.

### KIT040-F04: the frozen product list omits intermediates the kernels must match, and freezes the rounded grid

- **Evidence:** KIT-040:31-34 lists "filter mask, TMM factors, logCPM; grid scores, best k and percentage, W; the edgeR table, the RUVs-corrected logCPM, RUVIII `newY`; PCA and RLE inputs". The notebook shows what that list leaves out or leaves ambiguous:
  - **Two TMM factor sets.** One is computed before filtering on `dge_raw`, and `filterByExpr` uses it: probe `edgeR:::filterByExpr.DGEList` uses `norm.factors`; rmd:166-170. The other is computed after filtering (rmd:177-178). The ticket says only "TMM factors".
  - **Zero-variance mask and `gene_vars`** (rmd:173-174).
  - **The preliminary QL ranking** `tt_prelim` and `ordered_features_idx` (rmd:229-233), which define the control sets.
  - **The selected control gene set** `best_control_genes` (rmd:296), which RUVs and RUVIII both consume (rmd:412, 707). It is not listed at all.
  - **Per-member `cancor_plot$data` curves.**
  - **Grid members that failed.** These are dropped by `purrr::compact()` (rmd:281), and the ticket does not require recording them.
  - **Unrounded scores.** `df_summary` rounds to 4 decimal places (rmd:288-290), and "grid scores" does not say which form to freeze.
  - **`design_final` columns and the contrast matrix** (rmd:530-547).
  - **Whether `counts(set_ruvs)` equals the filtered counts** after `round()` (rmd:407, 513).
  - **PCA scores and every variance fraction**, and RLE box statistics. The ticket freezes only the "inputs".

  TRX-002 §2 required all of these. The runbook requires unrounded scores, visible failed members, the control IDs and the first winning tie.
- **Impact:** KIT-093 must show "every kernel output identical" (KIT-093:46) for the filter, grid, RUVs, RUVIII and PCA/RLE kernels. Without these intermediates, a mismatch in a final product cannot be localised. A kernel that selects a different but equally scoring control set, or tie, cannot be falsified. A winner chosen from rounded scores could pass.
- **Smallest repair:** extend the product list in KIT-040:31-34 to these named objects, each with its notebook object name and line, and freeze the grid from `grid_results` (unrounded) plus a per-percentage status row.
- **Closure evidence:** the reference rds schema listed in the outcome record; a harness test that each named product is present and of the right type and dimension.

### KIT040-F05: no reproducibility check or environment binding for a reference that downstream tickets compare against with `identical`

- **Evidence:**
  - Each corpus gets the reference "once" (KIT-040:46).
  - Nothing requires a second render to give identical frozen products, unlike TRX-002's "Oracle reproducibility" criterion.
  - Nothing requires recording package versions, `sessionInfo()` or BLAS in the reference. `conventions/R-bioinformatics.md:15-16` says "record parameters and `sessionInfo()` with results". TRX-A06 (the TRX series audit) flagged BLAS recorded late.
  - "Output is `reference-<corpus>.rds`, with its hashes" (KIT-040:36) does not say whether the hashes bind the inputs: corpus manifest sha256, notebook sha256, driver sha256.
- **Impact:** a downstream identity mismatch cannot be told apart from reference non-determinism, a package-version drift (F01) or an input change. The reference is not bound to its inputs or environment.
- **Smallest repair:**
  - Embed the following in each reference rds and in the outcome record: the corpus manifest sha256, the notebook copy sha256, the driver sha256, `sessionInfo()`, the BLAS/LAPACK paths and the versions of the statistical packages.
  - Add one repeat render on the 12-sample corpus (development class) whose frozen products must be `identical` to the first.
- **Closure evidence:** the bound fields present in the rds; the recorded result of the repeat render.

### KIT040-F06: the corpus cannot falsify "no positional design fallback", and it skips the annotation and count edge cases

- **Evidence:**
  - KIT-040:25-29 specify NB counts, groups, batches and DE genes only.
  - TRX-001 §2-3 required count columns shuffled relative to the design, with the permutation recorded, and DE checked by Run ID. It also required missing and duplicated `GeneName` values, all-zero and constant genes, an integer-valued-double valid case, and invalid fixtures that each differ from their parent in exactly one value.
  - The runbook's verification paragraph names "shuffled rows/columns, missing design, duplicated gene symbols, ties, all-zero/constant genes, invalid counts".
  - The notebook reorders by Run (rmd:150), falls back from `GeneName` to `GeneID` (rmd:143-147) and filters zero variance (rmd:173-174).
- **Impact:** if count columns are written in design order, a later kernel that attaches groups by position matches the reference exactly. The series contract "No positional design fallback" (KIT-040:42) then has no falsifier in its reference corpus. The annotation join and filter edge paths go unexercised.
- **Smallest repair:** add to the corpus spec:
  - count columns shuffled against `design_matrix.tab`, with the permutation in `truth.json`;
  - a few NA, empty and duplicated `GeneName` rows;
  - all-zero genes and one constant non-zero gene;
  - invalid fixtures that each differ from the valid parent in one value;
  - one integer-valued-double valid fixture.
- **Closure evidence:** the harness test asserts each property on the generated files.

### KIT040-F07: technical replicates are not specified, yet KIT-093 needs them for the selectable RUVIII `M`

- **Evidence:**
  - User decision 2 (KIT-040:20) makes RUVIII `M` selectable as `group` or `technical_replicate`.
  - KIT-093:48 requires "RUVIII with technical replicates runs and is recorded as a distinct product" on KIT-040's corpora.
  - KIT-040:25-26 mention "design.tsv with `replicates`", the biological replicate number. That column is the notebook's `replicates`, not a technical replicate.
  - Whether the `trx` recipe sets `technical_replicates = 2L` (`kit_omics_corpus.R:26`, `kitOmicsDesign` :48-49) is not stated. With one technical replicate per sample, `replicate.matrix` on technical IDs is the identity and the alternative `M` is degenerate.
- **Impact:** KIT-093's third acceptance criterion may have no valid input, or a degenerate one.
- **Smallest repair:** state that the `trx` recipe keeps two technical replicates per biological sample, which gives 24 = 6 x 2 x 2, 12 = 3 x 2 x 2 and 16 = 4 x 2 x 2. State that `design.tsv` carries `tech_rep_group`, and that `truth.json` records the technical structure.
- **Closure evidence:** the manifest's `technical_replicates = 2`; a harness test on `tech_rep_group` cardinality.

### KIT040-F08: corpus parameters and non-degeneracy observables are not declared

- **Evidence:**
  - KIT-040:26 says "negative-binomial whole-number counts with library-size factors, two groups, batches and differentially expressed genes". It gives no mean distribution, dispersion (or its trend), library-size range, DE fraction or fold-change, batch effect size, or low-count fraction.
  - The edge corpora's gene count is not given.
  - No observable shows that the reference exercises the notebook's branches:
    - `filterByExpr` removes some genes;
    - the grid has several distinct composite scores and no total failure;
    - `best_k` above 1 on at least one corpus;
    - a DE count strictly between 0 and the gene count.
  - If every grid member fails, `which.max` on an empty list errors (rmd:294-295).
- **Impact:** with no declared parameters, the corpus's realism and the reference's branch coverage cannot be checked. A degenerate reference (k = 1 everywhere, no DE, nothing filtered) would satisfy KIT-040's acceptance and give KIT-093's identity tests little power.
- **Smallest repair:** declare the recipe parameters in the ticket and write them to `truth.json`. Add an acceptance line listing the branch observables above, recorded from each reference with the result stated whatever it is. A degenerate result is reported, not tuned away after the outcome is seen.
- **Closure evidence:** `truth.json` holds the parameters; the outcome record holds a branch-coverage table per corpus.

### KIT040-F09: nothing in KIT-040 is named to refuse the invalid fixtures, and the notebook itself accepts fractional counts

- **Evidence:**
  - KIT-040:50 says "A harness test ... refuses the invalid fixtures". No refusing function or condition class is named, and KIT-040:65 forbids package code. The count policy's classed refusal belongs to KIT-093.
  - Probe on edgeR 4.10.1: `DGEList` accepts fractional counts silently, and errors on negative and NA counts.
  - The notebook then applies `round()` before RUVs (rmd:407), so a fractional fixture would render through to a reference with a silent correction.
- **Impact:** the Simulation Setup criterion cannot be executed as written. If the driver does not pre-check, a fractional corpus would produce a "reference".
- **Smallest repair:** say that the driver runs a preflight count check before any render: whole numbers by `x == round(x)`, no NA, no negative values, with integer-valued doubles accepted. It refuses with a named class, for example `kit040_invalid_counts`. The harness test asserts that refusal. KIT-093's package-level refusal stays separate.
- **Closure evidence:** the harness test fails before the change and passes after; one test per invalid fixture.

### KIT040-F10: the experiment budget contradicts itself and has no runtime ceiling

- **Evidence:**
  - KIT-040:46 says "Development-class runs (≤ 24 samples) unlimited; the reference on the screen and edge corpora once each". Every corpus here (12, 16 and 24 samples) is ≤ 24 samples, so each reference render is both unlimited and once-only.
  - No runtime ceiling or per-render limit is given. AGENTS.md (Verification and experiment budgets) requires the purpose, build/input bindings, runtime ceiling, allowed attempts and decision before an expensive run.
  - The "one bounded correction pass" is not tied to a failure class: instrumentation versus an unfavourable valid result.
- **Impact:** the attempt rule cannot be enforced. A failed reference render can be retried indefinitely as "development", or a valid but surprising reference can be re-rendered until it looks right.
- **Smallest repair:**
  - Development renders are unlimited only before the driver's bindings are frozen.
  - After that, there is one counted reference render per corpus, plus the F05 repeat, under a declared ceiling (for example 900 s each).
  - The single correction pass is allowed only for a demonstrated driver or instrumentation failure. Every attempt is preserved.
- **Closure evidence:** the outcome record lists each attempt with its bindings, elapsed time and classification.

### KIT040-F11: the lifecycle pair does not say how failure is forced, what counts as observed, or what kind of member it is

- **Evidence:**
  - KIT-040:37 and :55 say only "completion and a forced failure on the 12-sample corpus".
  - The notebook may not be edited, so the failure must come from the driver or the inputs. TRX-002 used `--fail-after ruvs-normalization`.
  - The protocol's lifecycle member is a harness member on the screen corpus (`protocol.md:34`).
  - KIT-068 recorded its lifecycle runs as "recorded, `exact = NA`, not protocol members" (plan line 141).
- **Impact:** the AGENTS.md lifecycle rule (instrumentation observes the work through completion and failure, not merely startup) cannot be checked. A failure forced before the render starts proves nothing about observing mid-render failure.
- **Smallest repair:** specify:
  - the failure injection, for example a `knitr` chunk hook set by the driver that raises an error after the `ruvs-normalization` chunk, which leaves the `.rmd` untouched;
  - the observables: terminal status, last completed chunk, partial products marked incomplete, non-zero exit;
  - that these are recorded uncounted runs, not `kit:lifecycle:*` members.
- **Closure evidence:** two lifecycle records in `kit-040/` carrying those fields.

### KIT040-F12: the TRX retirement has no named destination, file permission or working rollback

- **Evidence:**
  - KIT-040:38 says "Retire TRX-001 to 005 into the archive". No archive exists: the approved plan names `archived_tickets`, and no such directory exists anywhere. The series precedent for retirements is index entries with `retirement_reason` and `successor_owners` plus a supersession map (`00-KIT-PLAN.md:5, 265`; `../openai/KIT-SUPERSESSION.md`).
  - `files_to_modify` (KIT-040:15) omits `tickets/transcriptomics/*`, `tickets/transcriptomics/tickets-index.json`, `00-TRANSCRIPTOMICS-PLAN.md` and `tickets/kit/00-KIT-PLAN.md`.
  - The rollback note "Per-ticket `git revert`" (KIT-040:69) cannot restore them: `tickets/` is gitignored (`.gitignore:81`), and `git ls-files tickets` lists 0 files.
  - The plan's TRX paragraph is stale. It still says "TRX-005 ... stays pending as KIT-040's out-of-index precondition" (`00-KIT-PLAN.md:266`).
  - The approved plan says the `transcriptomics` worktree is "listed for the user". KIT-040:38 says only "left as it is". The worktree exists at `../MultiScholaR-transcriptomics` (`e7db197f`, 0 commits ahead).
- **Impact:** the third acceptance criterion has no definition of done. Ticket-state changes would fall outside the declared files and outside any rollback.
- **Smallest repair:**
  - Name the mechanism: set the five TRX index entries' status to `superseded`, with `retirement_reason` and `successor_owners` (KIT-040, KIT-093, KIT-099 as applicable).
  - Write a short supersession map that sends each dropped TRX requirement to an owner or an explicit drop. Examples: TRX-001's isolated library goes to F01, TRX-002's guarded wrapper to the kit harness, TRX-003 to KIT-099, TRX-005's preflight errors to KIT-093.
  - Correct plan line 266.
  - Add these files to `files_to_modify`, and give a file-copy rollback for the gitignored ticket files.
- **Closure evidence:** validator passes on both ticket directories; the map file exists; plan line 266 is consistent.

### KIT040-F13: the acceptance criteria do not cover most of the bounded outcome

- **Evidence:**
  - Acceptance (KIT-040:54-56) has three items: the reference rds, the lifecycle and the archive.
  - Nothing covers corpus generation itself: the manifests, `truth.json`, the three sizes and the invalid fixtures. Only the Simulation Setup test touches them.
  - The expected notebook sha256 is not written into the ticket. "Its sha256 recorded" can be met by recording any hash.
  - Nothing requires the reference to label the control genes "inferred", although the series contract says so (KIT-040:42).
  - Nothing labels the notebook's "raw" diagnostic as filtered, non-TMM logCPM (rmd:490).
  - Nothing labels `newY` as RUVIII with `M` from `group`, an unconfirmed replicate assumption (runbook, Scientific decisions, "Replicates").
- **Impact:** the ticket could close with a reference that is unlabelled, built from an edited notebook copy, or built on an incomplete corpus.
- **Smallest repair:** add acceptance lines for:
  - the corpus manifests and sizes;
  - the notebook copy sha256 equal to `44d5f643946212c18cee10b5ebf77b1e7ebbe8d22ca37b40bafd08b5a0da652a`;
  - provenance labels in the rds: `control_genes_inferred = TRUE`, `raw_logcpm = "filtered, no TMM"` and `ruviii_M = "group (notebook default; replicate assumption unconfirmed)"`.
- **Closure evidence:** the outcome record ticks each line, with the hashes.

### KIT040-F14: the k-cap boundary is not falsified, and the gene count of the edge corpora is not stated

- **Evidence:**
  - The notebook rule is `if (ncol(count_matrix) < 15) 3L else 4L` (rmd:242). The 12-sample and 16-sample corpora (KIT-040:28) cover both branches.
  - With 2 groups x 2 technical replicates, sample counts are multiples of 4, so a kernel written as `<= 15` or `< 16` gives the same caps on 12, 16 and 24 samples.
  - The edge corpora's gene count is not given.
- **Impact:** small. An off-by-one in KIT-093's port would pass identity against the reference.
- **Smallest repair:** note in the handoff that KIT-093 covers the boundary with a unit test at 14 and 15 samples, and state the edge corpora's gene count (for example 6,000, so `max(10, …)` never binds).
- **Closure evidence:** the handoff line and the edge-corpus dimensions in the manifests.

### KIT040-F15: where the reference lives and which output hashes count are left open

- **Evidence:**
  - KIT-040:36 puts full reference products in `docs/qa/.../kit-040/`. Precedent oracle rds files there are fingerprint stores of at most about 116 KB (`kit-053`, `kit-057`, `kit-086`). Full matrices of 6,000 x 24 doubles for several products and three corpora would be several MB of tracked binary.
  - "Output file hashes" (KIT-040:34) include files that are non-deterministic by construction:
    - `Sys.time()` in the summary text (rmd:759);
    - PDF and xlsx timestamps;
    - unseeded `geom_text_repel` placement (rmd:618);
    - the HTML outputs.
- **Impact:** a consumer may treat a varying hash as a reference mismatch, or fetch products from an unbound location.
- **Smallest repair:**
  - State where full products live (`docs/qa`, or under the proof root bound by sha256 with a fingerprint rds in `docs/qa`) and how KIT-093's tests reach them.
  - Class each output file's hash as "expected stable" (CSV and TSV data) or "recorded only" (figures, xlsx, html, summary text).
- **Closure evidence:** the outcome record's product table with location and stability class.

### KIT040-F16: plan and neighbour records still name KIT-040 where the approved plan moved them to KIT-099

- **Evidence:**
  - `00-KIT-PLAN.md:240` still gives KIT-088's dependencies as "KIT-040, KIT-087", and `:242` gives KIT-050's as including "KIT-040". The index (`tickets-index.json`) and the approved plan have KIT-099.
  - The KIT-088 description says "authored in full when KIT-040 closes" (KIT-088:5, :24).
  - Plan line 301 says "Full suite next due at KIT-040's close". The approved plan's checkpoints are KIT-090's close, then KIT-095 and KIT-099, and KIT-040:61 has no full suite.
  - The ordering KIT-090 → KIT-092 → KIT-040 is not in the graph: KIT-092 depends on KIT-089, not KIT-090, which is `in_progress`.
  - KIT-040's dependency on KIT-092 is ordering only (KIT-092 changes the verifier for partitioned steps; KIT-040 uses none). Its only shared file is `test-workflow-kit-harness.R`.
- **Impact:** small: confusion over when the full suite runs and when KIT-088 is authored. The R CMD check baseline "3 WARNINGs, 3 NOTEs" (KIT-040:61) could change at KIT-090's close.
- **Smallest repair:** correct plan rows 240, 242 and 301 and the KIT-088 stub text. Optionally add KIT-090 to KIT-040's dependencies, or note the ordering. Take the check baseline from the immediate predecessor's outcome.
- **Closure evidence:** the plan table matches the index; the validator passes.

### KIT040-F17: small reference inaccuracies

- **Evidence:**
  - The kind is called `trx` (KIT-040:25). The generator's guard admits only `lipid_screen` and `metab_screen` (`kit_omics_corpus.R:22`) and its CLI `switch` (:192) dispatches on the same names. The existing convention is `<omic>_screen`, and three sizes need either one kind with the per-group argument (CLI argument 3) or three bound corpus keys.
  - The approved plan's file list (`contrast_strings.tab`, `truth.json`) omits `design_matrix.tab`. The ticket adds it, correctly, because the notebook needs it.
  - The notebook line count (793) and path, the `design.tsv` `replicates` column (`kit_omics_corpus.R:142, 179`), `kitWriteOmicsManifest` (:83), the seed rule (:23), the k caps (rmd:242), `newY` (rmd:708), the existence of `test-workflow-kit-harness.R`, and TRX-004's checks being carried in KIT-093 were all checked and are accurate.
- **Impact:** small. The implementer has to choose a kind name and how sizes are bound.
- **Smallest repair:** name the kind or kinds, for example `trx_screen` with the per-group argument (6, 3, 4). Name the comparator for the lipid/metab byte-for-byte regression: the manifest sha256 bound by `bindings-generator.R` `bind_corpus` (:60-67) for `lipid-screen`, `metab-screen`, `lipid-no-standards` and `metab-no-standards`, plus the KIT-086 large corpora if they are regenerable.
- **Closure evidence:** the ticket text; a regression record with the hash comparisons.

## Verdict per dimension

1. **Completeness: not ready.** The outcome names the right deliverables but leaves too much out to execute or close:
   - a dependency the run needs and the authority to install it (F01);
   - the notebook's input contract (F02) and its render mechanics (F03);
   - the archive mechanism and files (F12);
   - the corpus and labelling acceptance criteria (F13).

   Compared with the superseded TRX-001/002/005, most fixture, oracle and reproducibility requirements were dropped without a recorded owner (F04, F05, F06, F12).
2. **Scientific rigor: sound in intent, not yet falsifiable.** The choice of the notebook, unchanged, as reference is sound. So are the "declare tolerances first" rule and the whole-number, inferred-control and no-positional-fallback contract. Gaps:
   - the reference lacks the intermediates and unrounded scores needed to falsify kernel identity (F04);
   - it is not bound to its environment or shown reproducible (F05);
   - the corpus cannot falsify positional design handling or show non-degeneracy (F06, F08);
   - technical replicates for the second `M` are unspecified (F07).

   Among the notebook traps: positional fallback is handled only in the contract wording (F02, F06). RUVIII `M` defaults correctly but carries no label (F13). Inferred controls are neither frozen nor labelled (F04, F13). `round()` before RUVs is safe only with a preflight (F09). The "raw" diagnostic carries no label (F13). The k caps are correct, with an untested boundary (F14).
3. **Implementation best practices: mostly proportionate, with three gaps.** Reusing `kit_omics_corpus.R`, the harness test file and no package code is proportionate and follows AGENTS.md. The gaps are a self-contradictory run budget with no ceiling (F10), a lifecycle with no defined failure injection or observables (F11), and a dependency install with no approval or receipt (F01). The conventions are respected where they apply. R-bioinformatics' `sessionInfo()` rule is unmet (F05).
4. **Codebase reference accuracy: mostly accurate.** Every cited file exists, and the notebook's line count, k caps, input formats and products match. Wrong or stale assumptions:
   - "reusing the design writer" yields groups that the notebook's hard-coded levels reject (F02);
   - the `trx` kind is outside the generator's guard and naming (F17);
   - the "archive" does not exist and `tickets/` is untracked, so `git revert` does not cover it (F12);
   - plan rows 240, 242, 266 and 301 and the KIT-088 text are stale (F16);
   - RUVSeq/EDASeq are assumed present (F01).

**Overall:** KIT-040 cannot start as written (F01 is blocking). After F01 to F13 are repaired in a ticket revision, it is executable within its stated scope.
