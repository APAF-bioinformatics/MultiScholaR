# KIT-203 outcome: whole-family Benjamini-Hochberg as the default flag of the proteomics clusterProfiler route

Ticket r1.4 (`tickets/kit/KIT-203.md`). Branch `kit-203-enrichment-family-bh`. Every run is listed in `runs.tsv`.

## Commits

| Commit | What |
|---|---|
| 82d9fc0e, 63b95f90 | G0 scripts (`gate-k203.R`, `capture-k203.R`); G0 records and the capture's instrument repair |
| 34cfcacd | fixtures `tests/testdata/enrichment-significance/` from the base capture (run 45329062402) |
| 1ee481a2, 5c18e57e | tests only (S1 to S4, S6 to S12, S14; 10.2's moved expectations); R1's parse repair |
| 31515bf8 | the fix (code, Rd, report templates, workbooks, NEWS, ledger row `VAL-<merge>`, KIT-173 `results.md`, `wording.md`) |
| 7b35ad65 | R1b records |
| 31ffbae6 | merge of `origin/v0.9-staging` at f6c86cfb (KIT-175, KIT-199 merged; replaces R2b). Its message names ea93d948, the staging head when the ticket was handed over; the merged head was f6c86cfb |
| a0e8cac6 | `mutate-k203.R` (X) |
| ee5b5c34, 3f737844, 198d5827 | R2 records; the correction pass (status text; the R3 A107-03 notice literal); S1's expectation corrected to the per-row column 5.3 specifies (m = 10 on each of 3 rows), made in the correction pass; S2's stand-in repaired (instrument, VAL-54); OD5 and X out dir (instrument) |
| c015d1e0, ecb224fe | `omicEnrichW1FrozenUi` re-frozen from R2c's head dump (OD5); R2c and X records |
| 056a151b | merge of `origin/v0.9-staging` at 6a397776 (KIT-194, KIT-185, KIT-197; no shared code file) |

## What was observed

- **G0 (run 45326277295, base 82d9fc0e, R/ = 9453e79b).** Environment equal to KIT-189's provenance (R 4.6.1;
  clusterProfiler 4.20.0, enrichit 0.2.0, qvalue 2.44.0, GO.db 3.23.1). 400 of 400 lists: `list_sha` equal, status ok,
  m = 1,249 on every list, m' 176 to 453; the shipped TSV's `significant` count equalled KIT-189's A-Q `calls` on every
  list (92 lists with a call, 12,151 calls), `p.adjust < 0.05` equalled A-BH (17 lists, 50 calls), the whole-family
  count A-BHf (4 lists, 13 calls). `g0/`. The capture script in the same run failed (plot-file stubs set after the
  project capture: instrument); rerun 45329062402 wrote the S3, S4 and S6 base fixtures (`g0b/`).
- **R1 (run 45330101688, 1ee481a2).** The new file stopped at parse (one closing parenthesis too many in S10:
  instrument, repaired in 5c18e57e). 10.2's moved expectations failed as recorded in `r1/failures.txt`.
- **R1b (run 45333590969, 5c18e57e; instrument rerun 1 of 3).** The new file: 9 expectations passed, 20 failed, 10
  errors; every test fails (`r1b/`). At its stated value: S1, S2, S6, S10, S11, S12 ("absent"), S3 (120 significant, no
  `p_adjust_family`), S4 ("argument unknown on base"), S14 (`significance_rule` NULL in the recorded call). **Not at the
  stated value:** S7 and S9 stop at the absent `ENRICHMENT_SIGNIFICANCE_RULES` (from which they read the expected
  texts) before they compare the notice or the paragraph, and S8 at the unknown `setEnrichmentPlotInputs` argument
  before it compares the plot; they fail on the tests-only commit, but the ticket's "old text", "old paragraph" and
  "plot by q-value" were not shown.
- **R2 (run 45335245494, merge 31ffbae6 against staging f6c86cfb, 18 files, check both sides).**
  - G1: match on every list (`r2/scripts/g1-summary.tsv`): the default's count equalled A-BHf on 400 of 400 (4 lists,
    13 calls), `p_adjust_family` and `family_size` equal to the formula and m on 400, and the `"qvalue"` run equalled
    A-Q on 400 (92 lists, 12,151 calls).
  - S5 (`replay-k203.R`): a record made on f6c86cfb replays on head with the q-value: 11 of 11 checks (declared names =
    base's ten, no `significance_rule`; values, service, contrasts, tables exact, plot inputs, file names and every TSV
    SHA-256 equal; no family column). KIT-173's T7(b) on its base record: 12 of 12; T7 file 21 tests, none failing.
  - S13 (`options-k203.R`): every one of 21 figures in the notes, report sentences and NEWS regenerated and found; the
    recorded figures: 8 random-list settings with the q-value estimate above 0.075, tested-list counts 88, 90, 78, 59,
    53, 42, 90, 125, 137, the q-value's pooled share on lists of 25 with parent terms 25.2 %, the mean FDP among lists
    with a call 46.5 % (q-value) and 3.3 % (default). The default's pooled share is 3.0 % (cell 10 without parent terms,
    11 of 361 = 3.05 %); the ticket up to r1.3 and the plan register (`AUDIT-KIT-203-plan-20261124.md:5`, `:131`,
    `:537`, `:545`, `:550`, `:554`) wrote 3.1 %, a double rounding; the coordinator ruled 3.0 % (r1.4).
  - Check: 0 errors, 3 warnings, 3 notes on both sides, identical.
  - **A result caused by the change:** the g:Profiler completed status text gained an empty line (`paste(..., sep =
    "\n")` turns the NULL of the clusterProfiler-only line into an empty element): `(b) in process: gprofiler`,
    `kernel_gprofiler`, `organism_filter`, `restore` failed on head, passed on base. S10 compared two head outputs and
    could not see it. The R3 A107-03 test asserted the old correction notice, which 5.4 replaces by design (not listed in
    10.2).
  - S1 compared the per-row `family_size` column with a scalar (all three rows were 10): an expectation corrected in
    the correction pass (audit K203-I04: under VAL-114 a changed literal is not an instrument repair). Instrument: S2's stand-in
    built a one-row ID from an empty vector (`paste0("GO:", integer(0))` is `"GO:"`) and errored before any assertion;
    OD5's dumps were not written (the out dir stayed relative and `test_file` runs in `tests/testthat`).
- **Correction pass (3f737844).** The Significance line is joined to the cutoff's line, so the g:Profiler text is
  byte-identical to base; the R3 A107-03 expectation now holds the 5.4 notice (flag sentence, then the correction
  sentence; type warning). In the same pass (198d5827): S1's expectation corrected to the per-row column 5.3 specifies
  (m = 10 on each of 3 rows); S2's stand-in, OD5's and X's out dirs repaired (instrument, VAL-54).
- **R2c (run 45341049018, 198d5827 against f6c86cfb).** New file 191 of 191. Against base, only these differ:
  `(a)` the UI capture (OD5) and the clusterProfiler captures `(b) clusterprofiler`, `kernel_clusterprofiler`, the
  KIT-113 W1 record test and the two W2 core tests, all failing on base since KIT-173 (section 12: not edited). Every
  other file `same` or `count-changed` (new tests passing). Check identical on both sides.
  - OD5 (`dumps/`): base dump equal to the base tree's stored block (precondition met; test (a) passes in the base
    shard); the two head dumps, separate processes, byte-identical (SHA-256 73f4acef...). Base-to-head diff
    (`dumps/base-head.diff`): the HTML digest, the input `significance_rule` (id, label, choices, default
    `bh_whole_family`), the cutoff's label ("Significance cutoff:") and help, the selector's help and two notes, the
    clusterProfiler table legend. Nothing else moved. The legend was not named in section 12's list; it is 5.4's static
    text. The block was replaced by the dump file's content (no hand edit; verified equal).
  - X (`x/`): S3 fails under (m1) (4) and (m3) (2), S1 under (m3) (3), S4 under (m4) (15), S14 under (m2) (1); the
    unmutated control passes every test.
- **F**: see section "F" below.
- **Audit repairs (run 45362501898, 27d9e5cd).** K203-I02: `test-general-enrichment-shared.R` on head, 139 passed, 0
  failed (the default run's two conversions asserted again). X2: S7, S8 and S9 each fail under m5, m6, m7; the control
  passes. 31.5 job-minutes, slightly above the declared "about 30", within X's row ceiling of 45.

## Acceptance criteria

| # | Met? | Evidence |
|---|---|---|
| A1 | met | G0: 400 of 400 lists equal to A-Q; 92 of 400 with a call (`g0/g0-summary.tsv`) |
| A2 | met | G1: default = A-BHf on 400 (4 lists), q-value = A-Q on 400 (`r2/scripts/g1-summary.tsv`) |
| A3 | met | S1, S2, S3 pass in R2c; S3 counts m independently (314 on the null draw, no parent terms) |
| A4 | met | S4 (TSV SHA-256, `@result`, plot inputs) and S6 (declared names and values without `results` ids, outcome names) pass in R2c; the identity hash itself is S5's |
| A5 | met | S5 11 of 11; T7(b) 12 of 12 |
| A6 | met | S7 to S12 pass in R2c; g:Profiler status byte-identical to base after the correction |
| A7 | met | S13 21 of 21 |
| A8 | met | R1b: every new test fails on the tests-only commit; S7, S8, S9 fail there at an absent object or argument before their stated comparison (above). X: every guard fails under its mutation. X2 (run 45362501898, audit K203-I03): S7 fails under m5 (base's notice), S8 under m6 (the plot never draws the family value), S9 under m7 (the LFQ template never states the family rule); so each reader test is shown able to fail at its comparison |
| A9 | met | F: differences only in the files of 10.2 and the section 12 captures; check identical |
| A10 | met | `dumps/`, the re-freeze commit |
| A11 | met | ledger `VAL-<merge>`, NEWS (no version change), KIT-173 `results.md`, `wording.md` |

## F

Run 45347178521, head 056a151b (the fix merged with staging 6a397776) against base 6a397776, every test file one
process each, 20 shards per side, R CMD check both sides; 1,166.3 job-minutes (ceiling 1,200). Observed
(`f/compare.tsv`): 456 files `same`; 3 `count-changed` with no failure (`test-go-propagation.R` 319 to 331,
`test-general-enrichment-shared.R` 136 to 139, `test-workflow-operation-registry.R` 1,093 to 1,095 passing
expectations); `test-enrichment-significance.R` only on head, 191 of 191; 2 `message-changed`
(`test-omic-enrich-module.R`, `test-omics-enrich-core.R`): only the five clusterProfiler captures that fail on base since
KIT-173 (`(b) clusterprofiler`, `kernel_clusterprofiler`, the KIT-113 W1 record test, the two W2 core tests; section 12:
not edited). No new failure, no timeout, no crash. The re-frozen UI capture `(a)` passes on head. The GUI-entry files
are `same`. Check: 0 errors, 3 warnings, 3 notes on both sides, identical (`f/check-compare.txt`).

## Recorded divergences

- `addEnrichmentFamilyAdjustment` reads `@geneSets` and reads and writes `@result` of enrichit's `enrichResult` directly
  (the write as `sortEnrichmentGeneIds` already does) (R-bioinformatics
  accessor rule): no accessor returns the size-filtered family (K203-P16).
- Moved expectation outside 10.2: `test-omic-enrich-module.R` R3 A107-03 (the notice literal), changed in the correction
  pass, not in the tests-only commit.
- Moved, not named in section 12 (audit K203-I06; `f/compare-failures.txt`): in the five clusterProfiler captures that fail
  since KIT-173, besides section 12's list, the completed status's "Significance:" line and the download summary's
  digest also move. Plan section 11's KIT-173 row gains: "and KIT-203's `significance_rule` (declared, kernel argument,
  outcome, notice, the completed status's Significance line, the download summary's digest, TSV columns)". In the UI
  capture the table legend moved (5.4's static text; recorded above).

## Hand-offs

**KIT-193** (audit K203-I05; the coordinator sends it):

> KIT-203 merged: the route's flag is `R/func_general_enrichment_orchestration.R:669`; under `significance_rule = "qvalue"` it is exactly `significant = qvalue < q_cutoff` (your `QvalueBh.flag`); under the default `"bh_whole_family"` it is `p_adjust_family < q_cutoff` with `p_adjust_family = p.adjust(p, "BH", n = m)` (your `callsWhole`, strict). (1) Fill row T4's digest on the merged `processEnrichments`. (2) The docstrings `Spec/QvalueBh.lean:15`, `:65`, `:136` change only through your 5.3 item 4 (no statement change): `:544` becomes `:669`, and 'the route's flag' becomes 'the route's flag under `significance_rule = "qvalue"`'. (3) Pair (ii): stored tables made before KIT-203 are q-value runs; a re-run must pass `significance_rule = "qvalue"`. The nesting `p_adjust_family < α ⇒ p.adjust < α ⇒ qvalue < α` is KIT-203 S3's check (optional statement 4e, yours).

## For the owner's wording pass (VAL-41; audit K203-I07; no text changed now)

- The tab notes show post-hoc figures without saying so: 5.8 and 8.8 points, 3.0 % and 82.1 %, 10.7 %. The marking lives in `source`, which no user sees.
- "On lists of 25 with parent terms the average was 10.7 %" means planted lists (cell 13). Suggested: "On lists of 25 built around one process, with parent terms, the average was 10.7 %".
- `R/func_general_enrichment_orchestration.R:171` and `man/processEnrichments.Rd:63`: "(in simulations it called terms on up to 42.9 % of random protein lists)". The figure is true. It omits "at a cutoff of 0.05", and 42.9 % is without parent terms; the shipped propagation gives 22.9 %. Suggested: "(in simulations at a cutoff of 0.05 it called terms on up to 42.9 % of random protein lists without parent terms and up to 22.9 % with them)", as NEWS says.

## Not checked

- Down lists, real (correlated) lists, cutoffs other than 0.05, organisms other than E. coli and yeast (stated limits).
- The g:Profiler service's own `fdr` family (no service call).
- S7, S8 and S9 failing at their stated values on the tests-only commit itself (they failed earlier, see R1b); X2 shows each fails at its comparison under a mutation restoring base's behaviour of that reader.
- That KIT-193's row T4 digest follows the new `processEnrichments` body (KIT-193 not merged; handed off, ticket 8 item 6 (b)).

## For the change report

**What was wrong, who was affected.** On the proteomics enrichment route for organisms outside g:Profiler's 13
(clusterProfiler `enricher`, UniProt GO annotations), a Gene Ontology term was called significant when its Storey
q-value was below the cutoff. enrichit computes that q-value over only the terms with at least one listed protein, so
the share of true null terms is estimated far below 1 on random lists and the q-value calls terms where there is no
signal. KIT-189 measured it: on simulated random lists at 0.05, 1.7 % to 42.9 % of lists had a significant term (6 of 18
settings clearly above 7.5 %; with the shipped parent terms, lists of 25 proteins: 22.9 % E. coli, 16.7 % its 2,000-protein
subset); up to 23.6 % of random lists had 10 or more significant terms. Every user of this route (any organism not in
g:Profiler's set) was affected; G0 confirmed the shipped route gives KIT-189's q-value calls list by list.

**Code.** Base 9453e79b: the flag `R/func_general_enrichment_orchestration.R:544` (`significant = .data$qvalue <
q_cutoff`); the calls `:443-450`, `:460-467`. Fix 31515bf8 (correction 3f737844): `enrichmentSignificanceShipped`
(`R/func_general_enrichment_orchestration.R:33`), `addEnrichmentFamilyAdjustment` (`:122`), the argument
`significance_rule` (`:190`), the call after each `enricher` (`:594`), the flag by rule (`:669`); the plot
`R/func_general_enrichment_engine_helpers.R:255`; the notice `R/utils_workflow_enrichment.R:914`; replay of older
records with the q-value `:1219` (lines at 31515bf8; on the merge 056a151b the notice is `:970` and the replay default
`:1279`, the other lines unchanged); declaration `R/utils_workflow_dia_da.R` (`diaEnrichmentDeclaredNames`); tab
`R/mod_omic_enrich.R`, `R/mod_omic_enrich_panels.R`; four report templates; five workbooks.

**What changed in results.** By default a term is significant when its Benjamini-Hochberg adjusted p-value over every
term with 10 to 500 measured proteins (terms with no listed protein entering with p = 1; `p_adjust_family`, with the
family size `family_size`) is below the cutoff. On KIT-189's lists of 25 with parent terms, lists with a call fall from
92 of 400 to 4 of 400 (calls 12,151 to 13). The q-value flag stays selectable (`significance_rule = "qvalue"`), computed
as before, byte-identical. p-values, `p.adjust` and `qvalue` columns are unchanged. Analyses recorded before the change
replay with the q-value, files and identity unchanged; their reports carry the measured sentence. Every clusterProfiler
run now records and is told the rule it applied.

**Evidence and limits.** KIT-189 (score run 45311934798; VAL-123): the whole-family rule controls in 36 of 36 deciding
settings (random lists 1.0 % to 2.5 %; planted mean FDP 0.1 % to 1.7 %, pooled at most 3.0 %), at a power cost of up
to 5.8 points with parent terms, 8.8 without. This ticket: G0 45326277295, G1 and S5 45335245494, X 45341049018.
Simulated on reviewed E. coli K-12 and S. cerevisiae proteomes, up lists of 25, 126 and 400, cutoff 0.05; not simulated:
down lists, correlated real lists, other cutoffs and organisms.

## Budget

Core-hours (job minutes x 4 / 60): G0 1.47, G0b 0.49, R1 1.84, R1b 0.64, R2 9.00, R2c 8.21, F 77.76, X2 2.10; total 101.51
against 83.33 expected and 235.67 ceiling (the correction pass used 8.21 of its 111.33; F 1,166.3 of 1,200 job-minutes).
Instrument reruns (VAL-54): G0b (capture), R1b (R1 parse), and in R2c the S2 and OD5 repairs (one each). S1's
expectation correction rode the single correction pass (R2c), not an instrument rerun.
