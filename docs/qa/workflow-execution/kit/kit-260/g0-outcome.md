# KIT-260 gate G0: outcome

Ticket r1.2, section 6.1. Branch `kit-260-g0` (script and records only; never merged), from `v0.9-staging`
16112917329a63ab062655f04db5bf1dd7df9729. Script `tools/kit260/g0.R` at e3e4fd25db2a74690d03b1803e8e210c2492d56b
(16112917 + the script). The predictions are in the script (`PREDICTIONS`), written to `predictions.tsv` before any case
ran.

- Run 45592957172 (request `kit260-g0-e3e4fd25-2`): every case, one scripts job, 18.0 job-minutes; groups d0 9 s,
  core 50 s, d6 312 s, d7 24 s, all exit 0. R 4.6.1, stringr 1.6.0 (ICU 70.1), limma 3.68.4.
- Run 45592901435 (request `kit260-g0-e3e4fd25`): the implementer pushed the blank request template by mistake; the
  plan job refused it in 11 s and no job ran. It is not an attempt at G0 and used no runner compute beyond the plan.
- Core-hours: 18.0 x 4 / 60 = 1.20 (G0 expected 2.0, ceiling 4.0). No instrument rerun used.
- Seeds: case k draws with `set.seed(260000 + k, "Mersenne-Twister", "Inversion", "Rejection")` before any package
  call: D1, D4, D5 260001 (one object); D2, D2n 260002; D3 260003; D6 260006; D7 260007. D0 and D0b draw nothing.
- C14 is not repeated (6.1): KIT-242's run 45584784285 observed it at 3eeaf100 (5 `w/o MS2: ` rows used as A1's
  standards and removed, the 2 `(d7)` rows kept, `itsd_counts_per_assay` A1 = 5, the integration input marking the 5),
  record `kit-242-g0:docs/qa/workflow-execution/kit/kit-242/outcome.md`. Between 3eeaf100 and 16112917 no file it
  exercised changed (`git diff --stat 3eeaf100 16112917 -- R/`: the ITSD methods, the omics declarations and the
  integration input are not in the list), so the record is reused as at the same base (INFERRED from the diff).

Tables in `g0-45592957172/`: `predictions.tsv`, `observations.tsv` (case, class, arm, item, value), `conditions.tsv`
(warnings, kept messages, logger lines), `d0-matches.tsv` (every pattern and alternative against each of the 20 names),
`d0-counts.tsv` (observed against predicted), `d0b-previews.tsv`, `d6-arms.tsv`, `d6-shift.tsv`, `groups.tsv`,
`provenance.txt`, the child logs and `console.log`.

## Observed against predicted, case by case

The synthetic objects are section 4's: two assays, 2 groups x 6 samples, 40 ordinary rows plus the named rows, id
`Alignment ID` integers (named rows are ids 41 to 60), the pattern taken from `get…ColumnDefaults()` at the commit.

| Case | Predicted (ticket 6.1) | Observed | Verdict |
|---|---|---|---|
| D0 | base MS-DIAL default matches 13 of the 20 names (the 4 intended, the 5 `w/o MS2: `, `w/o MS2:Glucose`, `w/o MS2:IS_Caffeine`, `CIS-ACONITIC ACID`, `HIS-LEU`); NA gives NA; (a) 6, (b) 4, (c) 5, (c) hyphen 6, (e) on (c) 6; other formats only their named labels | stringr on the runner: base 13 with exactly those names, per alternative `^IS_` 1, `_d[0-9]+$` 1, `ISTD` 1, `IS-` 3 (`CIS-ACONITIC ACID`, `HIS-LEU`, `IS-Leucine`), `w/o MS2:` 7; the missing name gives NA under every pattern. (a) 6, (b) 4, (c) 5 (+ `PC 15:0_18:1(d7)`), (c) with `[_-]d[0-9]+$` 6 (+ `Caffeine-d9`), (e) on (c) 6 (+ `w/o MS2:IS_Caffeine`), each with exactly the predicted names. Progenesis, XCMS, Compound Discoverer and LipidSearch defaults match `IS_Caffeine` and `Caffeine-ISTD` only; custom (NA) matches none. Lipid and metab strings identical per format; the base string equals the ticket's | as predicted (the Python-`re` counts are confirmed by ICU) |
| D0b | the MS-DIAL default for every metab layout and every lipid layout except LipidSearch, which gets `^IS_\|ISTD\|Internal\|_IS` | metab msdial, custom, malformed, progenesis, xcms, compound_discoverer: `^IS_\|_d[0-9]+$\|ISTD\|IS-\|w/o MS2:` (detected msdial, unknown, unknown, progenesis, xcms, compound_discoverer). Lipid the same six, plus lipidsearch (detected lipidsearch): `^IS_\|ISTD\|Internal\|_IS` | as predicted |
| D1 | metab, regex mode, median: `itsd_counts_per_assay` 13 per assay; those 13 and the NA row absent; `PC 15:0_18:1(d7)`, `Caffeine-d9`, `no MS2: Glucose`, `low score: Glucose`, `cis-Aconitic acid` present; console "Identified NA"; study-parameters 13 integer ids, first 5, "... +8 more"; record `itsd_applied` TRUE, no names | declaration (runSavedOperation) and worker (`omicsWorkerRun`, handoff `complete`, product rebuilt from its checkpoint) alike: 60 -> 46 rows per assay; absent: the 13 matched rows and the missing-annotation row; present: `cis-Aconitic acid`, `PC 15:0_18:1(d7)`, `Unknown`, `no MS2: Glucose`, `low score: Glucose`, `Caffeine-d9` and all 40 ordinary rows; args `itsd_counts_per_assay` {A1 13, A2 13}, features ids 41-51, 55, 58. Console "Identified NA ITSD features via regex." then "Using 13 ITSD features for normalization." and "Removed 13 ITSD features after normalization." Study-parameters: "A1 (13 features): 41, 42, 43, 44, 45, ... +8 more" (A2 the same). Record step: parameters `itsd_aggregation=median; itsd_feature_ids=NULL; method=ITSD`, outcomes rows 92, `itsd_applied` true; the record JSON contains neither `w/o MS2` nor `CIS-ACONITIC` | as predicted |
| D2 | lipid twin, annotation `Name`: as D1 | identical to D1 on every item (declaration and worker) | as predicted |
| D2n | lipid, name column `Metabolite name`: annotation NA; ITSD skips both assays; integration marks nothing | preview of a lipid table with MS-DIAL's header: detected msdial, id `Alignment ID`, annotation selected NULL (none), pattern the MS-DIAL default. Object with annotation NA: warning per assay "None of the specified ITSD pattern columns ('NA') found. Cannot identify ITSDs. Skipping normalization."; 60 -> 60 rows; args applied TRUE, no counts or features. Integration input (ITSD unticked): 0 rows marked per view, the model's drop 0; rule columns `Alignment ID` and "NA" | as predicted |
| D3 | LipidSearch object: normalisation matches nothing (class column); the integration input marks `IS_PC 15:0` (id column) | warning per assay "No ITSD features identified using regex '^IS_\|ISTD\|Internal\|_IS' in columns 'LipidClass'. Skipping normalization."; 42 -> 42 rows, `IS_PC 15:0` and `PC(15:0_18:1(d7))` kept. Integration input: 1 row marked per view, the model's drop 1 | as predicted in count. Which row was marked is not shown in the table (the script printed the class column, "PC" for every row); that it is `IS_PC 15:0` is INFERRED: it is the only id or class value any alternative can match (`PC(15:0_18:1(d7))` matches none, as D0's lipidsearch row shows for `(d7)` names) |
| D4 | D1's object, ITSD unticked: `is_internal_standard` TRUE on 13 per view, NA row FALSE; model drop 13 | 13 rows marked per view (the same 13 annotations as D1), the missing-annotation row FALSE; `integrationInternalStandards()` 13 per view, also with the column read back as text | as predicted |
| D5 | QC analysis, box empty: the 13 rows reported (per assay; 26 over two assays) | pattern used = the object's MS-DIAL default; "Found IS in column: Metabolite name"; 13 per assay (the 13 matched names), total 26 | as predicted |
| D6 | recorded, not judged; complete version, A against B: none one common constant (= minus the group contrast of log2(a_A / a_B)); cyclicloess and quantile identical to 1e-9; scale nonzero | see below | cyclicloess, quantile, scale as predicted; **none not observed** (no DA table) |
| D7 | where the pattern lives; replay reads the stored pattern, not today's default | module-CI MS-DIAL import with the box as proposed: `column_mapping$is_pattern` and `column_mapping.json` hold the MS-DIAL default; it is the only project file holding the text. With `getMetabolomicsColumnDefaults` mocked to `^KIT260_CHANGED_DEFAULT$` (instrument shown), `metabDesignObject` from the stored mapping gives the slot the stored default. A project with the stored default (ITSD 13 per assay, log2, cyclic loess) recorded and replayed under the mocked default: record valid; no record step holds a pattern parameter and the record file does not contain the text; replay preconditions `initial`, `metab_raw_data_s4` taken from the stored project; replayed states equal the history, no mismatch, ITSD again "Using 13 ITSD features" | as predicted: the pattern lives in `column_mapping.json` and the object slot, not in the record; replay reads the stored state. OD-3 needs no marker for OD-1 |

### D6 in detail (recorded; it cannot fail G0)

Object: section 4's r1.0 rows (54 per assay), +1 log2 in group B on the 7 ordinary matched rows; missing version: 17
cells of those 7 rows and 28 of the other 47 removed per assay, lowest first. Arms: A default pattern (11 used and
removed, 43 kept), B the 4 intended by id with the 7 removed beforehand (43 kept), C the 4 intended (50 kept). All
arms: `itsd_count` 11 / 4 / 4, rows kept 43 / 43 / 50, as designed.

Complete version, A against B, RUV off (43 common features per assay):
- cyclicloess: max |ΔlogFC| 1.8e-14 (A1), 1.4e-14 (A2). Quantile: 1.1e-14, 1.1e-14. Identical to 1e-9, as derived.
  No cyclic-loess difference above 1e-9 appeared, so R1-N05's interpolation cause did not arise.
- scale: max 0.125 (A1), 0.096 (A2), median 0.031, 0.037; no common constant (max deviation 0.125, 0.095). Nonzero, as
  predicted.
- none: **no DA table on any arm, either version.** The DA step failed with "Can't join `x$metabolite_id` with
  `y$metabolite_id` due to incompatible types (`x` character, `y` integer)". Only the method "none" fails; the other
  methods produce tables on the same objects. The realised predicted constant (minus the group contrast of
  log2(a_A / a_B)) is -0.716 (A1) and -0.665 (A2), recorded but not compared with anything.
- calls at q < 0.05: 0 on both sides in every RUV-off route.

Complete version, RUV Automatic with cyclicloess: k and controls A1 1/8 and A2 5/5 in arms A and B alike; max
|ΔlogFC| 5.8e-12 (A1), 7.8e-11 (A2); calls A2 2 up 1 down on both sides. C: k 2/5 (A1), 1/7 (A2).

Missing-value version (not predicted): A against B, cyclicloess max 0.0012 (A1), 0.0079 (A2): a small residual,
consistent with section 5 (5); quantile 1.4e-14, 1.1e-14 (exact, as section 5 derives for quantile with missing
values); scale 0.121, 0.074. RUV on: in arms A and B, **assay A1 has no RUV k or controls (NA) and no DA rows (0)**:
the assay is dropped, as KIT-242's C4 found, so those rows are not "equal to RUV off" (R1-N04's concern); A2 k 1,
controls 7. Arm C: A1 k 2/5, A2 2/10.

B against C (the deletion alone), complete, RUV off: median |ΔlogFC| 0.080 / 0.101 (cyclicloess), 0.061 / 0.141
(quantile), 0.092 / 0.258 (scale); C finds 7 up calls per assay (INFERRED to be the 7 planted rows, which C keeps; the
table holds counts, not names) and 0 to 9 down calls; B finds none. With RUV on, C's k and controls differ from B's (A1 2/5 against 1/8; A2 1/7 against 5/5).

## Verdict

G0's stop rule (0) is **not triggered**: D0, D0b, D1, D2, D2n, D3, D4, D5 and D7 are as predicted (D3's marked row
identified by inference, above). D6 is recorded; three of its four complete-version A-against-B predictions are
observed as predicted and the fourth (none) is not observable at this commit because the DA step refuses the objects.
The Python-`re` option counts are confirmed under stringr/ICU, so OD-1's string table stands as written.

## Observed versus inferred

Observed: every value in the tables above. Inferred: that C14's record holds at 16112917 (from the diff); D3's marked
row being `IS_PC 15:0`; that D6's "none" failure would also occur with the real import's id type (the import reads
`Alignment ID` as numeric; not checked whether double or integer).

## Adjacent discovery (not in scope; recorded only)

With a numeric id column (MS-DIAL's `Alignment ID`) and between-sample method "none", the metabolomics DA step fails
("Can't join `x$metabolite_id` <character> with `y$metabolite_id` <integer>"); the other methods do not. For the
coordinator to route.

## Not checked

- D6's "none" route (no DA table); the size of the factor's shift under None therefore remains derived, not observed.
- The integration model's own view summary (`internal_standards_dropped`): D4 counts with `integrationInternalStandards()`,
  the function the summary sums (`R/utils_workflow_integration_model.R:699-702`); the model was not run.
- D3's marked row by name (count only).
- The managed tab's log and the report text (not in G0's cases).
- G-R (optional; not run, by instruction).
