# KIT-281 L1 declaration (blind lifecycle, once; committed before the send)

Authority: ticket KIT-281 r1.1 section 7, owner's acceptance of the ceiling 74.67 core-hours (VAL-371), the script audit AUDIT-KIT-281-scripts-20261126 (S01 to S11 and Re-check 1, RC1-1 to RC1-6), the coordinator's rulings (non-finite estimate = Stop; q-method reported, not a stop; self-test prediction repair; L1 sent once after it). This is an instrument lifecycle: it decides nothing about calibration. No counted job, no G3 table is used for a counted range, no scoring of counted data.

## Additions after the script audit (stated as such)
- **E11** (RC1-5): `sim-k281.R` on the counted range r = 401 (`K211_CELLS=A`) with no `g3-rule.tsv` in the tree must refuse through the last-net writer before `load_all`: exit 1, `instrument-stop.txt` containing "a counted job needs the G3 table", no `metrics.tsv`, no part file (no draw made). Change to `l1-k281.R` only.
- **E4** is the scorer self-test with the S01/S07 cases and the repaired prediction (instrument defect 1, `smoke/outcome.md`).
- G2 plans 200-, 171- and 77-replicate jobs from L1's seconds (S05); E1' records a missing unit (S04); calls name their arguments (S08).

## Binding
- `package_commit`: the commit of branch `kit-281` carrying this declaration (named in the request); `git diff --stat a5de1c2eb <package_commit> -- R/ src/ inst/ tests/ DESCRIPTION NAMESPACE` empty; `threshold_values_in_filter` absent from `R/` (KIT-280 not in the tree); G0: the three filter files equal 894b9cc8's (`l1/g0.txt`, checked again by L1).
- `base_commit` "none", `tests` [], `blas` "reference", `extra_packages` [], one `scripts` entry `Rscript docs/qa/workflow-execution/kit/kit-281/l1-k281.R l1-out`, `out` `l1-out`, `out_pattern` `"\\.(tsv|json|txt)$"`, `timeout_minutes` 55.
- SHA-256 of the files that matter: `sim-k281.R` 10fa5bfda5cb05c3ee84def2fd528b64a2bb84b4c51ce242731f5eb0b8b4a236 (unchanged since the smoke send); `score-k281.R` 0f838f47e3eb21bf4fe389e6e776c0639fba751f73880d0e589e4b3ddf6d8628; `l1-k281.R` 7d5d4cb295c4ea3968d33e93ea737cd662ad1efead1accfc7beb7e31fdb88c32; `selftest/make-selftest.py` 11b7891239df593d2988792da12ca4b7213bf3626d15d14cf5d9b1a351d9c49c; `selftest/expected.tsv` f846f88b3ef7958c9739bc0e820c60233217b68c2b673afd22f5a382d1104dd1; `selftest/cases.tsv` ea2b665542d51b0f3aaa324ae1f683a10f37d863a9e040c50d3e70abcd044d90. (L1 records the hashes it sees in `hashes.tsv` and `provenance-l1.json`.)
- Seeds: r 9201, 9202 (lifecycle), 201-205 (E1'), 9201-9204 (E6), 401 (E11, refused before any draw).

## Ceiling and ledger
`timeout_minutes` 55 = 3.67 core-hours at the ceiling (expected about 20 min = 1.33). Running ceiling before L1: 0.94 used + 3.67 + 69.00 planned = 73.61 <= 74.67. Boundary rule: a limit passes iff its computed double is <= the double of the written limit; no tolerance.

## Predictions (fixed now)
G0 1/1; **E11 1/1**; E2, E3, E7, E10 24/24 each (12 replicate-masks per range), E10 maximum absolute deviation 0; E8 2/2; E1' 120/120 (the DA-dependent flag may read FALSE: KIT-201/KIT-221 move p-values; flagged only); **E4: `SELFTEST: 2007 checks, 0 differ`** (units 2007, passed 2007); E5 2/2; E6 1/1; G2 3/3 (planned jobs of 200, 171, 77 replicates each <= 75 min; `g2-plan.tsv`); G3rule 24/24; G3 1/1 with `g3-count.tsv` (12 rows) and `g3-rule.tsv` (4 rows, one R in 200 to 571). All checks pass at the first attempt.

## Output (public)
`checks.tsv`, `timing.tsv`, `hashes.tsv`, `provenance-l1.json`, `progress-l1.tsv`, `g2-plan.tsv`, `g3-count.tsv`, `g3-rule.tsv` (and `l1-stop.txt` on an uncaught stop). Blind: nothing derived from a p-value, a call or a fold change leaves the rows except the one count of rule G3 (K281-P08).

## Failure handling
A difference in E4 or any mechanics failure with no result produced is an instrument defect of L1: hand back, do not repair or re-send without the coordinator's ruling (VAL-54: at most three re-sends, nothing about values, margins, rules or seeds changes).
