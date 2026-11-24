# KIT-191 outcome: the integration enrichment's tests as a list of options

## Resume here (2026-11-24, after the implementation audit)

Register `tickets/kit/audits/AUDIT-KIT-191-impl-20261124.md`: ready after repairs, one major (I01). The repairs are
written and committed but UNRUN (the one correction pass was used by R2c; the verification run is the owner's, pool
X114): a tests-only commit (d021463a) then the code commit (c85aa863). Repaired: I01 (the floor printed as the exact
fraction "1/1,001" by one function, `integrationEnrichmentFloorFraction()`, the form for the owner to confirm, O1; the
`_none` sentence when p_min > 0.05; k_min in whole numbers; U10 gains the round-up cases 1/10,001 and 2/1,001 and the
whole-number case B = 1,000, d = 10, m = 273 giving exactly 60), I02 (the limits clause names the kind's terms of row
E), I03 (the declared test read first; a value that is no offered test gives no figure, never camera's), I04 (exact
reads of `test` in the tab's summary and the record's ranking and Methods), I05 (the two refusals name the permutation
test), I06 (record slips), I07 (T10's score column exact; X's m2 also replays the branch-point store's old entries).
Not repaired: I08 to I13 (notes); the O3 wording items, below, wait for the owner. Staging not merged.

The one run needed: R2's 13 files, base against head with check and its scripts, on the merge with `v0.9-staging`;
ceiling 420 job-minutes; plus X's repaired script if the owner wants m2's old-entry guard observed.

### V: the one verification run allowed by VAL-115 (declared before sending)

- Purpose: verify the implementation audit's repairs (I01 to I07) and VAL-115's changes (O1 to O3) on the merge with
  staging. Commit 8fa64739 (merge of 7d90a4c8 into the branch at 67da0230); base 7d90a4c8 (the staging head merged).
  Scope: R2's 13 files base against head, one shard per side, R CMD check both sides; scripts U3, U14, base-head (U7 (b)
  and U8's store, base 7d90a4c8), gate-off, and X's repaired mutations script with the base (m2 on the old entries).
  Ceiling 420 job-minutes; no further correction without the owner (VAL-115 (1)); VAL-54 for mechanics only.
- What each new or changed test must show:
  - U10: the floor as "1/1,001", "9/1,001", "1/10,001" (k_min 4), "2/1,001" (k_min 4), "11/1,001" with k_min exactly 60
    at m = 273, and the `_none` sentence with "51/1,001"; the permutation note names "10,000" when chosen and "1,000" by
    default; every offered test's note carries the paired-samples caution; the lipid strength sentence reads "simulated
    metabolite data (no lipid data were simulated), on factors".
  - U9: the permutation's Methods sentence followed by its floor clause with "1/1,001"; the `_none` clause with
    "51/1,001"; the limits clause for a gene-and-protein gsea entry names "features that vary together with them in up to
    60.1 %"; a KIT-183-made block entry (test = 12, declared gsea) prints 34.4 % and an unknown test prints nothing;
    the power rows' `applies_to` carries "(expected, not measured".
  - U6: both refusals name the permutation test; a pre-ranked product's summary is not read as correlated.
  - U12: the refusal reads "about 121 seconds" and "longer than the 120 seconds allowed".
  - T10: the score column's exact text; T10/T11 strength with "metabolite data".
  - Unchanged tests: same as base or as R2c; U14 120/120 and 58/58 (permutation notes read with "1,000"); U3, U7 (b),
    gate-off pass; X: head passes, m1 fails U7, m2 fails U8 and the old entries' replay, m3 fails U5.
- A failure that is a result stops the ticket and goes to the owner.

### O3 wording items, held for the owner (file:line at the repair commit; the auditor's proposed text)

| Item | Where | Proposed |
|---|---|---|
| I05 refusal texts | `R/utils_integration_enrichment.R:85-91` | applied as the register gives them (the owner may still edit) |
| I10 "about 2.0 minutes ... longer than the 2 minutes allowed" | `R/utils_integration_enrichment.R:94-98` (`integration_enrichment_test_too_large`) | wording for the owner; no text proposed beyond noting the reading |
| I12 permutation notes: "1,000 random reorderings" also when 10,000 is chosen | `R/utils_integration_gui.R:700-701`, `:713-714` | "with what 1,000 (or the chosen 10,000) random reorderings" |
| I12 camera and gsea notes lack the independent-samples assumption | `R/utils_integration_gui.R:657`, `:669`, `:725`, `:736` | add the sentence the permutation notes carry ("samples taken in pairs, repeatedly from one subject, or in batches" break the assumption of independent samples) |
| I12 lipid wording "simulated metabolite (no lipid data were simulated) data" | `R/utils_integration_gui.R:765` (`enrich_strength_kind_lipids`) | owner's wording |
| I12 camera notes attach the C4b figure to "terms made only of the most variable features"; C4 on weaker factors reaches 10.9 % (genes and proteins) and 5.4 % (metabolites) | `R/utils_integration_gui.R:657-668`, `:669-679` | owner's wording |
| I12 `applies_to` of G, H, I: "every test finds fewer, the permutation test most" is inferred | `R/utils_integration_enrichment.R:1020-1022` | add "(expected, not measured)" |
| I12 camera metabolite note "0.4 to 4.5 %" covers C0 to C2; C3 reaches 4.6 % | `R/utils_integration_gui.R:669-679` | owner's wording |

Ticket `tickets/kit/KIT-191.md` revision 1.3 (1.2: the re-check's minors B01 to B06; 1.3: AUDIT-KIT-193 U1 and U2, by
the coordinator's instruction). Owner decisions: VAL-85 (2), VAL-97 (D1 to D9 as recommended). Branch
`kit-191-enrichment-test-options` from `v0.9-staging` 147adde8. Runs in `runs.tsv`.

## What was run, at which commit

| Step | Commit | What |
|---|---|---|
| R1 | ee5a33fb (tests only; R/ = 147adde8) | the five files, failing first |
| L1 | 53290941 (fix) | the permutation's time and memory (`lifecycle.md`) |
| R2 | 2f54469e (fix + L1 coefficients + U1/U2 corrections) | 13 files base against head, check both sides, scripts U3, U14, U7 (b) with U8's store, gate-off |
| R2c | 677fb251 (the one correction) | the same scope |
| X | 677fb251 | mutations m1, m2, m3 against U7, U8, U5 |
| B | run-k191-bundle 1141f52d (= 677fb251 + the bundle lane) | the new file on the installed bundle, ubuntu-22.04 and windows-2025; B24 ubuntu-24.04 |
| F | 677fb251 against 147adde8 | full suite with check, both sides |

## Observed

- **R1 (45284547817).** The new file: 14 passed, 10 failed, 11 errors; every new test fails, each because the
  behaviour is absent (`camera_floor`/`permutation` refused as not a test, `permutations` not a setting, texts and
  functions absent, the coefficients absent). KIT-183's restated T10, T11 and T15 fail at base by the texts and the
  RNG kinds this ticket changes. The other three files pass (unchanged).
- **L1 (45285106025).** Six cases completed, largest increase 561 MB (stop rule 1 GB not reached); the invalid call
  refused. Coefficients a = 5.12e-9, c = 2.13e-7, s = 5.48e-8 (an integer overflow in the script's coefficient
  arithmetic, an instrument defect, was worked from the measured table without a re-run; `lifecycle.md`).
- **R2 (45288809984).** One result caused by this change: KIT-183 T10 new-failure. Cause: the record block read an
  entry's test with `value$test`, which on a pre-ranked product partially matches its `tests` (the count), so the
  pre-ranked entry's strength sentence carried camera's figure. Two script instrument defects (VAL-54): U3 looked for
  the comparator at the script's top level (it sits inside `local({...})`), base-head loaded without the test helpers.
  U14 and gate-off passed. The new file 232 passed, 1 skipped (the cross-platform reference not yet recorded); its
  log gave the reference.
- **R2c (45294896497), the one correction pass.** 13 files: 11 `same`, 1 `count-changed` (KIT-183's file, its
  restated tests) and 1 `only-on-head` (the new file: 237 passed, 0 failed, 0 skipped); no new failure; KIT-183's file 180 passed on head (178 on base: the restated tests' extra
  expectations). Check: 0 ERROR, 3 WARNING, 3 NOTE on both, nothing only on head. Scripts:
  - U3: t within 1e-10; p-values of the product and the committed `simulate-k183.R` comparator identical on all 90
    sets (n = 6, 12, 20); no draw listed; no draw left the factor unchanged in any of the three fixtures (none at 6 samples either,
    where about 1.4 were expected), so U3 did not exercise the unchanged-draw rule; U2's 4-sample case tests it against
    the loops (`u3/`).
  - U14: 120 of 120 evidence cells and 58 of 58 text figures regenerated equal (`u14/`).
  - U7 (b) and U8's store: camera's and gsea's terms identical (SHA-256 equal) base against head on T1's fixture and
    the corpus store; the store made at base replays at head with each entry's test, terms identical; the camera
    product's identity equal in every component but `software` (`u7b/`).
  - gate-off: pages identical to base (`gate-off/`).
- **X (45294927461).** Head passes; m1 fails U7, m2 fails U8, m3 fails U5 (`x/`).
- **B (45294928896), B24 (45297583402).** The new file on the installed bundle: 237 passed, 0 failed, 0 skipped on
  ubuntu-22.04, windows-2025 and ubuntu-24.04; the cross-platform test within 1e-6 of the Linux reference. Bundles
  built and discarded inside their jobs (the lane uploads nothing).
- **F.** See "Full suite" below.

## The coordinator's two findings (AUDIT-KIT-193 U1, U2), during implementation

- **U1, the BH-floor sentence.** As shipped at 53290941 it read "at least k_min must have a p-value of 0.05 k_min / m
  or less before any can be called". Checked by hand on the register's counterexample: m = 100, B = 1,000, d = 0,
  p_min = 1/1001, k_min = 2; four terms at 2/1001 and 96 at 1; BH's step-up index is 4 because 2/1001 = 0.001998 is at
  most 0.05 x 4/100 = 0.002, so four are called, and no p-value is at or below 0.05 x 2/100 = 0.001. The sentence is
  false. True form (BH calls k* terms with p_(k*) <= 0.05 k*/m and p_(k*) >= p_min, so k* >= m p_min / 0.05): "none can
  be called at a 5 % false discovery rate unless at least k_min are called together". Shipped at 2f54469e (text,
  formatter, roxygen); U10 asserts the new sentence and the counterexample (BH calls 4, none at or below 0.001, 4 >= 2);
  the 40-of-2,000 example holds in the new form. Ticket 5.2 step 8, 5.4 and History (revision 1.3) corrected.
- **U2, reflection ties.** What the implementation at 53290941 did: draws leaving f unchanged were counted explicitly;
  every other draw was compared as |t*| >= |t| in floating point, so a draw reflecting a symmetric factor about its
  mean (t* = -t exactly) was counted or not by rounding. Now (2f54469e) a draw is counted when |t*| >= |t| - 1e-9
  (`INTEGRATION_ENRICHMENT_PERMUTATION_TIE`). Why this option: the rounding error of t is of order 1e-13 at the
  package's sizes, so no exact tie is lost, whatever its kind (reflections, or any other symmetry); a draw within 1e-9
  that is not an exact tie is counted too, which only makes p larger, so (1 + b)/(1 + B) stays valid (it never
  under-counts). Enumerating reflections exactly would need a floating comparison of 2 mean(f) - f with f[perm] and
  covers one kind of tie only. The canonical text gains `ties=abs_within_1e-9`; U3's rule (differences only on unchanged
  draws or within 1e-9) is the same tolerance, and U3 found no listed draw. Test: U2's symmetric fixture f = (1, 2, 3,
  4), 1,000 draws: reversals present, every p >= (1 + d + r)/1001, p identical to the loops' reference.
  `min_p` = (1 + d)/(1 + B) stays a valid lower bound (reflections can only raise p).

## Acceptance criteria

| # | Criterion | State | Evidence |
|---|---|---|---|
| A1 | camera_floor and permutation equal direct computations; permutation equals the simulation's comparator | met | U1, U2 (R2c); U3 `u3/` |
| A2 | RNG independent of the session and restored; seed in the identity | met | U4, U5 (R2c) |
| A3 | resolver, refusals, STRING; KIT-184 combinations | met (KIT-184 not merged: 8.2 is its own step 0) | U6 (R2c) |
| A4 | camera and gsea unchanged | met | U7 (a) (R2c); `u7b/` identical SHA-256, identity but software; X m1; F |
| A5 | keys, entries, replay of old and new entries | met | U5/U8 (R2c); `u7b/replay.tsv`; X m2 |
| A6 | record: evidence sheet, Methods (no format code), limits, strength, several tests | met | U9 (R2c) |
| A7 | tab: four choices, camera, the kind's notes and evidence, seed, few-samples, floor | met (seed visibility asserted on the panel's condition, not in a browser) | U10 (R2c) |
| A8 | workbook, Rd, starter | met | U11 (R2c); Rd and starter by hand |
| A9 | guard and L1 | met | U12 (R2c); `lifecycle.md` |
| A10 | corpus positive control | met | U13 (R2c) |
| A11 | every figure equals the regenerated table | met | U14 `u14/` |
| A12 | failing first; guards fail on mutations | met (m2 is caught by the replay of the head's new entries in U5/U8; the branch-point store's old entries were replayed unmutated in U7 (b), `u7b/replay.tsv`; the mutant logs do not name the failing expectation. The repaired X script also replays the old entries under m2: unrun) | R1; X `x/` |
| A13 | cross-platform | met | B, B24, `cross-platform.tsv` |
| A14 | full suite base against head with check | met | F 45300767142 |
| A15 | ledger row, NEWS under the unreleased heading, no version change, wording list | met | `docs/methods-decision-ledger.md` (VAL-<merge>), `NEWS.md`, `wording.md` |

## Full suite

F (45300767142, 677fb251 against 147adde8, 20 shards per side, check both): 456 files; 454 `same`, 1 `only-on-head`
(the new file, 237 passed), 1 `count-changed` (KIT-183's file: T15 renamed and restated, its counts moved from the old
name to the new); 0 `new-failure`. Every failing expectation on head fails identically on base (stored laptop oracles,
the KIT-156 record). R CMD check: 0 ERROR, 3 WARNING, 3 NOTE on both, nothing only on head. A14 met.

Core-hours: 1,382.91 job-minutes = 92.19 core-hours, against 88.25 expected and 231.00 ceiling (the correction pass R2c
used 8.06 of its 28.00).

## Inferred, not observed

- The cost prediction over-states every measured case (`lifecycle.md`); a user's computer may be slower or faster.
- D9's fix is inferred to leave default-session results unchanged from U7 (b) (gsea terms identical base against
  head) and T15 (value identical at seed 11); whether gsea at base differed under non-default kinds was not observed
  (U4 errored at base on the permutation call before its gsea part).

## Departures from the ticket's text (recorded)

- KIT-183's T10 and T11 restated besides T15: their assertions named texts this ticket changes by design (the score
  column, the limits paragraph, the strength sentence per test and kind, the seed label, the notes and few-samples per
  kind). The record block's test read exactly (`value[["test", exact = TRUE]]`): a latent KIT-183 defect this ticket's
  strength sentence exposed (R2).
- Notes' power figures follow U14's rule (one decimal below 10 %): "2.0 to 2.2 %", "1.0 %", "1.0 to 1.3 %" where 5.4's
  indicative text had "2 %", "1 %", "1 %". camera's and gsea's metabolite notes carry the metabolite figures (KIT-183's
  text with the kind's numbers), so no kind is shown another kind's figures.
- Record sheet named `integration_enrichment_tests` (the record's prefix; ticket: `enrichment_tests`), written only
  when a KEGG, Reactome or GO entry exists (`INTEGRATION_RECORD_EVIDENCE_SHEET`; one line in the workbook writer,
  `R/utils_workflow_analysis_record.R`).
- The scripted message is given for every scripted permutation run (VAL-97 (4) read literally), through `rlang::inform`.
- `integrationSettings()$enrichment` names the choices in its `default` text (its columns are fixed by KIT-154's test).
- U2's comparison of `effect` is within 1e-12 (the product computes it from sums), not `identical()`.
- The symmetric-factor case of U2 was added after R1 (finding U2); on base it fails because the function is absent.

## Not checked

- The app in a browser (the seed and reorderings panels' visibility; the tab's responsiveness for up to 120 s).
- Gene Ontology (KIT-181 not merged); KIT-184's combinations (not merged; 8.2 is its step 0).
- Power of `camera_floor` with removed factors; the permutation below 6 samples, on lipid data, on real data.
- Whether fgsea draws from R's RNG (F1's premise): U4 shows gsea invariant at head; base behaviour under other kinds was
  not reached in R1.
- The SHA-256 of the bundle platforms' tables (the lane uploads nothing).
