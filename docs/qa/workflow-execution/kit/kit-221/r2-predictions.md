# KIT-221 R2: purpose, binding, price, predictions (written before the request is sent)

**Purpose.** Fix commit against the tests-only side: the 45 files of the regression scope plus R CMD check on both
sides, and the head script job (G1 and G3 head digests). Decides stop rules 3 to 5.

**Binding.** base_commit = e80e558d68709b75556b02a6b8e383652dfa3a30 (the tests-only side: nothing under `R/` differs
from base 176f9d3ff; the new test file is the same on both sides). package_commit = the tip of `kit-221` when the
request is written (4654fd600 plus NEWS and records; nothing under `R/` after 4654fd600). 3 shards per side, check
on both sides, `sva`, `blas: reference`, `test-gui-entry-proteomics.R` timeout 50 min, one script (g1.R, head).

**Price (coordinator's ruling: pro rata from the ticket's row, 45 files against the ticket's 40).** The ticket's R2
row: 8 jobs, ceiling 630 min = 42.0 core-hours, expected 28.0 (about 36 to 40 files). Scaled by 45/36: expected
35.0 core-hours; ceiling 787.5 min = 52.5 core-hours. Plus the head script job: expected 1.0, ceiling 40 min = 2.67.
R2 total: expected 36.0, ceiling 55.2 core-hours. Instrument rerun 1 (g1 base digests, run 45715339444): expected 1.0.

## Predictions

- `test-kit221-da-several-comparisons.R`: base 43 passed / 38 failed / 0 errors (R1's observation); head 124 passed /
  0 failed / 0 errors / 0 skipped (W1 10, W1b 10, W2 14, W3 17, W4 24, W5 18 with its guard, W7 14, W8 15, W9 2).
  Verdict `fixed`.
- Each of the other 44 files: verdict `same` (per-file counts and masked failure messages equal on both sides),
  including `test-workflow-saved-da.R` failing the same 171 expectations on both sides (stored references recorded on
  the laptop) and `test-kit218-da-wrapper-tail.R` 80 passed / 2 skipped on both sides.
- R CMD check: nothing only on head; `check-compare.txt` identical.
- g1.R head: exit 0; every digest of `g1-head.tsv` equal to `g1-base.tsv` row for row (cases limma, limpa, ruv,
  kernel, alias; same fields).
