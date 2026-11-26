# KIT-260 second correction (VAL-276): predictions, written before the run

Run: the new file plus the same 130 files, base staging abfc99ac, head = kit-260 after 35496462 and these records
(the merge of abfc99ac included), R CMD check both sides.

- `test-kit260-itsd-default-pattern.R`: only on head; every test passes, I05 included (the older record replays with no
  mismatch under the marker's earlier token) and its mirror (a new record replays exactly; its ITSD product carries the
  rows used).
- `test-omic-norm-module.R`: the 20 keys VAL-275 confirmed pass on head; **the `provenance` keys of the lipid and
  metabolomics manual and automatic captures fail on head** (the ITSD step records the fixed marker `itsd_outcome`):
  up to 14 expectations, not approved; they are listed for the owner by file and not re-pinned.
- Any other file that pins the ITSD step's recorded parameters or step provenance (for example analysis-record or
  session captures) may fail on head for the same reason: listed, not re-pinned.
- `test-omics-norm-core.R`: `count-changed` (one more passing expectation), as in R2b.
- `test-integration-gui.R`, the two pins, `test-workflow-omics-import.R`: `same` (VAL-275 allowed them; base and head
  each pass their own expectations).
- Every other file: `same`; the files failing on staging fail identically.
- Check: findings identical on both sides.

If any expectation outside VAL-275's list moves, the work stops there and the list goes to the coordinator.
