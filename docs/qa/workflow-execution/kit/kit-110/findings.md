# KIT-110 findings

Append only. Ids are F110-NN. A fix that does not change the methodology is made here, with a test that fails first. Anything else is parked.

| Id | Where | What | Status |
|---|---|---|---|
| (instrumentation) | `kit-105/lifecycle/lifecycle.R` | The lead's `ListReporter` version of the lifecycle harness fails before a lane starts: the reporter had no current file. Found by W2. | Repaired by W2 (`start_file`/`end_file`); committed `be0a2601`, with a correction note in `kit-105/entry-suite.md`. Not a product defect. |
| (test harness note) | `test-trx-norm-module.R`, the managed cases | Under `omicsManagedChainMocks` the in-test worker runs in the test process. A managed run after a refused run in the same session then fails with "native reference roots contain a live payload". W1 gives the managed refusal its own session. A real worker process was not checked. | Recorded, not claimed as a product defect. |
