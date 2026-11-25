# KIT-252: the failure sentences, DRAFTS for the owner's wording pass (VAL-41, row QC-8)

Each sentence follows "Error applying <step>: " in the error notice of a step run on the managed worker (every omic:
the proteomics QC sub-steps, normalisation, correlation filter and DA publication; the lipidomics, metabolomics and
transcriptomics tab actions). The proteomics normalisation tab shows it a second time after "Error in normalization: "
(finding F1, out of scope). It is shown only when the step's worker did not write an error of its own; a worker's own
error is shown as it is (GUI audit A111-03). The internal status is written to the application's log as
"managed action <operation> ended: <status>".

None of these words is final. The code holds them in one place, `workflowManagedFailureText`
(`R/utils_workflow_execution.R`); the tests pin them in `tests/testthat/test-kit252-managed-failure-text.R`
(`kit252Text`) and in the five declared edits (`test-prot-norm-core.R`, `test-workflow-execution.R`,
`test-omic-da-module.R`, `test-workflow-omics-application-reads.R`). A change of words changes those places only.

| Id | Shown for (internal status) | What is true of the cause (code at v0.9-staging 16112917) | Draft |
|---|---|---|---|
| WP-1 | `memory_guard` | the monitor stopped the step when the application and the step together reached 15.5 GB (`R/utils_workflow_execution.R:675`, `:1608`, `:1694`), or the kernel killed the worker at its cap, the smaller of 16 GB less what the application holds (and a 1 GB reserve) and the computer's free memory less 3 GB (`:686-696`, `:1994`); nothing of the step is committed | "the step needed more memory than the application allows for one step (at most 16 GB, less when the computer has less free memory), so it was stopped before it finished. Its result was not saved. Closing other programs may help; if it stops again, this data set is too large for this step on this computer." |
| WP-2 | `host_headroom` | the computer's free memory fell below 2 GB (`:676`, `:1609`, `:1695`) | "the computer's free memory fell below 2 GB, so the step was stopped to protect the computer. Its result was not saved. Close other programs and press the button again." |
| WP-3 | `worker_failed`, no error of its own, worker `Result` "timeout" | the worker ran past the step's time limit: 1,800 s from the server's policy (`:677`) unless the step has its own (`workflowActionWorkerDeadline`); `<N>` is that limit in minutes (30 for every step today). The same data reach it again | "the step reached its time limit of <N> minutes and was stopped. Its result was not saved. Running it again on the same data will stop at the same limit." |
| WP-4 | `worker_failed`, no error of its own, otherwise | the worker process ended with no result: killed (KIT-233 run L3w), crashed, or killed for memory under the plain launcher (no monitor there) | "the step's process stopped before it finished (for example, the computer ran out of memory or the process was closed). Its result was not saved. Press the button again; if it stops again, the application's log has the details." |
| WP-5 | `admission_refused`, `monitor_start_timeout` | the step's memory cap would be below 768 MB (`:1951`), the worker could not be launched (`:1974`), or the monitor did not start in 30 s (`:1943`); nothing ran | "the step could not be started (too little free memory, or its process could not start). Nothing was run. Close other programs and press the button again." |
| WP-6 | `publication_pending` | the result was accepted and committed (`:2072-2073`); its files were not written (`:2078-2083`); pressing the same button retries the publication without recomputing (`:919-931`) | "the result was saved, but its tables and figures could not be written to the project folder. Press the button again to write them; the step is not run again." |
| WP-7 | `stale`, `project_replaced` | the head moved, or another project was opened, while the step ran (`:2063`, `:1923`, `:2030`) | "the table or the project changed while the step ran, so its result was not used. Press the button again." |
| WP-8 | `monitor_lost`, `invalid_result`, `idle`, the six guard reasons, any other status | the monitor stopped, the result was over 8 MB or unreadable, no job was found, or a guard check failed (systemd only); nothing of the step is committed | "the step's result could not be received (an internal check failed). Its result was not saved. Press the button again; if it happens again, the application's log file has the details." |
| WP-9 | `cleanup_failed`, `duplicate` | the worker or monitor unit could not be emptied (`:2196-2199`), which can follow a committed and published result (`:2083`); or a completion of a request already complete (`:2052`) | "the step's process could not be closed cleanly. Its result may have been saved; check this tab before running the step again, and restart the application if this happens again." |

"Its result was not saved." is said of the action that failed. In a chain (the normalisation of every omic), the
actions delivered before it keep their commits (test G3: the lipidomics head stays at the normalised state when the
RUV action ends `memory_guard`).

Out of scope, for the same pass (ticket findings F1, F2): the second notice of the proteomics normalisation tab with
the generic label "proteomics step"; the controller's own messages "project was replaced", "result was superseded",
"browser acknowledgement timed out" and the condition messages shown through `fail(action, error)` (for example
"retry pending publication with its original parameters first").
