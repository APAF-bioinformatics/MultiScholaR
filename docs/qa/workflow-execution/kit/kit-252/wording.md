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
| WP-1 | `memory_guard` | in the application, the kernel stopped the worker at its own memory cap, the smaller of 16 GB less what the application holds (and a 1 GB reserve) and the computer's free memory less 3 GB (`:686-696`, `:1478`, `:2041`; systemd only). The application's own stop at 15.5 GB for the application and its step ends the whole application under systemd (`:1655`, `:1714-1718`), and the plain launcher has no cap. Nothing of the step is committed | "the step needed more memory than the application allows for one step (at most 16 GB, less when the computer has less free memory), so it was stopped before it finished. Its result was not saved. Closing other programs may help; if it stops again, this data set is too large for this step on this computer." |
| WP-2 | `host_headroom` | not shown in the application: under systemd a stop below 2 GB free ends the whole application (`:1656`, `:1714-1718`; the shell then says 'The analysis engine stopped'); the plain launcher has no monitor. Kept for completeness of the mapping. | "the computer's free memory fell below 2 GB, so the step was stopped to protect the computer. Its result was not saved. Close other programs and press the button again." |
| WP-3 | `worker_failed`, no error of its own, worker `Result` "timeout" | the worker ran past the step's time limit: 1,800 s from the server's policy (`:677`) unless the step has its own (`workflowActionWorkerDeadline`); `<N>` is that limit in minutes (30 for every step today). The same data reach it again | "the step reached its time limit of <N> minutes and was stopped. Its result was not saved. Running it again on the same data will stop at the same limit." |
| WP-4 | `worker_failed`, no error of its own, otherwise | the worker process ended with no result: killed (KIT-233 run L3w), crashed, or killed for memory under the plain launcher (no monitor there) | "the step's process stopped before it finished (for example, the computer ran out of memory or the process was closed). Its result was not saved. Press the button again; if it stops again, the application's log has the details." |
| WP-5 | `admission_refused`, `monitor_start_timeout` | the step's memory cap would be below 768 MB (`:1951`), the worker could not be launched (`:1974`), or the monitor did not start in 30 s (`:1943`); nothing ran | "the step could not be started (too little free memory, or its process could not start). It did not run. Close other programs and press the button again." |
| WP-6 | `publication_pending` | the result was accepted and committed (`:2072-2073`); its files were not written (`:2078-2083`); pressing the same button retries the publication without recomputing (`:919-931`) | "the result was saved, but its tables and figures could not be written to the project folder. Press the button again to write them; the step is not run again." |
| WP-7 | `stale`, `project_replaced` | the head moved, or another project was opened, while the step ran (`:2063`, `:1923`, `:2030`) | "the table or the project changed while the step ran, so its result was not used. Press the button again." |
| WP-8 | `monitor_lost`, `invalid_result`, `idle`, the six guard reasons, any other status | the monitor stopped, the result was over 8 MB or unreadable, no job was found, or a guard check failed (systemd only); nothing of the step is committed. The guard reasons and `monitor_lost` end the application under systemd in the same way; in the application WP-8 is met for `invalid_result` (a result over 8 MB or not valid) and `idle`. | "the step's result could not be received (an internal check failed). Its result was not saved. Press the button again; if it happens again, the application's log has the details." |
| WP-9 | `cleanup_failed`, `duplicate` | the worker or monitor unit could not be emptied (`:2196-2199`), which can follow a committed and published result (`:2083`); or a completion of a request already complete (`:2052`) | "the step's process could not be closed cleanly. Its result may have been saved; check this tab before running the step again, and restart the application if this happens again." |

"Its result was not saved." is said of the action that failed. In a chain (the normalisation of every omic), the
actions delivered before it keep their commits (test G3: the lipidomics head stays at the normalised state when the
RUV action ends `memory_guard`).

Line numbers in the cause cells are at 16112917, except WP-1's, WP-2's and WP-8's added sentence, which are
AUDIT-KIT-252-impl I01's at the branch head (the base's lines + 45 before the step's observer, + 47 after it).

"The application's log" (WP-4, WP-8) is the Log Output panel at the bottom of the window (in memory, the last 1,000
lines; not the logs folder). For a failed step it holds 'managed action <operation> ended: <status>' and the notice,
not more (AUDIT-KIT-252-impl I03).

## Every place the sentences appear (AUDIT-KIT-252-impl I04)

Besides "Error applying <step>: " the same sentence is shown after these frames, so the owner words each sentence to
hold in every one of them (no code change by KIT-252):

| Frame | Where | Note |
|---|---|---|
| "Error in normalization: " (second notice, proteomics) | `R/mod_prot_norm_workflow_helpers.R:447-460` | F1 |
| "Error in DA analysis: " (second notice, proteomics DA runner) | `R/mod_prot_da_handlers_run_analysis.R:54-61`, `:125-128` | |
| "Analysis error: " (omics DA) | `R/mod_omic_da.R:425`, `:446-450` | |
| "Warning: Could not save results to disk: " (omics DA, the publication action failing after the tab's rows arrived) | `R/utils_workflow_omics_actions.R:610`, `:628` | with WP-6 it reads "Could not save results to disk: the result was saved, ..." |
| "Error: " in the final-QC PCA plot (lipid/metab) | `R/mod_omic_norm_panels.R:1344` | a display: "Press the button again" names a button the display does not have; it asks again when next drawn (`R/utils_workflow_omics_actions.R:131-139`) |
| the volcano page's error (lipid/metab DA) | `R/utils_workflow_omics_actions.R:164` | as above |
| "ERROR: " in the tab's log; "Error: " notice under the RUV design gate | `R/mod_omic_norm.R:471-475` | |

## Open for the owner's wording pass (QC-8; none changed now)

- WP-1 (I06): "this data set is too large for this step on this computer": when the cap was the application's 16 GB
  limit (`source = "slice"`, written to `admitted.rds`), a larger computer gives the same cap.
- WP-3 (I06): "will stop at the same limit" is a prediction (a busy computer can make the difference); no setting can
  raise the limit from the window, so the sentence offers no way forward.
- WP-5 (I06): "Close other programs" does not help a launch failure, or a cap made small by the application's own
  memory.
- WP-6 (I06): "the result was saved" against "could not be written" is a fine distinction for non-specialists; the
  retry works only with the same settings (other settings are refused with "retry pending publication with its
  original parameters first"). Whether a second press of the omics DA button reaches the retry was not established.
- WP-7 (I06): it does not say the data are fine; `stale` also arises from the user's own Reset, QC Revert or Save
  Design (KIT-253 G0), where pressing again may not be wanted; the poll's `project_replaced` means the project's store
  failed its identity check (replaced, removed or corrupt), where pressing again cannot work.
- WP-9 (I05): "check this tab": after `cleanup_failed` following a committed result the step's outputs never show
  that result (its delivery did not run); the place that shows whether it was kept is the step history (the current
  state), not the step's own outputs.

Out of scope, for the same pass (ticket findings F1, F2): the second notice of the proteomics normalisation tab with
the generic label "proteomics step"; the controller's own messages "project was replaced", "result was superseded",
"browser acknowledgement timed out" and the condition messages shown through `fail(action, error)` (for example
"retry pending publication with its original parameters first"). Of these, the memory refusal a user is most likely to
meet is "insufficient host headroom" (`R/utils_workflow_execution.R:1866-1867` at the branch head: free memory under
4 GB at any press, both launchers; AUDIT-KIT-252-impl I08). A user who runs low on memory during a step under systemd
does not see WP-1 or WP-2: the application closes ("The analysis engine stopped (exit ...)"; I01, I09).
