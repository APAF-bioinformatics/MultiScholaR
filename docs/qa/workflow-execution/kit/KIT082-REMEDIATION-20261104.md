# KIT-082 remediation, 2026-11-04

User: "fix literally every issue before another run. write a task list against each one and do it systematically, then rerun after all the appropriate tests". The audit [KIT082-AUDIT-20261104.md](KIT082-AUDIT-20261104.md) stays as written; closure is tracked here. The user's instruction accepts the monitor trade-off in A-01 and the monitor's place in this ticket (A-12).

| Task | Finding | Work | Closure evidence | State |
|---|---|---|---|---|
| T01 | A-01 | One read decides; a test on real starting and stopping units; the old function shown failing on it | Stress result; execution tests; members | code and tests done: on 30 real units, aborts second-look 30, the original check-then-read 11, one-read 0 ([stress](kit-082/monitor-window-stress.txt)); members pending |
| T02 | A-02 | Wait for "Found N internal standards" with N above 0; fail on an error | Negative member (bad pattern); members | open |
| T03 | A-03 | Non-blank image by pixels; a data row in every table; download with data rows | Negative member (nothing significant); members | open |
| T04 | A-04 | A managed tracking-figure failure warns and finishes, as in process | Test on both omics | done: tests pass on both omics |
| T05 | A-05 | The ITSD page read guarded again | Test on both omics | done: tests pass on both omics |
| T06 | A-06 | The lipid download reads the tab's results at each click | Test with two runs | done: tests pass, two runs and before any run |
| T07 | A-07 | A failed tab read clears the tab's results | Test on both omics | done: test passes |
| T08 | A-08 | A dropped display asks again; a failed one retries on the next invalidation; retrying a publication keeps the new callbacks | Tests | done: tests pass: dropped, own failure once then retry, retry resumes, other work queues |
| T09 | A-09 | A failure branch in each open-ended wait | Negative members; driver review | open |
| T10 | A-10 | The final driver on `kit_081`, so control and candidate compare like for like | Members on `kit_081` | open |
| T11 | A-11 | KIT-084 declared in its ticket; builds recorded in the trial record; lipid trial and a managed member at size, after KIT-082 closes | Declaration; trial lines | done: declared in KIT-084 revision 1.1; trial record carries provenance; its runs wait for KIT-082 to close |
| T12 | A-12 | Monitor and other files named in the ticket; `kit-gui-082` on the final build | Ticket; member | open |
| T13 | A-13 | Plan, memory, links, allowlist, DA size tool and figures, evidence committed | Links and status checked | done: ticket, plan, memory, links, allowlist and tool done; figures 0.37 -> 0.08 MiB and 0.59 -> 0.12 MiB; evidence commits with this change |
| T14 | A-14 | The quiet-host wording corrected | Script comment; report | done: script comment corrected; reported to the user |
| T15 | Reviewer | The probe also bounds rows kept by a row read | Test | done: cells kept by each store read below the whole table's, in the journey probe |
| T16 | Reviewer | Probe comment corrected; the download's reads shown as the one exception | Test | done: download shown reading contrast products; nothing else does |
| T17 | Reviewer | Download writer: a zero-row piece; column types across pieces | Test | done: test passes: empty piece, type across pieces, no rows |
| T18 | Reviewer | ITSD pages with a multi-column annotation equal the whole table | Test, repaired if it differs | done: unreachable: `annotation_id_column` is validated as one string (R/func_lipid_s4_assay_data_class.R:89, R/func_metab_s4_objects.R:110), so the multi-column branch does not run through these sites |
| T19 | Reviewer | Trace every remaining `read*WorkflowState` site and the export's source read; assert the tracking `metrics` argument | Site list; test | done: every remaining site traced behind a managed branch or on the in-process path; the journey probe is the runtime check; metrics forwarding test passes |
