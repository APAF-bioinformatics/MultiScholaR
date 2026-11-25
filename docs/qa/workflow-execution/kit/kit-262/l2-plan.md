# KIT-262 L2: the M2 lifecycle (plan, committed before the run)

Purpose: prove the instrument observes the intended work through completion and failure, blind (VAL-11): nothing result-bearing leaves it (row files under tempdir(); only checks.tsv, timing.tsv, hashes.tsv and the console counts are written). Binding: the tip of branch kit-262 holding the scripts (base 15d86bd4). Ceiling: ticket row L2, timeout 40 script-minutes, expected 25 job-minutes = 1.67 core-hours, ceiling 45 job-minutes = 3.00 core-hours. Attempts: once; instrument defects under VAL-54 (at most three further runs).

Scripts, in order: (1) e5b-k262.R (E5b); (2) score-k262.R selftest (E5 M2 part, E6 on constructed inputs); (3) m2-k262.R, K262_LIFECYCLE=1, all 198 cells at r = 91, E1' on every data set, the E8 cells also as metabolomics at r = 91 (E1, E1', E2, E3, E4, E7, E8, COMPLETION, E6 on the lifecycle's own rows); (4) m2-k262.R on the ten E8 cells at r = 92 to 95 (E8, with E1'); (5) the failure check: cell 1, r = 91 and 92, K262_FAULT=1 (the first task signals after its data are generated; the error row must be written and be the only one).
Local, before the run: e0-k262.py (E0) and e5hand-k262.py (hand values) pass.
Pass: every row of every checks.tsv has status pass. Decides only whether the instrument may be used; G2 recomputes the counted jobs' plan from timing.tsv.
