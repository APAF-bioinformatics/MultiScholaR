# KIT-262 M3: the continuation rule and the timeouts (written with the scripts; declared again with each counted job)

Applies the audit's M11 (AUDIT-KIT-262-m1m2, "Timing") and the lessons of KIT-275's RH and RJ checks. It binds the scripts as drafted here; the counted declaration repeats it with the job lists.

## 1. What a unit is, and when it is complete

A unit is one data set: a cell and a replicate (tag `cNN-rR`, `NN` the cell index 1 to 10, `R` the replicate), with all its arms (the P journey; the C1 Manual journey where C1's choice differs; the checks). `m3-k262.R` appends ONE line to the job's `summary.tsv` when the unit ends, and nothing before: a unit is complete when its line is there, whole. The line carries status `ok` (a result), `error` (a journey of the product stopped, class and message kept, a trace file written: an errored data set, never a denominator, never re-run) or `instrument` (the guard fired, the script failed, or the worker returned nothing: a stop for the scorer). A unit that was running when the job was killed has no line, or a short last line; the scorer drops a short line (`score-m3-k262.R`, `readSource`).

## 2. The rule

1. **Missing units run once, by tag.** After a counted job ends (completed or timed out), `Rscript docs/qa/workflow-execution/kit/kit-262/score-m3-k262.R missing <plan.tsv> <the job's source dir(s)>` (counts only; it reads cell, r and purpose, never a value) prints the planned units that no source holds a complete line for, as `K262_M3_ONLY_UNITS=<tags>`. If the list is empty there is no continuation. Otherwise ONE continuation request per job, declared before it is sent (the tag list, the plan and the timeout below), runs `m3-k262.R` with exactly `K262_M3_ONLY_UNITS=<tags>`, the same cells and replicate settings and the same package commit, into its own output directory. `m3-k262.R` stops when a tag is outside the job's plan, when the setting is empty, or when it is used in the lifecycle.
2. **No partial re-run.** A unit with a line, whether `ok` or `error`, is never run again. An unfinished unit is discarded whole, not repaired from what it had computed (it computed nothing that was written). The continuation runs no unit that any source holds complete. A continuation that is itself interrupted is not continued again: its missing units are reported "did not complete" and their cell is "not determined (did not complete)".
3. **Each data set is scored once.** The scorer is run once, on the union of the jobs and their continuations, each in its own source directory. A unit complete in two sources, or twice in one, stops the scorer (`m3_unit_in_two_sources`, `m3_duplicate_key`); a row whose replicate is outside the planned set, a cell outside the ten, a seed that is not the declared formula, or a header that is not the row format also stop it.
4. **Nothing is re-run for what it showed.** An unfavourable or inconclusive valid result is recorded and never re-run (VAL-54). An error row is a recorded outcome of the product (VAL-242: a failed computation stops), not an instrument failure.
5. **Instrument rows stop the scorer** and are handled by VAL-54 (a script or test mechanics failure with no result: fixed and re-run only under its four conditions, at most three, in the instrument reserve). The continuation is not a way round that: it never runs a unit that has an instrument row.

## 3. Timeouts

Every `timeout_minutes` is written in the request, for every script, and is at least two times that script's plan:

- Counted job: plan = (L3's seconds per data set of the job's cells, summed over its units) x 1.93 / 4 workers, in minutes; `timeout_minutes` = 2 x plan, rounded up to 5. 1.93 is the largest ratio of a counted job's sweep seconds to its lifecycle's seen in M2 (AUDIT-KIT-262-m1m2 Q3, M11; the re-measurement's jobs were 0.79 to 1.42 of a plan made at 1.0). Jobs are cut so that the plan is at most 150 minutes (the cell and replicate ranges are split by seeds, the same R); the timeout is then at most 300 minutes. If the cut would need more jobs than the ceiling allows, the plan goes to the coordinator (G2); R is never reduced.
- Continuation: plan from the same seconds for the tags listed; timeout 2 x plan, at least 30 minutes.
- Lifecycle scripts and scorer runs: the figures in `m3-lifecycle-declaration.md`; the scorer, the blind scorer and the self-tests carry their own explicit timeouts.

## 4. Where the continuation is paid

The continuation, if any, is charged to the correction reserve (ticket section 11: 13.33 core-hours, 4.56 used by M2's continuation, 8.77 left on paper) and counted in the running ceiling line against the ticket's 250.32 (VAL-390), not to M3's own line. If the reserve is not enough, the cells it would have filled stay "did not complete" and the coordinator decides; nothing is shortened.
