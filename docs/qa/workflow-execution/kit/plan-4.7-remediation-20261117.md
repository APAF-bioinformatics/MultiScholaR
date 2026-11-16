# Plan 4.7 audit: remediation record (2026-11-17)

The register is [plan-4.7-audit-20261117.md](plan-4.7-audit-20261117.md) and is left as written. Authority for the repairs: user, 2026-11-17, "1 yes"; on integration, "i need to source the workbook and verify it first so let me do that in the next few days - integration comes last anyway. GUI FOR THE INTEGRATION offered the same way as the others but import data would source data processed through MS and likely from the available data in the project folder structure that is currently in place for a project." The repairs are ticket and plan text (plan 4.8); no package code, no run.

| Finding | State | Where | Closure evidence |
|---|---|---|---|
| B-01 | **open** | KIT-134 "The comparison across machines"; plan 4.8 working rule | The laptop's exact comparison is the gate again. Closes when one full run on a Linux runner has been read and the comparison rule is written. |
| B-02 | repaired in text | KIT-134 P1, P5 | Lockfile from the verified library; verification from the installed artifact; the paper's numbers on the bundle. Closes at KIT-134's acceptance. |
| B-03 | repaired | KIT-134 new; KIT-050 revision 1.4; KIT-052 revision; KIT-132, KIT-133 revisions; index | KIT-132 and KIT-133 depend on KIT-134; KIT-050 and KIT-052 are the final qualification. Inventory: no cycle, no missing dependency, 132 tickets. |
| B-04 | **open, with the user** | KIT-060 revision 1.2 | Form decided (a GUI like the other omics, importing processed data from the project's folders). The authority is the workbook the user will source and verify; asked again at the 7 December checkpoint. The test corpus is specified in KIT-060 and built in KIT-062. |
| B-05 | **open, with the user** | KIT-134 "What the release build needs" | The pipeline runs without the two packages; the release candidate needs the repositories public or a token by 14 December. |
| M-01 | repaired | KIT-125 revision 1.1; KIT-126, KIT-127, KIT-129 dependencies | KIT-125 at kernel level beside KIT-126 to 128; it gates KIT-129. |
| M-02 | repaired | KIT-125 revision 1.1 | 12 replicates or more, distinct splits, one dataset one unit, a batch-aligned arm, the limit stated. |
| M-03 | repaired in text | KIT-125 revision 1.1; KIT-133 revision 1.1 | The survey comes first. Closes when it is done. |
| M-04 | repaired in text | KIT-124 revision 1.4; KIT-129 revision 1.1; KIT-126 revision 1.1 | Both readings of rule 4 at the gate; the verdict-to-option rule proposed for the user's confirmation; the warning's own check; refusal = no RUV with a notice. Closes at the gate. |
| M-05 | repaired in text | KIT-126 revision 1.1; KIT-061 revision; KIT-060 revision 1.2; KIT-129 revision 1.1 | The W0 table owns the consumers; the integration input is named by role; the templates' wording is KIT-129's. |
| M-06 | repaired in text | KIT-129 revision 1.1 (30 hours, to be declared in full at authoring) | Sequence with store comparisons, one harness revision, one final re-record. |
| M-07 | repaired in text | KIT-088 revision; KIT-133 revision 1.1 | The claim restated; measurement windows. Closes when the windows are dated at authoring. |
| M-08 | repaired | KIT-130 revision 1.1; the running spike was told | The resolve-only step is a piece of the outcome table. |
| M-09 | repaired | KIT-131 revision 1.1; KIT-134 P6; KIT-052 revision | Acceptance and fallback written. |
| M-10 | repaired | plan 4.8 row: 30 November, 7 December, 14 December | |
| M-11 | repaired in text | KIT-124 revision 1.4 | Closes in KIT-124's outcome. |
| m-01 | repaired | front matter and index agree for every ticket of the plan | inventory: `status_consistent` true for all. |
| m-02 | part | KIT-130, KIT-131, KIT-134 carry the sections | KIT-125 to 129, 132, 133 are stubs; the sections are added when each is authored in full. |
| m-03 | part | plan 4.8 row | KIT-091 and KIT-067 placed after the submission. KIT-114's status was not read: open. |
| m-04 | open | | With B-01. |
| m-05 | repaired | KIT-132 revision 1.1 | |
| m-06 | open | | When the plan is next edited. |

The independent read (`plan-4.7-independent-review-20261117.md`) had not reported when this record was written; its findings are added below when it does.

## Added by KIT-131's W0 read (2026-11-17, checked against the code)

| ID | Finding | Evidence | State |
|---|---|---|---|
| B-06 | **No DIA route runs on Windows today.** The precursor rollup, a step of DIA and DIA_LIMPA, is refused outside Linux x86_64 with 80-bit long double, and needs `fastduckdb`. The paper claims Windows. | `R/utils_workflow_dia_relational_rollup.R:11-13`; `src/Makevars` (`-Werror`, `-fPIC`, `uname`). | **open, with the user**: whether Windows must include the DIA routes at submission. The discovery run says first whether the native code builds there. |
| B-07 | **The shipped application would not use the worker at all.** Nothing in `R/` calls `launchWorkflowApplication`; `run_app()` runs every heavy step in process, on Linux too. The memory figures of the series are measured through the harness's launcher. | `grep launchWorkflowApplication R/`: its definition only. | repaired in text: KIT-134 revision 1.1, P3 launches through it. Closes at KIT-134's lifecycle test. |
| M-12 | The on-disk state owner compares `normalizePath()` output with `/`-joined paths: 39 calls, 9 containment checks. On Windows the worker's store would be refused. | `R/utils_workflow_saved_entities.R:266-267`, `:742-745` (W0 design, section 5). | repaired in text: KIT-131 revision 1.2, D4, bounded. |
| B-05 (sharpened) | `fastduckdb` is required by the DIA routes, not only preferred. CI and the bundle pipeline cannot verify a DIA lane without it. | as B-06. | **open, with the user**: the token or the public repositories are needed when KIT-134's P2 starts. |
