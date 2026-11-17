# Plan 4.7 audit: remediation record (2026-11-17)

The register is [plan-4.7-audit-20261117.md](plan-4.7-audit-20261117.md) and is left as written. Authority for the repairs: user, 2026-11-17, "1 yes"; on integration, "i need to source the workbook and verify it first so let me do that in the next few days - integration comes last anyway. GUI FOR THE INTEGRATION offered the same way as the others but import data would source data processed through MS and likely from the available data in the project folder structure that is currently in place for a project." The repairs are ticket and plan text (plan 4.8); no package code, no run.

| Finding | State | Where | Closure evidence |
|---|---|---|---|
| B-01 | **closed 2026-11-17** (run 44054974678, `kit-134/b01-ci-run.md`) | KIT-134 "The comparison across machines"; plan 4.8 working rule | The laptop's exact comparison is the gate again. Closes when one full run on a Linux runner has been read and the comparison rule is written. |
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

## Added by the independent read (`plan-4.7-independent-review-20261117.md`; it judged the text at 04:45, after plan 4.8)

| ID | Finding, in one line | State |
|---|---|---|
| IR-B1 | No DIA route on Windows, no owner (= B-06). | **open, with the user** |
| IR-B2 | Two dates for the token; a bundle built with it would publish the internal packages as an artifact. | repaired in text: KIT-134 revision 1.2 (no artifact, no cache while internal; one date) |
| IR-B3 | Managed execution needs Linux 6.12 or later, user systemd, cgroup v2 (checked: `R/utils_workflow_execution.R:1019-1020`); the bundle targets older hosts and would refuse at the first heavy step. | repaired in text: KIT-134 revision 1.2 (start-up probe, fallback order); KIT-131 off the cut list |
| IR-B4 | KIT-125 had no definition of a contradiction. | repaired: KIT-125 revision 1.2 |
| IR-B5 | Integration cannot be authored today; KIT-060 held behind KIT-124; KIT-133 behind KIT-063. | part: KIT-060 can start (revision 1.3); KIT-133 decoupled (revision 1.2); KIT-061 to 063 wait for the workbook: **open, with the user** |
| IR-M1 | P1's proof was a rebuild on the laptop. | repaired: KIT-134 revision 1.2 |
| IR-M2 | The tolerance was to be written after the differences were seen. | repaired: KIT-134 revision 1.2 (two parts fixed now) |
| IR-M3 | The workflow does not pin R; KIT-129's exit asked for an Actions suite before a rule exists. | part: KIT-129 revision 1.2; pinning R in `module-ci.yml` is open (with B-01) |
| IR-M4 | W from complete proteins fails at high missingness (KIT-118); the 50% cell is a stress cell. | repaired in text: KIT-126 revision 1.2 (shown at the gate, or a stated limit) |
| IR-M5 | Consumers missing from KIT-126's list. | repaired: KIT-126 revision 1.2 |
| IR-M6 | The likely winning arms take k from `sva::num.sv` (checked: `kit-124/sim/s1_*.R:90`); `sva` is not a dependency. | repaired in text: KIT-129 revision 1.2; the choice is **with the user at the gate** |
| IR-M7 | KIT-126, 127, 129 fixed the design before the gate. | repaired: KIT-126 revision 1.2 (re-read against the gate) |
| IR-M8 | No block in `glmQLFit`; no limit; needless dependency. | repaired: KIT-128 revision 1.1 |
| IR-M9, IR-M10 | MOFA2 input: no claim limit, four meanings of "corrected", hash identity with a Python model. | repaired in text: KIT-060 revision 1.3 |
| IR-M11 | KIT-133 had no consequence for a bad result. | repaired: KIT-133 revision 1.2 |
| IR-M12 | KIT-131's exit named a lane that cannot run on Windows. | repaired: KIT-131 revision 1.3 |
| IR-M13 | The spike did not install MultiScholaR on Windows. | repaired: the running spike was told (KIT-130 revision 1.1, piece 2) |
| IR-M14 | The laptop's gate between KIT-126 and the re-record. | repaired: KIT-126 revision 1.2 |
| IR-M15 | No ticket holds the laptop's calendar; the confirming studies run under load with a wall ceiling. | **open**: the measurement windows are dated at KIT-088's authoring; the paused time is given back to the simulators' ceiling (recorded in `scripts-k124/paused-seconds.txt`) |
| IR-M16 | Estimates known to be wrong; KIT-050's release case needs a corpus nothing provides. | **open**: re-estimated as each stub is authored; KIT-050's case is decided at the final gate |

The review's nine minors are in its own file and are not repaired here.
