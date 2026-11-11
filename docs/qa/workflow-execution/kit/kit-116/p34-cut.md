# KIT-116, P-34: the cut point, fixed before any new member (2026-11-11)

**Measure:** the controller's anon (`controller_anon`, the last column of `memory-trace.tsv`) over the journey's `rollup` step (`journey-phases.tsv` start to complete), its minimum. The step comes before DA, so R5's volcano page at DA does not reach it.

**Calibration:** every earlier `m` member with both files (read-only).

| Member | Build | `rollup` min (GiB) | max | end |
|---|---|---|---|---|
| 101m | kit_101 | 0.745 | 0.874 | 0.745 |
| 102m | kit_102 | 0.743 | 0.863 | 0.743 |
| 103m | kit_103 | 0.709 | 0.843 | 0.709 |
| 104m | kit_104 | 0.733 | 0.824 | 0.733 |
| 105m | kit_105 | 0.735 | 0.825 | 0.735 |
| 106m | kit_106 | 0.732 | 0.861 | 0.732 |
| 107m | kit_107 | 0.733 | 0.823 | 0.733 |
| 108m | kit_108 | 0.729 | 0.862 | 0.729 |
| 110m | kit_110 | 0.734 | 0.825 | 0.734 |
| **111m** | kit_111 | **0.873** | 0.964 | 0.874 |
| 111n | kit_111 | 0.731 | 0.821 | 0.731 |
| **112m** | kit_112 | **0.872** | 0.962 | 0.872 |
| 113m | kit_113 | 0.732 | 0.861 | 0.732 |

The two groups do not overlap: released 0.709 to 0.745 GiB (11 members), not released 0.872 to 0.873 (2 members), a gap of 0.13 GiB.

**Cut point: 0.80 GiB.** A new member whose `rollup` minimum is at or above 0.80 is counted as "not released".

**Base rate so far:** 2 of 13 overall, 2 of 4 on builds with KIT-111's code (111m, 111n, 112m, 113m).

**Decision (declared in KIT-116):**
- 0 of 6 new members at or above the cut: P-34 is closed as intermittent and rare.
- 1 or more: diagnose from that member's trace, and run nothing new without the user's word.
