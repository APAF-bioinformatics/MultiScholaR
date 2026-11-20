# KIT-151 c, task 2 (SCI-05): pinning the maths library

Jobs x conditions: 36; CPU models seen: AMD EPYC 7763 64-Core Processor (6 jobs); AMD EPYC 9V45 96-Core Processor (2 jobs, AVX-512); AMD EPYC 9V74 80-Core Processor (4 jobs, AVX-512).
Instrumentation: every job and condition passed its checks.

| condition | arms | arms with one class | most classes in an arm | same-CPU pairs not identical | cross-CPU pairs not identical | cross-CPU iterations differ | largest cross-CPU abs dZ (raw) | lowest cross-CPU matched abs r |
|---|---|---|---|---|---|---|---|---|
| U | 8 | 0 | 3 | 24/176 | 304/352 | 0 | 5.83e-13 |       1 |
| P | 8 | 8 | 1 | 0/176 | 0/352 | 0 |    0 |       1 |
| PN | 8 | 8 | 1 | 0/176 | 0/352 | 0 |    0 |       1 |

- U: agree to rounding, not bit-identical
- P: identical across the CPU models seen
- PN: identical across the CPU models seen

History check (instrumentation, not a gate): unpinned runs on the CPU model of a committed run; see pin-history.tsv. Rows with a same-CPU job: 7; iterations equal to history: 7.
