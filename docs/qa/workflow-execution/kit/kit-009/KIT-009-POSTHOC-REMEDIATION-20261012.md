# KIT-009 post-hoc remediation (2026-10-12)

Register: [KIT-009-POSTHOC-AUDIT-20261012.md](KIT-009-POSTHOC-AUDIT-20261012.md). User: "fix findings".

| ID | Repair | Closure |
|---|---|---|
| PH-01 | `tickets/kit/KIT-009.md` carries a claim-limit section: the partitioner's admission is refuse-only at this ticket; the batches are the cursor's delivery bound, the same as the loop's; the pair proved mechanism and exactness, not a reduction. | Ticket text (rev 1.1, section added). |
| PH-02, PH-07 | `tickets/kit/KIT-014.md` carries the hand-off: the family's peak is the ordered output write and the finishing pass (about 0.9 GiB above the kernel plateau in process) plus the slice floor; a reduction needs bounded output passes, not kernel partitioning; the `grouped` arm's fate is decided there. | KIT-014 ticket section. |
| PH-03 to PH-06, PH-08 | Accepted without change (recorded). | — |

No code changed; the completion and the control (`kit_009`) stand.

## Addendum 2026-10-13: the attribution probe's second mode (PH-03, and KIT-010 pre-execution PE-07)

[memory-attribution-probe.json](memory-attribution-probe.json) (both runs in one process; the second starts from the first's retained 1.78 GiB):

| Run | Engine limit | Start | Kernel-phase peak (500 to 540 s) | Output-phase peak (about 210 s) |
|---|---|---|---|---|
| normal cap | 1 GB | 0.58 GiB | 1.87 GiB | 2.48 GiB |
| reduced cap | 268 MB | 1.78 GiB | 2.55 GiB | 2.49 GiB |

The kernel phase rises 1.29 GiB above its start at the 1 GB engine limit and 0.77 GiB at 268 MB: the engine's budget is used during the batches, about 0.5 GiB of it. The output phase peaks at the same reading in both runs, so the step's peak is the ordered output write whatever the engine limit, as the KIT-009 audit found. PH-03 is refined, not changed.
