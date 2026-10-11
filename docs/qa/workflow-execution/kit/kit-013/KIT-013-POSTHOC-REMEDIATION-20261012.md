# KIT-013 post-hoc remediation (2026-10-12)

Closes findings of [KIT-013-POSTHOC-AUDIT-20261012.md](KIT-013-POSTHOC-AUDIT-20261012.md) (user: "fix findings"). The audit register is unchanged; this record tracks closure. No science, threshold or product surface changed.

| ID | Action | Closure evidence | State |
|---|---|---|---|
| PX-01 | Reproduced in `test-workflow-native-session.R` ("audit PX-01: a QC-time configuration update does not fire the design sync observer"): bus registered on the session, owner adopted, the design initial-state reactive and sync observer registered with counters, then `updateConfigParameter` run through the bus environment. Result: neither the initial-state reactive nor the sync observer runs again after the configuration write (counts unchanged across the flush); the sync observer is bound to the data reference and the Reset input (`bindEvent`, `R/mod_prot_design_builder_server_helpers.R:130`), and only Reset rebuilds the state. **No code change**: the plausible defect does not occur. | `test-workflow-native-session.R` 488 pass; the one failure is the pre-existing process-loss publication artefact. | closed, verified |
| PX-02 | `open-probe.R` rerun at `16b1f056`; `open-probe.json` `head_commit` now names a tree containing the audited code; result unchanged (replaced folder removed after each Open; 266 → 123 entities). | `open-probe.json` | closed |
| PX-03 | `open-probe.R` extended: after the second Open, a native session export from the opened owner (`runDIANativeSessionExport` with the session's metadata, extra references and resolved dependencies) and a reopen of the produced archive. Result: archive 184.1 MB written in 5.8 s, reopened at `correlation_filtered` with the same 16 state names. Probe-side note: `resolveDIAStateValue` returns the metadata list itself. | `open-probe.json` `save_after_open` | closed |
| PX-04 | No change (application never reaches the state; scripts unchanged). | — | accepted |
| PX-05 | No change (decision recorded in the item's note). | — | accepted |
| PX-06 | `MEMORY.md` index hook updated to `kit_011`. | memory index | closed |
| PX-07 | No change. | — | accepted |
| PX-08 | No change (claim limits already record them). | — | accepted |

Failed attempts recorded: the first two probe reruns aborted in the new export step ("native session metadata contains legacy history fields") because the probe read the metadata reference as `…$data`, which is NULL; fixed in the probe only.
