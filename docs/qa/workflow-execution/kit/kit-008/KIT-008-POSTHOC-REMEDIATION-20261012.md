# KIT-008 post-hoc remediation (2026-10-12)

Closes findings of [KIT-008-POSTHOC-AUDIT-20261012.md](KIT-008-POSTHOC-AUDIT-20261012.md) (user: "fix please"). The audit register is unchanged. No science, threshold, member or decision changed.

| ID | Action | Closure evidence | State |
|---|---|---|---|
| PH-01 | The coordinator passes `api` to `kitGuiMemberDecision`, so every `kit-gui` launch decision recomputes both members' state fingerprints from their immutable stores under the current canonical form (the branch correction 2 exercised on the real stores of `kit-gui-008`). Bindings regenerated (the coordinator is a hashed harness file). | `workflow_kit_coordinator.R`; `bindings.json` harness hashes; the KIT-008 decision is unchanged (`gui_member_observed`). | closed |
| PH-02 | Outcome test count corrected to 227. | `outcome.md` | closed |
| PH-03 to PH-08 | No change, as the audit proposed. | — | accepted |
