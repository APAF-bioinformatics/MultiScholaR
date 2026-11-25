# KIT-251 wording for the owner's pass (VAL-41)

Drafts only; the owner sets the final wording.

| Where | Draft (as on the branch) | Note |
|---|---|---|
| `NEWS.md`, unreleased section, `### Fixes` | Proteomics, Protein Replicate Filter: one click is handled once. Before, the installed application showed 'A protein replicate filter action is already running' at every Apply although the step ran once, and a refusal was shown twice (KIT-251). | AUDIT-KIT-251-impl-20261125 I05: "the installed application" in place of "the application" (the warning came only from the managed worker; in an R session it never appeared). "A refusal was shown twice" is true on both launchers. |

No new user-facing text in the application (one false warning stops appearing).
