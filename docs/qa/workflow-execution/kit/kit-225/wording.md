# KIT-225: user-facing texts (draft; the owner sets the final wording)

Drafts under VAL-41, from the ticket's section 7. Each is shown in the proteomics Normalisation tab after
"Reset to Pre-Normalization". Code: `R/mod_prot_norm_session_helpers.R`.

| Where | When | Draft text |
|---|---|---|
| Notice (warning, 5 s), `buildProtNormResetNotificationMessage` | Reset returned to the table the normalisation started from | "Normalisation reset: the current table is again the one the normalisation started from (<step>)." where <step> is the Quality Control tab's title of the state (`DIA_QC_STEP_TITLES`, e.g. "Duplicate Removal", "Protein Replicate Filter", "Peptide Count Filter"); a state with no title is named "the peptide quality-control steps", the fallback `diaPreNormalisationStepName` uses today (KIT-223 T4 (a), OD-1b may change both) |
| Notice (warning, 5 s), same function | No normalisation result was current (no normalisation-tab state in the active history), and with no state manager | "Normalisation reset: no normalisation result was current, so the Quality Control steps were left as they are." |
| Error notice (10 s), `revertProtNormStateManagerToPreNormalization`, shown by `handleProtNormResetError` after its prefix "Error resetting normalization: " | The table the normalisation started from is in no state of the active history (class `prot_norm_reset_source_unknown`; nothing changed) | "the table this normalisation started from can no longer be found in this session (for example, the design was saved again after the normalisation, or the project was saved without that table), so nothing was reset. In the Quality Control tab press Revert until the step you want is current, then press "Run Normalization & RUV" again." |
| NEWS.md, Fixes | release notes | see `NEWS.md` (KIT-225 entry) |

Replaced: "Normalization has been reset to pre-normalization state (<internal state name>)" and
"Normalization has been reset (no valid previous state found)".
