# DRAFT: protocol amendment v1.11, for the user

Status: DRAFT, not in force. On the laptop day, after the user's go, the text between the rules below is appended to
`docs/qa/workflow-execution/kit/protocol.md` as section 13 and the version line is updated. The header still reads 1.9
although section 12 made it 1.10; the same edit corrects that. Then the bindings are regenerated once (step 4 of
`re-record-plan.md`).

Items 1 to 5 are KIT-129 section 4 as the user approved them (U5 and U11, 2026-11-18). Item 6 is the user's decision on
KIT-052 of 2026-11-19 (option A: "shorten the worker's request-folder names in the approved harness revision (protocol
v1.11)"). It is new text and needs the user's read before it is in force.

---

## 13. Amendment v1.11 (KIT-129: the RUV design presets and the single re-record; KIT-052: the request-folder names)

**Authority.**
- The user, 2026-11-18, U5 and U11 on KIT-129: "yes", after reading items 1 to 5. The CORE pin on
  `run_workflow_checkpoint_proof.R` is lifted for the one edit of item 3.
- The user, 2026-11-19, on KIT-052: option A, item 6.

No rule, threshold, guard, gap tolerance, ceiling, workload or member kind changes.

1. **RUV design presets.**
   - The six step maps' RUV and DA rows, `kitJourneySteps` and the driver follow the presets. The bindings carry
     `ruv_design` per journey: a preset id, or `default` for the tab's pre-selection.
   - A binding without `ruv_design` reads as before:
     - `automatic`, `manual` and `covariates` are Legacy;
     - `correction` is the legacy correction;
     - `skip` is None.
   - The driver reproduces such a binding by selecting the Legacy preset and the bound mode. With `skip` that runs
     today's route without RUV.
   - A build without the panel runs only the Legacy presets. A binding for any other preset is refused there.
   - Decisions recorded before v1.11 stand under their own protocol hash.
2. **Intermediate builds.** Intermediate builds are judged store against store with the build before. They are
   recorded, not counted, and none is an oracle. In-process oracles (7.1) are recorded once, on the final build.
3. **`run_workflow_checkpoint_proof.R`.**
   - The pin is lifted for one edit of `checkCheckpointSummaryFiles`. That edit:
     - checks the integration files of KIT-061's contract under `integration/proteomics/`, validated by the package's
       reader;
     - keeps the RUV table names per KIT-126's W0 table.
   - The new hash is bound with the kit harness files (`core_pin_lifted`).
   - The other CORE-pinned files stay as they are.
4. **Transcriptomics digests off the laptop.**
   - Off the laptop, the nine transcriptomics digest expectations are judged by `.github/scripts/save-trx-products.R`
     and `packaging/verify/compare-results.R` under KIT-134's rule (2.0.4: relative 1e-6, absolute 1e-9, discrete
     exact). There, `MULTISCHOLAR_TRX_DIGEST_RULE=values` is set and a differing digest is a note.
   - On the laptop the exact digest is the gate.
   - The corpus manifest digest covers counts and design, not `truth.json`. The four files the products read are
     compared exactly with `kit-040/bindings.json` on every machine.
5. **Transcriptomics' independent check (v1.8).**
   - The notebook describes the legacy design on library counts.
   - For the default route, the independent check is a direct call of `RUVSeq::RUVs` and edgeR on the journey's saved
     inputs.
   - The notebook reference stays the check of the Legacy preset, on a corpus without declared technical replicates
     (question U6).
6. **Request-folder names (KIT-052 option A).**
   - The package writes shorter folder names. Old name to new name:

     | Folder | Before | From v1.11 |
     |---|---|---|
     | Session folder | `.dia-native-<14>` | `.n-<8>` |
     | Store | `.workflow-entities` | `.we` |
     | Request and staging root | `.workflow-export` | `.wx` |
     | Worker request | `wfx-<32>` | `w-<16>` |

   - The native archive's store members carry the store folder's name.
   - A folder that exists keeps its name: projects and archives written before open as before, and a user's folders
     are never renamed.
   - The harness and the GUI-entry recorder read both layouts:
     - `workflow_kit.R`: `kitStoreDirectory`, `kitRequestDirectory`, `KIT_*_RE`;
     - the coordinator, the checkpoint driver, `workflow_core.R` and `workflow_perf106.R`;
     - `helper-gui-entry.R` and `helper-workflow-bulk-rollup.R`.
   - Stores from builds before and after the change compare store against store. Fingerprints are of content; the
     paths are read through these helpers.
   - The GUI-entry oracles' `<store entries>` keys and the frozen keys of the omic module tests name the folders. They
     change with this item and are re-recorded in the single re-record, attributed to it.
7. **Effect.**
   - The changed harness files' hashes, the six step maps' hashes and the protocol hash change. The bindings are
     regenerated once.
   - The stored references move as listed in `harness-v1.11/inventory.md` section 5. Each difference is attributed to
     one mechanism (`re-record-plan.md`), as in KIT-118 and KIT-123.
   - Every rule, threshold, guard, gap tolerance, ceiling, workload, limit, member kind and selector is unchanged.
   - Decisions recorded before v1.11 stand under their own protocol hash and remain readable:
     - their receipts name the earlier harness hashes;
     - their stores use the earlier folder names, which the v1.11 readers still accept;
     - a binding without `ruv_design` keeps its meaning.

---

## What the amendment does not change

- Member kinds, selectors, the 16 GiB cap, the 15.5 GiB stop, the 12 GiB start floor, the 2 GiB host floor and the
  gap rules (250 ms; 500 ms for journeys).
- The 0.95 and 1.02 decision rules.
- The run ledger's authority and the correction-rerun rule.
- The KIT-005 oracles, the CORE populations and the other two CORE-pinned files (`run_workflow_checkpoint_proof.sh`,
  `workflow_perf106_repair.R`).

## Bindings to regenerate (once, after the protocol text is in force)

- `harness.protocol_sha256`.
- The six `step_map*_sha256`.
- `harness.files`: the kit harness files, now with `run_workflow_checkpoint_proof.R`, and `core_pin_lifted`.
- The B0-B5 builds (install, build identity, exact-sum identity).
- The ticket keys for the B-stage members and the lifecycle pair, with `ruv_design` per journey.
- The in-process oracles recorded on B5.

## How v1.10 evidence stays readable

- Receipts carry their own bindings and protocol hash. The v1.11 readers accept the earlier folder names, so earlier
  member directories and stores are still fingerprinted and compared.
- `kitRuvDesign` gives every earlier binding its earlier meaning.
- The earlier `protocol.md` and `bindings.json` stay in git history.
