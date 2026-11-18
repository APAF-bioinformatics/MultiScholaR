# KIT-124 scored tables

Outputs of `../aggregate.R`, one folder per scoring run; each folder's `k124-aggregate.meta.txt` names the script hash, inputs (sha256) and command.

| Folder | Content | State |
|---|---|---|
| s1-core, s2-core, s3-core, s4-core | Main cells C0 to C7, 1,000 runs each | Final for the main cells |
| s34-stress-partial | Extra-hard cells of S3 and S4 finished by 2026-11-18 04:30 | **Superseded 2026-11-18 by `final-20261118`** (kept, not deleted; its 483 extra-hard rows are byte-identical in the final tables) |
| final-20261118 | The final scoring, one folder per simulator (`s1` to `s4`): main cells C0 to C7 and every extra-hard cell (S1, S2, S4: S1 to S14; S3: S1 to S5, S7 to S14), 500 runs per extra-hard cell. Laptop rows for the cells finished there; hosted-runner run 44201244793 (label `hard`) for S1 S9 to S14, S2 S1 to S14, S3 S14. With `inventory.tsv` (where each cell ran, rows and statuses), `jobs.tsv` (jobs, attempts, artifact digests), `option-arms-extra-hard.tsv` (the option-set designs side by side, written by `option-arms.R` from the folder's tables), `not-committed.sha256`, and `run.log` per simulator | Final, scored 2026-11-18. Main-cell rows byte-identical to s1-core to s4-core, which stay the reference for the main cells. No cell mixes platforms; no cell exists on both, so the runner's cells rest on `../gha-proof-20261118/`. Runner figures of arm 3 are platform-dependent (stated there) |
| newarms-gha | Hosted-runner run 44200924581: arms 1, 2, 22, 23, 24 (S4: 1, 2, 22, 23), C0 to C7, 1,000 runs, package tree `0b6fcb6d`. With `compare/` (arms 1 and 2 against the laptop's rows under the fixed rule), `jobs.tsv` (jobs, attempts, artifact digests) and `new-arms-beside-reference.tsv` (arms 22 to 24 beside arms 8 and 5 of s*-core) | Scored 2026-11-18, all eight cells of every simulator. In S3 C0, C1, C5, C6, C7 the gate as first declared was missed on bias and one rmse (8 of 16,000 rows; every value the rules use identical); the user accepted the exception (decision M19) and the five cells were scored afterwards; `compare/gate-summary.tsv` holds both readings. C6 of S3 and S4 is void on this tree |
| r0a-gha | Hosted-runner run 44203261310: the technical-duplicate cell C6 on kit-128's tree (`ace1f4de`), S3 arms 1, 2, 7 and S4 arms 1, 2, 7, 15, 16, 1,000 runs. With `jobs.tsv` | Scored 2026-11-18. No laptop rows exist on that tree: not compared with the laptop |

Copied from the session scratch folder on 2026-11-18 so they are durable. KIT-144's figures read from here.

The hosted-runner folders were scored from the runs' artifacts; the per-run tables are not committed (sha256 in `not-committed.sha256` and in the sidecars, artifact digests in `jobs.tsv`). Reading: `../outcome.md`, the sections of 2026-11-18 (the last one for `final-20261118`).

`final-20261118`: S1's input was the laptop's `s1-runs.tsv` and log followed by the runner's rows for S9 to S14 under one header (S1 has no per-cell files); S2 and S3 were read from the laptop's and the runner's per-cell files side by side in one `cells/` folder; S4 from the laptop's folder alone. The sidecars list every input with its sha256 (scratch paths for the assembled copies; `not-committed.sha256` repeats the runner's and the assembled files by relative path).
