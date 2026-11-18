# KIT-124 scored tables

Outputs of `../aggregate.R`, one folder per scoring run; each folder's `k124-aggregate.meta.txt` names the script hash, inputs (sha256) and command.

| Folder | Content | State |
|---|---|---|
| s1-core, s2-core, s3-core, s4-core | Main cells C0 to C7, 1,000 runs each | Final for the main cells |
| s34-stress-partial | Extra-hard cells of S3 and S4 finished by 2026-11-18 04:30 | Partial; replaced by the final scoring when every simulator's extra-hard cells are done |
| newarms-gha | Hosted-runner run 44200924581: arms 1, 2, 22, 23, 24 (S4: 1, 2, 22, 23), C0 to C7, 1,000 runs, package tree `0b6fcb6d`. With `compare/` (arms 1 and 2 against the laptop's rows under the fixed rule), `jobs.tsv` (jobs, attempts, artifact digests) and `new-arms-beside-reference.tsv` (arms 22 to 24 beside arms 8 and 5 of s*-core) | Scored 2026-11-18, all eight cells of every simulator. In S3 C0, C1, C5, C6, C7 the gate as first declared was missed on bias and one rmse (8 of 16,000 rows; every value the rules use identical); the user accepted the exception (decision M19) and the five cells were scored afterwards; `compare/gate-summary.tsv` holds both readings. C6 of S3 and S4 is void on this tree |
| r0a-gha | Hosted-runner run 44203261310: the technical-duplicate cell C6 on kit-128's tree (`ace1f4de`), S3 arms 1, 2, 7 and S4 arms 1, 2, 7, 15, 16, 1,000 runs. With `jobs.tsv` | Scored 2026-11-18. No laptop rows exist on that tree: not compared with the laptop |

Copied from the session scratch folder on 2026-11-18 so they are durable. KIT-144's figures read from here.

The hosted-runner folders were scored from the runs' artifacts; the per-run tables are not committed (sha256 in `not-committed.sha256` and in the sidecars, artifact digests in `jobs.tsv`). Reading: `../outcome.md`, last section.
