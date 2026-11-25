# KIT-237: Windows journey through the installed application (records)

Run branch `run-k237-journey` (from `v0.9-staging` dd769671; never merged to staging). Ticket text: KIT-237 r1.2.

## State

Built, not run. The stage `journey` of `.github/workflows/distributable.yml` and the files of this folder exist on the
run branch. No run has been sent: WL waits for KIT-233's lifecycle run L to be recorded "instrument proven" (at
2026-11-25 L was in its VAL-54 instrument reruns: `kit-233-journey` 4793c9e8 and 0260ac73 fix runs 45523435786 and
45525410012; no `runs.tsv` or `outcome.md` on that branch yet). WJ' is not this item's.

## What is here

| File | What |
|---|---|
| `journey-win.R` | the Windows wrapper: sources KIT-233's `journey.R` (at the request's `journey_commit`), rebinds four of its names (listed in `k237Rebind`), launches the installed `entry.R` with the installed `R\bin\Rscript.exe`, stops by process tree, and runs WL (WL-0 to WL-5, one journey per further route, the chooser once) |
| `journey-route-lfq.R` | the label-free (MaxQuant) route on `prot_lfq_maxquant` (5.2): steps, QC table, named files per step, the S9 download |
| `journey-route-omics.R` | the three further Windows routes (OD-1 (b), VAL-225): `lipid_canonical`, `metab_lc`, `trx_screen` (5.2a) |
| `known.tsv` | known Windows items (2.2), expected findings (KIT-241's download; manifest export names no code writes) |
| `runs.tsv` | one row per run (none yet) |

## The request (one new file per run under `requests/`, never changed after it is added)

```json
{
  "run": "WL",
  "commit_under_test": "<40-hex: the v0.9-staging commit just merged into run-k237-journey (SC-03 (a))>",
  "journey_commit": "<40-hex: kit-233-journey commit whose L is recorded instrument proven (at or after 628545be)>",
  "route_commit": "<40-hex: the run-k237-journey commit holding these files>",
  "route_files": ["journey-route-lfq.R", "journey-route-omics.R"],
  "routes": ["prot_lfq", "lipid_canonical", "metab_lc", "trx_screen"],
  "mode": "test",
  "platforms": ["windows"],
  "e2e": true,
  "gui_entry": false,
  "deadline_minutes": 150,
  "inputs": {
    "prot_lfq": {"dir": "prot_lfq_maxquant", "files": {
      "proteinGroups.txt": "4d94c48aecdd6d3d533fe8e21464d37befa3ea2638fbf622acb0a02cf0d2b04e",
      "proteins.fasta": "beb40a412d8850ad734fefee5747a9f54ce5b94097cfc4414e25a358634a0b32",
      "design.tsv": "fed1f2da738f786f9a0ad38d0bd81b2c3ea306dd19e645e03d092122806411ec"}},
    "lipid_canonical": {"dir": "lipid_canonical", "files": {
      "lipidsearch_lcms_pos.txt": "c23154c95596cc055c5acfe9c14de4f9e2ca0d3b4821455f398595b957f8b0e4",
      "lipidsearch_lcms_neg.txt": "14b100d399ffe23d381de9e5da4953a8e13cc07bb2916ad00abf3722cbd4d15f",
      "design.tsv": "46a038d8823cd81c81f33cc61a086308001b985ce49753d495a9077a341f4310"}},
    "metab_lc": {"dir": "metab_lc", "files": {
      "lcms_pos_features.tsv": "77ea1cd00bf4292bf3112c64b8466bbd371a2971fe66ef351e5318c447ec21e5",
      "lcms_neg_features.tsv": "f32cc0a3e60eed918f5fd42127eeecc43964f3b53dd248e2e40efa57f91cbf04",
      "design.tsv": "46a038d8823cd81c81f33cc61a086308001b985ce49753d495a9077a341f4310"}},
    "trx_screen": {"dir": "trx_screen", "files": {
      "counts.tsv": "47f5b6ce4b388f38431ce82825b6712f80155e77eea45245ae8830a065a06a8d",
      "design.tsv": "9437ea87a5a0175a52500e684c250f6f14c92c1f83ee4deb46c3d1652dd8dadc"}}
  }
}
```

Before the request is written (SC-03 (a)): merge the current `v0.9-staging` into `run-k237-journey` (a merge commit
on the run branch only) and set `commit_under_test` to that staging commit; the build's package commit is then that
commit, and the early check after "The package commit" step (SC-03 (b)) stops a mismatch in seconds. While WL runs,
nothing else is dispatched on the branch (SC-10: any other stage's dispatch would still cancel it).

SHA-256 from `tests/testdata/e2e/manifest.json` (prot_lfq :172-174, lipid_canonical :365-368, metab_lc :251-254,
trx_screen :407-409). Dispatch: `gh workflow run distributable.yml --ref run-k237-journey -f stage=journey
-f platform=windows -f journey_request=<name>`, one at a time (stop rule 7).

WL as declared in the ticket (r1.2): purpose "instrument proven on Windows"; binding the commit under test, the
KIT-233 commit and the route commit named in the request; journey step `timeout-minutes` 160, internal deadline 150
minutes; WL ceiling 220 job-minutes (14.67 core-hours), expected 102; at most three VAL-54 instrument reruns.

## r1.3 (script check before WL, SC-01 to SC-08)

- SC-01: lipid and metab QC press intensity filter, duplicate detection and resolution, the internal-standard analysis,
  then Finalize, at rendered defaults; `sanitize_names` left at its rendered TRUE. The ITSD analysis is expected to
  refuse on these fixtures (no standards; `known.tsv` K12, KIT-242); the refusal is recorded and QC goes on.
- SC-02: user-mode chooser import on every route; WL-0u (blocking) per route, the chooser selecting every input.
- SC-03: early package-commit check in the stage; staging merged into the branch before the request.
- SC-04: the omic's own stepper as the DOM signal (state script rebound per session); exact texts elsewhere.
- SC-05: named screen outputs (import status; DA status, summary statistics, results table); N2 on metab_lc; WL-4 (a)
  planted on metab_lc. The trx route names only its DA status (its results tab was not read).
- SC-06: WL-4 (c) requires a non-empty deleted-file name.
- SC-07: WL-2, WL-3 and the four user-mode sessions run before the three added journeys.
- SC-08: `LongPathsEnabled` restored after the driver (finally).

## WL: declaration and predictions (written before the run, 2026-11-25)

- **Purpose:** the lifecycle on the installed Windows application, one journey per omic shipped on Windows (OD-1 (b),
  VAL-225); decides "instrument proven on Windows" (WL-0, WL-0u, WL-1 to WL-4 without a fail) and WL-5.
- **Binding:** commit under test `7327607f` (`v0.9-staging` head, merged into this branch as c4829951; its product
  tree equals e0747fc1, KIT-233's fourth L, which differs from it only in `docs/methods-decision-ledger.md`; the bundle
  is built from this branch's tree, whose package paths equal `7327607f`'s, so the early check and the journey step's
  check bind it); journey driver `76e4b519` (`kit-233-journey`, instrument proven by run 45601782272, VAL-234); route
  files at the commit that adds this declaration; mode test (WL; the user mode of VAL-216 is WJ''s, with WL-0u
  proving the chooser); e2e lanes on, GUI-entry lanes off; `deadline_minutes` 150.
- **Machine:** a GitHub-hosted `windows-2025` runner (Windows Server 2025 Datacenter, 4 cores), Windows' long-path
  support switched off during the journey (`LongPathsEnabled` 0, restored afterwards), the application's default
  project folder.
- **Bounds:** journey step `timeout-minutes` 160, internal deadline 150 minutes, hard wait per step 10 minutes; WL
  ceiling 220 job-minutes (14.67 core-hours), expected 102 (6.80); at most three VAL-54 instrument reruns, each
  recorded here and in `runs.tsv`; a result is not re-run.
- **Predictions:**
  1. WL-0 passes in full (commit checks; Chrome; launch and attach; launcher `process`; firewall; LongPathsEnabled 0;
     the application's `.libPaths()` inside the bundle; no test-features line).
  2. WL-1 on `prot_lfq`: S1 to S10.5 complete with at least one heavy step run by a worker. Known findings, not new:
     S9.d's download not a zip that opens (K1, KIT-241 not merged), three kind-3 placeholder export names (K2).
     N1 to N3 equal.
  3. WL-4 (a), (b), (c) pass on LFQ; WL-4 (a) on `metab_lc` passes.
  4. WL-2 and WL-3 pass.
  5. WL-0u: uncertain (the chooser has never been driven on Windows); predicted to pass on all four routes if the
     chooser's first volume is the drive holding the default project folder.
  6. The three added journeys reach S10.5. Known: on `lipid_canonical` and `metab_lc` "Analyze Internal Standards"
     refuses with "No internal standards found" (K12, KIT-242) and QC goes on; placeholder export names (K3, K4).
     `trx_screen`'s sample-key check is equal (KIT-222 merged; a shared key would be a finding). `metab_lc` N2 equal.
  7. WL-5 passes on LFQ (KIT-233's L-5 passed on the same signals); on the three added routes uncertain (first check of
     their DOM signals).
  8. The whole journey finishes inside the 150-minute deadline.

## Observed

Nothing has run yet.

## Checked locally (no R)

- The plan job's guard, run with its own awk on the edited file: 4 uploads, 1 in snapshot, 2 guarded
  `parity-results`, 1 guarded `journey-results` (windows job), upload jobs "linux: windows: windows: ".
- The plan step's and the guard's bash parse (`bash -n`); the workflow parses as YAML.
- Fixture SHA-256 in the working tree equal the manifest's for `prot_lfq_maxquant`.

## Not checked

- No R, PowerShell or browser code of this folder or the stage has run; every selector, completion signal and named
  file pattern of the three added routes is by reading only.
- Whether `processx` starts the installed `Rscript.exe` with the wrapper file, and whether `AppDriver` attaches, on
  Windows (WL-0).
- Whether `taskkill /T /F` plus `kill_tree()` makes the port refuse within 10 s (WL-3).
- Whether git's zip archive keeps the committed bytes of the fixtures on the runner (`inputs.tsv` decides).
- The DOM completion signals of the lipid, metab and trx routes (WL-5 is their first check).
