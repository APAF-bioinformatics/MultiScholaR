# KIT-237: Windows journey through the installed application (records)

Run branch `run-k237-journey` (from `v0.9-staging` dd769671; never merged to staging). Ticket text: KIT-237 r1.2.

## State

**WL: instrument proven on Windows** (instrument rerun 2, run 45637098269, 2026-11-25): WL-0, WL-0u, WL-1 to WL-4
60 pass, 0 fail, 1 recorded; WL-5 passed on all four routes. Earlier: run 45605975875 and instrument rerun 1
45619531044, not proven, every fail in the instrument's mechanics (below). WJ' waits for the release commit.

Earlier state: built, not run. The stage `journey` of `.github/workflows/distributable.yml` and the files of this folder exist on the
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

## WL, run 45605975875 (2026-11-25 11:36 to 12:57 UTC; 81 job-minutes, 5.40 core-hours)

Records: `wl-45605975875/` (the tables; screenshots and per-step JSON stay in the artifact `journey-windows`, 30 days).
Binding as declared: commit under test 7327607f (early check and full check passed; the planted wrong commit stopped),
driver 76e4b519, routes 47ef001f. Machine: `windows-2025`, R 4.6.1 ucrt, Chrome 157.0.8103.79, 4 cores (AMD64 Family
26), LongPathsEnabled 0 read by the installed R, project root `C:/Users/runneradmin/Documents` (30 characters).
Withheld by the check: `WL1-study_parameters.txt` and `workflow-diff.txt` (forbidden string: package names).

**Verdict (`verdict.txt`): instrument not proven on Windows (49 pass, 11 fail, 1 recorded); WL-5 passed.**

Observed:
- WL-0 passes in full (launch 3.1 s, launcher `process`, attach, 40 programs firewalled, LongPathsEnabled 0, the
  application's library the bundle's only, no test-features line).
- WL-1 on `prot_lfq`: S1 to S10.5 reached; workers seen at every QC, normalisation and DA step; N1 to N4 equal; every
  named file found. S9 refused offline: "Enrichment needs protein annotation and this project has none" (prediction 2
  said a test-mode substitute: wrong; KIT-233's L met the same refusal, its expected outcome). S9.d: 0 bytes, "Unable
  request data from server" (K1, KIT-241). Three placeholder export names missing (K2).
- WL-1 on `trx_screen`: S1 to S10.5 reached; K222 equal; N3 equal; one kind-3 finding at S8 from a wrong named-file
  pattern (the DA writes `results/transcriptomics/da_genes/de_results_<contrast>.csv`).
- WL-1 on `lipid_canonical` and `metab_lc`: S5.4 refused "No internal standards found" (K12, as predicted); S6.1
  refused "A step is still running" (the tab's entry QC), then the lane's values were applied (S6.1r) and completed;
  S7 completed with both signals already true at the press (0 s); S8.1 timed out at 600 s although "Session loaded
  successfully!" appeared at 0.5 s (the digest signal named a state the load never commits); the route stopped there.
- WL-2 and WL-3 pass (port refused; S2 crashed from the browser session). WL-4 (a), (b) pass on LFQ.
- WL-5 passes on all four routes (LFQ 9 agree, 10 DOM only; lipid and metab 9/1; trx 7/6).
- e2e lanes: lipid 88/0, metab 190/0, trx 16/0, TMT-LFQ 162/0.
- Paths met (RA-01): the `timeout` outcome (S8.1, with kind-1 and kind-5 findings), "not reachable as a user" (the
  chooser). Not met: the internal deadline (41 minutes of 150 used), a progress modal over 60 s.

The 11 fails and their cause (each in the instrument's mechanics; no product result is judged by them):
1. WL-4 (c): the plant deleted one copy of `post_norm_pca.png` (`results/proteomics/protein_qc/`) while the pattern
   still matched other copies, so no file was missing: the plant, not the check, was wrong.
2-5. WL-0u on all four routes: the chooser opened and listed drive C: (screenshot), but journey.R's
   `k233ChooseFile` looks for folders among `.sF-file` rows; shinyFiles 0.9.3 (the bundle's, lock file) marks folders
   `.sF-directory` (`shinyFiles.js:1138`) and opens them on their own double-click handler (`:1155`). `k233ChooseFile`
   had never been run before this (KIT-233's L-u was not sent).
6-11. `lipid_canonical` and `metab_lc`: S8.1's digest signal (state `loaded_for_de`) is never committed by the load
   (its page signal came at 0.5 s), so S8.1 timed out and the route stopped: "reached the last action", N3, N2 and
   WL-4 (a) on `metab_lc` could not be judged.

Instrument rerun 1 (VAL-54; no criterion, expected value or stop rule changed; ceiling unchanged): `journey-win.R`
deletes every copy the plant names; `k237ChooseFile` replaces `k233ChooseFile` (folders and files, folders by their own
double-click, the listed names kept on a miss); S8.1 and S7 on the omics routes read the page's own message, without a
digest signal; S7 and trx S8 name the files these routes write; a "still running" refusal is pressed again with the
same settings after 30 s (the lane's values only for a default refused on the fixture). With the result audit of
KIT-233's L (RA-02, RA-03, RA-05): a terminal screenshot at the moment of a refusal; the DA table by an exact two-sided
match of the contrast's group names, the candidates and the path read recorded.

Predictions for rerun 1 (written before it): WL-0 as before; WL-0u passes on all four routes; WL-1 reaches S10.5 on
all four routes, with S6.1 on lipid and metab completing at the rendered defaults after the second press (or refused by
a default: then the lane's values, recorded); WL-4 (a) on LFQ and `metab_lc`, (b), (c) pass; N2 on `metab_lc` equal;
WL-2, WL-3 pass; WL-5 passes on all four; known findings only (K1, K2 to K4, K12); the internal deadline and a progress
modal over 60 s still not met (RA-01: not planted; owed before WJ', ticket section 20).

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

## WL instrument rerun 1, run 45619531044 (2026-11-25 13:02 to 14:39 UTC; 97 job-minutes, 6.47 core-hours)

Records: `wl-45619531044/`. Binding as run 1, routes at becc4930. **Verdict: instrument not proven on Windows (52
pass, 8 fail, 1 recorded); WL-5 not passed (lipid and metab, S6.1w).**

Observed: WL-0, WL-2, WL-3, WL-4 (a), (b), (c) on LFQ pass (the plant now deletes every copy; the missing-file finding
names `post_norm_pca`). LFQ, lipid and trx reach S10.5; N1 to N4 equal on LFQ (the N2 table, by the exact two-sided
match, `results/proteomics/da_proteins/da_proteins_groupKO-groupWT_long_annot.tsv`, 1 candidate: RA-03, RA-05); K222
equal on trx (6 keys, 6 biological samples; KIT-222 merged). 32 terminal screenshots taken (RA-02). e2e lanes 88, 190,
16, 162 passed, none failed. Known findings only, plus one observation for the coordinator: on lipid and metab the first
press of "Run Normalization Pipeline" is refused "A step is still running" (the tab's own entry QC) while the button is
shown enabled; the second press, 30 s later with the same settings, completes (runs 1 and rerun 1 alike).

The 8 fails and their cause (instrument mechanics):
1-4. WL-0u: the chooser now reaches the input folder, but "Select" stays disabled: shinyFiles selects a row only on a
   click with `event.button === 0` (`shinyFiles.js:77-78`), and jQuery's `trigger('click')` carries no button.
5, 7. WL-5 on lipid and metab at S6.1w: the digest state is committed at 27.7 s, before the pipeline ends ("Normalization
   pipeline complete!" appeared during S6.2), so S6.2 was pressed early. journey.R made its own S6.1 page-only for the
   same reason (its fix 2).
6, 8 (and N2, WL-4 (a) on metab): `metab_lc` S8.3 timed out: its 6 metabolites have 0 significant, so the results
   table's default view (significant only) has no rows, and the step's signal waited for rows.

Instrument rerun 2 (VAL-54; no criterion, expected value or stop rule changed): the chooser dispatches native mouse
events with button 0; lipid and metab S6.1 use the page's completion message only; S8.3 waits for the summary
statistics and checks the table after (non-empty required when the summary counts significant features); WL-4 (a) on
metab sets "Show:" to "All Results" (a user's choice) before emptying the table.

Predictions for rerun 2 (written before it): as rerun 1, plus WL-0u passing on all four routes, WL-5 passing on all
four, `metab_lc` reaching S10.5 with N2 equal (0 and 0) and WL-4 (a) on metab passing.

## WL instrument rerun 2, run 45637098269 (2026-11-25 14:43 to 15:43 UTC; 61 job-minutes, 4.07 core-hours)

Records: `wl-45637098269/`. Binding: commit under test 7327607f (early and full checks passed; the planted wrong commit
stopped), driver 76e4b519, routes 8f3b08ee. Machine as run 1. Journey step 23.7 minutes (internal deadline 150 not
reached); `LongPathsEnabled` restored to 1 after the journey (job log).

**Verdict (`verdict.txt`): instrument proven on Windows (WL-0, WL-0u, WL-1 to WL-4: 60 pass, 0 fail, 1 recorded); WL-5
passed.**

Observed (each against its prediction):
- WL-0 passes in full (prediction 1: met).
- WL-0u passes on all four routes: the chooser clicked through the folders to every input in user mode (prediction:
  met after the two instrument fixes).
- WL-1: all four routes reach S10.5 with workers seen (LFQ 245 s of steps, lipid 280, metab 275, trx 139). N1 to N4
  equal on LFQ (394 proteins, 10 samples; 18 significant; the report's sentence 18); N3 equal on all four; N2 on
  `metab_lc` equal (0 and 0, from `de_negmode_…` + `de_posmode_metabolites_groupKO-groupWT_long_annot.tsv`); K222 equal
  on trx (6 keys, 6 biological samples).
- WL-4 (a) on LFQ and on `metab_lc`, (b), (c) pass. WL-2, WL-3 pass.
- WL-5 passes on all four routes.
- e2e lanes: lipidomics 88/0, metabolomics 190/0, transcriptomics 16/0, TMT-LFQ 162/0 (laptop's 153; the lane gained
  tests since).
- Findings, all known except the last: S9 refused offline, "Enrichment needs protein annotation…" (as KIT-233's L; my
  prediction 2 of a test-mode substitute was wrong); S9.d download 0 bytes, not a zip (K1, KIT-241); placeholder export
  names (K2 to K4); "No internal standards found" on lipid and metab (K12, KIT-242); **new, for the coordinator:** on
  lipid and metab the first press of "Run Normalization Pipeline" after opening the tab is refused "A step is still
  running…" while the button is shown enabled (the tab's own entry QC); the same press 30 s later completes (all three
  runs). WL3w (observation): a worker ended during the longest LFQ QC step; the step refused "managed action ended:
  worker_failed".

Inferred: the instrument (KIT-233's driver 76e4b519 under `journey-win.R` at 8f3b08ee) observes, on the installed
Windows application at 7327607f, completion, refusal, crash, timeout (rerun 1's S8.1/S8.3), "not reachable" (runs 1 and
rerun 1's chooser), planted silent outputs, and the digest/page agreement on each of the four routes, through their
completion.

## Core-hours

WL 81 + rerun 1 97 + rerun 2 61 = 239 job-minutes = 15.93 core-hours (expected for WL 102 minutes; ceiling WL 220 plus
instrument reruns within 1,410; two of three instrument reruns used).

## Not checked

- The internal deadline and a progress modal over 60 s were not met by any WL run (RA-01; owed before WJ', ticket 20).
- That the firewall rules block outbound traffic from the application's R (shown only by their existence).
- User mode beyond S2 (WJ''s; VAL-216), and WJ' itself (release commit).
- Screenshots and per-step JSON are in the run artifacts only (30 days), not committed.

