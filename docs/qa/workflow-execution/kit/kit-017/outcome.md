# KIT-017: DIA stage closure

Status 2026-10-15: **complete; the closure pair measured, both members complete and verified; decision `rejected` on the per-step guard (the DA step 1.129 ×, the whole slice 1.077 ×) with the bound rule's required reduction met (resident end 0.927 ×).** The stage closes as measured-not-accepted on the pair rule, with the shortfall stated below, and as accepted on its standing claim against the production-default baselines. Ticket revision 1.1: control `kit_018` (`core/browser-smoke.yXPrcl`, `6a7304e7`), candidate `kit_017` (`core/browser-smoke.c1cAdn`, `7efb6ccd` = the `kit_024` source plus the screen-run and counted-run repairs), rule `resident_end`. Four counted journey members (two from the ledger, two added by the user on 2026-10-15); four screen lifecycle runs and one earlier counted attempt, all archived.

## What the closure measures and why

The pair is a protocol v1.4 ticket journey pair: each build one journey member as users run the GUI (no launch record) on the 300-sample corpus, in one session, the candidate's states verified against the KIT-005 oracles. The bound rule is `resident_end`: the one memory mechanism between the two builds is KIT-019's DA publication through the worker, whose claim is end-of-journey resident memory (candidate ≤ 0.95 × control), with the whole-slice anonymous peak not worse than 1.02 × and the per-step guard from the first saved state. The stage's whole-slice reduction is reported against the production-default baselines (KIT-005 `default-1` on `kit_001`, KIT-011's `kit_006` control; both guard readings, so lower bounds), as protocol section 4 says; no mechanism between `kit_018` and `kit_024` claims a whole-slice reduction, so the pair does not judge one ([audit](KIT-017-PREEXEC-AUDIT-20261015.md) PE-02).

## Before the counted launch

- [Pre-execution audit](KIT-017-PREEXEC-AUDIT-20261015.md) (ten findings) and its [remediation record](KIT-017-PREEXEC-REMEDIATION-20261015.md), which also carries the three screen lifecycle attempts. PE-01: the bindings committed with KIT-034 hashed the journey driver at its pre-KIT-034 content (regenerated). PE-03: neither build had been observed through a managed journey past design, so the pair ran on the screen corpus first (the harness pair branch takes a bound corpus key, the KIT-018 pattern).
- Screen attempt 1 failed both members and exposed two plain code paths only a live browser reaches: the server-side paged tables (KIT-024) resolved the state owner inside Shiny's request handler, which carries no reactive domain (SR-01), and the application monitor died reading the worker registration file in the window after the controller unlinked it, taking the application with it (SR-02). Attempt 2, on the repaired candidate, hung at the same point: the served page was a named data frame, so its rows serialised as JSON objects where the widget reads arrays, and DataTables' native error dialog blocked the renderer (SR-04). The managed driver now records and dismisses native dialogs and fails the next command by name (SR-05). Attempt 3 completed and verified on both builds (control 561 s, candidate 563 s; [decision](screen-lifecycle-attempt3-decision.rds)).
- Installs: `kit_017` first as `core/browser-smoke.UpIfdA` (`9293655c`), then `core/browser-smoke.sOtQzu` (`c9635938`), then `core/browser-smoke.c1cAdn` (`7efb6ccd`, the counted candidate); all recorded, none counted. The control build keeps SR-02's race; a recurrence there is a demonstrated instrumentation failure (one correction rerun of that member).
- Predeclared before the outcome was visible (remediation record, attempt 3): the DA publication worker is a second R process inside the measured slice (about 0.5 GiB on the screen corpus, where the DA step read 1.28 × and the whole slice 1.14 × while resident end read 0.825), so the DA step, which is the journey's whole-slice peak on 300 samples, and possibly the whole slice may exceed the frozen 1.02 bound. The rule is not changed; a rejection on that ground closes the stage as measured-not-accepted with the shortfall stated.

## Result

**First counted attempt** (12:04–15:09; summary [kit-journey-017-attempt1-stage.json](kit-journey-017-attempt1-stage.json); runtime `kit-journey-017.attempt1`): the control completed and verified after correction 1 (verifier only: the KIT-016 rendered volcano page reported as an additional product, [control-verifier-correction.R](control-verifier-correction.R)); the candidate failed at the DA step because the managed step handed the DA publication the 300-sample filtering progress as its completion context, above the 4 MiB bound (remediation record CP-02; repaired as `workflowManagedCompletionContext`, the peptide completion role only). Stage `invalid`; the user added two members ("amend and go"); the repair was observed live on the screen pair (attempt 4) before the rerun.

**Counted rerun** (19:44–22:55; tracked summary [kit-journey-017-stage.json](kit-journey-017-stage.json); 300-sample corpus, no launch record, 14,400 s ceiling, one session, control first):

| Member | Build | Outcome | Elapsed | Slice peak | Anonymous peak | Resident end | Max gap | Verified |
|---|---|---|---|---|---|---|---|---|
| control-1 | `kit_018` | complete, 18 steps | 4,730 s | 7.20 GiB | 5.67 GiB | 2.95 GiB | 110 ms | every state equal to the KIT-005 oracles; rendered volcano page reported as additional |
| candidate-1 | `kit_017` (`c1cAdn`) | complete, 18 steps | 4,698 s | 7.20 GiB | 6.11 GiB | 2.73 GiB | 190 ms | every state equal to the KIT-005 oracles; rendered volcano page reported as additional |

Decision (`kitJourneyTicketDecision`, rule `resident_end`, protocol v1.4): **`rejected`, `step_not_worse`.**

- Required reduction on resident end: **met**, 0.927 (2.95 → 2.73 GiB; KIT-005's noise on this metric is 0.983, so the reduction is four times the noise). This is the KIT-019 mechanism's claim: the application process no longer holds Glimma's namespace tree and the publication pages at the end of the journey.
- Whole-slice anonymous peak: **not within 1.02**, 1.077 (5.67 → 6.11 GiB). On the control the journey's peak is the P14 rollup (5.67 GiB); on the candidate it is the DA step (6.11 GiB).
- Per-step guard: **fails on the DA step only**, 1.129 (5.41 → 6.11 GiB, +0.70 GiB). The other guarded steps: q-value 0.973, rollup 0.972, IQ 0.969, normalisation 1.014, correlation 1.013 are within 1.02; session export reads 1.034 (2.07 → 2.14 GiB, +0.07 GiB), above the bound but inside the 3 % noise band KIT-005 recorded (0.971 on the whole-slice metric), and the decision names the DA step because the guard reports the first failing rule, not each step. Import and design (1.034, 1.040) are reported, not judged.

What the shortfall is: the DA publication worker is a second R process inside the measured slice. It loads the package namespace, Glimma and the hydrated DA result to build the pages, about 0.70 GiB on 300 samples (0.50 GiB on the screen corpus), while the application side sheds about 0.44 GiB. The mechanism moves memory from the resident application to a transient worker; the slice, which is what a user's machine must hold, peaks higher at that moment. Predeclared before the result was visible (remediation record, attempt 3): this reading closes the stage as measured-not-accepted on the pair rule; the rule is not changed. The follow-on is a mechanism, not a protocol change: a lighter publication worker that builds the pages from the saved tables without hydrating the whole DA result and without the application's namespace (KIT-036, authored from this outcome).

Time (reported, never judged): the two members ran 4,730 and 4,698 s; the normalisation and RUV sweep is 45 % of each (2,124 and 2,114 s).

**Stage claim against the production-default baselines** (protocol section 4; KIT-005 `default-1` and KIT-011's `kit_006` control both ended at the memory guard in imputation at 15.40 GiB anonymous after 8 of 18 steps): the final stage-1 build completes the 300-sample journey at 6.11 GiB anonymous, 0.40 × the guard reading, which is a lower bound of the in-memory default's true peak, and ends at 2.73 GiB resident. Per step against the guarded steps of the in-memory default (q-value 6.31, rollup 10.24, imputation 14.3 GiB): 4.85, 5.51, 4.32 GiB. This is the stage's standing whole-slice evidence, first recorded by KIT-011 and unchanged in kind by this pair.

**Partial evidence towards the DIA-IQ journey case** (PERF-128 CASE-PERF-JOURNEY-001, owner KIT-050): two complete, verified 300-sample journeys on the synthetic corpus within the 16 GiB cap; the case itself needs a non-synthetic corpus and is not claimed here.

## Regression (standing rule: the whole suite once per ticket, classified at the previous commit)

| Commit | Suite (343 files) | `R CMD check` (tarball from a clean worktree under `/home`, proof library on `R_LIBS`, suggests not forced) |
|---|---|---|
| `c439ce88` (binding, harness only) | 29,627 pass, 0 fail, 1 error, 139 skip | 3 warnings, 3 notes ([logs](check-logs/)) |
| `9293655c` (round-one repairs) | 29,631 pass, 0 fail, 1 error, 139 skip | 3 warnings, 3 notes |
| `c9635938` (round-two repairs, the candidate's source) | 29,632 pass, 0 fail, 1 error, 140 skip | 3 warnings, 3 notes |

The one error in every run is `test-workflow-perf101-diagnostic.R` "request source snapshots are self-contained and detect missing or changed bytes", which passes alone (20 pass) and in KIT-034's run: the regression script builds the tarball in the same worktree while that test snapshots source bytes. An artefact of the runner, not of the code; recorded, the runner left as is. KIT-002's check record had 4 warnings and 3 notes; the code/documentation-mismatch warning is gone and no warning or note is new.

## Kit completeness checklist

[kit-checklist.md](kit-checklist.md): every kit function exists with no DIA logic in its body (a naming debt recorded: the storage primitives keep the `DIA` prefix), all 16 data steps are registered declarations, `partitionByPlex` is exempt until KIT-020, the three kernel adapters are named, the quarantine state and the overtaken-tests record are stated.

## Hand-off

Stage 2 and stage 3 inputs are in the checklist's section G; the technical-replicate declaration (KIT-035), the lighter publication worker (KIT-036) and the analysis-record series (KIT-025 to KIT-029, KIT-033) follow this closure. The control stays `kit_018` (protocol section 2: the latest accepted build; this pair was rejected on the per-step guard, so `kit_017` is not promoted by it). Recommended order: KIT-036 next, measured as a pair against `kit_018` with its candidate the lighter worker on the `kit_017` source, so the stage ends on an accepted build before KIT-035 and the analysis-record series; the user decides.
