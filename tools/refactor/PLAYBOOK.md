# Refactor Playbook

This playbook turns the Wave 1 lessons into the default protocol for later waves.

The prioritized stabilization queue lives in
[GOD_MODULE_STABILIZATION_BACKLOG.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md:1).

## Core Rules

- Move code by exact source extraction, not by hand-rewriting function bodies.
- Stage every wave outside `R/` first.
- Treat `DESCRIPTION` `Collate:` as part of the refactor output.
- Assume tests, dev scripts, workbooks, and reports may contain filename-coupled references.
- Resolve compatibility aliases explicitly instead of carrying duplicate implementations.
- Enforce file-size budgets instead of treating “small enough” informally.

## Codex Skill Support

An installed Codex skill now exists for this workflow at:

- [god-module-stabilization](</home/doktersmol/.codex/skills/god-module-stabilization/SKILL.md>)

Use it in two modes:

- `survey`: identify oversized files, classify god modules, and generate/update handovers before structural work
- `stabilize`: run the staged stabilization protocol, including characterization, seam introduction, staged extraction, and wave apply helpers

The skill is an orchestrator around the existing repo tooling in `tools/refactor/`.
It does not replace `extract_blocks.R`, `verify_refactor.R`, or `check_wave_apply.R`.

## Refactor Guards

- Do not introduce new live `R/` helper files or rewritten helper functions during exploration.
- If a split is needed, extract exact code blocks into staging first, review there, and only then apply the reviewed result to live `R/`.
- Do not rename local symbols into a different style during stabilization. Follow the dominant naming convention already used in the touched area.
- For proteomics import and adjacent modules, new helper/function names must stay in the existing camelCase style unless a reference file explicitly says otherwise.
- If a test exposes a production bug, a direct live bugfix is acceptable, but structural refactors still must go through the staging-first extraction path above.

## File Size Budget

- Default target: `150-500` LOC.
- Acceptable upper band: `501-800` LOC.
- Soft cap: `801-1000` LOC.
- Above `1000` LOC: treat as a failed split unless there is a documented reason.
- Below `150` LOC: acceptable only for true leaf files, not artificial fragmentation.

Operational rule:

- Split by stable responsibility first.
- If a file still lands above `800-1000` LOC after a wave, queue a follow-up split immediately.
- Helper, render, result, plotting, IO, and compatibility files should usually stay under `500`.
- Orchestrators and tightly coupled S4 class/method groupings may occasionally land in the `500-800` band.

## Wave Protocol

1. Pre-wave audit

- Run the filename-coupling audit before touching the wave.
- Sweep `tests/`, `dev/`, `inst/`, `docs/`, and `Workbooks/` for direct references to `R/*.R`.
- Convert file-based checks to symbol-based checks where possible.

Commands:

```bash
Rscript tools/refactor/audit_refactor_coupling.R
Rscript tools/refactor/audit_refactor_coupling.R --output tools/refactor/AUDIT-filename-coupling.md
Rscript tools/refactor/audit_file_sizes.R --output tools/refactor/AUDIT-file-sizes.md
Rscript tools/refactor/audit_file_sizes.R --manifest tools/refactor/<manifest>.yml --output tools/refactor/AUDIT-<wave>-file-sizes.md
```

2. Manifest design

- Keep the manifest as the human-authored source of truth.
- Prefer `symbol` selectors first, `setMethod` / `setClass` / `setGeneric` second, and `anchor_range` only as a fallback.
- Mark real collisions as `manual_merge`.
- Document compatibility expectations in the manifest `note` field.

3. Pre-apply verification

- Verify selector coverage before generating files.
- Fail fast on ambiguous or missing selectors.

Command:

```bash
Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/<manifest>.yml
```

4. Stage generation

- Generate target files into a staging directory first.
- Review staged files before touching live package sources.
- Resolve every `manual_merge` in staging.

Commands:

```bash
Rscript tools/refactor/extract_blocks.R \
  --manifest tools/refactor/<manifest>.yml \
  --write-targets \
  --output-root tools/refactor/staging/<wave> \
  --emit-collate collate-<wave>.txt
```

5. Apply to live `R/`

- Only apply after staged manual merges are reviewed.
- If the extractor rewrites sources mechanically, copy the reviewed staged targets over the generated live targets afterward.
- Remove any leftover duplicate symbols from source files if they were part of a manual merge.

6. Post-apply verification

- Check that all target files exist.
- Check that moved selectors no longer exist in their source files.
- Check for duplicate top-level symbols within the wave scope.
- Parse the wave scope.

Command:

```bash
Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/<manifest>.yml
```

7. Package integration

- Update `DESCRIPTION` `Collate:` with the new target files in the correct order.
- Run `devtools::document()` in a dependency-complete environment.
- Run the tests most relevant to the wave.

## Touchpoints To Audit Every Wave

- `DESCRIPTION` `Collate:`
- `NAMESPACE` regeneration via `devtools::document()`
- `tests/`
- `dev/`
- `inst/reports/`
- `inst/workbooks/`
- `docs/`
- wrapper files carrying old naming conventions such as `de_*` vs `da_*`
- legacy source files left as placeholders after extraction

## Operational Lessons From Wave 1

- `manual_merge` alone is not enough. It identifies collisions, but it does not remove the old block from the source file.
- Staging-first prevented a live duplicate symbol from surviving into the package.
- Filename-coupled tests and scripts are a recurring tax. They must be audited before each wave, not discovered mid-apply.
- Load order is part of correctness in this package. Collate updates are not optional cleanup.
- Empty placeholder source files are acceptable as migration breadcrumbs, but they should be an explicit choice, not an accidental outcome.

## Current Baseline Findings

- The current filename-coupling audit is recorded in
  [AUDIT-filename-coupling.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/AUDIT-filename-coupling.md:1).
- The current file-size baseline should be tracked in `AUDIT-file-sizes.md`, and each wave should generate its own scoped file-size audit.
- The main remaining filename-coupling risk surface is `dev/`, especially the lipid dev scripts:
  - [dev/test_lipid_app.R](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_app.R:6)
  - [dev/test_lipid_core.R](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_core.R:7)
  - [dev/test_lipid_de.R](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_de.R:6)
  - [dev/test_lipid_qc.R](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_qc.R:6)
  - [dev/test_lipid_s4.R](/home/doktersmol/Documents/MultiScholaR/dev/test_lipid_s4.R:7)
- Some of those scripts already point at old names such as `func_lipid_de.R` and `mod_lipid_de.R`.
- Treat them as optional manual-dev cleanup, not as a mandatory pre-wave blocker for package stabilization or for separate-worktree lipid/metabolite pilot lanes.
- Current docs references are low-risk and informational, but they should still be updated once waves are stable.
- Wave 1 did not fully meet the size budget and should be treated as a structural first pass, not the final DA decomposition.

## Recommended Order For Later Waves

1. Proteomics vertical slices with the best regression surface.
2. Metabolomics and lipidomics in mirrored layouts.
   If run in parallel, do so only from separate git worktrees with isolated
   `.god-module-stabilization/` state.
3. Shared omics helpers.
4. Large general utility monoliths last, especially `func_general_filemgmt.R`.
