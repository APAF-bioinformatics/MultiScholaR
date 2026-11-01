#!/usr/bin/env bash
# KIT-072 Simulation Setup: every owner-or-double branch, `owned` sub-branch and double tolerance in the lipid/metab
# modules and the omics declaration/session helpers. Run before and after; the after-list keeps only refusals.
cd "$(git rev-parse --show-toplevel)"
grep -n -E 'inherits\([^)]*"DIAStateAdapter"\)|usesDIAReferenceStates|\bowned\b|isTRUE\(committed\)|!is\.null\(.*\$\.context\)' \
  R/mod_lipid_*.R R/mod_metab_*.R R/mod_lipidomics.R R/mod_metabolomics.R \
  R/utils_workflow_lipid_declarations.R R/utils_workflow_metab_declarations.R | grep -v '^\s*#'
echo "--- shared session helpers (omics) ---"
awk '/^[A-Za-z.][A-Za-z0-9_.]* <- function/{f=""} /^(bindWorkflowOwnerTypeAtDesign|workflowSummaryFinalStateName|ingestWorkflowOmicsLegacySession|resolveWorkflowOmicsOwnedSession) <- function/{f=$1} f && /inherits\([^)]*"DIAStateAdapter"\)/{print FILENAME":"FNR": ["f"] "$0}' R/utils_workflow_dia_sessions.R
