#!/usr/bin/env bash
# KIT-074 inventory: the DIA-named kit functions the omics code calls (rule a) or a later stage-3b ticket names (rule b),
# with their neutral names and the callers of each name by omic; then every DIA-named call left in the omics files.
# Usage: inventory.sh <out.tsv>   (run from the repository root). Exit 0 always; the residual list is printed.
set -euo pipefail
out=${1:?out.tsv}
lipid=(R/mod_lipid_*.R R/mod_lipidomics.R R/utils_workflow_lipid_declarations.R)
metab=(R/mod_metab_*.R R/mod_metabolomics.R R/utils_workflow_metab_declarations.R)
core=(R/utils_workflow_omics_declarations.R)
# The record's omics branches: the step replay and the DA re-run keyed by the omics bindings (KIT-073).
record_branches=$(awk '/else if \(!is.null\(omicsDeclarationBindings\(\)/{on=1} /if \(!is.null\(omics_binding\)\)/{on=1}
    on{print} on && /^ *} else \{/{on=0}' R/utils_workflow_analysis_record.R)
rows='saveWorkflowStateValue	saveDIAStateValue	a
recordWorkflowInputs	recordDIAWorkflowInputs	a
clearWorkflowStepProvenance	clearDIAStepProvenance	a
collectWorkflowNativeReferences	collectDIANativeReferences	a
runWorkflowOwnedSessionExport	runDIAOwnedSessionExport	a
createWorkflowSessionOwner	createDIASessionOwner	a
saveWorkflowOwnedValue	saveDIAOwnedValue	kit045
publishWorkflowStateHead	publishDIAStateHead	kit045
withWorkflowStateScope	withDIAStateScope	kit045
selectWorkflowStateReference	selectDIAStateReference	kit045
closeWorkflowStateOwner	closeDIAStateAdapter	kit045
readWorkflowStateGuiMetadata	readDIAStateGuiMetadata	b:KIT-075
getWorkflowOwnedPayload	getDIAOwnedPayload	b:KIT-075
setWorkflowOwnedPayload	setDIAOwnedPayload	b:KIT-075
getWorkflowOwnedInput	getDIAOwnedInput	b:KIT-075
resolveWorkflowOwnedValue	resolveDIAOwnedValue	b:KIT-075
getWorkflowReactiveReference	getDIAReactiveReference	b:KIT-075
assertWorkflowExportReady	assertDIAExportReady	b:KIT-078
publishWorkflowActionOutputs	publishDIAActionOutputs	b:KIT-078'
count() { local name=$1; shift; { grep -ohE "\b$name\b" "$@" 2>/dev/null || true; } | wc -l; }
calls_in() { local name=$1; printf '%s\n' "$record_branches" | { grep -oE "\b$name\b" || true; } | wc -l; }
others() { local name=$1; { grep -lE "\b$name\b" R/*.R 2>/dev/null || true; } | grep -vxF -f <(printf '%s\n' "${lipid[@]}" "${metab[@]}" "${core[@]}") | wc -l; }
defined() { { grep -lE "^$1 *<- *function" R/*.R || true; } | head -1; }
printf 'neutral\tdia_alias\trule\tdefined_in\tlipid_dia\tmetab_dia\tcore_dia\trecord_omics_dia\tlipid_neutral\tmetab_neutral\tcore_neutral\tother_R_files_dia\n' > "$out"
while IFS=$'\t' read -r neutral dia rule; do
    def=$(defined "$neutral"); [ -n "$def" ] || def="$(defined "$dia") (as $dia)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$neutral" "$dia" "$rule" "$def" \
        "$(count "$dia" "${lipid[@]}")" "$(count "$dia" "${metab[@]}")" "$(count "$dia" "${core[@]}")" "$(calls_in "$dia")" \
        "$(count "$neutral" "${lipid[@]}")" "$(count "$neutral" "${metab[@]}")" "$(count "$neutral" "${core[@]}")" "$(others "$dia")" >> "$out"
done <<< "$rows"
echo "residual DIA-named calls in the omics files (modalDialog, *Diagnostics and the MS-DIAL importers are not kit names):"
{ grep -nE '\b[A-Za-z_.]*(DIA|Dia)[A-Za-z_.]*\b *\(' "${lipid[@]}" "${metab[@]}" "${core[@]}" || true; } | grep -vE 'modalDialog|Diagnostics|MSDIAL' || echo "  none"
echo "record omics branches:"; printf '%s\n' "$record_branches" | { grep -nE '\b[A-Za-z_.]*(DIA|Dia)[A-Za-z_.]*\b *\(' || echo "  none"; }
