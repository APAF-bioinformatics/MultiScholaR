#!/usr/bin/env bash
# KIT-159 section 5.2: the pre-push allow-list and hash check of an evidence folder, without R.
#
# Usage: check-evidence.sh <folder>
#   1. every file matches outputs.sha256 (sha256sum -c), and outputs.sha256 lists exactly the files;
#   2. every file name is on the allow-list (the fixed names of section 5.7; ledger-*.tsv and
#      members-*.tsv; scored.tsv.gz; ticket/ paths only when listed in allowed-paths.txt, which the
#      job writes from the declaration's expanded outputs);
#   3. the fixed tables' header lines are exactly their columns;
#   4. no file over 5 MB, no archive or binary other than scored.tsv.gz;
#   5. a hand-off folder (ticket/handoff/) matches its handoff.sha256.
# The column lists repeat msrval's EVIDENCE_COLUMNS (licence.R); a change there changes them here.
set -u
dir="${1:?usage: check-evidence.sh <folder>}"
fail() { echo "check-evidence: $*" >&2; exit 1; }
cd "$dir" || fail "no folder $dir"
[ -f outputs.sha256 ] || fail "outputs.sha256 is missing"
sha256sum --quiet --strict -c outputs.sha256 >/dev/null 2>&1 || fail "a file differs from outputs.sha256"
listed="$(sed -E 's/^[0-9a-f]{64}  //' outputs.sha256 | LC_ALL=C sort)"
present="$(find . -type f ! -name outputs.sha256 | sed 's#^\./##' | LC_ALL=C sort)"
[ "$listed" = "$present" ] || fail "outputs.sha256 does not list exactly the folder's files"

join() { local IFS=$'\t'; echo "$*"; }
expected_header() {
    case "$1" in
        cells.tsv) join unit run arm exit stop_reason step class message_sha256 fallback \
            options_at_start ruv_k ruv_controls seconds chunk ;;
        calls.tsv) join unit arm contrast call_statistic fallback_count ;;
        metrics.tsv) join ticket dataset truth_id unit contrast arm assay view class call_rule metric \
            statistic_class value method form scope ci_low ci_high achieved_level n_units unit_kind B \
            seed na_share status primary family blinded ;;
        comparisons.tsv) join arm comparator truth_id metric view call_rule diff method scope ci_low \
            ci_high margin conf outcome no_worse outcome_family n_refused regime flags member \
            component not_sensitive amended headline blinded ;;
        coverage.tsv) join truth_id unit arm contrast class keys quantified testable called_own \
            called_bh completeness ;;
        attribution.tsv) join role kind licence redistribution citation ;;
        determinism.tsv) join table sha256_job1 sha256_job2 identical sha256_10_job1 sha256_10_job2 \
            max_abs_diff max_rel_diff cause ;;
        ticket-coverage.tsv) join regime cell method class form g level simulations bound coverage \
            mc_se miss critical pass ;;
        ticket-coverage-cells.tsv) join regime cell parameter value ;;
        ledger-*.tsv) join step action input_id value source location ;;
        parity-da.tsv) join route comparison status digest sha256_10 ;;
        uniprot.tsv) join taxon accessions rows release release_source sha256 ;;
        parity-*.tsv) join action input_id driver_value browser_value equal excluded_class ;;
        *) echo "" ;;
    esac
}

while IFS= read -r f; do
    [ -n "$f" ] || continue
    size="$(stat -c %s "$f")"
    [ "$size" -le 5242880 ] || fail "$f is over 5 MB"
    case "$f" in
        scored.tsv.gz) continue ;;
        ticket/handoff/*) continue ;;
        *.zip|*.gz|*.tgz|*.tar|*.bz2|*.xz|*.7z|*.rar|*.rds|*.RDS|*.rdata|*.RData|*.so|*.o|*.a|*.dll)
            fail "$f is an archive or binary" ;;
    esac
    base="${f##*/}"
    want="$(expected_header "$f")"
    [ -n "$want" ] || want="$(expected_header "$base")"
    if [ -n "$want" ] && [ "$f" = "$base" ]; then
        got="$(head -n 1 "$f")"
        [ "$got" = "$want" ] || fail "$f does not have its columns"
        continue
    fi
    case "$f" in
        provenance.json|results.md|tables.sha256|not-committed.sha256|inputs.tsv|library.tsv|\
        coverage-table.tsv|coverage-detail.tsv|coverage-run.txt|coverage-tasks.tsv|\
        harness-tests.tsv|harness-failures.txt|ticket-tests.tsv|ticket-failures.txt|\
        code-hashes.tsv|fsim.tsv|cost-probe.tsv|allowed-paths.txt|pre-outputs.tsv|members-*.tsv)
            continue ;;
    esac
    if [ -f allowed-paths.txt ] && grep -qxF -- "$f" allowed-paths.txt; then continue; fi
    fail "$f is not on the allow-list"
done <<< "$present"

if [ -d ticket/handoff ]; then
    [ -f ticket/handoff/handoff.sha256 ] || fail "the hand-off folder has no handoff.sha256"
    (cd ticket/handoff && sha256sum --quiet --strict -c handoff.sha256 >/dev/null 2>&1) \
        || fail "the hand-off folder differs from handoff.sha256"
fi
echo "check-evidence: ok"
