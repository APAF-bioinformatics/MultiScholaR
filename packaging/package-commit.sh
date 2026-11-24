# KIT-176: the package commit of a bundle build. Sourced by the bundle lane's rule step
# (.github/workflows/distributable.yml) and by its dry test of the rule (T3), so both run these lines.
# Defines PKG_PATHS and two functions; does nothing else.
#
# PKG_PATHS: what `R CMD INSTALL` of the checkout reads for code and data (packaging/r/install-library.R).
# NEWS.md, LICENSE and README.md are also copied by R and are left out on purpose (text only).
#
# package_tree_checks <expected head>: in the current directory, stops unless the clone has its full history,
#   HEAD is <expected head>, and no tracked file under PKG_PATHS is changed and no untracked file is there.
# package_commit_rule <head> <line>: prints "<commit> <true|false>": the merge base of <head> and <line> and true
#   when the package source there equals <head>'s; otherwise <head>'s commit and false. Stops on a missing line,
#   no common history, several merge bases or a git error.
# Every stop is an explicit `exit 1` with a message on standard error: the rule runs inside a command
# substitution, where `bash -e` is not inherited.

PKG_PATHS="DESCRIPTION NAMESPACE R inst src man data configure configure.win configure.ucrt demo exec .Rinstignore"

package_tree_checks() {
    local expected="$1" shallow head dirty
    shallow=$(git rev-parse --is-shallow-repository) || { echo "package_tree_checks: git rev-parse failed" >&2; exit 1; }
    if [ "$shallow" != false ]; then
        echo "package_tree_checks: the clone is shallow (is-shallow-repository: $shallow); the package commit needs the full history" >&2
        exit 1
    fi
    head=$(git rev-parse HEAD) || { echo "package_tree_checks: git rev-parse HEAD failed" >&2; exit 1; }
    if [ "$head" != "$expected" ]; then
        echo "package_tree_checks: HEAD is $head, not the expected $expected" >&2
        exit 1
    fi
    # shellcheck disable=SC2086
    dirty=$(git status --porcelain --untracked-files=all -- $PKG_PATHS) || { echo "package_tree_checks: git status failed" >&2; exit 1; }
    if [ -n "$dirty" ]; then
        echo "package_tree_checks: the package source differs from $head (git status of the package paths):" >&2
        echo "$dirty" >&2
        exit 1
    fi
}

package_commit_rule() {
    local head="$1" line="$2" B bases n s commit
    B=$(git merge-base "$head" "$line")
    if ! [[ "$B" =~ ^[0-9a-f]{40}$ ]]; then
        echo "package_commit_rule: no merge base of $head and $line (got '$B')" >&2
        exit 1
    fi
    bases=$(git merge-base --all "$head" "$line") || { echo "package_commit_rule: git merge-base --all failed" >&2; exit 1; }
    n=$(printf '%s\n' "$bases" | wc -l | tr -d ' ')
    if [ "$n" != 1 ]; then
        echo "package_commit_rule: $head and $line have $n merge bases (criss-cross history)" >&2
        exit 1
    fi
    s=0
    # shellcheck disable=SC2086
    git diff --quiet "$B" "$head" -- $PKG_PATHS || s=$?
    if [ "$s" = 0 ]; then
        echo "$B true"
    elif [ "$s" = 1 ]; then
        commit=$(git rev-parse "$head") || { echo "package_commit_rule: git rev-parse $head failed" >&2; exit 1; }
        echo "$commit false"
    else
        echo "package_commit_rule: git diff failed (status $s)" >&2
        exit 1
    fi
}
