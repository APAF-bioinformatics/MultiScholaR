#!/usr/bin/env bash
# KIT-134 (the Linux library path): host-check.R with the installed application's R, on this
# machine and in containers of other distributions (network none, the application mounted
# read-only at a path with a space, an unprivileged user without a home of its own).
#
# Usage: linux-host-check.sh <installed application folder (holds app/resources)> <image>...
# Images are pulled by the caller. Prints one table; exits 1 when a host fails.
set -uo pipefail

installed="$1"
shift
checks_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
resources="$installed/app/resources"
summary="$(mktemp)"
status=0

record() { # label, exit status, log
    local verdict=passed
    if [ "$2" -ne 0 ]; then verdict=FAILED; status=1; fi
    printf '%s\t%s\t%s\t%s\n' "$1" "$verdict" \
        "$(grep -c '^ok  ' "$3")/$(grep -c -E '^(ok  |FAIL) ' "$3")" \
        "$(sed -n 's/^for the record, bash under the path exported before the repair: //p' "$3")" >>"$summary"
}

log="$(mktemp)"
label="this runner"
"$resources/runtime/bin/Rscript" --vanilla "$checks_dir/host-check.R" "$resources/runtime" "$resources/entry.R" "$label" 2>&1 | tee "$log"
record "$label" "${PIPESTATUS[0]}" "$log"

for image in "$@"; do
    log="$(mktemp)"
    docker run --rm --network none --user "$(id -u):$(id -g)" --env HOME=/tmp \
        --volume "$installed:/opt/msr installed:ro" --volume "$checks_dir:/checks:ro" "$image" \
        "/opt/msr installed/app/resources/runtime/bin/Rscript" --vanilla /checks/host-check.R \
        "/opt/msr installed/app/resources/runtime" "/opt/msr installed/app/resources/entry.R" "$image" 2>&1 | tee "$log"
    record "$image" "${PIPESTATUS[0]}" "$log"
done

echo "== hosts"
printf 'host\tresult\tchecks passed\tbash under the old exported path\n'
cat "$summary"
exit "$status"
