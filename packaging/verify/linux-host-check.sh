#!/usr/bin/env bash
# KIT-134 (the Linux library path): host-check.R with the installed application's R, on this
# machine and in containers of other distributions (network none, the application mounted
# read-only at a path with a space, an unprivileged user without a home of its own).
#
# Usage: linux-host-check.sh <installed application folder (holds app/resources)> <phase> <image>...
#        linux-host-check.sh --no-which <installed application folder> <image>...
# Images are pulled (and, for the first form, given /usr/bin/which) by the caller. <phase> labels
# the rows (for example "before" and "after" a change to the installed copy). The first form prints
# one table and the loaded-library rows (LIBS); it exits 1 when a host fails. The second form starts
# the installed launcher in images that have no /usr/bin/which (decision M25): it passes when the
# launcher stops with exit status 69 and its message naming the requirement.
set -uo pipefail

if [ "$1" = --no-which ]; then
    installed="$2"
    shift 2
    status=0
    for image in "$@"; do
        out="$(docker run --rm --network none --user "$(id -u):$(id -g)" --env HOME=/tmp \
            --volume "$installed:/opt/msr installed:ro" "$image" \
            sh -c 'if [ -e /usr/bin/which ]; then echo "IMAGE HAS /usr/bin/which"; exit 99; fi; exec "/opt/msr installed/app/resources/runtime/bin/Rscript" --vanilla -e "cat(\"R STARTED\n\")"' 2>&1)"
        code=$?
        printf '%s\n' "$out" | sed "s|^|[$image, no which] |"
        if [ "$code" = 69 ] && printf '%s' "$out" | grep -q 'no /usr/bin/which' && ! printf '%s' "$out" | grep -q 'R STARTED'; then
            printf 'NOWHICH\t%s\tpassed\texit %s, the message, R not started\n' "$image" "$code"
        else
            printf 'NOWHICH\t%s\tFAILED\texit %s\n' "$image" "$code"
            status=1
        fi
    done
    exit "$status"
fi

installed="$1"
phase="$2"
shift 2
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
label="this runner ($phase)"
"$resources/runtime/bin/Rscript" --vanilla "$checks_dir/host-check.R" "$resources/runtime" "$resources/entry.R" "$label" 2>&1 | tee "$log"
record "$label" "${PIPESTATUS[0]}" "$log"

for image in "$@"; do
    log="$(mktemp)"
    label="$image ($phase)"
    docker run --rm --network none --user "$(id -u):$(id -g)" --env HOME=/tmp \
        --volume "$installed:/opt/msr installed:ro" --volume "$checks_dir:/checks:ro" "$image" \
        "/opt/msr installed/app/resources/runtime/bin/Rscript" --vanilla /checks/host-check.R \
        "/opt/msr installed/app/resources/runtime" "/opt/msr installed/app/resources/entry.R" "$label" 2>&1 | tee "$log"
    record "$label" "${PIPESTATUS[0]}" "$log"
done

echo "== hosts ($phase)"
printf 'host\tresult\tchecks passed\tbash under the old exported path\n'
cat "$summary"
exit "$status"
