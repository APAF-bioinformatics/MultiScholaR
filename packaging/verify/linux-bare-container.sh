#!/usr/bin/env bash
# KIT-134: run one check script with the bundled R inside a bare ubuntu:22.04 container. The image
# has no R, no Python and none of R's shared libraries; the bundle is mounted read-only at a path
# of the caller's choosing; with network "none" the process cannot reach anything.
# From the KIT-130 spike (tools/dist-spike/checks/linux-in-container.sh).
#
# Usage: linux-bare-container.sh <none|bridge> <bundle on the host> <mount point> <check>.R [arguments...]
# The check is taken from this folder (packaging/verify) and receives the mount point first.
set -euo pipefail

network="$1"
bundle="$2"
mount_point="$3"
check="$4"
shift 4
checks_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

docker run --rm \
    --network "$network" \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp \
    --volume "$bundle:$mount_point:ro" \
    --volume "$checks_dir:/checks:ro" \
    ubuntu:22.04 \
    "$mount_point/bin/Rscript" --vanilla "/checks/$check" "$mount_point" "$@"
