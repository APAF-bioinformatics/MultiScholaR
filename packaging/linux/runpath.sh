#!/usr/bin/env bash
# KIT-134, build time only: write the run path to the bundle's syslib/ into shared objects, relative
# to each object's own place ($ORIGIN), in front of whatever run path it had. DT_RUNPATH serves only
# the object that carries it, so every object that needs a bundled library gets one. An object that
# already has the entry is left untouched, so a second pass changes nothing.
#
# Usage: runpath.sh <syslib directory (need not exist yet)> <file or directory to scan>...
# Run twice per build: by packaging/r/install-library.R on duckdb, duckplyr and fastduckdb at the
# build prefix before MultiScholaR is installed (its install records the sha256 of their binaries,
# src/build_dia_rsum.R, and the DIA rollup and the direct writer refuse binaries changed after
# that), and by assemble.sh on the whole bundle. Never on a library a running R has loaded:
# patchelf rewrites the file in place.
set -euo pipefail

syslib="$1"
shift
patched=0
kept=0
while IFS= read -r -d '' object; do
    old="$(patchelf --print-rpath "$object" 2>/dev/null)" || continue # not a dynamic ELF object
    new="\$ORIGIN/$(realpath -m --relative-to="$(dirname "$object")" "$syslib")"
    case ":$old:" in *":$new:"*) kept=$((kept + 1)); continue ;; esac
    chmod u+w "$object"
    patchelf --set-rpath "$new${old:+:$old}" "$object" || { echo "patchelf could not write the run path of $object" >&2; exit 1; }
    patched=$((patched + 1))
done < <(find "$@" -type f \( -name '*.so' -o -name '*.so.*' -o -path '*/bin/exec/R' \) -print0)
echo "== run path to syslib written into $patched objects ($kept had it already)"
