#!/usr/bin/env bash
# KIT-134, build time only: the steps that make the Linux bundle, one phase per workflow step
# (.github/workflows/distributable.yml). On the ubuntu-22.04 runner they run on the runner itself;
# on any other runner they run in an ubuntu:22.04 container (MSR_BUILD_CONTAINER=1, as root, the
# runner's temporary folder and the checkout mounted at their own paths), so the bundle is built
# against Ubuntu 22.04 wherever it is then installed and tested.
#
# Usage: build.sh <tools|r|plan|sysreqs|rust|library|versions|python|pandoc|assemble>
# Environment: R_VERSION, RUNNER_TEMP, BUNDLE, TOOLS_LIB, OUT, RLIB, PBS_RELEASE, PBS_PYTHON,
# PANDOC_VERSION, GITHUB_STEP_SUMMARY; for library: MSR_INTERNAL and the read token in GITHUB_PAT
# (packaging/r/install-library.R reads it once and removes it from its environment).
set -euo pipefail

phase="$1"
if [ -d "$HOME/.cargo/bin" ]; then PATH="$HOME/.cargo/bin:$PATH"; fi
r_prefix="/opt/R/${R_VERSION}"
rscript="$r_prefix/bin/Rscript"

case "$phase" in
tools)
    if [ "${MSR_BUILD_CONTAINER:-}" = 1 ]; then
        # What the hosted runner's image already has and a bare ubuntu:22.04 lacks.
        apt-get update -qq
        apt-get install -y --no-install-recommends sudo ca-certificates curl jq git build-essential gfortran cmake \
            pkg-config patchelf file unzip xz-utils locales tzdata
        git config --global --add safe.directory '*'
    else
        command -v patchelf >/dev/null || { sudo apt-get update -qq && sudo apt-get install -y --no-install-recommends patchelf; }
    fi
    patchelf --version
    echo "runner's own R: $(command -v R || echo none); Python: $(command -v python3 || echo none)"
    ;;
r)
    curl -fsSL -o "$RUNNER_TEMP/r.deb" "https://cdn.posit.co/r/ubuntu-2204/pkgs/r-${R_VERSION}_1_amd64.deb"
    sudo apt-get update -qq
    sudo apt-get install -y --no-install-recommends "$RUNNER_TEMP/r.deb"
    sudo chown -R "$(id -u):$(id -g)" "$r_prefix"
    mkdir -p "$TOOLS_LIB" "$OUT"
    "$rscript" --vanilla -e 'install.packages("jsonlite", lib = commandArgs(TRUE)[[1]], repos = "https://packagemanager.posit.co/cran/latest")' "$TOOLS_LIB"
    ;;
plan)
    R_LIBS="$TOOLS_LIB" "$rscript" --vanilla packaging/r/install-library.R plan . "$OUT"
    ;;
sysreqs)
    query() { # repository, then package names; prints the apt packages
        local repo="$1"; shift
        local extra=""; [ "$repo" = bioconductor ] && extra="&bioc_version=3.23"
        local url="https://packagemanager.posit.co/__api__/repos/$repo/sysreqs?all=false&distribution=ubuntu&release=22.04$extra$(printf '&pkgname=%s' "$@")"
        if ! curl -fsS "$url" | jq -r '.requirements[].requirements.packages[]?' 2>/dev/null; then
            # an unknown name fails the whole query: ask one at a time and say which were skipped
            if [ "$#" -gt 1 ]; then for p in "$@"; do query "$repo" "$p"; done
            else echo "no system requirements known for $repo/$1" >&2; fi
        fi
    }
    : > "$RUNNER_TEMP/apt.txt"
    for pair in CRAN:cran Bioconductor:bioconductor; do
        mapfile -t names < <(awk -F'\t' -v r="${pair%%:*}" 'NR > 1 && $2 == r { print $1 }' "$OUT/plan-packages.tsv")
        for ((i = 0; i < ${#names[@]}; i += 60)); do query "${pair##*:}" "${names[@]:i:60}" >> "$RUNNER_TEMP/apt.txt"; done
    done
    sort -u "$RUNNER_TEMP/apt.txt" | tee "$OUT/apt-packages.txt"
    # shellcheck disable=SC2046
    sudo apt-get install -y --no-install-recommends $(sort -u "$RUNNER_TEMP/apt.txt")
    ;;
rust)
    # The internal packages need rustc 1.88 or later.
    if ! command -v rustup >/dev/null; then
        curl --proto '=https' --tlsv1.2 -fsSL https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain none
        PATH="$HOME/.cargo/bin:$PATH"
    fi
    rustup toolchain install stable --profile minimal --no-self-update
    rustup default stable
    rustc --version; cargo --version
    ;;
library)
    MSR_BUNDLE_RUNPATH=1 R_LIBS="$TOOLS_LIB" "$rscript" --vanilla packaging/r/install-library.R install . "$OUT" \
        2>&1 | tee "$OUT/install-library.log"
    grep '^INTERNAL ' "$OUT/install-library.log" || true
    echo "installed as binary packages: $(grep -c 'installing \*binary\* package' "$OUT/install-library.log" || true)"
    echo "built from source:"; grep 'installing \*source\* package' "$OUT/install-library.log" || true
    ;;
versions)
    R_LIBS="$TOOLS_LIB" "$rscript" --vanilla packaging/r/version-table.R . "$RLIB" \
        "$OUT/bundle-sources.tsv" "$OUT/install-library.log" "$OUT" 2>&1 | tee "$OUT/version-table.log"
    cat "$OUT/version-table.md" >> "$GITHUB_STEP_SUMMARY"
    ;;
python)
    mkdir -p "$BUNDLE"
    curl -fsSL -o "$RUNNER_TEMP/python.tar.gz" \
        "https://github.com/astral-sh/python-build-standalone/releases/download/${PBS_RELEASE}/cpython-${PBS_PYTHON}%2B${PBS_RELEASE}-x86_64-unknown-linux-gnu-install_only.tar.gz"
    tar -xzf "$RUNNER_TEMP/python.tar.gz" -C "$BUNDLE"
    PYTHONNOUSERSITE=1 "$BUNDLE/python/bin/python3" -m pip install --no-cache-dir --no-warn-script-location \
        --require-hashes --only-binary=:all: --no-deps -r packaging/python/requirements.txt
    "$BUNDLE/python/bin/python3" -m pip list | tee "$OUT/python-packages.txt"
    ;;
pandoc)
    curl -fsSL -o "$RUNNER_TEMP/pandoc.tar.gz" \
        "https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-linux-amd64.tar.gz"
    mkdir -p "$RUNNER_TEMP/pandoc" "$BUNDLE/pandoc"
    tar -xzf "$RUNNER_TEMP/pandoc.tar.gz" -C "$RUNNER_TEMP/pandoc"
    install -m 0755 "$RUNNER_TEMP/pandoc/pandoc-${PANDOC_VERSION}/bin/pandoc" "$BUNDLE/pandoc/pandoc"
    "$BUNDLE/pandoc/pandoc" --version | head -n 1
    ;;
diamond)
    # KIT-167: DIAMOND from the pin file (inst/extdata/diamond-pin.tsv): bytes and SHA-256 of the release asset, then
    # the executable's SHA-256; a mismatch stops the build.
    pin_row="$(awk -F'\t' '$1 == "linux"' inst/extdata/diamond-pin.tsv)"
    pin_url="$(printf '%s' "$pin_row" | cut -f3)"; pin_bytes="$(printf '%s' "$pin_row" | cut -f4)"
    pin_sha="$(printf '%s' "$pin_row" | cut -f5)"; pin_exe_sha="$(printf '%s' "$pin_row" | cut -f7)"
    curl -fsSL -o "$RUNNER_TEMP/diamond.tar.gz" "$pin_url"
    [ "$(stat -c %s "$RUNNER_TEMP/diamond.tar.gz")" = "$pin_bytes" ] || { echo "DIAMOND archive size differs from the pin" >&2; exit 1; }
    echo "$pin_sha  $RUNNER_TEMP/diamond.tar.gz" | sha256sum -c - || { echo "DIAMOND archive SHA-256 differs from the pin" >&2; exit 1; }
    mkdir -p "$RUNNER_TEMP/diamond" "$BUNDLE/diamond"
    tar -xzf "$RUNNER_TEMP/diamond.tar.gz" -C "$RUNNER_TEMP/diamond"
    install -m 0755 "$RUNNER_TEMP/diamond/diamond" "$BUNDLE/diamond/diamond"
    echo "$pin_exe_sha  $BUNDLE/diamond/diamond" | sha256sum -c - || { echo "DIAMOND executable SHA-256 differs from the pin" >&2; exit 1; }
    "$BUNDLE/diamond/diamond" version
    ;;
assemble)
    MSR_LDPATHS_ORIGINAL="$OUT/ldpaths.original" bash packaging/linux/assemble.sh "$r_prefix/lib/R" "$BUNDLE" 2>&1 | tee "$OUT/assemble.log"
    sudo rm -rf "$r_prefix"
    {
        echo "### Linux bundle, size (MiB)"
        echo '```'
        for part in "" /R /R/library /python /syslib /pandoc /diamond; do du -sm "$BUNDLE$part"; done
        echo "files: $(find "$BUNDLE" -type f | wc -l)"
        echo "largest R packages:"; du -sm "$BUNDLE"/R/library/* | sort -rn | sed -n 1,12p
        echo '```'
    } | tee -a "$GITHUB_STEP_SUMMARY" "$OUT/sizes.txt"
    ;;
*)
    echo "build.sh: unknown phase $phase" >&2
    exit 2
    ;;
esac
