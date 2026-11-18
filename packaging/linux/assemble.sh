#!/usr/bin/env bash
# KIT-134 P2, build time only: turn an R installed at its build prefix, with the bundle's library
# installed into it, into a folder that runs from anywhere. From the KIT-130 spike
# (tools/dist-spike/build/linux-assemble.sh on branch kit-130-dist-spike).
#
# Usage: assemble.sh <installed R home, e.g. /opt/R/4.6.1/lib/R> <bundle directory>
# The bundle directory may already hold python/ and pandoc/. The result is
#   <bundle>/R       the R home, its front-end script patched to locate itself
#   <bundle>/syslib  every shared library R and the packages need, except the C library; R's binary,
#                    its modules, every package library and the libraries of syslib itself carry a
#                    run path to it relative to their own place ($ORIGIN), so nothing exports
#                    LD_LIBRARY_PATH to the host programs R starts (needs patchelf at build time)
#   <bundle>/bin     the launcher (packaging/linux/Rscript)
set -euo pipefail

r_home_src="$1"
bundle="$2"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_prefix="$(dirname "$(dirname "$r_home_src")")" # /opt/R/4.6.1

mkdir -p "$bundle/bin" "$bundle/syslib"
cp -a "$r_home_src" "$bundle/R"

# 1. The front-end script carries the build prefix as a literal. Make it compute its own home,
#    and make every other mention of the prefix point nowhere, so a copy of R that happens to
#    sit at the build prefix on the user's machine can never be picked up.
front="$bundle/R/bin/R"
echo "== build-prefix literals in the front-end script before patching"
grep -n "$build_prefix" "$front" || true
sed -i \
    -e 's|^R_HOME_DIR=.*$|R_HOME_DIR="$(cd "$(dirname "$(readlink -f "$0")")/.." \&\& pwd -P)"|' \
    -e "s|$r_home_src|\${R_HOME_DIR}|g" \
    -e "s|$build_prefix|/nonexistent-r-build-prefix|g" \
    "$front"
if grep -q "$build_prefix" "$front"; then
    echo "the front-end script still names the build prefix" >&2
    exit 1
fi
echo "== other files under etc/ and bin/ that name the build prefix (recorded, not patched)"
grep -rl "$build_prefix" "$bundle/R/etc" "$bundle/R/bin" || true

# 2. Shared libraries: everything the R binary, its modules and the package libraries resolve to
#    outside the bundle, except glibc itself (the one thing taken from the host).
glibc='^(ld-linux-x86-64\.so\.2|libc\.so\.6|libm\.so\.6|libdl\.so\.2|libpthread\.so\.0|librt\.so\.1|libresolv\.so\.2|libutil\.so\.1|libanl\.so\.1|libmvec\.so\.1|libBrokenLocale\.so\.1|libthread_db\.so\.1|libnss_.*)$'
ldd_out="$(mktemp)"
find "$bundle/R" -type f \( -name '*.so' -o -name '*.so.*' -o -path '*/bin/exec/R' \) -print0 |
    LD_LIBRARY_PATH="$bundle/R/lib" xargs -0 -n 40 ldd >"$ldd_out" 2>/dev/null || true
echo "== unresolved at build time"
grep 'not found' "$ldd_out" | sort -u || true
# Each library's Ubuntu package, for the third-party notices (packaging/r/notices.R).
printf 'file\tpackage\tversion\tsource_package\n' >"$bundle/syslib/PACKAGES.tsv"
awk '/=> \//{print $3}' "$ldd_out" | sort -u | while read -r lib; do
    case "$lib" in "$bundle"/*) continue ;; esac
    if basename "$lib" | grep -Eq "$glibc"; then continue; fi
    cp -L "$lib" "$bundle/syslib/"
    real="$(readlink -f "$lib")"
    owner="$(dpkg-query -S "$real" 2>/dev/null || dpkg-query -S "${real#/usr}" 2>/dev/null || true)"
    package="$(printf '%s' "$owner" | head -n 1 | cut -d: -f1)"
    if [ -n "$package" ]; then
        printf '%s\t%s\t%s\t%s\n' "$(basename "$lib")" "$package" \
            "$(dpkg-query -W -f='${Version}' "$package" 2>/dev/null)" \
            "$(dpkg-query -W -f='${source:Package}' "$package" 2>/dev/null)" >>"$bundle/syslib/PACKAGES.tsv"
    else
        printf '%s\tunknown\t\t\n' "$(basename "$lib")" >>"$bundle/syslib/PACKAGES.tsv"
    fi
done
echo "== bundled shared libraries: $(find "$bundle/syslib" -type f | wc -l)"
ls "$bundle/syslib"

# 3. The run path (runpath.sh). Each object finds syslib/ from its own location; every object gets
#    one, the bundled libraries included (they need each other). duckdb, duckplyr and fastduckdb
#    have theirs since before MultiScholaR was installed and are not touched again. The launcher
#    then sets no LD_LIBRARY_PATH.
bash "$here/runpath.sh" "$bundle/syslib" "$bundle/R" "$bundle/syslib"
#    R's etc/ldpaths (sourced by the front-end script for every R, the systemd workers' included)
#    puts the host's library folders after R's own lib folder in the LD_LIBRARY_PATH it exports, and
#    that path is searched before any run path: on Debian-family hosts the host's libreadline,
#    libcurl and libgfortran then won over the bundled ones (outcome.md, defect 2). Decision M26
#    (the user, 2026-11-18): R's default keeps only its own lib folder (libR, libRblas, libRlapack),
#    so the bundled copies are taken on every host. Nothing from syslib/ is exported, so the host
#    programs R starts still load the host's libraries. A user-set R_LD_LIBRARY_PATH is removed by
#    the launcher. The original is kept outside the bundle (MSR_LDPATHS_ORIGINAL) for the
#    before-and-after measurement of the host check.
echo "== etc/ldpaths: R's default library path before"
grep -n 'R_LD_LIBRARY_PATH=' "$bundle/R/etc/ldpaths" || true
if [ -n "${MSR_LDPATHS_ORIGINAL:-}" ]; then cp "$bundle/R/etc/ldpaths" "$MSR_LDPATHS_ORIGINAL"; fi
sed -i 's|^\(: \${R_LD_LIBRARY_PATH=\${R_HOME}/lib\):[^}]*}|\1}|' "$bundle/R/etc/ldpaths"
echo "== etc/ldpaths: after (only R's own lib folder)"
grep -n 'R_LD_LIBRARY_PATH=' "$bundle/R/etc/ldpaths" || true
if grep -E '^: \$\{R_LD_LIBRARY_PATH=' "$bundle/R/etc/ldpaths" | grep -vqE '^: \$\{R_LD_LIBRARY_PATH=\$\{R_HOME\}/lib\}$'; then
    echo "etc/ldpaths still names a host folder in R's default library path" >&2
    exit 1
fi
grep -qE '^: \$\{R_LD_LIBRARY_PATH=\$\{R_HOME\}/lib\}$' "$bundle/R/etc/ldpaths" ||
    { echo "etc/ldpaths has no default R_LD_LIBRARY_PATH line in the expected form" >&2; exit 1; }
#    Checked here as the user's machine will resolve them: with no path but R's own lib folder (what
#    R's front-end script sets), nothing but the C library may resolve outside the bundle.
find "$bundle/R" "$bundle/syslib" -type f \( -name '*.so' -o -name '*.so.*' -o -path '*/bin/exec/R' \) -print0 |
    LD_LIBRARY_PATH="$bundle/R/lib" xargs -0 -n 40 ldd >"$ldd_out" 2>/dev/null || true
outside="$(awk '/=> \//{print $1 "\t" $3}' "$ldd_out" | sort -u | BUNDLE="$bundle/" GLIBC="$glibc" awk -F'\t' \
    'index($2, ENVIRON["BUNDLE"]) != 1 { n = split($2, parts, "/"); if (parts[n] !~ ENVIRON["GLIBC"]) print }')"
echo "== resolved outside the bundle after patching (must be none)"
if [ -n "$outside" ]; then
    printf '%s\n' "$outside" >&2
    exit 1
fi
echo "== unresolved after patching"
grep 'not found' "$ldd_out" | sort -u || true

# 4. The launcher.
install -m 0755 "$here/Rscript" "$bundle/bin/Rscript"
