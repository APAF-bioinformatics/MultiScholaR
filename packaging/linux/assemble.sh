#!/usr/bin/env bash
# KIT-134 P2, build time only: turn an R installed at its build prefix, with the bundle's library
# installed into it, into a folder that runs from anywhere. From the KIT-130 spike
# (tools/dist-spike/build/linux-assemble.sh on branch kit-130-dist-spike).
#
# Usage: assemble.sh <installed R home, e.g. /opt/R/4.6.1/lib/R> <bundle directory>
# The bundle directory may already hold python/ and pandoc/. The result is
#   <bundle>/R       the R home, its front-end script patched to locate itself
#   <bundle>/syslib  every shared library R and the packages need, except the C library
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
awk '/=> \//{print $3}' "$ldd_out" | sort -u | while read -r lib; do
    case "$lib" in "$bundle"/*) continue ;; esac
    if basename "$lib" | grep -Eq "$glibc"; then continue; fi
    cp -L "$lib" "$bundle/syslib/"
done
echo "== bundled shared libraries: $(find "$bundle/syslib" -type f | wc -l)"
ls "$bundle/syslib"

# 3. The launcher.
install -m 0755 "$here/Rscript" "$bundle/bin/Rscript"
