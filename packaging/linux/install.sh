#!/bin/sh
# Install MultiScholaR for the current user (KIT-134 P4, D2): no administrator rights, no FUSE,
# nothing downloaded. Run it from the unpacked archive:
#     ./install.sh
# It copies the application to ~/.local/share/multischolar/<version>, adds a menu entry and an
# icon, and puts uninstall.sh beside the application. Needs glibc 2.35 or later (Ubuntu 22.04,
# Debian 12, Fedora 36 and later; D1) and the program /usr/bin/which (M25).
set -eu

here="$(cd "$(dirname "$0")" && pwd -P)"
version="$(sed -n 's/^ *"release": *"\([^"]*\)".*$/\1/p' "$here/VERSION.json" | head -n 1)"
[ -n "$version" ] || { echo "install.sh: cannot read the version from $here/VERSION.json" >&2; exit 1; }
if [ "$(id -u)" -eq 0 ]; then
    echo "install.sh: run this as yourself, not as root: MultiScholaR installs for one user" >&2
    exit 1
fi

glibc="$(getconf GNU_LIBC_VERSION 2>/dev/null | awk '{print $2}')"
case "$glibc" in
    2.3[5-9]* | 2.[4-9][0-9]* | [3-9].*) ;;
    *)
        echo "install.sh: this system's C library is ${glibc:-unknown}; MultiScholaR needs glibc 2.35 or later" >&2
        exit 1
        ;;
esac
# Decision M25: /usr/bin/which is needed by the bundled R (the launcher refuses to start without it).
if [ ! -x /usr/bin/which ]; then
    echo "install.sh: MultiScholaR needs the program /usr/bin/which, which this system does not have." >&2
    echo "  Install it first (Fedora, RHEL: sudo dnf install which; Arch Linux: sudo pacman -S which;" >&2
    echo "  Debian, Ubuntu: sudo apt install debianutils; openSUSE: sudo zypper install which)." >&2
    exit 1
fi

data="${XDG_DATA_HOME:-$HOME/.local/share}"
dest="$data/multischolar/$version"
echo "Installing MultiScholaR $version into $dest"
rm -rf "$dest.partial"
mkdir -p "$dest.partial"
cp -a "$here/app" "$dest.partial/app"
cp "$here/VERSION.json" "$here/THIRD-PARTY-NOTICES.txt" "$here/multischolar.png" "$dest.partial/"
install -m 0755 "$here/uninstall.sh" "$dest.partial/uninstall.sh"

# The launcher the menu entry starts. Electron's sandbox needs unprivileged user namespaces, which
# Ubuntu 23.10 and later restrict through AppArmor for programs installed per user; only there is
# the sandbox turned off (the window shows only the application's own local server, with context
# isolation and no Node in the page).
cat >"$dest.partial/multischolar" <<'LAUNCHER'
#!/bin/sh
here="$(cd "$(dirname "$(readlink -f "$0")")" && pwd -P)"
if [ "$(cat /proc/sys/kernel/apparmor_restrict_unprivileged_userns 2>/dev/null)" = "1" ]; then
    exec "$here/app/multischolar" --no-sandbox "$@"
fi
exec "$here/app/multischolar" "$@"
LAUNCHER
chmod 0755 "$dest.partial/multischolar"

rm -rf "$dest"
mv "$dest.partial" "$dest"

mkdir -p "$data/applications" "$data/icons/hicolor/512x512/apps"
cp "$here/multischolar.png" "$data/icons/hicolor/512x512/apps/multischolar.png"
sed -e "s|@EXEC@|$dest/multischolar|g" -e "s|@VERSION@|$version|g" -e "s|@ICON@|$data/icons/hicolor/512x512/apps/multischolar.png|g" \
    "$here/multischolar.desktop" >"$data/applications/multischolar.desktop"
chmod 0644 "$data/applications/multischolar.desktop"
if command -v update-desktop-database >/dev/null 2>&1; then
    update-desktop-database "$data/applications" >/dev/null 2>&1 || true
fi

echo "MultiScholaR $version is installed. Start it from the applications menu, or run:"
echo "    $dest/multischolar"
echo "To remove it: $dest/uninstall.sh"
