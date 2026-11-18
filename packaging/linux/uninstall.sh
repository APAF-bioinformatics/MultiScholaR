#!/bin/sh
# Remove this user's MultiScholaR installation (KIT-134 P4): the application of this version, the
# menu entry and the icon. The logs (~/.config/MultiScholaR/logs) and the run folders
# (~/.local/state/multischolar) are left; delete them yourself if you no longer need them.
set -eu

here="$(cd "$(dirname "$(readlink -f "$0")")" && pwd -P)"
data="${XDG_DATA_HOME:-$HOME/.local/share}"
case "$here" in
    "$data"/multischolar/*) ;;
    *) echo "uninstall.sh: run the copy inside the installation ($data/multischolar/<version>)" >&2; exit 1 ;;
esac

rm -rf "$here"
rmdir "$data/multischolar" 2>/dev/null || true
if [ ! -d "$data/multischolar" ] || [ -z "$(ls -A "$data/multischolar" 2>/dev/null)" ]; then
    rm -f "$data/applications/multischolar.desktop" "$data/icons/hicolor/512x512/apps/multischolar.png"
    if command -v update-desktop-database >/dev/null 2>&1; then
        update-desktop-database "$data/applications" >/dev/null 2>&1 || true
    fi
fi
echo "MultiScholaR removed from $here"
