#!/bin/bash
# run-tests: write a blank request file .github/run-tests/requests/<name>.json (every field empty or at its default;
# nothing is copied from another request). Fill it in, commit it, push run-tests: the push starts one run for it.
#
# Usage: .github/run-tests/new-request.sh <name>
#   <name>  short slug naming this request, e.g. kit129-figures-303217ac (a-z, 0-9, '.', '_', '-'; up to 80)
# Run from the root of a run-tests worktree. A name already used on the branch (now or in its history) is refused:
# one file per request, never edit or reuse another request's file.
set -euo pipefail
name="${1:-}"
if ! [[ "$name" =~ ^[a-z0-9][a-z0-9._-]{0,79}$ ]]; then
    echo "usage: $0 <name>   (a short slug: a-z, 0-9, '.', '_', '-'; up to 80 characters)" >&2; exit 1
fi
root="$(git rev-parse --show-toplevel)"
rel=".github/run-tests/requests/$name.json"
f="$root/$rel"
if [ -e "$f" ]; then echo "$rel exists: choose a new name (never edit another request's file)" >&2; exit 1; fi
git -C "$root" fetch --quiet origin run-tests 2>/dev/null || echo "warning: could not fetch origin run-tests; history check is local only" >&2
if [ -n "$(git -C "$root" log --all --oneline -1 -- "$rel")" ]; then
    echo "$rel was used before on the branch's history: choose a new name" >&2; exit 1
fi
mkdir -p "$(dirname "$f")"
cat > "$f" <<'JSON'
{
  "requested_by": "",
  "label": "",
  "package_commit": "",
  "base_commit": "",
  "tools_ref": "",
  "tests": [],
  "shards": 8,
  "file_timeout_minutes": 45,
  "file_timeouts": {},
  "extra_packages": [],
  "env": {},
  "check": false,
  "scripts": [],
  "fetch_refs": [],
  "blas": "reference",
  "os": "ubuntu-22.04"
}
JSON
echo "$rel written (blank). Fill in requested_by, label, package_commit, base_commit (40 digits, or \"none\"), tests; then:"
echo "  git add $rel && git commit -F <message file> && git push origin run-tests"
