#!/usr/bin/env bash
set -euo pipefail
repo=/tmp/perf108-fast-duckdb
evidence=/home/doktersmol/Projects/APAF/MultiScholaR/docs/qa/workflow-execution/perf-108/packaging-repair
[[ -z "$(git -C "$repo" status --porcelain)" && ! -e "$evidence/published-repair-release.json" ]]
python3 - "$evidence" <<'PY'
import hashlib,json,sys
from pathlib import Path
p=Path(sys.argv[1]);x=json.loads((p/'remediation.json').read_text())
assert x['status']=='verification_passed_publication_pending'
assert [r['success'] for r in x['runs']]==[True,False,True]
for key in ['installer','guide','clean_install','corrected_check','generic_tests']:
    v=x[key];assert hashlib.sha256(Path(v['path']).read_bytes()).hexdigest()==v['sha256'],key
PY
gh api --method GET repos/APAF-bioinformatics/fast_duckdb/releases/tags/v0.2.0 > "$evidence/pre-repair-release.json"
python3 - "$evidence" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);old=json.loads((p.parent/'closeout/published-release.json').read_text());now=json.loads((p/'pre-repair-release.json').read_text())
assert now['body'].strip()==old['body'].strip(), 'release notes changed independently'
assets={x['name']:x['digest'] for x in now['assets']}
for x in old['assets']: assert assets[x['name']]==x['digest']
for name in ['install-fastduckdb-0.2.0.sh','install-fastduckdb-0.2.0.md','packaging-validation.json']:assert name not in assets
PY
git -C "$repo" push origin feat/scoped-native-query
gh release upload v0.2.0 --repo APAF-bioinformatics/fast_duckdb \
    "$repo/scripts/install-fastduckdb-0.2.0.sh" "$repo/docs/install-fastduckdb-0.2.0.md" \
    "$evidence/packaging-validation.json"
gh release edit v0.2.0 --repo APAF-bioinformatics/fast_duckdb --notes-file "$evidence/release-notes-updated.md"
gh api --method GET repos/APAF-bioinformatics/fast_duckdb/releases/tags/v0.2.0 > "$evidence/published-repair-release.json"
