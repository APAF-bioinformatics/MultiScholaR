#!/usr/bin/env bash
set -euo pipefail
evidence=/home/doktersmol/Projects/APAF/MultiScholaR/docs/qa/workflow-execution/perf-108/native-null-repair
[[ ! -e "$evidence/published-release.json" ]]
gh api --method GET repos/APAF-bioinformatics/fast_duckdb/releases/tags/v0.2.0 > "$evidence/pre-update-release.json"
python3 - "$evidence" <<'PY'
from pathlib import Path
import hashlib,json,sys
p=Path(sys.argv[1]);r=json.loads((p/'remediation.json').read_text())
assert r['status']=='closed' and r['installed_assertions']==408
for key in ['source_freeze','patch','fixed_build','isolated_checks','installed_tests','rollback','compatibility']:
    x=r[key];assert hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['sha256'],key
old=json.loads((p.parent/'packaging-repair/published-repair-release.json').read_text())
current=json.loads((p/'pre-update-release.json').read_text())
assert old['body'].strip()==current['body'].strip()
assert all(x['name']!='native-provider-repair-validation.json' for x in current['assets'])
PY
gh release upload v0.2.0 --repo APAF-bioinformatics/fast_duckdb "$evidence/native-provider-repair-validation.json"
gh release edit v0.2.0 --repo APAF-bioinformatics/fast_duckdb --notes-file "$evidence/release-notes-updated.md"
gh api --method GET repos/APAF-bioinformatics/fast_duckdb/releases/tags/v0.2.0 > "$evidence/published-release.json"
