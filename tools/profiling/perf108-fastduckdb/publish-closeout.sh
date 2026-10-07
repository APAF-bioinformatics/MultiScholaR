#!/usr/bin/env bash
set -euo pipefail
repo=/tmp/perf108-fast-duckdb
evidence=/home/doktersmol/Projects/APAF/MultiScholaR/docs/qa/workflow-execution/perf-108/closeout
archive=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/perf-108/closeout/browser-smoke.sopzDb/fastduckdb_0.2.0.tar.gz
[[ -z "$(git -C "$repo" status --porcelain)" ]]
python3 - "$repo" "$evidence" "$archive" <<'PY'
import hashlib,json,subprocess,sys
from pathlib import Path
repo,evidence,archive=map(Path,sys.argv[1:])
q=json.loads((evidence/'release-qualification.json').read_text())
exception=json.loads((evidence/'timing-exception.json').read_text())
assert exception['status']=='authorized'
assert hashlib.sha256((evidence/'timing-exception.json').read_bytes()).hexdigest()==q['exception_sha256']
assert hashlib.sha256(archive.read_bytes()).hexdigest()==q['source_archive']['sha256']
assert subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()==q['source_commit']
for name,expected in q['source_files'].items():
    actual=subprocess.check_output(['git','-C',str(repo),'show','HEAD:fastduckdb/'+name])
    assert hashlib.sha256(actual).hexdigest()==expected,name
for name,expected in q['tested_artifacts'].items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==expected,name
PY
git -C "$repo" tag -a v0.2.0 -m 'fastduckdb 0.2.0: qualified native query and bounded I/O interfaces'
git -C "$repo" push --atomic origin feat/scoped-native-query refs/tags/v0.2.0
gh release create v0.2.0 --repo APAF-bioinformatics/fast_duckdb --verify-tag \
    --title 'fastduckdb 0.2.0: native queries and bounded data transfer' \
    --notes-file "$evidence/release-notes.md" "$archive" "$evidence/release-qualification.json"
gh api --method GET repos/APAF-bioinformatics/fast_duckdb/releases/tags/v0.2.0 > "$evidence/published-release.json"
