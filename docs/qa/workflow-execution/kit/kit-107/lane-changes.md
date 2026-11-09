# KIT-107: entry-lane changes

**None.** No entry oracle was re-frozen in KIT-107.

## F-16 (the user, 2026-11-07: "f16 pass to kit-107"): the prot_dia record diff (D6)

`prot_dia` was recorded on the merged tree at `2b1e02c8` (`MULTISCHOLAR_GUI_ENTRY_RECORD_DIR`, alone, through `kit-105/lifecycle/lifecycle.R`) and JSON-diffed against `tests/testdata/gui-entry/prot_dia.json`. The result is [lifecycle/prot_dia-record.diff.txt](lifecycle/prot_dia-record.diff.txt): 6 differences, all header fields. The oracle file carries `frozen`, `source_sha`, `load_mode` and `r_version`, and the raw record carries `omic` and `seconds`. **No step value differs.** As D6 predicted, F-16 changed no capture: the analysis record and `study_parameters.txt` are captured only as `<present>`, and provenance adds no store entity.

W1 (b)'s attributed re-freezes (F107-01, F107-03, F107-05) are unit captures, not entry oracles: [w3/w1-refreeze.md](w3/w1-refreeze.md).
