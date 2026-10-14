# Installed-runtime run (KIT-034, quarantine Q-04 and Q-14)

The tests that need the installed package or a browser carry skip conditions (`tests/testthat/helper-installed-runtime.R`: `skipIfNoInstalledRuntime()`, `skipIfNoBrowserRuntime()`) and were run once against the installed candidate library, `core/browser-smoke.70tjXp` (kit_024, `28075aa1`), on the measurement host, 2026-10-14. Command shape: `MULTISCHOLAR_INSTALLED_RUNTIME=1 R_LIBS=<library> Rscript -e "library(MultiScholaR); testthat::test_file(<file>, package = 'MultiScholaR')"`, with `CHROMOTE_CHROME=/usr/bin/brave NOT_CRAN=true` for the browser files (the harness's own browser choice; chromote's discovery does not find Brave) and `MULTISCHOLAR_DIA008_LIBRARY=<library>` for the passive diagnostic wrapper.

| File | Result against the installed library | Target test |
|---|---|---|
| test-perf-003.R | 173 pass, 0 fail, 0 error, 0 skip | "registered dispatch freezes inputs ..." ran and passed (its job fixture gained the KIT-008 `monitor_borrowed` flag) |
| test-workflow-native-session.R | 541 pass, 0 fail, 0 error, 1 skip (unrelated) | "native publication recovers after process loss" ran and passed |
| test-workflow-saved-recovery.R | 108 pass, 0 fail, 0 error, 0 skip | "SIGKILL at storage boundaries", "concurrent writers" ran and passed |
| test-workflow-ruv-pearson-diagnostic.R | 326 pass, 0 fail, 0 error, 4 skips (unrelated: qualified builds, retained plots) | "design assignment flushes buffered values" ran and passed with the browser set |
| test-workflow-gui-passive-diagnostic.R | 15 pass, 0 fail, 1 error, 0 skip | "the actual diagnostic browser wrapper starts in passive mode" ran and passed (its library requirement is now a skip condition); the one error is "recorded release comparison preserves exact corresponding output bytes", a pinned-release comparison that only runs with the library variable set and expects a different build: outside this ticket, unchanged (it skips under `load_all`) |

Under `devtools::load_all` these tests skip with the reason naming this record.
