# KIT-259 G0: triage of KIT-247 G0-4's "Problem while converting geom to grob"

**Verdict: an artefact of the test launch, not an application defect a user meets.** The application's process
started by shinytest2 (`AppDriver$new` with a function app, as `tests/testthat/helper-e2e-browser.R:66-148` does)
runs at `options(warn = 2)`. There, the warning ggplot2 raises for the label-free fixture's Pearson histogram
("Removed 1 row containing missing values or values outside the scale range (`geom_bar()`)") becomes an error, which
ggplot2 wraps as "Problem while converting geom to grob". Under the default `warn = 0` (testServer, a plain R
process, the packaged entry's `Rscript --vanilla`) the same figure is drawn with that warning. A real, smaller figure
defect underlies it and is a user's: the histogram's fixed y limits drop every 0.001 bin holding more than 4 pairs
(here one bar of 5 of the 20 pairs), silently. Recorded as a candidate (plan section 13), not a ticket.

Run 45610170354 (one scripts job, ubuntu-22.04, R 4.6.1, ggplot2 4.0.3, ragg 1.5.2, rlang 1.3.0, cli 3.6.6, shiny
1.14.0; shinytest2 0.5.1 by the lockfile), at branch `kit-259-g0` 9ea11ba2 = v0.9-staging 63143aa3 plus
`g0-diag.R` and `k259-probe.R` (no package code). The drawing code is byte-identical to 16112917, the lane's
binding. Predictions were written in `g0-diag.R` before the run. Tables: `g0-45610170354/`.

## Observed

| Part | What | Observed |
|---|---|---|
| 1 | testServer lanes `prot_lfq` and `prot_dia` to `norm`, spies on each figure save | every save ok (16/16 each: the pre-normalisation figures are drawn three times, at tab entry, at the run's `redraw_pre` and again; post 4); both heads `normalized`; label-free completion notice shown; `warn` inside each save = 1 |
| 1 | data of the correlation figures | label-free: 20 pairs, pre 0.96345 to 0.97358, post 0.96342 to 0.97373, fullest 0.001 bin 5 (pre and post), 1 bar censored (ymax NA); DIA-NN: 12 pairs, fullest bin 3 (pre) and 4 (post), none censored; no non-finite value, none above 1 |
| 2 | plain process, warn 0: PCA, RLE, density, correlation (the app's producer) and `plotPearson`'s native plot, by ggsave (its png device: `png_dev(..., res = dpi, units = "in")`, ragg) and `png(type = "cairo")` | 20/20 renders ok on each fixture; label-free correlation and native plot: warning "Removed 1 row containing missing values or values outside the scale range (`geom_bar()`)" |
| 2 | the same correlation figures under `warn = 2` | label-free pre and post: error "Problem while converting geom to grob. Error occurred in the 1st layer. Caused by error: (converted from warning) Removed 1 row ... (`geom_bar()`)", call `ggplot2::geom_histogram(breaks = breaks)`; DIA-NN: ok |
| 3 | probe inside shinytest2's application process, after `pkgload::load_all` and before `run_app(test_mode = TRUE)` (app_top) | `getOption("warn")` = **2** (the driver's own process: 0); no global calling handlers; no DISPLAY, cairo TRUE; label-free pre and post correlation: the same chain as the browser's (outer "Problem while converting geom to grob", parent "(converted from warning) Removed 1 row ..."); DIA-NN pre and post: ok |
| 3 | probe inside an observer under `withProgress` (app_observer) | not observed: `AppDriver$new` failed, "Chrome debugging port not open after 10 seconds", after the app_top probe had run |

## Observed against predicted

20 of 26 held. Not held (2): `p1.*.saved_ok` predicted 8/8 and observed 16/16 (my count missed the redraws; every
save was ok, the substance held). Not observed (4): the three app_observer rows and `p3.app_observer.env.warn` (Chrome
did not open its port). The warn = 2 prediction held in the app's process (app_top); the app_top failure prediction
held. Not re-run under VAL-54: the app_top probe already ran in the launch chain that failed and answered the
launch-mode question; the observer probe would add a reactive context only.

## Inferred

- The parent error is the ggplot2 removed-rows warning, raised in `GeomBar`'s missing-value handling of the 1st layer
  (`geom_histogram`) and promoted to an error by `warn = 2`; the bar is censored because `scale_y_continuous(limits =
  c(0, 4))` (`R/utils_workflow_dia_products.R:84`; native `R/func_prot_s4_qc_correlation_methods.R:68`) puts the count
  5 out of range. Not the device, a font, the display or the ggplot2 version.
- Nothing in the package sets `warn` (no `options(warn` in `R/`, `inst/`, `packaging/`, `tools/` at 63143aa3); the
  driver's process loaded the same package with `warn = 0`. Which layer of the launch chain (shinytest2 0.5.1, callr,
  the runner's profile) sets 2 was not established.
- A user is not affected by the error unless they run R with `warn >= 2` themselves. They are affected by the
  censored bar: any proteomics within-group correlation figure with more than 4 pairs in one 0.001 bin loses that bar
  without a word (legacy: `main:R/func_prot_s4_objects.R:2661`; lipid and metabolomics have the limits commented
  out, `R/func_lipid_s4_correlation_methods.R:287`, `R/func_metab_s4_qc_methods.R:768`).

## Not checked

- The probe inside a Shiny observer (Chrome port failure); which launcher layer sets `warn = 2`.
- The worker-drawn figures (`computeWorkflowDiaQcFiguresAction`) on a managed application: by reading they use the
  same producer and ggsave at the worker's default `warn`.
- Real data; TMT; the pass-through route of DIA with limpa.
