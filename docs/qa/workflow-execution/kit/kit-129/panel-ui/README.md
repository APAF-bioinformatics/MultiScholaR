# KIT-129 panel UI: the RUV design section in a browser, before and after the layout pass

Branch `kit-129-panel-ui`, cut from `kit-129-panel` 5b34c4e9. Presentation only: the section's texts, numbers,
options, ids, input values, default selection, gate and server behaviour are those of 5b34c4e9.

## How the images were produced

`screenshots.R` (this folder), run once on the unedited tree (`before/`) and once on the branch's tip (`after/`):

    R_LIBS=<the test library> Rscript screenshots.R <worktree> <output dir> <before|after>

- Browser: headless Chromium (`~/.cache/ms-playwright/chromium_headless_shell-1223`), driven by `shinytest2` /
  `chromote`. Nothing was installed.
- Page: a harness, not the full app. It is the app's own shell (`app_ui()`: shinydashboard, the app's styles and
  scripts, dark mode as the app opens) with the home tab's content replaced by the real `mod_omic_norm_ui()` of
  one omic, gate on (`options(MultiScholaR.ruv_design_panel = TRUE)`). The server runs the section's own observers
  only (`registerOmicNormRuvDesign`), so the choices of the k rule follow the selected design and the DIA with
  limpa table replaces the classic one as in the app.
- Widths 1440 and 1024 px. Cases: proteomics classic, proteomics DIA with limpa, transcriptomics, lipidomics
  (metabolomics shares lipidomics' table and layout). States: `1-legacy-default` (the default selection, today's
  controls shown), `2-recommended-estimated`, `3-supplied-separation-pending` (the pending k-rule notice),
  `4-firstpass-number` (the number input), `5-none` and `6-recommended-light` (light mode) for classic proteomics.
  `_page` files are the whole page; the others are the section (after) or the options column that held it
  (before).
- The files here are the screenshots reduced to 64 colours to keep the folder small.

## What the harness does not show

- The tab inside the omic's workflow page (stepper and step tabs above it); the content width is the app's.
- Any other server output of the tab: the plots, `ruv_design_summary` and the k diagnostic slots are empty.
- Tooltips: they are `title` attributes and appear on hover only, which a screenshot does not capture.
- A real session's sequence (import to normalisation); the inputs were set through Shiny's input bindings.
