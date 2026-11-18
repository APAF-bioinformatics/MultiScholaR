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

## Second pass (2026-11-18): wording, three layout changes, the figures now measured

`after/` was produced again on the branch's tip after the second pass (commits "wording for bench scientists",
"three layout changes", "figures now measured"); `before/` is unchanged (the tree before the first layout pass). The
images of the first pass are in the history (cd353645). The complete old/new wording is in `wording-review.md`.

- Run from a directory that is not the package's (shinytest2 otherwise looks for a loaded development package):
  `cd <scratch> && R_LIBS=<the test library> Rscript <worktree>/docs/.../screenshots.R <worktree> <output dir> after`.
- New state `7-supplied-rules-open`: the k rules' explanation opened from the keyboard (focus on the summary, Enter;
  the script prints that the summary took the focus and that the disclosure is open). In every other state it is
  closed, as the tab opens.
- `3-supplied-separation-pending` keeps its file name, but nothing is pending any more: the separation rule's figure
  for supplied controls is the second line of that row.

Height of the section in px (first pass -> second pass; the first pass showed the k rules' explanation open, the
second shows it closed, and `7-...` open with supplied controls selected):

| Case | State | 1440 | 1024 |
|---|---|---|---|
| Proteomics classic | 2 recommended | 1289 -> 1493 | 1674 -> 1961 |
| Proteomics classic | 1 earlier versions (its controls shown) | 1780 -> 2234 | 2017 -> 2702 |
| Proteomics classic | 3 supplied controls | 1365 -> 1493 (open: 1966) | 1864 -> 1961 (open: 2676) |
| Proteomics classic | 5 none | 842 -> 1296 | 1079 -> 1764 |
| DIA with limpa | 2 recommended | 1425 -> 1743 | 1823 -> 2299 |
| Transcriptomics | 2 recommended | 1233 -> 1548 | 1625 -> 2025 |
| Lipidomics | 2 recommended | 1269 -> 1571 | 1625 -> 2064 |

The section is taller, not shorter: each design's description is now visible text (it was a hover title), the terms
are explained in the section (what unwanted variation is, how to read the figures), supplied controls and controls
from the data carry a second figure line, and lipid and metab's technical-replicate row carries a figure. The
disclosure takes about 470 px (1440) and 715 px (1024) out of the states that offer a k rule. In the first pass
"none" and "earlier versions" hid the whole k block; the explanatory lines above it are now longer.
