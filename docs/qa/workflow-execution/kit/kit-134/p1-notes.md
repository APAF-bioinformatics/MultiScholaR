# KIT-134 P1: notes to the lockfile

Written 2026-11-17 at commit 015d6577 (branch `performance-v2`). Read-only work: nothing was installed, upgraded or removed. Files: `lockfile-source.tsv` (one row per package), `multischolar.lock.json` (the lockfile), this note. Not committed.

## What was read

The per-file runner sets `R_LIBS` to the proof library and nothing else, so R 4.6.1 (2026-06-24, `x86_64-pc-linux-gnu`) searches, in this order:

1. `/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/core/browser-smoke.y5ordS/library` (21 packages, MultiScholaR itself among them)
2. `/home/doktersmol/R/x86_64-pc-linux-gnu-library/4.6` (557)
3. `/usr/lib/R/library` (30: base and recommended)

`utils::installed.packages()` over those paths, first hit per package. The closure starts from DESCRIPTION's Imports, Depends, LinkingTo and Suggests and recurses over Imports, Depends and LinkingTo of the installed packages. Bioconductor is 3.23 (BiocVersion 3.23.1), still the current release on bioconductor.org today.

## Counts

| | Hard closure | Reached only through Suggests | Total |
|---|---|---|---|
| CRAN | 158 | 190 | 348 |
| Bioconductor 3.23 | 2 (limma, qvalue) | 44 | 46 |
| GitHub, with a recorded sha | 1 (RUVIIIC) | 1 (Glimma) | 2 |
| Local, no sha recorded | 0 | 2 (fastduckdb, fastlimpadpc) | 2 |
| Base R 4.6.1 | 10 | 1 | 11 |
| Not installed | 0 | 5 | 5 |
| **Total** | **171** | **243** | **414** |

By library: 20 from the proof library (all Suggests-only), 367 from the user library, 22 from R's own (11 base, 11 recommended). One package is shadowed: limpa 1.4.2 in the proof library hides limpa 1.4.0 in the user library; the lockfile has 1.4.2.

## Five Suggests are not installed, so the suite has never run with them

| Package | Where it would come from | Windows binary |
|---|---|---|
| MOFA2 | Bioconductor 3.23 (1.22.1 now) | yes |
| GSEABase | Bioconductor 3.23 (1.74.0 now) | yes |
| dynamicTreeCut | CRAN (1.63-1 now) | yes |
| RhpcBLASctl | CRAN (0.23-42 now) | yes |
| GlimmaV2 | nowhere under that name (see below) | no |

They are in the lockfile with `Installed: false` and no version. Their own dependencies are not in the closure, because the closure is walked through installed packages; MOFA2 in particular would add its Python side (basilisk or reticulate, and mofapy2). If the bundle is to carry MOFA2, its version and its dependencies have no verified counterpart on the laptop.

`GlimmaV2` in Suggests names a package that does not exist under that name: the repository `APAF-bioinformatics/GlimmaV2` installs a package called `Glimma` (2.15.2, sha `743af094`), which is what the library holds and which takes the place of Bioconductor's Glimma (2.22.1 in 3.23). Recorded, not changed.

## What cannot be reproduced from a repository

**No commit recorded (2).** `fastduckdb` 0.2.0 (built 2026-10-06) and `fastlimpadpc` 0.1.0 (built 2026-10-21) were installed from local source trees: their installed DESCRIPTION files carry no `RemoteSha`, no `Repository` and no `Packaged` field. DESCRIPTION's Remotes pins them at

- `APAF-bioinformatics/fast_duckdb/fastduckdb@d78511504fbdb89fbf34bc26b22b3ff1524b0555`
- `APAF-bioinformatics/fast_limpa_dpc@d4f9cc1dcfeb0183c15bd08872dbe5cf3c2b09a2`

but nothing in the library proves that the installed builds are those commits. The lockfile records the installed version and the pin from DESCRIPTION as two separate fields and does not claim they are the same. Both repositories are internal (the CI token cannot read them, commit 5534aa97), and both need Rust at build time. I did not look for the local clones to compare; that is the check that would close this.

**GitHub by sha (2), reproducible only from GitHub.** `RUVIIIC` 1.0.19 from `cran/RUVIIIC` at `c534151f6b2265b36294ff630f4dd13c2e86a9fa` (hard dependency; no longer in CRAN's current index, so the mirror or CRAN's Archive are the only sources; it compiles C++). `Glimma` 2.15.2 from `APAF-bioinformatics/GlimmaV2` at `743af0944957ce2515b17b38d329506edd06ca5a` (no compiled code). Both were installed from `HEAD`; the sha is what pins them. DESCRIPTION's Remotes names neither sha.

**Installed version is no longer the repository's current version (108 of 394).** 97 from CRAN (40 in the hard closure) and 11 from Bioconductor (limma in the hard closure). 57 of the 108 are compiled (22 in the hard closure: arrow, bit64, curl, data.table, duckdb, igraph, iq, lattice, limma, Matrix, openssl, openxlsx, parallelly, qs2, RcppParallel, readxl, ruv, stringfish, stringi, writexl, xfun, zip).

- CRAN keeps old sources in `src/contrib/Archive`, so the 97 can be rebuilt at the locked version, from source only. Binaries of old versions are not on CRAN; a dated snapshot repository would have them, and choosing one is a P2 decision.
- Bioconductor does not keep superseded versions inside a release. The 11 (biomaRt, Biostrings, clusterProfiler, edgeR, enrichplot, ggtree, limma, S4Arrays, S4Vectors, SparseArray, UniProt.ws) can only be had at the locked version from `git.bioconductor.org` by commit; the lockfile carries each package's `GitLastCommit` (seven characters, as Bioconductor records it). 3.23 is still the live release, so more of the 46 can move before P2 runs.
- Several gaps are large and touch stored numbers or formats: arrow 24.0.0 against 26.0.1, duckdb 1.5.5 against 2.0.0, qs2 0.2.2 against 0.4.0, future 1.70.0 against 1.82.0, curl 7.1.0 against 8.1.0, writexl 1.5.4 against 2.0.1, RcppParallel 5.1.11-2 against 6.2.1, limma 3.68.4 against 3.68.5, edgeR 4.10.1 against 4.10.7. P2's rule ("same version, or the difference named") has to be applied to these one by one.

| Package | Installed (locked) | Repository now | Source | Closure | Compiled | Windows binary now |
| --- | --- | --- | --- | --- | --- | --- |
| aplot | 0.3.1 | 0.3.2 | CRAN | suggests-only | no | 0.3.2 |
| arrow | 24.0.0 | 26.0.1 | CRAN | hard | yes | 26.0.1 |
| beepr | 2.0 | 2.0.1 | CRAN | suggests-only | no | 2.0.1 |
| biomaRt | 2.68.0 | 2.68.6 | Bioconductor | suggests-only | no | 2.68.6 |
| Biostrings | 2.80.1 | 2.80.2 | Bioconductor | suggests-only | yes | 2.80.2 |
| bit64 | 4.8.2 | 4.8.6 | CRAN | hard | yes | 4.8.6 |
| bitops | 1.0-9 | 1.1-0 | CRAN | suggests-only | yes | 1.1-0 |
| bslib | 0.11.0 | 0.12.0 | CRAN | hard | no | 0.12.0 |
| caTools | 1.18.3 | 1.18.4 | CRAN | suggests-only | yes | 1.18.4 |
| cluster | 2.1.8.2 | 2.1.8.3 | CRAN | suggests-only | yes | 2.1.8.3 |
| clusterProfiler | 4.20.0 | 4.20.3 | Bioconductor | suggests-only | no | 4.20.3 |
| colorspace | 2.1-2 | 2.1-4 | CRAN | suggests-only | yes | 2.1-4 |
| configr | 0.3.5 | 0.3.6 | CRAN | hard | no | 0.3.6 |
| covr | 3.6.5 | 4.0.0 | CRAN | suggests-only | yes | 3.6.5 |
| curl | 7.1.0 | 8.1.0 | CRAN | hard | yes | 8.1.0 |
| data.table | 1.18.4 | 1.18.6.1 | CRAN | hard | yes | 1.18.6.1 |
| Deriv | 4.2.0 | 4.3.6 | CRAN | suggests-only | no | 4.3.6 |
| diffobj | 0.3.6 | 0.3.10 | CRAN | suggests-only | yes | 0.3.10 |
| duckdb | 1.5.5 | 2.0.0 | CRAN | hard | yes | 2.0.0 |
| duckplyr | 1.2.1 | 1.2.2 | CRAN | suggests-only | no | 1.2.2 |
| edgeR | 4.10.1 | 4.10.7 | Bioconductor | suggests-only | yes | 4.10.7 |
| enrichit | 0.2.0 | 0.2.7 | CRAN | suggests-only | yes | 0.2.7 |
| enrichplot | 1.32.0 | 1.32.1 | Bioconductor | suggests-only | no | 1.31.2 |
| flextable | 0.10.0 | 0.10.4 | CRAN | suggests-only | no | 0.10.4 |
| Formula | 1.2-5 | 1.2-6 | CRAN | suggests-only | no | 1.2-6 |
| future | 1.70.0 | 1.82.0 | CRAN | hard | no | 1.82.0 |
| future.apply | 1.20.2 | 1.21.0 | CRAN | suggests-only | no | 1.21.0 |
| gargle | 1.6.1 | 1.7.0 | CRAN | suggests-only | no | 1.7.0 |
| gdtools | 0.5.1 | 0.6.0 | CRAN | suggests-only | yes | 0.6.0 |
| gert | 2.3.1 | 2.5.0 | CRAN | suggests-only | yes | 2.5.0 |
| GGally | 2.4.0 | 2.5.0 | CRAN | suggests-only | no | 2.5.0 |
| ggfun | 0.2.1 | 0.2.2 | CRAN | suggests-only | no | 0.2.2 |
| ggplot2 | 4.0.3 | 4.0.4 | CRAN | hard | no | 4.0.4 |
| ggplotify | 0.1.3 | 0.1.4 | CRAN | hard | no | 0.1.4 |
| ggsci | 5.1.0 | 5.2.0 | CRAN | suggests-only | no | 5.2.0 |
| ggstats | 0.13.0 | 0.14.0 | CRAN | suggests-only | no | 0.14.0 |
| ggtangle | 0.1.2 | 0.1.5 | CRAN | suggests-only | no | 0.1.5 |
| ggtree | 4.2.0 | 4.2.1 | Bioconductor | suggests-only | no | 4.2.1 |
| gh | 1.6.0 | 1.6.1 | CRAN | hard | no | 1.6.1 |
| globals | 0.19.1 | 0.20.0 | CRAN | hard | no | 0.20.0 |
| googledrive | 2.1.2 | 2.2.0 | CRAN | suggests-only | no | 2.2.0 |
| googlesheets4 | 1.1.2 | 1.1.3 | CRAN | suggests-only | no | 1.1.3 |
| graphlayouts | 1.2.4 | 1.2.5 | CRAN | suggests-only | yes | 1.2.5 |
| gson | 0.2.0 | 0.2.2 | CRAN | suggests-only | no | 0.2.2 |
| gt | 1.3.0 | 1.4.0 | CRAN | suggests-only | no | 1.4.0 |
| httr | 1.4.8 | 1.4.9 | CRAN | hard | no | 1.4.9 |
| httr2 | 1.2.3 | 1.3.0 | CRAN | hard | no | 1.3.0 |
| igraph | 2.3.3 | 2.3.4 | CRAN | hard | yes | 2.3.4 |
| iq | 2.0.1 | 2.3.0 | CRAN | hard | yes | 2.3.0 |
| KernSmooth | 2.23-26 | 2.23-27 | CRAN | suggests-only | yes | 2.23-27 |
| knitr | 1.51 | 1.52 | CRAN | hard | no | 1.52 |
| lambda.r | 1.2.4 | 1.2.6 | CRAN | suggests-only | no | 1.2.6 |
| lattice | 0.22-9 | 0.23-1 | CRAN | hard | yes | 0.23-1 |
| limma | 3.68.4 | 3.68.5 | Bioconductor | hard | yes | 3.68.5 |
| listenv | 1.0.0 | 1.1.0 | CRAN | hard | no | 1.1.0 |
| litedown | 0.9 | 1.0 | CRAN | suggests-only | no | 1.0 |
| lme4 | 2.0-1 | 2.1-0 | CRAN | suggests-only | yes | 2.1-0 |
| logger | 0.4.2 | 0.4.3 | CRAN | hard | no | 0.4.3 |
| markdown | 2.0 | 2.1 | CRAN | suggests-only | no | 2.1 |
| MASS | 7.3-65 | 7.3-66 | CRAN | suggests-only | yes | 7.3-66 |
| Matrix | 1.7-5 | 1.7-6 | CRAN | hard | yes | 1.7-6 |
| mgcv | 1.9-4 | 1.9-5 | CRAN | suggests-only | yes | 1.9-5 |
| nlme | 3.1-169 | 3.1-172 | CRAN | suggests-only | yes | 3.1-172 |
| NLP | 0.3-2 | 0.3-3 | CRAN | suggests-only | no | 0.3-3 |
| nnet | 7.3-20 | 7.3-21 | CRAN | suggests-only | yes | 7.3-21 |
| officer | 0.7.5 | 0.7.8 | CRAN | suggests-only | no | 0.7.8 |
| openssl | 2.4.2 | 2.5.0 | CRAN | hard | yes | 2.5.0 |
| openxlsx | 4.2.8.1 | 4.2.9 | CRAN | hard | yes | 4.2.9 |
| otel | 0.2.0 | 0.3.0 | CRAN | hard | no | 0.3.0 |
| pak | 0.10.0 | 0.11.2 | CRAN | suggests-only | yes | 0.11.2 |
| parallelly | 1.48.0 | 1.50.0 | CRAN | hard | yes | 1.50.0 |
| plotly | 4.12.0 | 4.12.1 | CRAN | hard | no | 4.12.1 |
| qs2 | 0.2.2 | 0.4.0 | CRAN | hard | yes | 0.4.0 |
| RcppArmadillo | 15.4.0-1 | 15.6.2-1 | CRAN | suggests-only | yes | 15.6.2-1 |
| RcppParallel | 5.1.11-2 | 6.2.1 | CRAN | hard | yes | 6.2.1 |
| RCurl | 1.98-1.19 | 1.98-1.20 | CRAN | suggests-only | yes | 1.98-1.20 |
| readxl | 1.5.0 | 1.6.0 | CRAN | hard | yes | 1.6.0 |
| reticulate | 1.46.0 | 1.47.0 | CRAN | suggests-only | yes | 1.47.0 |
| rmarkdown | 2.31 | 2.32 | CRAN | hard | no | 2.32 |
| roxygen2 | 8.0.0 | 8.1.1 | CRAN | suggests-only | yes | 8.1.1 |
| rstatix | 1.0.0 | 1.1.0 | CRAN | suggests-only | no | 1.1.0 |
| ruv | 0.9.7.1 | 0.9.7.2 | CRAN | hard | yes | 0.9.7.2 |
| S4Arrays | 1.12.0 | 1.12.1 | Bioconductor | suggests-only | yes | 1.12.1 |
| S4Vectors | 0.50.1 | 0.50.3 | Bioconductor | suggests-only | yes | 0.50.3 |
| segmented | 2.2-1 | 2.2-2 | CRAN | suggests-only | no | 2.2-2 |
| selectr | 0.6-0 | 0.8-0 | CRAN | hard | no | 0.8-0 |
| shinyFiles | 0.9.3 | 0.9.4 | CRAN | hard | no | 0.9.4 |
| shinytest2 | 0.5.1 | 0.6.0 | CRAN | suggests-only | yes | 0.6.0 |
| slam | 0.1-55 | 0.1-56 | CRAN | suggests-only | yes | 0.1-56 |
| sp | 2.2-1 | 2.2-4 | CRAN | suggests-only | yes | 2.2-4 |
| SparseArray | 1.12.2 | 1.12.3 | Bioconductor | suggests-only | yes | 1.12.3 |
| stringfish | 0.19.0 | 0.19.2 | CRAN | hard | yes | 0.19.2 |
| stringi | 1.8.7 | 1.8.9 | CRAN | hard | yes | 1.8.9 |
| survival | 3.8-6 | 3.8-12 | CRAN | suggests-only | yes | 3.8-12 |
| systemfonts | 1.3.2 | 1.3.3 | CRAN | suggests-only | yes | 1.3.2 |
| tidydr | 0.0.6 | 0.0.7 | CRAN | suggests-only | no | 0.0.7 |
| tidytree | 0.4.8 | 0.4.9 | CRAN | suggests-only | no | 0.4.9 |
| tinytex | 0.60 | 0.61 | CRAN | hard | no | 0.61 |
| tm | 0.7-18 | 0.7-20 | CRAN | suggests-only | yes | 0.7-20 |
| UniProt.ws | 2.52.0 | 2.52.1 | Bioconductor | suggests-only | no | 2.52.1 |
| usethis | 3.2.1 | 3.2.2 | CRAN | suggests-only | no | 3.2.2 |
| V8 | 8.2.0 | 8.2.1 | CRAN | suggests-only | yes | 8.2.1 |
| writexl | 1.5.4 | 2.0.1 | CRAN | hard | yes | 2.0.1 |
| xfun | 0.60 | 0.62 | CRAN | hard | yes | 0.62 |
| XML | 3.99-0.23 | 3.99-0.25 | CRAN | suggests-only | yes | 3.99-0.25 |
| yulab.utils | 0.2.4 | 0.2.5 | CRAN | hard | no | 0.2.5 |
| zip | 3.0.0 | 3.0.2 | CRAN | hard | yes | 3.0.2 |
| zoo | 1.8-15 | 1.9-1 | CRAN | suggests-only | yes | 1.9-1 |

**Annotation data (1).** `GO.db` 3.23.1 is a Bioconductor annotation package (source only, no compilation, current).

## Compiled packages and what the runners must provide

176 of the 398 installed non-base packages contain compiled code (78 in the hard closure, 98 Suggests-only): 153 CRAN, 20 Bioconductor, RUVIIIC, fastduckdb and fastlimpadpc. `lockfile-source.tsv` marks each (`needs_compilation`, `has_compiled_libs`). The two internal packages have no `NeedsCompilation` field; they are counted from their `libs/` directory. Every compiled package was built on the laptop (`x86_64-pc-linux-gnu`, R 4.6.1), none taken as a repository binary; the ticket names the laptop's `-march=native` as the B-01 concern (the flag itself was not checked here).

The 46 compiled packages that declare SystemRequirements:

| Package | Version | Source | Closure | SystemRequirements |
| --- | --- | --- | --- | --- |
| arrow | 24.0.0 | CRAN | hard | C++20; for AWS S3 support on Linux, libcurl and openssl (optional); cmake >= 3.26 (build-time only, and only for full source build) |
| BiocParallel | 1.46.0 | Bioconductor | suggests-only | C++11 |
| curl | 7.1.0 | CRAN | hard | libcurl (>= 7.73): libcurl-devel (rpm) or libcurl4-openssl-dev (deb) |
| duckdb | 1.5.5 | CRAN | hard | xz (for building from source) |
| fastduckdb | 0.2.0 | local | suggests-only | Cargo (Rust's package manager), rustc >= 1.88.0, xz |
| fastlimpadpc | 0.1.0 | local | suggests-only | Cargo (Rust's package manager), rustc >= 1.65.0, xz |
| fs | 2.1.0 | CRAN | hard | libuv: libuv-devel (rpm) or libuv1-dev (deb). Alternatively to build the vendored libuv 'cmake' is required. GNU make. |
| gdtools | 0.5.1 | CRAN | suggests-only | cairo, freetype2, fontconfig |
| gert | 2.3.1 | CRAN | suggests-only | libgit2 (>= 1.0): libgit2-devel (rpm) or libgit2-dev (deb) |
| ggiraph | 0.9.6 | CRAN | suggests-only | libpng |
| git2r | 0.36.2 | CRAN | suggests-only | libgit2 (>= 1.0) and its dependencies libssh2 and openssl header and libraries. |
| haven | 2.5.5 | CRAN | suggests-only | GNU make, zlib: zlib1g-dev (deb), zlib-devel (rpm) |
| httpuv | 1.6.17 | CRAN | hard | GNU make, zlib |
| igraph | 2.3.3 | CRAN | hard | libxml2 (optional), glpk (>= 4.57, optional) |
| jpeg | 0.1-11 | CRAN | suggests-only | libjpeg |
| locfit | 1.5-9.12 | CRAN | suggests-only | USE_C17 |
| lubridate | 1.9.5 | CRAN | hard | A system with zoneinfo data (e.g. /usr/share/zoneinfo). On Windows the zoneinfo included with R is used. |
| magick | 2.9.1 | CRAN | hard | ImageMagick++: ImageMagick-c++-devel (rpm) or libmagick++-dev (deb) |
| minqa | 1.2.8 | CRAN | suggests-only | GNU make |
| nanoarrow | 0.9.0 | CRAN | suggests-only | libzstd (optional), liblz4 (optional) |
| nloptr | 2.2.1 | CRAN | suggests-only | cmake (>= 3.2.0) which is used only on Linux or macOS systems when no system build of nlopt (>= 2.7.0) can be found. |
| openssl | 2.4.2 | CRAN | hard | OpenSSL >= 1.0.2 |
| png | 0.1-9 | CRAN | hard | libpng |
| qs2 | 0.2.2 | CRAN | hard | GNU make, C++17 |
| ragg | 1.5.2 | CRAN | suggests-only | freetype2, libpng, libtiff, libjpeg, libwebp, libwebpmux |
| RcppParallel | 5.1.11-2 | CRAN | hard | GNU make, Intel TBB, Windows: cmd.exe and cscript.exe, Solaris: g++ is required |
| RcppTOML | 0.2.3 | CRAN | hard | A C++17 compiler |
| RCurl | 1.98-1.19 | CRAN | suggests-only | GNU make, libcurl |
| reticulate | 1.46.0 | CRAN | suggests-only | Python (>= 2.7.0) |
| rgl | 1.3.36 | CRAN | suggests-only | OpenGL and GLU Library (Required for display in R. See "Installing OpenGL support" in README.md. Not needed if only browser displays using rglwidget() are wanted.), zlib (optional), libpng (>=1.2.9, optional), FreeType (optional), pandoc (>=1.14, needed for vignettes), LaTeX needed for optional latex3d() examples. |
| Rhtslib | 3.8.0 | Bioconductor | suggests-only | libbz2 & liblzma & libcurl (with header files), GNU make |
| Rsamtools | 2.28.0 | Bioconductor | suggests-only | GNU make |
| RUVIIIC | 1.0.19 | GitHub | hard | C++11 |
| sass | 0.4.10 | CRAN | hard | GNU make |
| seqinr | 4.2-44 | CRAN | suggests-only | zlib headers and library. |
| stringfish | 0.19.0 | CRAN | hard | GNU make, C++17 |
| stringi | 1.8.7 | CRAN | hard | ICU4C (>= 61, optional) |
| svglite | 2.2.2 | CRAN | suggests-only | libpng |
| systemfonts | 1.3.2 | CRAN | suggests-only | fontconfig, freetype2 |
| textshaping | 1.0.5 | CRAN | suggests-only | freetype2, harfbuzz, fribidi |
| timechange | 0.4.0 | CRAN | hard | A system with zoneinfo data (e.g. /usr/share/zoneinfo). On Windows the zoneinfo included with R is used. |
| V8 | 8.2.0 | CRAN | suggests-only | On Linux you can build against libv8-dev (Debian) or v8-devel (Fedora). We also provide static libv8 binaries for most platforms, see the README for details. |
| websocket | 1.4.4 | CRAN | suggests-only | GNU make, OpenSSL >= 1.0.2 |
| writexl | 1.5.4 | CRAN | hard | zlib |
| XML | 3.99-0.23 | CRAN | suggests-only | libxml2 (>= 2.6.3) |
| xml2 | 1.6.0 | CRAN | hard | libxml2: libxml2-dev (deb), libxml2-devel (rpm) |

Read as a list for a Linux runner building from source: a C++17/C++20 toolchain and GNU make, cmake, Cargo and rustc (>= 1.88 for fastduckdb, >= 1.65 for fastlimpadpc) and xz, libcurl, OpenSSL, libxml2, zlib, libbz2 and liblzma (Rhtslib), libuv, libgit2, ImageMagick++ (magick, a hard dependency), libpng, libjpeg, libtiff, libwebp, freetype2, fontconfig, harfbuzz, fribidi, cairo, ICU, glpk (igraph, optional), TBB (RcppParallel), V8 or its static download, OpenGL/GLU (rgl, optional), zoneinfo. On Windows, CRAN and Bioconductor binaries carry these libraries statically; only what is built from source needs Rtools 4.6 and, for the two internal packages, Rust with the GNU target.

Packages without compiled code whose SystemRequirements matter at run time (the bundle must carry or do without them):

| Package | Closure | SystemRequirements |
| --- | --- | --- |
| chromote | suggests-only | Google Chrome or other Chromium-based browser. chromium: chromium (rpm) or chromium-browser (deb) |
| clipr | hard | xclip (https://github.com/astrand/xclip) or xsel (http://www.vergenet.net/~conrad/software/xsel/) for accessing the X11 clipboard, or wl-clipboard (https://github.com/bugaevc/wl-clipboard) for systems using Wayland. |
| credentials | suggests-only | git (optional) |
| GetoptLong | suggests-only | Perl, Getopt::Long |
| GGally | suggests-only | openssl |
| gitcreds | hard | git |
| knitr | hard | Package vignettes based on R Markdown v2 or reStructuredText require Pandoc (http://pandoc.org). The function rst2pdf() requires rst2pdf (https://github.com/rst2pdf/rst2pdf). |
| pkgdown | suggests-only | pandoc (>= 2.10.1) |
| remotes | hard | Subversion for install_svn, git for install_git |
| reprex | suggests-only | pandoc (>= 2.0) - https://pandoc.org/ |
| rmarkdown | hard | pandoc (>= 1.14) - http://pandoc.org |

Pandoc (rmarkdown, knitr) is already in P3. Chrome or Chromium (chromote, for webshot2 and shinytest2) is a test and export dependency, not bundled by any part of the ticket today.

## Windows

`available.packages()` with the `win.binary` contrib URL works from Linux under R 4.6.1 for CRAN (`bin/windows/contrib/4.6`) and for Bioconductor software (`packages/3.23/bioc/bin/windows/contrib/4.6`). Bioconductor's annotation, experiment and workflow repositories return 404 for a Windows index: they are source-only by design.

Of the 394 CRAN and Bioconductor packages in the lockfile:

- 285 have a Windows binary at exactly the locked version.
- 108 have a Windows binary at another version: 105 newer than the lock (the packages of the table above, less enrichplot, and less covr and systemfonts, whose Windows binaries are still at the locked version) and 3 older than the lock because the Windows build lags the source (DOSE 4.5.0, enrichplot 1.31.2, KEGGREST 1.51.1; none compiled). 55 of the 108 are compiled (22 in the hard closure), so matching the lock on Windows means building them from source under Rtools, or taking a snapshot repository's binaries, or accepting the newer binary and naming the difference.
- 1 has no Windows binary: GO.db (annotation, installs from source without a compiler).

No Windows binary exists anywhere for the four non-repository packages: RUVIIIC (C++, hard dependency, must build under Rtools), fastduckdb and fastlimpadpc (Rust; whether they build and hold their results under Rtools is the open question of KIT-130 revision 1.1), Glimma from GlimmaV2 (no compiled code). The query says nothing about whether a binary works, only that it is listed.

## Python

No `mofapy2` is installed: not for the system Python (3.14.9, `/usr/bin/python3`), not in a virtualenv, and there is no reticulate or basilisk environment directory. reticulate 1.46.0 is installed; MOFA2 is not. The lockfile's Python section is a placeholder that says mofapy2 is pinned in P2. There is no verified Python state on the laptop to export, which also means the MOFA route has not run against the verified library.

## For P1's proof (revision 1.2)

The proof is `installed.packages()` against the lockfile as a table, then P5 on a runner. The lockfile was generated from that same call, so the table is trivially equal today; it becomes evidence only if it is re-run when P2 starts, since the user library (path 2) is a live library that anything can upgrade in the meantime. The lockfile is not committed yet.
