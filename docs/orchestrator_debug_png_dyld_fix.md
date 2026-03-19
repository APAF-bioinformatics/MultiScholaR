# Orchestrator Report: Resolving `png` Package Load Error (dyld @rpath)

**Date**: 2026-03-19
**Subject**: Fix for R `png` package failing to load dynamic library (`libpng16.16.dylib`) during `dyn.load`.

## Issue Description
The user reported an error when building and loading the `png` R package from source under macOS:
```
Error: package or namespace load failed for ‘png’ in dyn.load(file, DLLpath = DLLpath, ...):
 unable to load shared object '.../png.so':
  dlopen(.../png.so, 0x0006): Library not loaded: @rpath/libpng16.16.dylib
  Reason: tried: '/Library/Frameworks/R.framework/Resources/lib/libpng16.16.dylib' (no such file)
```
The error was identified as a consequence of having Anaconda in the `PATH` during the package build. The build tools link against Anaconda's `libpng` correctly, but the resulting `.so` object lacks the `@rpath` to Anaconda's `lib` directory. Conda builds macOS libraries with an `@rpath` dependency by design, so `dyld` needs to know *where* this `@rpath` points at loaded time.

## Root Cause Analysis
- `CONDA_PREFIX=/Users/ignatiuspang/anaconda3` was set and Conda paths were listed first in `PATH`.
- R calls `libpng-config --cflags` and `libpng-config --ldflags` which might resolve to Conda's version or macOS/homebrew's version depending on path precedence.
- If it links to Conda's `libpng16.16.dylib` (`-L/Users/ignatiuspang/anaconda3/lib -lpng16`), the dynamic linker `dyld` at runtime attempts to resolve `@rpath/libpng16.16.dylib`.
- Because the standard macOS R compile chain does not append the conda `lib` directory to the linked runpath bindings (`-Wl,-rpath,...`), the runtime dynamic load fails.

## Implemented Fix
To ensure that R can load any Conda-provided shared libraries dynamically at runtime without manually updating `DYLD_LIBRARY_PATH` or `DYLD_FALLBACK_LIBRARY_PATH` (which is frowned upon due to SIP in modern macOS):
- We appended the explicit rpath argument to the compiler's Linker flags (`LDFLAGS`).
- Modifying `~/.R/Makevars` natively instructs R to embed the search path in all built `.so` objects.

**Modified File:** `~/.R/Makevars`
```make
FC=/opt/homebrew/bin/gfortran
F77=/opt/homebrew/bin/gfortran
FLIBS=-L/opt/homebrew/lib
LDFLAGS += -Wl,-rpath,/Users/ignatiuspang/anaconda3/lib
```

## Validation
- Ran `Rscript -e 'install.packages("png", type="source", repos="https://cloud.r-project.org")'`
- Output confirmed the injection of `-Wl,-rpath,/Users/ignatiuspang/anaconda3/lib` into the linking command.
- The `png` package installed correctly and the `library(png)` call passed during the temporary validation stage of R's standard installation protocol.

---

<!-- APAF Bioinformatics | orchestrator_debug_png_dyld_fix.md | Approved | 2026-03-19 -->
