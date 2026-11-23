#!/usr/bin/env bash
# KIT-157 (branch run-k157-rc only): build fastlimpadpc at the PREVIOUS pin (fast_limpa_dpc 87a90b89, -C target-cpu=native)
# into $RUNNER_TEMP/k157-old, beside the bundle's library (which holds this checkout's pin), so the two builds are compared
# in the same job on the same CPU. Runs inside the workflow's library install step (the only step that holds the read
# token), on the runner or in the build container like packaging/linux/build.sh. The token is read once by old-limpa.R
# and removed from its environment before R CMD INSTALL; the build output stays under $RUNNER_TEMP and is not printed.
set -euo pipefail
if [ -d "$HOME/.cargo/bin" ]; then PATH="$HOME/.cargo/bin:$PATH"; fi
"/opt/R/${R_VERSION}/bin/Rscript" --vanilla packaging/verify/k157/old-limpa.R "$RUNNER_TEMP/k157-old"
