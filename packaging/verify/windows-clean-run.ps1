# KIT-134: run one R script with the bundled R on Windows, in the environment the shell gives R
# (packaging/electron/main.js, rEnvironment): the runner's own R, Rtools and Python off the PATH,
# every variable that could steer R or Python removed, the library paths inside the bundle.
# From the KIT-130 spike (tools/dist-spike/checks/windows-clean-run.ps1).
#
# Usage: windows-clean-run.ps1 -Bundle <bundle directory> -Script <file.R> [-Arguments ...]
# The script receives the bundle directory first, then the arguments.
param(
    [Parameter(Mandatory)][string]$Bundle,
    [Parameter(Mandatory)][string]$Script,
    [string[]]$Arguments = @(),
    # The script takes only the arguments (run-lanes.R, test-tools.R), not the bundle directory first.
    [switch]$NoBundleArgument
)

Get-ChildItem Env: |
    Where-Object { $_.Name -match '^(PATH|RHOME|R_.*|RTOOLS.*|RETICULATE_.*|PYTHON.*|CONDA.*|VIRTUAL_ENV|BASILISK_.*|LD_LIBRARY_PATH)$' } |
    ForEach-Object { Remove-Item "Env:$($_.Name)" }
$env:PATH = "$env:SystemRoot\System32;$env:SystemRoot;$env:SystemRoot\System32\WindowsPowerShell\v1.0"
$env:R_LIBS_USER = Join-Path $Bundle 'R\no-user-library'
$env:R_LIBS_SITE = Join-Path $Bundle 'R\no-site-library'
# R for Windows puts Rtools' folders at the front of its own PATH; this one does not exist.
$env:RTOOLS45_HOME = Join-Path $Bundle 'no-rtools'
$env:RETICULATE_PYTHON = Join-Path $Bundle 'python\python.exe'
$env:PYTHONNOUSERSITE = '1'
$env:RSTUDIO_PANDOC = Join-Path $Bundle 'pandoc'
$env:MSR_DIAMOND = Join-Path $Bundle 'diamond\diamond.exe'
# The saved-entity store must not take WSL's bash or other GNU tools found on a host (F6).
$env:MULTISCHOLAR_DISABLE_TOOLS = 'bash,dd,sha256sum,sync,tar,cp'
$env:MSR_BUNDLE = $Bundle

if ($NoBundleArgument) {
    & (Join-Path $Bundle 'R\bin\Rscript.exe') --vanilla $Script @Arguments
} else {
    & (Join-Path $Bundle 'R\bin\Rscript.exe') --vanilla $Script $Bundle @Arguments
}
exit $LASTEXITCODE
