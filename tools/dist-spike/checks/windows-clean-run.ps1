# KIT-130 packaging spike: run one check script with the bundled R on Windows, with the runner's
# own R and Python off the PATH and every variable that could steer R or Python removed. The
# variables set here are the ones the launcher sets (main.js, rEnvironment).
#
# Usage: windows-clean-run.ps1 -Bundle <bundle directory> -Check <check>.R [-CheckArguments ...]
param(
    [Parameter(Mandatory)][string]$Bundle,
    [Parameter(Mandatory)][string]$Check,
    [string[]]$CheckArguments = @()
)

Get-ChildItem Env: |
    Where-Object { $_.Name -match '^(RHOME|R_.*|RTOOLS.*|RETICULATE_.*|PYTHON.*|CONDA.*|VIRTUAL_ENV|BASILISK_.*)$' } |
    ForEach-Object { Remove-Item "Env:$($_.Name)" }
$env:PATH = "$env:SystemRoot\System32;$env:SystemRoot"
$env:R_LIBS_USER = Join-Path $Bundle 'R\no-user-library'
$env:R_LIBS_SITE = Join-Path $Bundle 'R\no-site-library'
$env:RETICULATE_PYTHON = Join-Path $Bundle 'python\python.exe'
$env:PYTHONNOUSERSITE = '1'
$env:MSR_BUNDLE = $Bundle

$checkScript = Join-Path $PSScriptRoot $Check
& (Join-Path $Bundle 'R\bin\Rscript.exe') --vanilla $checkScript $Bundle @CheckArguments
exit $LASTEXITCODE
