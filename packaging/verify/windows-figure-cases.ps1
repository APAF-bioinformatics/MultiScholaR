# KIT-134 / KIT-052 (Windows figure paths): run figure-paths.R once per case with the installed
# application's R in the shell's environment (windows-clean-run.ps1). Build-time verification.
#
# Usage: windows-figure-cases.ps1 -Runtime <installed runtime> -Checkout <checkout> -TestTools <library>
#            -Out <evidence folder> -Cases <name>=<TMPDIR length>:<item>[,<item>...] ...
# A TMPDIR length of 0 keeps the runner's temporary folder; otherwise TMPDIR is a new folder under
# RUNNER_TEMP whose path has that many characters, and every fixture makes its project under it
# (withr::local_tempdir), which sets the project root's length. Prints the lines the cases' tables
# are made from; exits 1 when a case's script fails.
param(
    [Parameter(Mandatory)][string]$Runtime,
    [Parameter(Mandatory)][string]$Checkout,
    [Parameter(Mandatory)][string]$TestTools,
    [Parameter(Mandatory)][string]$Out,
    [Parameter(Mandatory)][string[]]$Cases
)

$status = 0
$runner = Join-Path $Checkout 'packaging\verify\windows-clean-run.ps1'
$script = Join-Path $Checkout 'packaging\verify\figure-paths.R'
$env:MSR_PROBE_DIR = Join-Path $env:RUNNER_TEMP 'dp'
foreach ($case in $Cases) {
    $name, $rest = $case -split '=', 2
    $length, $list = $rest -split ':', 2
    $items = $list -split ','
    Remove-Item Env:TMPDIR -ErrorAction SilentlyContinue
    if ([int]$length -gt 0 -and [int]$length -lt ($env:RUNNER_TEMP.Length + 3)) {
        # Shorter than a folder under RUNNER_TEMP can be: a folder at the root of RUNNER_TEMP's drive (D:\o is 4).
        $dir = (Split-Path -Qualifier $env:RUNNER_TEMP) + '\' + ('o' * ([int]$length - 3))
        New-Item -ItemType Directory -Force $dir | Out-Null
        $env:TMPDIR = $dir
        "case ${name}: TMPDIR $dir ($($dir.Length) characters)"
    } elseif ([int]$length -gt 0) {
        $pad = [int]$length - $env:RUNNER_TEMP.Length - 1
        $dir = Join-Path $env:RUNNER_TEMP ('t' + ('x' * ($pad - 1)))
        New-Item -ItemType Directory -Force $dir | Out-Null
        $env:TMPDIR = $dir
        "case ${name}: TMPDIR $dir ($($dir.Length) characters)"
    } else {
        "case ${name}: the runner's temporary folder (TEMP $env:TEMP, $($env:TEMP.Length) characters)"
    }
    $log = Join-Path $Out "figures-$name.log"
    $clock = [Diagnostics.Stopwatch]::StartNew()
    & $runner -Bundle $Runtime -Script $script -NoBundleArgument -Arguments (@($Checkout, $TestTools, (Join-Path $Out "figures-$name"), $name) + $items) *> $log
    $code = $LASTEXITCODE
    "case ${name}: exit $code after $([math]::Round($clock.Elapsed.TotalSeconds)) s"
    Get-Content $log | Where-Object { $_ -match '^(CASE|TEST|ITEM|FIGSUM|FIGURES|DEVICES|DEVROW|FIGROW|FIGPAIR|FIGPAIRSUM|RAISED|\s+expectation_)' }
    if ($code -ne 0) {
        $status = 1
        "case ${name}: the end of its log:"
        Get-Content $log -Tail 30
    }
}
Remove-Item Env:TMPDIR -ErrorAction SilentlyContinue
exit $status
