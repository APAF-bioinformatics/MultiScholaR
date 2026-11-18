# KIT-134, build time only, Windows. From the KIT-130 spike's Windows steps
# (.github/workflows/dist-spike.yml on branch kit-130-dist-spike).
#
#   assemble.ps1 -Action r -Bundle <bundle> -RVersion 4.6.1
#       R from the CRAN installer, silently, into <bundle>\R (no registry entries, no icons), and
#       Rtools45 present for the build (it is not part of the bundle).
#   assemble.ps1 -Action pathlength -Folder <application folder> [-Limit 170]
#       Fails if any path under the folder, relative to it, is longer than the limit: 170 leaves
#       room for a 20-character user name under the 260-character limit once installed under
#       %LOCALAPPDATA%\Programs\<product> (revision 2.0, 2.0.6, P4).
param(
    [Parameter(Mandatory)][ValidateSet('r', 'pathlength')][string]$Action,
    [string]$Bundle,
    [string]$RVersion = '4.6.1',
    [string]$Folder,
    [int]$Limit = 170
)
$ErrorActionPreference = 'Stop'

if ($Action -eq 'r') {
    $installer = Join-Path $env:RUNNER_TEMP "R-$RVersion-win.exe"
    curl.exe -fsSL -o $installer "https://cran.r-project.org/bin/windows/base/old/$RVersion/R-$RVersion-win.exe"
    if ($LASTEXITCODE -ne 0) { throw "could not download R $RVersion" }
    $process = Start-Process -FilePath $installer -Wait -PassThru -ArgumentList @(
        '/VERYSILENT', '/SUPPRESSMSGBOXES', '/NORESTART', '/NOICONS',
        "/DIR=$Bundle\R",
        '/MERGETASKS=!recordversion,!associate,!desktopicon,!quicklaunchicon'
    )
    if ($process.ExitCode -ne 0) { throw "the R installer exited with $($process.ExitCode)" }
    & "$Bundle\R\bin\Rscript.exe" --vanilla -e "cat(R.version.string, R.home(), sep = '\n')"
    if ($LASTEXITCODE -ne 0) { throw "the bundled R does not start" }

    if (-not (Test-Path 'C:\rtools45\usr\bin\make.exe')) {
        $rtools = Join-Path $env:RUNNER_TEMP 'rtools45.exe'
        curl.exe -fsSL -o $rtools 'https://github.com/r-hub/rtools45/releases/download/latest/rtools45.exe'
        if ($LASTEXITCODE -ne 0) { throw 'could not download Rtools45' }
        $process = Start-Process -FilePath $rtools -Wait -PassThru -ArgumentList '/VERYSILENT', '/SUPPRESSMSGBOXES', '/NORESTART'
        if ($process.ExitCode -ne 0) { throw "the Rtools installer exited with $($process.ExitCode)" }
    }
    "Rtools45 (build time only): $((Get-Item 'C:\rtools45\usr\bin\make.exe').FullName)"
    exit 0
}

$root = (Resolve-Path -LiteralPath $Folder).Path.TrimEnd('\')
$lengths = Get-ChildItem -LiteralPath $root -Recurse -File -Force |
    ForEach-Object { [pscustomobject]@{ Length = $_.FullName.Length - $root.Length - 1; Path = $_.FullName.Substring($root.Length + 1) } } |
    Sort-Object Length -Descending
"files: $($lengths.Count); longest path under the folder: $($lengths[0].Length) characters (limit $Limit)"
$lengths | Select-Object -First 10 | ForEach-Object { "{0,4}  {1}" -f $_.Length, $_.Path }
$over = @($lengths | Where-Object { $_.Length -gt $Limit })
if ($over.Count -gt 0) { throw "$($over.Count) paths under the folder are longer than $Limit characters" }
