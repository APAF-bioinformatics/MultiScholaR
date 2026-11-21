# KIT-134 / KIT-131 (branch run-k134-rc): record the scratch folders the package makes under a folder (the TMP of the
# worker lane), so that where a worker staged the session archive and drew its figures is observed, not inferred.
# Polls every 100 ms until <StopFile> exists; writes one line per folder the first time it is seen:
#   WATCH <seconds since start> <depth 0|1> <path> <path length>
# Depth 0 = made directly in the watched folder (a worker's workflowShortScratchRoot() = TMP); depth 1 = inside a
# folder of it (an R session's own tempdir(), RtmpXXXX). Only names of the package's scratch folders are recorded:
# archive-stage-*, wfx-*, dia-qc-figures-*, dia-qc-composite-*, native-export-*.
param(
    [Parameter(Mandatory)][string]$Folder,
    [Parameter(Mandatory)][string]$Log,
    [Parameter(Mandatory)][string]$StopFile
)
$seen = @{}
$clock = [Diagnostics.Stopwatch]::StartNew()
$pattern = '^(archive-stage-|wfx-|dia-qc-figures-|dia-qc-composite-|native-export-)'
"WATCH START $Folder" | Out-File -FilePath $Log -Encoding utf8
while (-not (Test-Path -LiteralPath $StopFile)) {
    foreach ($depth in 0, 1) {
        $items = if ($depth -eq 0) { Get-ChildItem -LiteralPath $Folder -Directory -Force -ErrorAction SilentlyContinue }
                 else { Get-ChildItem -LiteralPath $Folder -Directory -Force -ErrorAction SilentlyContinue |
                        ForEach-Object { Get-ChildItem -LiteralPath $_.FullName -Directory -Force -ErrorAction SilentlyContinue } }
        foreach ($item in $items) {
            if ($item.Name -match $pattern -and -not $seen.ContainsKey($item.FullName)) {
                $seen[$item.FullName] = $true
                "WATCH {0:N1} {1} {2} {3}" -f $clock.Elapsed.TotalSeconds, $depth, $item.FullName, $item.FullName.Length |
                    Out-File -FilePath $Log -Append -Encoding utf8
            }
        }
    }
    Start-Sleep -Milliseconds 100
}
"WATCH STOP after $([math]::Round($clock.Elapsed.TotalSeconds)) s; $($seen.Count) folders seen" | Out-File -FilePath $Log -Append -Encoding utf8
