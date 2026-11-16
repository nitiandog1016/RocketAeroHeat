# 模块职责：仅对已逐字证明相同的开关控制副本建立硬链接，保留两套路径与内容。
$ErrorActionPreference = 'Stop'
$taskRoot = (Resolve-Path -LiteralPath 'E:\RocketAeroHeat\outputs\v246_validation').Path
$taskPrefix = $taskRoot.TrimEnd('\') + '\'
$names = @('state_000000.bin','state_000001.bin','state_000005.bin','coarse_initial.bin','coarse_final.bin','history.csv','restart_flow.csv')
$records = @()
foreach ($name in $names) {
    $sourcePath = (Resolve-Path -LiteralPath (Join-Path $taskRoot ('control_observe01\' + $name))).Path
    $targetPath = (Resolve-Path -LiteralPath (Join-Path $taskRoot ('control_off01\' + $name))).Path
    if (-not $sourcePath.StartsWith($taskPrefix,[StringComparison]::OrdinalIgnoreCase) -or -not $targetPath.StartsWith($taskPrefix,[StringComparison]::OrdinalIgnoreCase)) { throw 'Outside expected workspace' }
    $sourceHash = (Get-FileHash -LiteralPath $sourcePath -Algorithm SHA256).Hash
    $targetHash = (Get-FileHash -LiteralPath $targetPath -Algorithm SHA256).Hash
    $length = (Get-Item -LiteralPath $sourcePath).Length
    if ($sourceHash -ne $targetHash -or $length -ne (Get-Item -LiteralPath $targetPath).Length) { throw 'Distinct data must be preserved' }
    $temporaryPath = $targetPath + '.verified-link'
    if (Test-Path -LiteralPath $temporaryPath) { throw 'Unexpected pending link' }
    New-Item -ItemType HardLink -Path $temporaryPath -Target $sourcePath | Out-Null
    if ((Get-FileHash -LiteralPath $temporaryPath -Algorithm SHA256).Hash -ne $targetHash) { throw 'Link verification failed' }
    Remove-Item -LiteralPath $targetPath -Force
    Move-Item -LiteralPath $temporaryPath -Destination $targetPath
    if ((Get-FileHash -LiteralPath $targetPath -Algorithm SHA256).Hash -ne $targetHash) { throw 'Final data verification failed' }
    $records += [ordered]@{ source=$sourcePath; target=$targetPath; sha256=$sourceHash.ToLowerInvariant(); duplicateLogicalBytes=$length; bothPathsPreserved=$true }
}
$report = [ordered]@{ passed=$true; uniqueDataDeleted=$false; duplicatedContentBytes=($records | Measure-Object -Property duplicateLogicalBytes -Sum).Sum; actualAllocatedBytesMeasured=$false; files=$records }
$json = $report | ConvertTo-Json -Depth 6
[IO.File]::WriteAllText((Join-Path $taskRoot 'deduplication.json'),$json+[Environment]::NewLine,[Text.UTF8Encoding]::new($false))
Write-Output ('Verified duplicate content linked: ' + $report.duplicatedContentBytes + ' bytes; both paths preserved.')
