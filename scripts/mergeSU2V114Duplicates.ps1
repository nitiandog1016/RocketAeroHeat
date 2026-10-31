# 只合并本版逐字相同的输出；全部原路径、独有状态及历史原清单保留。
$ErrorActionPreference = 'Stop'
$taskRoot = (Resolve-Path -LiteralPath 'E:\RocketAeroHeat').Path
$engineRoot = (Resolve-Path -LiteralPath (Join-Path $taskRoot 'outputs/embedded_workbench/96284313-167b-4afc-a198-39d8a5cc2385')).Path
$freshRoot = (Resolve-Path -LiteralPath (Join-Path $taskRoot 'outputs/v114_validation/area01/corrected_input')).Path
$recordPath = Join-Path $taskRoot 'outputs/v114_validation/storage_receipt.json'
if (Test-Path -LiteralPath $recordPath) { throw '已有维护回执，不得重复。' }
$pairs = @(
    @{ target=(Join-Path $engineRoot 'restart_flow.csv'); retained=(Join-Path $engineRoot 'restart_flow_002999.csv') },
    @{ target=(Join-Path $engineRoot 'surface_flow.vtk'); retained=(Join-Path $engineRoot 'surface_flow_002999.vtk') },
    @{ target=(Join-Path $freshRoot 'mesh.su2'); retained=(Join-Path $engineRoot 'mesh.su2') },
    @{ target=(Join-Path $freshRoot 'case.cfg'); retained=(Join-Path $engineRoot 'case.cfg') }
)
foreach ($pair in $pairs) {
    foreach ($path in @($pair.target,$pair.retained)) {
        $resolved = (Resolve-Path -LiteralPath $path).Path
        if ([IO.Path]::GetDirectoryName($resolved) -notin @($engineRoot,$freshRoot)) { throw '超出已核对目录。' }
        if (-not $resolved.StartsWith($taskRoot + '\',[StringComparison]::OrdinalIgnoreCase)) { throw '超出工程。' }
    }
    $h = (Get-FileHash -LiteralPath $pair.target -Algorithm SHA256).Hash.ToLower()
    if ($h -ne (Get-FileHash -LiteralPath $pair.retained -Algorithm SHA256).Hash.ToLower()) { throw '原字节不一致。' }
    $pair.sha256=$h; $pair.bytes=(Get-Item -LiteralPath $pair.target).Length
}
$freeBefore=(Get-PSDrive E).Free
foreach ($pair in $pairs) {
    Remove-Item -LiteralPath $pair.target
    try { New-Item -ItemType HardLink -Path $pair.target -Value $pair.retained | Out-Null }
    catch { Copy-Item -LiteralPath $pair.retained -Destination $pair.target; throw }
    if ((Get-FileHash -LiteralPath $pair.target -Algorithm SHA256).Hash.ToLower() -ne $pair.sha256) { throw '后验失败。' }
}
$freeAfter=(Get-PSDrive E).Free
[ordered]@{ action='merge-byte-identical-files-with-hardlinks';pairs=$pairs;freeBytesBefore=$freeBefore;freeBytesAfter=$freeAfter;observedFreeIncrease=($freeAfter-$freeBefore);allPathsPreserved=$true;uniqueScientificStatesDeleted=0 } | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $recordPath -Encoding utf8
Get-Content -LiteralPath $recordPath -Raw -Encoding UTF8
