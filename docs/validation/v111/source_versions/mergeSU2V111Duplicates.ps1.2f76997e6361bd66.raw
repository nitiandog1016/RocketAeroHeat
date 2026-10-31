# 模块用途：仅将v111原生具名/基础文件的逐字重复副本合并硬链接；独有结果不删除。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath('E:\RocketAeroHeat\outputs\v111_validation')
$receiptPath = Join-Path $taskRoot 'duplicate_merge.json'
if (Test-Path -LiteralPath $receiptPath) { throw '已有清理回执，禁止重复执行。' }
$pairs = @()
foreach ($case in @('first_probe01','control01','cfl01')) {
    $iteration = if ($case -eq 'first_probe01') { '000000' } else { '000599' }
    foreach ($item in @(@('restart_flow','csv'),@('surface_flow','vtk'))) {
        $pairs += [pscustomobject]@{ target=(Join-Path $taskRoot "$case/$($item[0])_$iteration.$($item[1])"); retained=(Join-Path $taskRoot "$case/$($item[0]).$($item[1])") }
    }
}
foreach ($name in @('writer_default','writer_named')) {
    $pairs += [pscustomobject]@{ target=(Join-Path $taskRoot "native01/$name/mesh.su2"); retained=(Join-Path $taskRoot 'control01/mesh.su2') }
}
# 所有目标先逐项解析、限制目录并核对SHA；不把跨壳字符串用作删除指令。
$checked = @()
foreach ($pair in $pairs) {
    $target = (Resolve-Path -LiteralPath $pair.target).Path
    $retained = (Resolve-Path -LiteralPath $pair.retained).Path
    foreach ($path in @($target,$retained)) {
        if (-not $path.StartsWith($taskRoot + [IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)) { throw '路径越界' }
    }
    if ($target -eq $retained) { throw '不能替换保留文件自身。' }
    $hash = (Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash
    if ($hash -ne (Get-FileHash -LiteralPath $retained -Algorithm SHA256).Hash) { throw '重复文件SHA不一致。' }
    $checked += [pscustomobject]@{ target=$target; retained=$retained; sha256=$hash.ToLower(); bytes=(Get-Item -LiteralPath $target).Length }
}
$freeBefore = (Get-PSDrive E).Free
foreach ($pair in $checked) {
    Remove-Item -LiteralPath $pair.target
    try { New-Item -ItemType HardLink -Path $pair.target -Value $pair.retained | Out-Null }
    catch { Copy-Item -LiteralPath $pair.retained -Destination $pair.target; throw }
    if ((Get-FileHash -LiteralPath $pair.target -Algorithm SHA256).Hash.ToLower() -ne $pair.sha256) { throw '后验SHA不一致。' }
}
$freeAfter = (Get-PSDrive E).Free
$receipt = [pscustomobject]@{ scope=$taskRoot; action='replace-byte-identical-duplicates-with-hardlinks'; pairs=$checked; freeBytesBefore=$freeBefore; freeBytesAfter=$freeAfter; freeBytesChange=($freeAfter-$freeBefore); uniqueScientificResultsDeleted=0; allOriginalPathsPreserved=$true }
$receipt | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $receiptPath -Encoding utf8
$receipt | Select-Object action,freeBytesChange,uniqueScientificResultsDeleted
