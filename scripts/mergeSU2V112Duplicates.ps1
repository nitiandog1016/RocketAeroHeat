# 模块用途：仅合并v112已核实的逐字重复输出，保留全部路径和独有数值状态。
$ErrorActionPreference = 'Stop'
$taskRoot = (Resolve-Path -LiteralPath 'E:\RocketAeroHeat').Path
$planPath = Join-Path $taskRoot 'outputs/v112_validation/storage_plan.json'
$receiptPath = Join-Path $taskRoot 'outputs/v112_validation/storage_receipt.json'
if (Test-Path -LiteralPath $receiptPath) { throw '已有回执，不得重复清理。' }
$plan = Get-Content -LiteralPath $planPath -Raw -Encoding UTF8 | ConvertFrom-Json
$allowed = @('outputs/v112_validation/mg01',
    'outputs/embedded_workbench/d7a55815-b446-4345-b076-cbac806f8b3e',
    'outputs/embedded_workbench/477a440e-6d6f-4179-b90e-66e82d6f2540') | ForEach-Object { (Resolve-Path -LiteralPath (Join-Path $taskRoot $_)).Path }
# 原生PowerShell逐项解析绝对目录并核SHA，删除仅限确有保留副本的单文件。
foreach ($pair in $plan.pairs) {
    foreach ($path in @($pair.target,$pair.retained)) {
        $resolved = (Resolve-Path -LiteralPath $path).Path
        if ([IO.Path]::GetDirectoryName($resolved) -notin $allowed) { throw '文件超出本轮明确范围。' }
        if ((Get-FileHash -LiteralPath $resolved -Algorithm SHA256).Hash.ToLower() -ne $pair.sha256) { throw 'SHA不符，停止。' }
    }
    if ($pair.target -eq $pair.retained) { throw '保留文件不能删除自身。' }
}
$freeBefore = (Get-PSDrive E).Free
foreach ($pair in $plan.pairs) {
    Remove-Item -LiteralPath $pair.target
    try { New-Item -ItemType HardLink -Path $pair.target -Value $pair.retained | Out-Null }
    catch { Copy-Item -LiteralPath $pair.retained -Destination $pair.target; throw }
    if ((Get-FileHash -LiteralPath $pair.target -Algorithm SHA256).Hash.ToLower() -ne $pair.sha256) { throw '后验SHA不符。' }
}
$freeAfter = (Get-PSDrive E).Free
$receipt = [ordered]@{ action='merge-byte-identical-files-with-hardlinks'; pairs=$plan.pairs; freeBytesBefore=$freeBefore; freeBytesAfter=$freeAfter; observedFreeIncrease=($freeAfter-$freeBefore); allPathsPreserved=$true; uniqueScientificStatesDeleted=0 }
$receipt | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $receiptPath -Encoding utf8
$receipt | Select-Object action,observedFreeIncrease
