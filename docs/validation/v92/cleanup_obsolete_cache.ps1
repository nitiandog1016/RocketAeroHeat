# 模块一：本次固定清理Q39派生系统；先验全部源码、输入、目标SHA和路径。
$ErrorActionPreference='Stop'
$taskRoot='E:\RocketAeroHeat'
$taskFolder=Join-Path $taskRoot 'docs/validation/v92'
$taskPlanPath=Join-Path $taskFolder 'obsolete_cache_plan.json'
$taskReceipt=Join-Path $taskFolder 'obsolete_cache_execution.json'
$taskJournal=Join-Path $taskFolder 'obsolete_cache_operations.jsonl'
$taskCachePrefix=Join-Path $taskRoot 'outputs/v92_validation/continuous01/step01/current'
$taskCachePrefix=[IO.Path]::GetFullPath($taskCachePrefix)+'\'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '本次已有执行记录，不得重复。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB仍在运行，拒绝清理。' }
$data=Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($data.kind -ne 'v92-verified-obsolete-Q39-system-cache' -or $data.workspace -ne $taskRoot -or @($data.deleteFiles).Count -ne 6) { throw '固定计划不符。' }
function Get-CheckedPath([string]$Relative) {
    if ([IO.Path]::IsPathRooted($Relative)) { throw '必须是项目相对路径。' }
    $path=[IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
    if (-not $path.StartsWith($taskRoot+'\',[StringComparison]::OrdinalIgnoreCase)) { throw '路径超出项目。' }
    $item=Get-Item -LiteralPath $path
    if (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw '不处理重解析文件。' }
    $parent=[IO.DirectoryInfo]::new([IO.Path]::GetDirectoryName($path))
    while ($null -ne $parent -and $parent.FullName.Length -ge $taskRoot.Length) {
        if (($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw '不处理重解析目录。' }
        $parent=$parent.Parent
    }
    return $path
}
function Assert-Record($record) {
    $path=Get-CheckedPath $record.path
    $item=Get-Item -LiteralPath $path
    if ($item.PSIsContainer -or $item.Length -ne $record.bytes -or (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $record.sha256) { throw "内容不符：$($record.path)" }
}
$seen=@{}
foreach ($entry in $data.deleteFiles) {
    $path=Get-CheckedPath $entry.path
    if (-not $path.StartsWith($taskCachePrefix,[StringComparison]::OrdinalIgnoreCase) -or [IO.Path]::GetDirectoryName($path)+'\' -ne $taskCachePrefix -or [IO.Path]::GetFileName($path) -notmatch '^current_peak_(cache|operator|packed|factor|scale|metric)_evidence\.mat$') { throw '不在固定六项派生缓存内。' }
    if ($seen.ContainsKey($path)) { throw '重复目标。' }
    $seen[$path]=$true; Assert-Record $entry
}
foreach ($entry in @($data.protectedFiles)+@($data.sourceModules)) { Assert-Record $entry }
$marker=Join-Path $taskCachePrefix 'CACHE_PRUNED.json'
if (Test-Path -LiteralPath $marker) { throw '清理标记已存在。' }
# 模块二：只删除六个明确文件，不递归删除目录；每项操作即时记录。
$receipt=[ordered]@{status='executing';startedUTC=[DateTime]::UtcNow.ToString('o');planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant();freeBefore=[IO.DriveInfo]::new('E:\').AvailableFreeSpace;deleted=@()}
foreach ($entry in $data.deleteFiles) {
    Assert-Record $entry
    Remove-Item -LiteralPath (Get-CheckedPath $entry.path)
    $receipt.deleted += $entry
    $entry | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
@{status='obsolete-Q39-derived-cache-pruned-after-all-v92-audits';originalPlan='docs/validation/v92/obsolete_cache_plan.json';rebuildFunction=$data.rebuildFunction;rebuildRequiresNewManifest=$true;historicalMATFileSHAReproductionClaimed=$false;note='输入、方向和结果已保留。旧清单不可直接恢复；历史复查应在新目录按对应源码重建新系统，不保证历史MAT文件SHA相同。下一次实际推进必须从Q40重建。'} |
    ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $marker -Encoding UTF8
# 模块三：删除后再次核对保留依赖与源码，记录真实磁盘可用空间变化。
foreach ($entry in @($data.protectedFiles)+@($data.sourceModules)) { Assert-Record $entry }
$receipt.status='completed-six-obsolete-Q39-caches-current-Q40-and-evidence-retained'
$receipt.protectedAndSourceSHAVerifiedAfter=$true
$receipt.completedUTC=[DateTime]::UtcNow.ToString('o')
$receipt.freeAfter=[IO.DriveInfo]::new('E:\').AvailableFreeSpace
$receipt.actualFreeIncreaseBytes=$receipt.freeAfter-$receipt.freeBefore
$receipt.deletedBytes=($data.deleteFiles | Measure-Object -Property bytes -Sum).Sum
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CACHE_CLEANUP_COMPLETED files=$($receipt.deleted.Count) bytes=$($receipt.deletedBytes) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
