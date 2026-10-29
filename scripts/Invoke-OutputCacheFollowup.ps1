# 模块一：固定十一项旧缓存；核对原路径、全部输入和来源后才允许删除。
$ErrorActionPreference='Stop'
$taskRoot=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskFolder=Join-Path $taskRoot 'docs/validation/storage_20261029/followup'
$taskPlanPath=Join-Path $taskFolder 'plan.json'
$taskReceipt=Join-Path $taskFolder 'execution.json'
$taskJournal=Join-Path $taskFolder 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，拒绝重复清理。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB仍在运行，拒绝清理。' }
$taskData=Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskData.kind -ne 'five-legacy-factors-and-six-obsolete-Q40-preparation-caches' -or $taskData.workspace -ne $taskRoot -or @($taskData.deleteFiles).Count -ne 11) { throw '固定清单不符。' }
$taskAllowed=@(
 'outputs/v83_validation/patchhierarchyreal01/patch_hierarchy_factor_evidence.mat',
 'outputs/v83_validation/blockfillreal01/block_fill1_factor_evidence.mat',
 'outputs/v83_validation/overlapFreal01/restricted_overlap_factor_evidence.mat',
 'outputs/v83_validation/floworderfactor01/flow_ordered_factor_evidence.mat',
 'outputs/v83_validation/spatialreductionreal01/spatial_coarse_factor_evidence.mat',
 'outputs/v93_validation/continuous01/step01/current/current_peak_cache_evidence.mat',
 'outputs/v93_validation/continuous01/step01/current/current_peak_operator_evidence.mat',
 'outputs/v93_validation/continuous01/step01/current/current_peak_packed_evidence.mat',
 'outputs/v93_validation/continuous01/step01/current/current_peak_factor_evidence.mat',
 'outputs/v93_validation/continuous01/step01/current/current_peak_scale_evidence.mat',
 'outputs/v93_validation/continuous01/step01/current/current_peak_metric_evidence.mat'
)
function Get-CheckedPath([string]$Relative) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw '只接受工程相对路径。' }
 $path=[IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 if (-not $path.StartsWith($taskRoot+'\',[StringComparison]::OrdinalIgnoreCase)) { throw '路径超出工程。' }
 $info=[IO.FileInfo]::new($path)
 if ($info.Exists -and (($info.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0)) { throw '拒绝重解析文件。' }
 $parent=$info.Directory
 while ($null -ne $parent -and $parent.FullName.Length -ge $taskRoot.Length) {
  if (($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw '拒绝重解析目录。' }
  $parent=$parent.Parent
 }
 return $path
}
function Assert-Record($Entry) {
 $path=Get-CheckedPath $Entry.path
 $item=Get-Item -LiteralPath $path
 if ($item.PSIsContainer -or $item.Length -ne $Entry.bytes -or (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $Entry.sha256) { throw "文件身份不符：$($Entry.path)" }
}
$taskSeen=@{}
foreach ($entry in $taskData.deleteFiles) {
 if ($taskAllowed -notcontains $entry.path -or $taskSeen.ContainsKey($entry.path)) { throw '清理项超范围或重复。' }
 $taskSeen[$entry.path]=$true
 Assert-Record $entry
}
foreach ($entry in $taskData.protectedFiles) {
 if ($taskSeen.ContainsKey($entry.path)) { throw '清理项与保留项冲突。' }
 Assert-Record $entry
}
foreach ($entry in $taskData.sourceFiles) {
 Assert-Record $entry
 if ($null -ne $entry.archivedFile) { Assert-Record $entry.archivedFile }
}
Assert-Record $taskData.inventoryBefore
foreach ($group in $taskData.groups) {
 if (Test-Path -LiteralPath (Get-CheckedPath ($group.folder+'/CACHE_PRUNED.json'))) { throw '已有清理标记。' }
}
Write-Output 'PREFLIGHT_PASSED eleven-explicit-obsolete-cache-files'

# 模块二：使用原生PowerShell逐文件删除；不递归删除目录，即时登记每项操作。
$receipt=[ordered]@{status='executing';startedUTC=[DateTime]::UtcNow.ToString('o');planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant();freeBefore=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace;deleted=@()}
foreach ($entry in $taskData.deleteFiles) {
 Assert-Record $entry
 $path=Get-CheckedPath $entry.path
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw '删除后目标仍存在。' }
 $receipt.deleted += $entry
 $entry | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
foreach ($group in $taskData.groups) {
 @{status='obsolete-derived-cache-pruned';plan='docs/validation/storage_20261029/followup/plan.json';removedArtifacts=$group.removedArtifacts;rebuildScope=$group.rebuildScope;historicalResultsAndInputsRetained=$true;rebuildInNewDirectory=$true;rebuildExecuted=$false;historicalMATFileSHAReproductionClaimed=$false;note='保留原结果、输入与构造源码。旧矩阵审计需要新目录重建；不得将旧清单或旧方向资格套到新文件。下一次推进从Q41新建系统。'} |
  ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Get-CheckedPath ($group.folder+'/CACHE_PRUNED.json')) -Encoding UTF8
}

# 模块三：再次核对保留文件，按磁盘空闲实际变化记录回收量。
foreach ($entry in $taskData.protectedFiles) { Assert-Record $entry }
foreach ($entry in $taskData.sourceFiles) {
 Assert-Record $entry
 if ($null -ne $entry.archivedFile) { Assert-Record $entry.archivedFile }
}
$receipt.status='completed-eleven-explicit-obsolete-caches-removed'
$receipt.completedUTC=[DateTime]::UtcNow.ToString('o')
$receipt.freeAfter=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes=$receipt.freeAfter-$receipt.freeBefore
$receipt.deletedBytes=($taskData.deleteFiles | Measure-Object -Property bytes -Sum).Sum
$receipt.protectedAndSourceVerifiedAfter=$true
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CLEANUP_COMPLETED files=$($receipt.deleted.Count) bytes=$($receipt.deletedBytes) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
