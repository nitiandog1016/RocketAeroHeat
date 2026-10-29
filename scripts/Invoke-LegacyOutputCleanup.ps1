# 模块一：仅清理固定两组九个过时构造缓存，逐项确认绝对路径及全部保留依赖。
$ErrorActionPreference='Stop'
$taskRoot=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskFolder=Join-Path $taskRoot 'docs/validation/storage_20261029'
$taskPlanPath=Join-Path $taskFolder 'plan.json'
$taskReceipt=Join-Path $taskFolder 'execution.json'
$taskJournal=Join-Path $taskFolder 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，不重复删除。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB仍在运行，拒绝清理。' }
$taskData=Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskData.kind -ne 'v83-two-obsolete-construction-cache-groups' -or $taskData.workspace -ne $taskRoot -or @($taskData.deleteFiles).Count -ne 9) { throw '固定计划不符。' }
$taskAllowed=@(
 'outputs/v83_validation/currentparent01/current_operator_evidence.mat',
 'outputs/v83_validation/currentparent01/current_packed_evidence.mat',
 'outputs/v83_validation/currentparent01/current_factor_evidence.mat',
 'outputs/v83_validation/currentparent01/current_image_scale_evidence.mat',
 'outputs/v83_validation/peakweightreal01/current_peak_cache_evidence.mat',
 'outputs/v83_validation/peakweightreal01/current_peak_operator_evidence.mat',
 'outputs/v83_validation/peakweightreal01/current_peak_packed_evidence.mat',
 'outputs/v83_validation/peakweightreal01/current_peak_factor_evidence.mat',
 'outputs/v83_validation/peakweightreal01/current_peak_scale_evidence.mat'
)
function Get-CheckedPath([string]$Relative) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw '必须是工程相对路径。' }
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
 if ($taskAllowed -notcontains $entry.path -or $taskSeen.ContainsKey($entry.path)) { throw '删除项超出固定九项或重复。' }
 $taskSeen[$entry.path]=$true
 Assert-Record $entry
}
foreach ($entry in $taskData.protectedFiles) {
 if ($taskSeen.ContainsKey($entry.path)) { throw '删除项与保留项冲突。' }
 Assert-Record $entry
}
foreach ($entry in $taskData.sourceFiles) {
 Assert-Record $entry
 if ($null -ne $entry.archivedFile) { Assert-Record $entry.archivedFile }
}
Assert-Record $taskData.inventoryBefore
foreach ($group in $taskData.groups) {
 if (Test-Path -LiteralPath (Get-CheckedPath ($group.folder+'/CACHE_PRUNED.json'))) { throw '清理标记已存在。' }
}
Write-Output 'PREFLIGHT_PASSED nine-explicit-cache-files'
# 模块二：只删除已核对的明确文件，不递归删除目录；立即记录每项操作。
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
 @{status='obsolete-v83-construction-caches-pruned';plan='docs/validation/storage_20261029/plan.json';removedArtifacts=$group.removedArtifacts;historicalResultsAndInputsRetained=$true;rebuildInNewDirectory=$true;rebuildExecuted=$false;historicalMATFileSHAReproductionClaimed=$false;note='构造模块的源码与输入已保留。历史方向和模型仅用于追溯；复查需新目录重建。本标记不允许将旧缓存身份套到新系统。'} |
  ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Get-CheckedPath ($group.folder+'/CACHE_PRUNED.json')) -Encoding UTF8
}
# 模块三：再次核对保留依赖，并以磁盘空闲实际变化报告回收量。
foreach ($entry in $taskData.protectedFiles) { Assert-Record $entry }
foreach ($entry in $taskData.sourceFiles) { Assert-Record $entry }
$receipt.status='completed-nine-obsolete-caches-removed'
$receipt.completedUTC=[DateTime]::UtcNow.ToString('o')
$receipt.freeAfter=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes=$receipt.freeAfter-$receipt.freeBefore
$receipt.deletedBytes=($taskData.deleteFiles | Measure-Object -Property bytes -Sum).Sum
$receipt.protectedAndSourceVerifiedAfter=$true
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CLEANUP_COMPLETED files=$($receipt.deleted.Count) bytes=$($receipt.deletedBytes) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
