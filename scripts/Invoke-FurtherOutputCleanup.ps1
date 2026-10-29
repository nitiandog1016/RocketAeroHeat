# 模块一：限定三个旧研究中的四项缓存，复核全部保留输入和历史源码。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskFolder = Join-Path $taskRoot 'docs/validation/storage_20261029/factors'
$taskPlanPath = Join-Path $taskFolder 'plan.json'
$taskReceipt = Join-Path $taskFolder 'execution.json'
$taskJournal = Join-Path $taskFolder 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，拒绝重复清理。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB仍在运行，拒绝清理。' }
$taskData = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
$taskAllowed = @(
 'outputs/v83_validation/blockilureal01/steady_graph_packed_evidence.mat',
 'outputs/v83_validation/blockilureal01/block_ilu_factor_evidence.mat',
 'outputs/v83_validation/coupledpatchreal01/coupled_patch_factor_evidence.mat',
 'outputs/v83_validation/fullFstrengthreal02/full_F_strength_factor_evidence.mat'
)
if ($taskData.kind -ne 'four-obsolete-Q17-graph-and-factor-caches' -or $taskData.workspace -ne $taskRoot -or @($taskData.deleteFiles).Count -ne 4) { throw '固定清单不符。' }
function Get-CheckedPath([string]$Relative) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw '只接受工程相对路径。' }
 $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 if (-not $path.StartsWith($taskRoot + '\',[StringComparison]::OrdinalIgnoreCase)) { throw '路径超出工程。' }
 $info = [IO.FileInfo]::new($path)
 if ($info.Exists -and (($info.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0)) { throw '拒绝重解析文件。' }
 $parent = $info.Directory
 while ($null -ne $parent -and $parent.FullName.Length -ge $taskRoot.Length) {
  if (($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw '拒绝重解析目录。' }
  $parent = $parent.Parent
 }
 return $path
}
function Assert-Record($Entry) {
 $path = Get-CheckedPath $Entry.path
 $item = Get-Item -LiteralPath $path
 if ($item.PSIsContainer -or $item.Length -ne $Entry.bytes -or (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $Entry.sha256) { throw "文件身份不符：$($Entry.path)" }
}
$taskSeen = @{}
foreach ($entry in $taskData.deleteFiles) {
 if ($taskAllowed -notcontains $entry.path -or $taskSeen.ContainsKey($entry.path)) { throw '清理项超范围或重复。' }
 $taskSeen[$entry.path] = $true
 Assert-Record $entry
}
foreach ($entry in $taskData.protectedFiles) {
 if ($taskSeen.ContainsKey($entry.path)) { throw '删除与保留范围冲突。' }
 Assert-Record $entry
}
foreach ($entry in $taskData.sourceFiles) {
 Assert-Record $entry
 if ($null -ne $entry.archivedFile) { Assert-Record $entry.archivedFile }
}
Assert-Record $taskData.inventoryBefore
foreach ($group in $taskData.groups) {
 if (Test-Path -LiteralPath (Get-CheckedPath ($group.folder + '/CACHE_PRUNED.json'))) { throw '已有清理标记。' }
}
Write-Output 'PREFLIGHT_PASSED four-explicit-obsolete-cache-files'

# 模块二：原生PowerShell逐文件删除，绝不递归删除或移动计算目录。
$receipt = [ordered]@{status='executing';startedUTC=[DateTime]::UtcNow.ToString('o');planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant();freeBefore=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace;deleted=@()}
foreach ($entry in $taskData.deleteFiles) {
 Assert-Record $entry
 $path = Get-CheckedPath $entry.path
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw '删除后目标仍存在。' }
 $receipt.deleted += $entry
 $entry | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
foreach ($group in $taskData.groups) {
 @{status='obsolete-Q17-derived-cache-pruned';plan='docs/validation/storage_20261029/factors/plan.json';removedArtifacts=$group.removedArtifacts;originalProcessExitCode=$group.originalExitCode;constructorStudy=$group.constructorStudy;historicalResultsAndSourceStateRetained=$true;rebuildInNewDirectory=$true;rebuildExecuted=$false;historicalMATFileSHAReproductionClaimed=$false;note='旧矩阵/因子和依赖它们的历史审计须先在新目录重建。先从保留Q17的firstOperator重建packed，再构造相应因子；原失败不改记成功，不套用旧SHA或旧清单。当前Q41和v94流线研究保留。'} |
  ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Get-CheckedPath ($group.folder + '/CACHE_PRUNED.json')) -Encoding UTF8
}

# 模块三：删除后再次验SHA，磁盘空闲实际变化与路径长度分开记录。
foreach ($entry in $taskData.protectedFiles) { Assert-Record $entry }
foreach ($entry in $taskData.sourceFiles) {
 Assert-Record $entry
 if ($null -ne $entry.archivedFile) { Assert-Record $entry.archivedFile }
}
$receipt.status = 'completed-four-obsolete-Q17-caches-removed'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfter = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes = $receipt.freeAfter - $receipt.freeBefore
$receipt.deletedBytes = ($taskData.deleteFiles | Measure-Object -Property bytes -Sum).Sum
$receipt.protectedAndSourceVerifiedAfter = $true
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CLEANUP_COMPLETED files=$($receipt.deleted.Count) bytes=$($receipt.deletedBytes) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
