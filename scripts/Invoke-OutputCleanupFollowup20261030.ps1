# 模块一：复核固定范围、绝对路径与全部SHA，拒绝重解析点和运行中的MATLAB。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskFolder = Join-Path $taskRoot 'docs/validation/storage_20261030/retired_contexts_and_archives'
$taskPlanPath = Join-Path $taskFolder 'plan.json'
$taskReceipt = Join-Path $taskFolder 'execution.json'
$taskJournal = Join-Path $taskFolder 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，禁止重复删除。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB运行中，拒绝维护。' }
$taskData = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskData.workspace -ne $taskRoot -or $taskData.kind -ne 'three-retired-QR-contexts-and-eight-v83-source-archives') { throw '清单身份不符。' }
$taskCanonical = @(
 'outputs/v88_validation/continuous01/step03/direction_context.mat',
 'outputs/v90_validation/continuous01/step01/direction_context.mat',
 'outputs/v90_validation/recovery01/direction_context.mat',
 'outputs/v83_validation/RocketAeroHeat_v83_current_state_complete_replacements.zip',
 'outputs/v83_validation/RocketAeroHeat_v83_current_state_complete_replacements_attempt01.zip',
 'outputs/v83_validation/RocketAeroHeat_v83_current_state_research_source.zip',
 'outputs/v83_validation/RocketAeroHeat_v83_current_state_research_source_attempt01.zip',
 'outputs/v83_validation/RocketAeroHeat_v83_peak_complete_replacements.zip',
 'outputs/v83_validation/RocketAeroHeat_v83_peak_research_source.zip',
 'outputs/v83_validation/RocketAeroHeat_v83_riemann_complete_replacements.zip',
 'outputs/v83_validation/RocketAeroHeat_v83_riemann_research_source.zip'
)
function Get-CheckedPath([string]$Relative) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw '拒绝绝对输入路径。' }
 $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 if (-not $path.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw '目标超出项目。' }
 $item = Get-Item -LiteralPath $path -ErrorAction SilentlyContinue
 if ($item -and ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) { throw '拒绝重解析文件。' }
 $parent = [IO.DirectoryInfo]::new([IO.Path]::GetDirectoryName($path))
 while ($parent -and $parent.FullName.Length -ge $taskRoot.Length) {
  if ($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw '拒绝重解析目录。' }
  $parent = $parent.Parent
 }
 return $path
}
function Assert-Record($Entry) {
 $path = Get-CheckedPath $Entry.path
 $item = Get-Item -LiteralPath $path
 if ($item.PSIsContainer -or $item.Length -ne $Entry.bytes -or (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $Entry.sha256) { throw "文件身份不符：$($Entry.path)" }
}
if (@($taskData.groups).Count -ne 11 -or @($taskData.deleteFiles).Count -ne 14) { throw '范围数量不符。' }
$taskAllowed = @{}
foreach ($group in $taskData.groups) {
 if ($taskCanonical -notcontains $group.file.path) { throw '对象不在固定范围。' }
 foreach ($entry in $group.aliases) {
  $candidate = Get-CheckedPath $entry.path
  $canonical = Get-CheckedPath $group.file.path
  if ($entry.path -ne $group.file.path) {
   $relativeTail = $group.file.path.Substring('outputs/'.Length)
   $expected = Join-Path $taskRoot ('outputs/' + $relativeTail.Split('/')[0] + '/delivery/isolated/RocketAeroHeat/outputs/' + $relativeTail)
   if ($candidate -ne [IO.Path]::GetFullPath($expected)) { throw '别名不在已知隔离目录。' }
  }
  if ($taskAllowed.ContainsKey($entry.path) -or $entry.sha256 -ne $group.file.sha256) { throw '重复路径或字节身份不符。' }
  $taskAllowed[$entry.path] = $true
 }
 foreach ($marker in $group.markerPaths) { if (Test-Path -LiteralPath (Get-CheckedPath $marker)) { throw '已有清理标记。' } }
}
foreach ($entry in $taskData.deleteFiles) {
 if (-not $taskAllowed.ContainsKey($entry.path)) { throw '删除路径超范围。' }
 Assert-Record $entry
}
foreach ($entry in @($taskData.protectedFiles) + @($taskData.sourceFiles) + @($taskData.inventoryBefore) + @($taskData.zipMembers)) {
 if ($taskAllowed.ContainsKey($entry.path)) { throw '保留与删除冲突。' }
 Assert-Record $entry
}
Write-Output "PREFLIGHT_PASSED objects=11 paths=14"

# 模块二：逐文件删除；保留每个非factor变量的原压缩记录及每个ZIP成员的准确替代字节。
$receipt = [ordered]@{status='executing';startedUTC=[DateTime]::UtcNow.ToString('o');planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant();freeBefore=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace;deleted=@()}
foreach ($entry in $taskData.deleteFiles) {
 Assert-Record $entry
 $path = Get-CheckedPath $entry.path
 if (-not $path.StartsWith((Join-Path $taskRoot 'outputs') + '\', [StringComparison]::OrdinalIgnoreCase)) { throw '删除目标不在outputs。' }
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw '删除后仍存在。' }
 $receipt.deleted += $entry
 $entry | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
foreach ($marker in @($taskData.groups.markerPaths | Sort-Object -Unique)) {
 $groups = @($taskData.groups | Where-Object { $_.markerPaths -contains $marker })
 @{status='verified-retired-output-pruned';plan='docs/validation/storage_20261030/retired_contexts_and_archives/plan.json';groups=$groups;note='原QR上下文所有非factor压缩记录原字节均保留在mat_metadata；旧QR需在新目录重建。ZIP每个成员均有SHA一致的保留文件或差异原文。原状态、独有结果、失败和检查点未删。'} | ConvertTo-Json -Depth 15 | Set-Content -LiteralPath (Get-CheckedPath $marker) -Encoding UTF8
}

# 模块三：核对保留SHA并记录实测磁盘变化；全目录身份由Python独立审计。
foreach ($entry in @($taskData.protectedFiles) + @($taskData.sourceFiles)) { Assert-Record $entry }
$receipt.status = 'completed-retired-context-and-archive-cleanup'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfter = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes = $receipt.freeAfter - $receipt.freeBefore
$receipt.deletedUniqueBytes = $taskData.plannedUniqueBytes
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CLEANUP_COMPLETED paths=$($receipt.deleted.Count) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
