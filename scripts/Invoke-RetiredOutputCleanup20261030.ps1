# 模块一：在任何删除前核对工程绝对路径、固定清单、SHA及保留依赖；拒绝重解析点。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskFolder = Join-Path $taskRoot 'docs/validation/storage_20261030'
$taskPlanPath = Join-Path $taskFolder 'plan.json'
$taskReceipt = Join-Path $taskFolder 'execution.json'
$taskJournal = Join-Path $taskFolder 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，拒绝重复删除。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB正在运行，拒绝维护。' }
$taskData = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskData.workspace -ne $taskRoot -or $taskData.kind -ne 'verified-redundant-intermediates-and-two-retired-QR-factors') { throw '清单身份不符。' }
$taskCanonical = @(
 'outputs/v70_validation/newton_step_comparison/newton_step_trial_raw.mat',
 'outputs/v71_validation/rocket_failure_final/newton_diagnostic_raw.mat',
 'outputs/v71_validation/rocket_failure_replay/newton_diagnostic_raw.mat',
 'outputs/v72_validation/continuation_trial/continuation_raw.mat',
 'outputs/v72_validation/warm_pressure_trial/warm_pressure_raw.mat',
 'outputs/v73_validation/pseudo_cfl_001/pseudo_cfl_raw.mat',
 'outputs/v73_validation/pseudo_cfl_010/pseudo_cfl_raw.mat',
 'outputs/v73_validation/pseudo_cfl_100/pseudo_cfl_raw.mat',
 'outputs/v73_validation/fine_difference_1em09/fine_difference_raw.mat',
 'outputs/v73_validation/fine_difference_7em08/fine_difference_raw.mat',
 'outputs/v73_validation/warm_hlle_davis/warm_flux_raw.mat',
 'outputs/v73_validation/warm_hlle_einfeldt/warm_flux_raw.mat',
 'outputs/v74_validation/local_difference_2em06/local_difference_raw.mat',
 'outputs/v74_validation/local_difference_2em07/local_difference_raw.mat',
 'outputs/v75_validation/spatial_root_double/spatial_root_raw.mat',
 'outputs/v75_validation/spatial_root_two_double_euler/spatial_root_raw.mat',
 'outputs/v83_validation/localrocket02/local_coarse_construction_evidence.mat',
 'outputs/v83_validation/localunit02/local_coarse_construction_evidence.mat',
 'outputs/v91_validation/continuous01/step01/direction_context.mat',
 'outputs/v92_validation/continuous01/step01/direction_context.mat'
)
function Get-CheckedPath([string]$Relative) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw '拒绝绝对输入路径。' }
 $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 if (-not $path.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw '目标超出工程。' }
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
$taskAllowed = @{}
foreach ($group in $taskData.groups) {
 if ($taskCanonical -notcontains $group.file.path) { throw '删除对象不在固定范围。' }
 foreach ($entry in $group.aliases) {
  if ($taskAllowed.ContainsKey($entry.path) -or -not $entry.path.StartsWith('outputs/') -or $entry.sha256 -ne $group.file.sha256) { throw '硬链接范围或身份不符。' }
  $taskAllowed[$entry.path] = $true
 }
 foreach ($marker in $group.markerPaths) { if (Test-Path -LiteralPath (Get-CheckedPath $marker)) { throw '已有维护标记。' } }
}
foreach ($entry in $taskData.deleteFiles) {
 if (-not $taskAllowed.ContainsKey($entry.path)) { throw '删除路径超范围。' }
 Assert-Record $entry
}
if ($taskAllowed.Count -ne @($taskData.deleteFiles).Count) { throw '删除路径数量不符。' }
foreach ($entry in @($taskData.protectedFiles) + @($taskData.sourceFiles) + @($taskData.inventoryBefore)) {
 if ($taskAllowed.ContainsKey($entry.path)) { throw '保留与删除冲突。' }
 Assert-Record $entry
}
Write-Output "PREFLIGHT_PASSED objects=$(@($taskData.groups).Count) paths=$($taskAllowed.Count)"

# 模块二：逐文件删除，不递归删除目录，不改写共享文件；磁盘变化单独实测。
$receipt = [ordered]@{status='executing';startedUTC=[DateTime]::UtcNow.ToString('o');planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant();freeBefore=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace;deleted=@()}
foreach ($entry in $taskData.deleteFiles) {
 Assert-Record $entry
 $path = Get-CheckedPath $entry.path
 if (-not $path.StartsWith((Join-Path $taskRoot 'outputs') + '\',[StringComparison]::OrdinalIgnoreCase)) { throw '删除目标不在outputs。' }
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw '删除后仍存在。' }
 $receipt.deleted += $entry
 $entry | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
foreach ($group in $taskData.groups) {
 foreach ($marker in $group.markerPaths) {
  @{status='verified-old-output-pruned';plan='docs/validation/storage_20261030/plan.json';group=$group;note='重复MAT的全部原变量已在保留最终MAT中逐值核对；QR因子须在新目录按准确历史源码及新准备清单重建。原方向、温度行、失败和检查点保留，不保证重建历史MAT容器SHA。'} | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath (Get-CheckedPath $marker) -Encoding UTF8
 }
}

# 模块三：再次核对保留内容，保存实际执行回执；独立全目录审计另行运行。
foreach ($entry in @($taskData.protectedFiles) + @($taskData.sourceFiles)) { Assert-Record $entry }
$receipt.status = 'completed-verified-output-maintenance'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfter = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes = $receipt.freeAfter - $receipt.freeBefore
$receipt.deletedUniqueBytes = $taskData.plannedUniqueBytes
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CLEANUP_COMPLETED paths=$($receipt.deleted.Count) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
