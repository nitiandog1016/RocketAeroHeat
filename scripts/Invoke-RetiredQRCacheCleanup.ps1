# 模块一：复核固定的十一份旧QR缓存和全部硬链接，删除操作只使用原生PowerShell。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskFolder = Join-Path $taskRoot 'docs/validation/storage_20261029/retired_qr'
$taskPlanPath = Join-Path $taskFolder 'plan.json'
$taskReceipt = Join-Path $taskFolder 'execution.json'
$taskJournal = Join-Path $taskFolder 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，拒绝重复清理。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB仍在运行，拒绝清理。' }
$taskData = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
$taskCanonical = @('outputs/v85_validation/qr_real01/augmented_qr_factor.mat', 'outputs/v86_validation/comparison01/temperature_projection_system.mat')
foreach ($version in @(87,89)) { foreach ($step in 1..3) { $taskCanonical += ('outputs/v{0}_validation/continuous01/step{1:00}/direction_context.mat' -f $version,$step) } }
foreach ($step in 1..2) { $taskCanonical += ('outputs/v88_validation/continuous01/step{0:00}/direction_context.mat' -f $step) }
$taskCanonical += 'outputs/v88_validation/recovery01/direction_context.mat'
if ($taskData.kind -ne 'eleven-retired-QR-construction-cache-objects' -or $taskData.workspace -ne $taskRoot -or @($taskData.groups).Count -ne 11) { throw '固定范围不符。' }
function Get-CheckedPath([string]$Relative) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw '拒绝绝对输入路径。' }
 $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 if (-not $path.StartsWith($taskRoot + '\',[StringComparison]::OrdinalIgnoreCase)) { throw '目标超出工程。' }
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
$taskAllowed = @{}
foreach ($group in $taskData.groups) {
 if ($taskCanonical -notcontains $group.file.path) { throw '构造缓存不在固定范围。' }
 foreach ($entry in $group.aliases) {
  if ($taskAllowed.ContainsKey($entry.path) -or -not $entry.path.StartsWith('outputs/') -or $entry.sha256 -ne $group.file.sha256) { throw '别名重复或身份不符。' }
  $taskAllowed[$entry.path] = $true
 }
 foreach ($marker in $group.markerPaths) { if (Test-Path -LiteralPath (Get-CheckedPath $marker)) { throw '已有清理标记。' } }
}
foreach ($entry in $taskData.deleteFiles) {
 if (-not $taskAllowed.ContainsKey($entry.path)) { throw '删除清单超范围。' }
 Assert-Record $entry
}
if ($taskAllowed.Count -ne @($taskData.deleteFiles).Count) { throw '删除与别名数量不同。' }
foreach ($entry in @($taskData.protectedFiles) + @($taskData.sourceFiles) + @($taskData.inventoryBefore)) {
 if ($taskAllowed.ContainsKey($entry.path)) { throw '保留与删除范围冲突。' }
 Assert-Record $entry
}
Write-Output 'PREFLIGHT_PASSED eleven-retired-QR-cache-objects'

# 模块二：逐个删除已核实的文件路径，不递归删除目录，不改写历史共享文件。
$receipt = [ordered]@{status='executing';startedUTC=[DateTime]::UtcNow.ToString('o');planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant();freeBefore=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace;deleted=@()}
foreach ($entry in $taskData.deleteFiles) {
 Assert-Record $entry
 $path = Get-CheckedPath $entry.path
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw '删除后文件仍在。' }
 $receipt.deleted += $entry
 $entry | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
foreach ($group in $taskData.groups) {
 foreach ($marker in $group.markerPaths) {
  @{status='retired-QR-construction-cache-pruned';plan='docs/validation/storage_20261029/retired_qr/plan.json';originalArtifact=$group.file;aliases=$group.aliases;preparationFile=$group.preparationFile;rowsExactlyRetained=$group.rowsExactlyRetainedInActiveDirection;rebuildInNewDirectory=$true;rebuildExecuted=$false;historicalMATFileSHAReproductionClaimed=$false;note='上游准备系统此前已清理；从保留原状态、几何及匹配历史源码重建准备系统，再重建QR/温度行/投影，生成新清单。原结果、失败和检查点保留，不能直接使用旧审计入口。'} |
   ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Get-CheckedPath $marker) -Encoding UTF8
 }
}

# 模块三：复核保留文件并记录实际磁盘空闲变化，不将硬链接路径长度累计为实收。
foreach ($entry in @($taskData.protectedFiles) + @($taskData.sourceFiles)) { Assert-Record $entry }
$receipt.status = 'completed-eleven-retired-QR-cache-objects'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfter = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes = $receipt.freeAfter - $receipt.freeBefore
$receipt.deletedPathBytes = ($taskData.deleteFiles | Measure-Object -Property bytes -Sum).Sum
$receipt.deletedUniqueBytes = $taskData.plannedUniqueBytes
$receipt.protectedAndSourceVerifiedAfter = $true
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CLEANUP_COMPLETED paths=$($receipt.deleted.Count) uniqueBytes=$($receipt.deletedUniqueBytes) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
