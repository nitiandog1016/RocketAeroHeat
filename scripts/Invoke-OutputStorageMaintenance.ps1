# 模块一：逐条核对固定清单和SHA；所有目标必须位于当前工程，拒绝重解析点。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskFolder = Join-Path $taskRoot 'docs/validation/storage_20261029/text_and_derived'
$taskPlanPath = Join-Path $taskFolder 'plan.json'
$taskReceipt = Join-Path $taskFolder 'execution.json'
$taskJournal = Join-Path $taskFolder 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，拒绝重复执行。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB正在计算，拒绝维护输出。' }
$taskData = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskData.workspace -ne $taskRoot -or $taskData.kind -ne 'three-retired-derived-caches-and-transparent-NTFS-text-compression') { throw '清单身份不符。' }
$taskCanonical = @('outputs/v83_validation/spatialreductionreal01/spatial_reduction_assembly_evidence.mat', 'outputs/v83_validation/fullFreal01/full_F_graph_evidence.mat', 'outputs/v83_validation/momentumlinereal03/momentum_line_factor_evidence.mat')
function Get-CheckedPath([string]$Relative, [bool]$OutputOnly = $false) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw '拒绝绝对输入路径。' }
 $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 $base = if ($OutputOnly) { Join-Path $taskRoot 'outputs' } else { $taskRoot }
 if (-not $path.StartsWith($base + '\', [StringComparison]::OrdinalIgnoreCase)) { throw '目标超出指定目录。' }
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
 $info = Get-Item -LiteralPath $path
 if ($info.PSIsContainer -or $info.Length -ne $Entry.bytes -or (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $Entry.sha256) { throw "文件身份不符：$($Entry.path)" }
}
if (@($taskData.deleteFiles).Count -ne 3) { throw '删除范围不符。' }
foreach ($entry in $taskData.deleteFiles) {
 if ($taskCanonical -notcontains $entry.path) { throw '删除对象不在固定范围。' }
 $null = Get-CheckedPath $entry.path $true
 Assert-Record $entry
}
foreach ($entry in @($taskData.compressFiles) + @($taskData.protectedFiles) + @($taskData.sourceFiles) + @($taskData.inventoryBefore)) { Assert-Record $entry }
foreach ($entry in $taskData.compressFiles) {
 $null = Get-CheckedPath $entry.path $true
 if ([IO.Path]::GetExtension($entry.path).ToLowerInvariant() -notin @('.csv','.js','.json','.html','.txt','.log','.xml','.svg','.m','.py','.md','.css')) { throw '压缩范围超出文本文件。' }
}
foreach ($group in $taskData.groups) {
 $null = Get-CheckedPath $group.markerPath $true
 if (Test-Path -LiteralPath (Get-CheckedPath $group.markerPath)) { throw '已有本轮清理标记。' }
}
Write-Output "PREFLIGHT_PASSED delete=3 compress=$(@($taskData.compressFiles).Count)"
$receipt = [ordered]@{status='executing'; startedUTC=[DateTime]::UtcNow.ToString('o'); planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant(); freeBefore=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace; deleted=@(); compressedObjects=0}

# 模块二：只删除三项已退役构造，不递归删除目录，不改写历史共享MAT。
foreach ($entry in $taskData.deleteFiles) {
 Assert-Record $entry
 $path = Get-CheckedPath $entry.path $true
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw '删除后仍存在。' }
 $receipt.deleted += $entry
 @{action='deleted-retired-derived-cache';file=$entry} | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
foreach ($group in $taskData.groups) {
 @{status='retired-derived-cache-pruned';plan='docs/validation/storage_20261029/text_and_derived/plan.json';originalArtifact=$group.file;upstreamInputs=$group.upstreamInputs;rebuildEntry=$group.rebuildEntry;rebuildInNewDirectory=$true;rebuildExecuted=$false;historicalMATFileSHAReproductionClaimed=$false;note='原状态、数值方向、失败和检查点保留。上游旧缓存也已清理；须按匹配历史源码在新目录恢复完整构造及新来源绑定，不能直接运行旧审计。'} |
  ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Get-CheckedPath $group.markerPath $true) -Encoding UTF8
}

# 模块三：普通NTFS透明压缩，保留原路径和原字节；不用WOF，不压缩整个目录。
$taskCompactLog = Join-Path $taskFolder 'compact.log'
$taskItems = @($taskData.compressFiles)
for ($i = 0; $i -lt $taskItems.Count; $i += 20) {
 $last = [Math]::Min($i + 19, $taskItems.Count - 1)
 $batch = @('/C','/A','/I','/Q')
 foreach ($entry in $taskItems[$i..$last]) { $batch += Get-CheckedPath $entry.path $true }
 $result = & compact.exe @batch 2>&1
 $code = $LASTEXITCODE
 $result | Add-Content -LiteralPath $taskCompactLog -Encoding UTF8
 if ($code -ne 0) { throw "NTFS压缩失败，已执行部分记录在operations.jsonl和compact.log，退出码$code。" }
 foreach ($entry in $taskItems[$i..$last]) {
  $item = Get-Item -LiteralPath (Get-CheckedPath $entry.path $true)
  if (($item.Attributes -band [IO.FileAttributes]::Compressed) -eq 0) { throw "文件没有NTFS压缩标志：$($entry.path)" }
  @{action='NTFS-compressed-original-bytes-retained';path=$entry.path;sha256=$entry.sha256} | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
 }
 $receipt.compressedObjects = $last + 1
 if (($i % 100) -eq 0) { Write-Output "COMPRESSED $($last + 1)/$($taskItems.Count)" }
}

# 模块四：复核内容SHA并记录实际磁盘变化；独立Python核验随后检查全目录身份。
foreach ($entry in @($taskData.compressFiles) + @($taskData.protectedFiles) + @($taskData.sourceFiles)) { Assert-Record $entry }
$receipt.status = 'completed-storage-maintenance'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfter = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes = $receipt.freeAfter - $receipt.freeBefore
$receipt.deletedUniqueFileLengthBytes = $taskData.plannedDeletedUniqueBytes
$receipt | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "STORAGE_MAINTENANCE_COMPLETED actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
