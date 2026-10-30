$ErrorActionPreference = 'Stop'
# 模块一：只接受固定十轮旧研究；复核绝对路径、原文件、保留文件和全部来源。
$taskRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$taskDestination = Join-Path $taskRoot 'docs/validation/storage_20261030/nested_operators'
$taskPlanPath = Join-Path $taskDestination 'plan.json'
$taskReceiptPath = Join-Path $taskDestination 'execution.json'
$taskJournalPath = Join-Path $taskDestination 'operations.jsonl'
if (Test-Path -LiteralPath $taskReceiptPath) { throw 'Existing execution receipt.' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB is running.' }
$taskPlan = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskPlan.workspace -ne $taskRoot -or $taskPlan.kind -ne 'retired-nested-operator-cache-cleanup') { throw 'Plan identity mismatch.' }
$taskVersions = @(1,2,3,4,5,6,7,8,18,19)
$taskCanonical = @($taskVersions | ForEach-Object { 'outputs/v83_validation/krylov{0:D2}/wall_velocity_krylov_evidence.mat' -f $_ })
$taskCacheFields = @('diagonal','lower','upper','minusCoupling','plusCoupling','secondLower','secondUpper','wlsNormalFactor')
function Get-CheckedPath([string]$Relative) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw 'Absolute input path refused.' }
 $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 if (-not $path.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Path outside workspace.' }
 $item = Get-Item -LiteralPath $path -ErrorAction SilentlyContinue
 if ($item -and ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) { throw 'Reparse file refused.' }
 $parent = [IO.DirectoryInfo]::new([IO.Path]::GetDirectoryName($path))
 while ($parent -and $parent.FullName.Length -ge $taskRoot.Length) {
  if ($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) { throw 'Reparse directory refused.' }
  $parent = $parent.Parent
 }
 return $path
}
function Assert-Record($Entry) {
 $path = Get-CheckedPath $Entry.path
 $item = Get-Item -LiteralPath $path
 if ($item.PSIsContainer -or $item.Length -ne $Entry.bytes -or (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $Entry.sha256) { throw "File identity mismatch: $($Entry.path)" }
}
if (@($taskPlan.groups).Count -ne $taskCanonical.Count) { throw 'Target count mismatch.' }
$taskSeen = @{}
foreach ($group in $taskPlan.groups) {
 $relative = $group.original.path
 if ($taskCanonical -notcontains $relative -or $taskSeen.ContainsKey($relative)) { throw 'Not in canonical scope.' }
 $taskSeen[$relative] = $true
 $tail = $relative.Substring('outputs/'.Length)
 if ($group.retained.path -ne ('outputs/retained_operator_results_20261030/' + $tail)) { throw 'Retained path mismatch.' }
 if ($group.markerPath -ne (($relative.Substring(0, $relative.LastIndexOf('/'))) + '/OPERATOR_CACHE_PRUNED_20261030.json')) { throw 'Marker mismatch.' }
 if (Test-Path -LiteralPath (Get-CheckedPath $group.markerPath)) { throw 'Existing cleanup marker.' }
 foreach ($field in $group.resultScan.removed.PSObject.Properties.Name) {
  if (-not $field.StartsWith('result.firstOperator.') -or $taskCacheFields -notcontains $field.Substring('result.firstOperator.'.Length)) { throw 'Unexpected removed field.' }
 }
 Assert-Record $group.original
 Assert-Record $group.retained
}
foreach ($entry in @($taskPlan.protectedFiles) + @($taskPlan.sourceFiles)) {
 if ($taskSeen.ContainsKey($entry.path)) { throw 'Protected target conflict.' }
 Assert-Record $entry
}
Write-Output "PREFLIGHT_PASSED targets=$($taskPlan.groups.Count)"

# 模块二：再查保留副本后逐文件删除；科学字段与当前状态不在删除范围。
$receipt = [ordered]@{
 status='executing'; startedUTC=[DateTime]::UtcNow.ToString('o')
 planSHA256=(Get-FileHash -LiteralPath $taskPlanPath -Algorithm SHA256).Hash.ToLowerInvariant()
 freeBeforeDelete=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace; deleted=@()
 plannedNetContentBytesReclaimed=$taskPlan.netContentBytesReclaimed
 retainedScientificResultBytes=$taskPlan.retainedBytes
}
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceiptPath -Encoding UTF8
foreach ($group in $taskPlan.groups) {
 Assert-Record $group.retained
 Assert-Record $group.original
 $path = Get-CheckedPath $group.original.path
 if (-not $path.StartsWith((Join-Path $taskRoot 'outputs\v83_validation') + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Deletion outside fixed outputs directory.' }
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw 'File remains after deletion.' }
 @{
  status='nested-operator-caches-pruned-science-retained'
  plan='docs/validation/storage_20261030/nested_operators/plan.json'
  original=$group.original; retained=$group.retained
  removedFieldPaths=@($group.resultScan.removed.PSObject.Properties.Name)
  retainedScientificFields=$group.retainedScientificFields
  originalContainerSHAReproductionClaimed=$false
  note='Only retired preconditioner matrices were removed. All other raw MAT matrix elements remain exact. Rebuild caches with preserved historical sources and inputs in a new directory; old container paths and SHA bindings cannot be reused.'
 } | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath (Get-CheckedPath $group.markerPath) -Encoding UTF8
 $receipt.deleted += $group.original
 $group.original | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournalPath -Encoding UTF8
 $receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceiptPath -Encoding UTF8
}

# 模块三：独立复核另作回执；删除窗口包含已提前迁存的结果，不能当作净回收。
foreach ($entry in @($taskPlan.protectedFiles) + @($taskPlan.sourceFiles) + @($taskPlan.groups.retained)) { Assert-Record $entry }
$receipt.status = 'completed-nested-operator-cache-cleanup'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfterDelete = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.deleteWindowFreeIncreaseBytes = $receipt.freeAfterDelete - $receipt.freeBeforeDelete
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceiptPath -Encoding UTF8
Write-Output "CLEANUP_COMPLETED files=$($receipt.deleted.Count) netContentBytes=$($taskPlan.netContentBytesReclaimed)"
