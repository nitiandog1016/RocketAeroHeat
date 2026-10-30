$ErrorActionPreference = 'Stop'
# 模块一：固定八份旧研究，只接受本轮计划和已核验的保留结果。
$taskRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$taskDestination = Join-Path $taskRoot 'docs/validation/storage_20261031/graph_operators'
$taskPlanPath = Join-Path $taskDestination 'plan.json'
$taskReceiptPath = Join-Path $taskDestination 'execution.json'
$taskJournalPath = Join-Path $taskDestination 'operations.jsonl'
if (Test-Path -LiteralPath $taskReceiptPath) { throw 'Existing execution receipt.' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB is running.' }
$taskPlan = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskPlan.workspace -ne $taskRoot -or $taskPlan.kind -ne 'retired-nested-operator-cache-cleanup' -or $taskPlan.scope -ne 'retired-v83-graph-operators-20261031') { throw 'Plan identity mismatch.' }
$taskVersions = @(9,10,11,12,13,14,15,17)
$taskCanonical = @($taskVersions | ForEach-Object { 'outputs/v83_validation/krylov{0:D2}/wall_velocity_krylov_evidence.mat' -f $_ })
$taskCacheFields = @('diagonal','lower','upper','minusCoupling','plusCoupling','secondLower','secondUpper','wlsNormalFactor','lineInverse','lineUpperFactor')
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
 if ($group.retained.path -ne ('outputs/retained_graph_results_20261031/' + $tail)) { throw 'Retained path mismatch.' }
 if ($group.markerPath -ne (($relative.Substring(0, $relative.LastIndexOf('/'))) + '/OPERATOR_CACHE_PRUNED_20261031.json')) { throw 'Marker mismatch.' }
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

# 模块二：逐文件复核并删除原容器；不递归删除目录、不改写共享历史文件。
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
 if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB started during cleanup.' }
 Remove-Item -LiteralPath $path
 if (Test-Path -LiteralPath $path) { throw 'File remains after deletion.' }
 @{
  status='retired-graph-operator-caches-pruned-science-retained'
  plan='docs/validation/storage_20261031/graph_operators/plan.json'
  original=$group.original; retained=$group.retained
  removedFieldPaths=@($group.resultScan.removed.PSObject.Properties.Name)
  retainedScientificFields=$group.retainedScientificFields
  originalContainerSHAReproductionClaimed=$false
  note='Only retired preconditioner matrices were removed. All other raw MAT elements remain exact. Rebuild caches and upstream inputs in a new directory with accurate historical sources and new bindings. Do not reuse old container paths or SHA.'
 } | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath (Get-CheckedPath $group.markerPath) -Encoding UTF8
 $receipt.deleted += $group.original
 $group.original | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournalPath -Encoding UTF8
 $receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceiptPath -Encoding UTF8
}

# 模块三：再查保留结果；删除窗口包含提前迁存的数据，不把它当作净回收量。
foreach ($entry in @($taskPlan.protectedFiles) + @($taskPlan.sourceFiles) + @($taskPlan.groups.retained)) { Assert-Record $entry }
$receipt.status = 'completed-nested-operator-cache-cleanup'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfterDelete = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.deleteWindowFreeIncreaseBytes = $receipt.freeAfterDelete - $receipt.freeBeforeDelete
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceiptPath -Encoding UTF8
Write-Output "CLEANUP_COMPLETED files=$($receipt.deleted.Count) netContentBytes=$($taskPlan.netContentBytesReclaimed)"
