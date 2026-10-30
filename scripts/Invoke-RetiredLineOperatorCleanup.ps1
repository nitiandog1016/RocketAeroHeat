$ErrorActionPreference = 'Stop'
# 模块一：固定20份旧结果，复核迁存文件、全部硬链接和当前输入，不递归删除目录。
$taskRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$taskDestination = Join-Path $taskRoot 'docs/validation/storage_20261031/retired_lines'
$taskPlanPath = Join-Path $taskDestination 'plan.json'
$taskReceiptPath = Join-Path $taskDestination 'execution.json'
$taskJournalPath = Join-Path $taskDestination 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceiptPath) -or (Test-Path -LiteralPath $taskJournalPath)) { throw 'Existing execution record.' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB is running.' }
& 'D:\anacondapython\python.exe' (Join-Path $PSScriptRoot 'compactRetiredLineOperators.py') 'predelete'
if ($LASTEXITCODE -ne 0) { throw 'Independent predelete check failed.' }
$taskPlan = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskPlan.workspace -ne $taskRoot -or $taskPlan.kind -ne 'retired-line-operator-cache-cleanup-20261031') { throw 'Plan identity mismatch.' }
$taskCanonical = @()
foreach ($version in @(20,21,22,23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)) {
 $taskCanonical += 'outputs/v83_validation/krylov{0:D2}/wall_velocity_krylov_evidence.mat' -f $version
}
function Get-CheckedPath([string]$Relative, [bool]$OutputOnly = $false) {
 if ([IO.Path]::IsPathRooted($Relative)) { throw 'Absolute input path refused.' }
 $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
 $base = if ($OutputOnly) { Join-Path $taskRoot 'outputs' } else { $taskRoot }
 if (-not $path.StartsWith($base + '\', [StringComparison]::OrdinalIgnoreCase)) { throw 'Path outside intended directory.' }
 $info = [IO.FileInfo]::new($path)
 if ($info.Exists -and ($info.Attributes -band [IO.FileAttributes]::ReparsePoint)) { throw 'Reparse file refused.' }
 $parent = $info.Directory
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
if (@($taskPlan.groups).Count -ne 20) { throw 'Target count mismatch.' }
$taskSeen = @{}; $taskAliases = @{}
foreach ($group in $taskPlan.groups) {
 $relative = $group.original.path
 if ($taskCanonical -notcontains $relative -or $taskSeen.ContainsKey($relative)) { throw 'Not in canonical scope.' }
 $taskSeen[$relative] = $true
 if ($group.retained.path -ne ('outputs/retained_line_results_20261031/' + $relative.Substring(8))) { throw 'Retained path mismatch.' }
 $fields = @($group.resultScan.removed.PSObject.Properties.Name)
 if ($fields.Count -ne 3) { throw 'Unexpected removed field count.' }
 foreach ($field in $fields) {
  if ($field -notin @('result.firstOperator.diagonal','result.firstOperator.lower','result.firstOperator.upper')) { throw 'Unexpected removed field.' }
 }
 Assert-Record $group.original
 Assert-Record $group.retained
 foreach ($alias in $group.aliases) {
  $null = Get-CheckedPath $alias $true
  if ($taskAliases.ContainsKey($alias)) { throw 'Duplicate alias.' }
  $taskAliases[$alias] = $true
  $marker = [IO.Path]::GetDirectoryName($alias) + '/LINE_CACHE_PRUNED_20261031.json'
  if (Test-Path -LiteralPath (Get-CheckedPath $marker $true)) { throw 'Existing cleanup marker.' }
 }
}
foreach ($entry in @($taskPlan.protectedFiles) + @($taskPlan.sourceFiles)) {
 if ($taskAliases.ContainsKey($entry.path)) { throw 'Protected target conflict.' }
 Assert-Record $entry
}
if ($taskAliases.Count -ne $taskPlan.originalAliasPaths) { throw 'Alias count mismatch.' }
Write-Output "PREFLIGHT_PASSED objects=20 aliases=$($taskAliases.Count)"

# 模块二：迁存结果再次核对后逐条移除原容器路径；不原位改写共享历史文件。
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
 foreach ($alias in $group.aliases) {
  $entry = @{path=$alias; bytes=$group.original.bytes; sha256=$group.original.sha256}
  Assert-Record $entry
  $path = Get-CheckedPath $alias $true
  if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB started during cleanup.' }
  Remove-Item -LiteralPath $path
  if (Test-Path -LiteralPath $path) { throw 'File remains after deletion.' }
  $marker = [IO.Path]::GetDirectoryName($alias) + '/LINE_CACHE_PRUNED_20261031.json'
  @{
   status='retired-line-matrices-pruned-science-retained'
   plan='docs/validation/storage_20261031/retired_lines/plan.json'
   original=$group.original; originalAliasPath=$alias; retained=$group.retained
   removedFieldPaths=@($group.resultScan.removed.PSObject.Properties.Name)
   retainedScientificFields=$group.retainedScientificFields
   originalContainerSHAReproductionClaimed=$false
   note='All noncache MAT elements remain exact. Rebuild historical matrices and any retired upstream cache in a new directory with accurate historical sources and new bindings. Retained results are not the original container or its SHA.'
  } | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath (Get-CheckedPath $marker $true) -Encoding UTF8
  $receipt.deleted += $entry
  $entry | ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournalPath -Encoding UTF8
  $receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceiptPath -Encoding UTF8
 }
 Write-Output "REMOVED $($group.original.path)"
}

# 模块三：保留输入与源码再次核对；删除窗口含已迁存结果，不能当成净释放量。
foreach ($entry in @($taskPlan.protectedFiles) + @($taskPlan.sourceFiles) + @($taskPlan.groups.retained)) { Assert-Record $entry }
$receipt.status = 'completed-retired-line-cache-cleanup'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfterDelete = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.deleteWindowFreeIncreaseBytes = $receipt.freeAfterDelete - $receipt.freeBeforeDelete
$receipt | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $taskReceiptPath -Encoding UTF8
Write-Output "CLEANUP_COMPLETED objects=20 paths=$($receipt.deleted.Count) netContentBytes=$($taskPlan.netContentBytesReclaimed)"
