# 模块一：仅执行已核对的四组旧矩阵清单；全部输入、源码和目标逐项验SHA。
param([string]$Plan = 'docs/validation/storage_20261028/followup/plan.json')
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskOutputs = (Join-Path $taskRoot 'outputs') + [IO.Path]::DirectorySeparatorChar
$taskEvidence = Join-Path $taskRoot 'docs/validation/storage_20261028/followup'
$taskReceipt = Join-Path $taskEvidence 'execution.json'
$taskJournal = Join-Path $taskEvidence 'operations.jsonl'
if ((Test-Path -LiteralPath $taskReceipt) -or (Test-Path -LiteralPath $taskJournal)) { throw '已有执行记录，不允许重复清理。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB仍在运行，不能清理计算数据。' }
$taskPlanPath = [IO.Path]::GetFullPath((Join-Path $taskRoot $Plan))
if ($taskPlanPath -ne (Join-Path $taskEvidence 'plan.json')) { throw '本次仅接受固定审核清单。' }
$taskData = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskData.kind -ne 'historical-output-cache-cleanup-v2' -or $taskData.workspace -ne $taskRoot -or @($taskData.deleteFiles).Count -ne 24) { throw '计划类型、工作区或文件数不符。' }
function Get-CheckedPath([string]$Relative, [bool]$OutputOnly) {
    if ([IO.Path]::IsPathRooted($Relative)) { throw '清单必须使用项目相对路径。' }
    $path = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
    $prefix = $taskRoot + [IO.Path]::DirectorySeparatorChar
    if ($OutputOnly) { $prefix = $taskOutputs }
    if (-not $path.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) { throw "路径越界：$Relative" }
    $item = [IO.FileInfo]::new($path)
    if ($item.Exists -and (($item.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0)) { throw '不处理重解析文件。' }
    $parent = $item.Directory
    while ($null -ne $parent -and $parent.FullName.Length -ge $taskRoot.Length) {
        if (($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw '不处理重解析目录。' }
        $parent = $parent.Parent
    }
    return $path
}
function Assert-Record($Record, [bool]$OutputOnly) {
    $path = Get-CheckedPath $Record.path $OutputOnly
    $item = Get-Item -LiteralPath $path
    if ($item.PSIsContainer -or $item.Length -ne $Record.bytes) { throw "文件大小不符：$($Record.path)" }
    if ((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $Record.sha256) { throw "文件SHA不符：$($Record.path)" }
}
$taskAllowed = @(
    'outputs/v84_validation/q24system',
    'outputs/v84_validation/full02/baseline/system2',
    'outputs/v84_validation/full02/minimum/system2',
    'outputs/v91_validation/continuous01/step01/current'
)
$taskSeen = @{}
foreach ($entry in $taskData.deleteFiles) {
    $parentRel = $entry.path.Substring(0, $entry.path.LastIndexOf('/'))
    if ($taskAllowed -notcontains $parentRel -or $entry.path -notmatch '/current_peak_(cache|operator|packed|factor|scale|metric)_evidence\.mat$') { throw '删除项不在已审核的四组派生矩阵内。' }
    if ($taskSeen.ContainsKey($entry.path)) { throw '重复删除项。' }
    $taskSeen[$entry.path] = $true
    Assert-Record $entry $true
}
foreach ($entry in $taskData.protectedFiles) {
    if ($taskSeen.ContainsKey($entry.path)) { throw '保留项与删除项冲突。' }
    Assert-Record $entry $true
}
foreach ($entry in $taskData.sourceFiles) { Assert-Record $entry $false }
Assert-Record $taskData.inventoryBefore $false
foreach ($system in $taskData.cacheSystems) {
    $markerRel = ($system.manifest.path -replace 'peak_system.json$', 'CACHE_PRUNED.json')
    if (Test-Path -LiteralPath (Get-CheckedPath $markerRel $true)) { throw '清理标记已存在。' }
}
Write-Output "PREFLIGHT_COMPLETE deletedCandidates=$($taskData.deleteFiles.Count) protected=$($taskData.protectedFiles.Count)"
# 模块二：只删除24个明确文件，不递归删除或移动任何目录。
$receipt = [ordered]@{status='executing'; startedUTC=[DateTime]::UtcNow.ToString('o');
    planSHA256=(Get-FileHash -LiteralPath $taskPlanPath).Hash.ToLowerInvariant();
    freeBefore=[IO.DriveInfo]::new($taskRoot).AvailableFreeSpace; deleted=@(); protectedVerifiedAfter=$false}
foreach ($entry in $taskData.deleteFiles) {
    Assert-Record $entry $true
    $path = Get-CheckedPath $entry.path $true
    Remove-Item -LiteralPath $path
    if (Test-Path -LiteralPath $path) { throw '删除后文件仍存在。' }
    $receipt.deleted += $entry
    @{action='deleted-obsolete-derived-cache';path=$entry.path;bytes=$entry.bytes;sha256=$entry.sha256} |
        ConvertTo-Json -Compress | Add-Content -LiteralPath $taskJournal -Encoding UTF8
}
# 模块三：原清单仅为历史记录；复算必须用对应源码在新目录生成新清单。
foreach ($system in $taskData.cacheSystems) {
    $markerRel = ($system.manifest.path -replace 'peak_system.json$', 'CACHE_PRUNED.json')
    @{status='obsolete-derived-cache-pruned-inputs-and-matching-source-retained';
      originalManifest=$system.manifest; source=$system.source; geometry=$system.geometry;
      sourceCommit=$taskData.sourceCommit; sourceModules=$system.sourceModules;
      removedArtifacts=$system.artifacts; rebuildFunction=$system.rebuildFunction;
      rebuildRequiresNewManifest=$true; historicalMATFileSHAReproductionClaimed=$false;
      rebuildExecuted=$false;
      note='旧清单不能直接恢复。按记录的源码版本及历史覆盖文件，在新目录重建并使用新清单；不保证MAT文件历史SHA复现。'} |
        ConvertTo-Json -Depth 15 | Set-Content -LiteralPath (Get-CheckedPath $markerRel $true) -Encoding UTF8
}
# 模块四：删除后再次验证全部保留文件，记录真实空闲空间变化。
foreach ($entry in $taskData.protectedFiles) { Assert-Record $entry $true }
foreach ($entry in $taskData.sourceFiles) { Assert-Record $entry $false }
$receipt.protectedVerifiedAfter = $true
$receipt.status = 'completed-sha-verified-historical-cache-cleanup'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfter = [IO.DriveInfo]::new($taskRoot).AvailableFreeSpace
$receipt.actualFreeIncreaseBytes = $receipt.freeAfter - $receipt.freeBefore
$receipt.deletedBytes = ($taskData.deleteFiles | Measure-Object -Property bytes -Sum).Sum
$receipt.noNumericalComputationExecuted = $true
$receipt | ConvertTo-Json -Depth 20 | Set-Content -LiteralPath $taskReceipt -Encoding UTF8
Write-Output "CLEANUP_COMPLETED files=$($receipt.deleted.Count) deletedBytes=$($receipt.deletedBytes) actualFreeIncreaseBytes=$($receipt.actualFreeIncreaseBytes)"
