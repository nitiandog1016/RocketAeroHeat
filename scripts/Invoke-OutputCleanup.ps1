# 模块一：只接受已逐项核对的本次计划，先验证所有路径和保留依赖，再执行文件操作。
param([string]$Plan = 'docs/validation/storage_20261028/outputs_cleanup_plan.json')
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskOutput = [IO.Path]::GetFullPath((Join-Path $taskRoot 'outputs')) + [IO.Path]::DirectorySeparatorChar
$receiptPath = Join-Path $taskRoot 'docs/validation/storage_20261028/outputs_cleanup_execution.json'
if (Test-Path -LiteralPath $receiptPath) { throw '执行回执已存在，不允许重复执行。' }
$planData = Get-Content -LiteralPath (Join-Path $taskRoot $Plan) -Raw -Encoding UTF8 | ConvertFrom-Json
if ($planData.kind -ne 'explicit-output-cleanup-v1' -or $planData.workspace -ne $taskRoot -or $planData.scope -ne 'outputs') { throw '清理计划的工作区或范围不符。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB 正在运行，暂不清理计算文件。' }
function Get-OutputPath([string]$Relative) {
    $resolved = [IO.Path]::GetFullPath((Join-Path $taskRoot $Relative))
    if (-not $resolved.StartsWith($taskOutput, [StringComparison]::OrdinalIgnoreCase)) { throw "超出 outputs 范围：$Relative" }
    $parent = [IO.DirectoryInfo]::new([IO.Path]::GetDirectoryName($resolved))
    while ($null -ne $parent -and $parent.FullName.Length -ge ($taskOutput.Length - 1)) {
        if (($parent.Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw "路径含重解析目录：$Relative" }
        $parent = $parent.Parent
    }
    if ((Test-Path -LiteralPath $resolved) -and (((Get-Item -LiteralPath $resolved).Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0)) { throw "不处理重解析文件：$Relative" }
    return $resolved
}
function Assert-File([string]$Relative, [long]$Bytes, [string]$SHA) {
    $path = Get-OutputPath $Relative
    $item = Get-Item -LiteralPath $path
    if ($item.PSIsContainer -or $item.Length -ne $Bytes) { throw "文件大小不符：$Relative" }
    if ((Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() -ne $SHA) { throw "文件 SHA 不符：$Relative" }
}
$checked = @{}
foreach ($item in @($planData.protectedFiles) + @($planData.deleteFiles)) {
    Assert-File $item.path $item.bytes $item.sha256
    $checked[$item.path] = $item.sha256
}
foreach ($system in $planData.cacheSystems) {
    foreach ($module in $system.sourceModules) {
        $source = [IO.Path]::GetFullPath((Join-Path $taskRoot $module.path))
        if (-not $source.StartsWith((Join-Path $taskRoot 'internalFlow') + '\', [StringComparison]::OrdinalIgnoreCase)) { throw '源码路径越界。' }
        if ((Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash.ToLowerInvariant() -ne $module.sha256) { throw "构造模块改变：$($module.path)" }
    }
}
foreach ($group in $planData.deduplicatePhysicalFiles) {
    foreach ($relative in @($group.canonical) + @($group.paths)) {
        if (-not $checked.ContainsKey($relative)) { Assert-File $relative $group.bytes $group.sha256; $checked[$relative] = $group.sha256 }
        if ($checked[$relative] -ne $group.sha256) { throw '计划中重复路径 SHA 冲突。' }
    }
}
Write-Output "PREFLIGHT_COMPLETE verifiedPaths=$($checked.Count)"
# 模块二：对旧重复数据建立硬链接，保留所有路径和内容。单文件交换失败时恢复原文件。
$receipt = [ordered]@{status='executing'; startedUTC=[DateTime]::UtcNow.ToString('o'); planSHA256=(Get-FileHash -LiteralPath (Join-Path $taskRoot $Plan)).Hash.ToLowerInvariant(); freeBefore=[IO.DriveInfo]::new('E:\').AvailableFreeSpace; deleted=@(); deduplicated=@(); protectedVerifiedAfter=$false}
$journal = Join-Path $taskRoot 'outputs/outputs_cleanup_operations_20261028.jsonl'
if (Test-Path -LiteralPath $journal) { throw '操作日志已存在，请先核对是否曾执行。' }
$mergeCount = 0
foreach ($group in $planData.deduplicatePhysicalFiles) {
    $canonical = Get-OutputPath $group.canonical
    foreach ($relative in $group.paths) {
        $target = Get-OutputPath $relative
        $temporary = Get-OutputPath ($relative + '.rah-dedupe-new')
        $backup = Get-OutputPath ($relative + '.rah-dedupe-original')
        if ((Test-Path -LiteralPath $temporary) -or (Test-Path -LiteralPath $backup)) { throw '文件交换临时名已存在。' }
        $newReady = $false; $originalMoved = $false; $replaced = $false
        try {
            $null = New-Item -ItemType HardLink -Path $temporary -Value $canonical
            $newReady = $true
            Move-Item -LiteralPath $target -Destination $backup
            $originalMoved = $true
            Move-Item -LiteralPath $temporary -Destination $target
            $newReady = $false; $replaced = $true
            Assert-File $relative $group.bytes $group.sha256
            Remove-Item -LiteralPath $backup
            $originalMoved = $false
        } catch {
            if ($replaced -and $originalMoved) { Remove-Item -LiteralPath $target }
            if ($originalMoved) { Move-Item -LiteralPath $backup -Destination $target }
            if ($newReady) { Remove-Item -LiteralPath $temporary }
            throw
        }
        @{action='identical-content-hardlink';path=$relative;canonical=$group.canonical;sha256=$group.sha256} | ConvertTo-Json -Compress | Add-Content -LiteralPath $journal -Encoding UTF8
    }
    $receipt.deduplicated += $group
    $mergeCount++
    if (($mergeCount % 100) -eq 0) { Write-Output "DEDUPLICATED physicalFiles=$mergeCount" }
}
# 模块三：仅删除清单内的具体旧缓存/错误导出文件，不递归删除任何目录。
foreach ($item in $planData.deleteFiles) {
    Assert-File $item.path $item.bytes $item.sha256
    Remove-Item -LiteralPath (Get-OutputPath $item.path)
    if (Test-Path -LiteralPath (Get-OutputPath $item.path)) { throw '删除后文件仍存在。' }
    $receipt.deleted += $item
    @{action='deleted';path=$item.path;bytes=$item.bytes;sha256=$item.sha256} | ConvertTo-Json -Compress | Add-Content -LiteralPath $journal -Encoding UTF8
}
foreach ($system in $planData.cacheSystems) {
    $marker = Get-OutputPath (($system.manifest -replace 'peak_system.json$', '') + 'CACHE_PRUNED.json')
    if (Test-Path -LiteralPath $marker) { throw '旧缓存清理标记已存在。' }
    @{status='obsolete-derived-cache-pruned-inputs-retained'; source=$system.source; geometry=$system.geometry;
      originalManifest=$system.manifest; removedArtifacts=$system.artifacts;
      rebuildFunction=$system.rebuildFunction; rebuildRequiresNewManifest=$true;
      historicalMATFileSHAReproductionClaimed=$false;
      note='旧清单仅作历史记录，不能直接恢复；如需复算，请在新目录重建并使用新的准备清单。'} |
        ConvertTo-Json -Depth 15 | Set-Content -LiteralPath $marker -Encoding UTF8
}
# 模块四：复核当前检查点、原几何、当前方向与旧缓存输入，记录实际磁盘变化。
foreach ($item in $planData.protectedFiles) { Assert-File $item.path $item.bytes $item.sha256 }
$receipt.protectedVerifiedAfter = $true
$receipt.status = 'completed-sha-verified-output-cleanup'
$receipt.completedUTC = [DateTime]::UtcNow.ToString('o')
$receipt.freeAfter = [IO.DriveInfo]::new('E:\').AvailableFreeSpace
$receipt.actualFreeIncrease = $receipt.freeAfter - $receipt.freeBefore
$receipt.deletedFileBytes = ($planData.deleteFiles | Measure-Object -Property bytes -Sum).Sum
$receipt.deduplicatedPhysicalFileBytes = ($planData.deduplicatePhysicalFiles | Measure-Object -Property bytes -Sum).Sum
$receipt | ConvertTo-Json -Depth 30 | Set-Content -LiteralPath $receiptPath -Encoding UTF8
Write-Output "CLEANUP_COMPLETED deletedFiles=$($receipt.deleted.Count) deduplicatedPhysicalFiles=$mergeCount freeIncreaseBytes=$($receipt.actualFreeIncrease)"
