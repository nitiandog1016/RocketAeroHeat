# 模块一：固定十八份历史结果，必须逐份通过独立原MAT完整SHA恢复审核。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskOutput = [IO.Path]::GetFullPath((Join-Path $taskRoot 'outputs'))
$taskFolder = Join-Path $taskRoot 'docs/validation/storage_20261031/result_deduplication'
$taskPlanPath = Join-Path $taskFolder 'plan.json'
$taskPlan = Get-Content -LiteralPath $taskPlanPath -Raw -Encoding UTF8 | ConvertFrom-Json
$taskAudit = Get-Content -LiteralPath (Join-Path $taskFolder 'predelete_archive_audit.json') -Raw -Encoding UTF8 | ConvertFrom-Json
if ($taskAudit.status -ne 'passed-independent-exact-original-MAT-SHA' -or $taskAudit.verifiedOriginals -ne 18 -or @($taskPlan.groups).Count -ne 18) { throw '独立恢复验证未通过。' }
if ((Get-FileHash -LiteralPath $taskPlanPath -Algorithm SHA256).Hash.ToLowerInvariant() -ne $taskAudit.planSHA256) { throw '审核后清单变化。' }
if ($taskPlan.kind -ne 'lossless-retired-MAT-shared-store-extension-20261031' -or $taskPlan.netContentBytes -le 0) { throw '范围或空间收益不符。' }
if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw 'MATLAB运行期间不维护结果。' }
if (Test-Path -LiteralPath (Join-Path $taskFolder 'execution.json')) { throw '已有执行回执，不重复删除。' }

function Get-TaskPath([string]$relative) {
    $candidate = [IO.Path]::GetFullPath((Join-Path $taskRoot $relative))
    if (-not $candidate.StartsWith($taskRoot + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) { throw '路径越界。' }
    $parent = $candidate
    while ($parent.Length -ge $taskRoot.Length) {
        if (Test-Path -LiteralPath $parent) {
            if ((Get-Item -LiteralPath $parent -Force).Attributes -band [IO.FileAttributes]::ReparsePoint) { throw '拒绝重解析点。' }
        }
        $parent = Split-Path -Path $parent -Parent
    }
    return $candidate
}
function Assert-TaskRecord($record) {
    $file = Get-TaskPath $record.path
    if ((Get-Item -LiteralPath $file).Length -ne $record.bytes) { throw "长度变化：$file" }
    if ((Get-FileHash -LiteralPath $file -Algorithm SHA256).Hash.ToLowerInvariant() -ne $record.sha256) { throw "SHA变化：$file" }
}
function Assert-TaskOriginal($group) {
    $file = Get-TaskPath $group.original.path
    $allowed = @('retained_operator_results_20261030', 'retained_graph_results_20261031') | ForEach-Object {
        [IO.Path]::GetFullPath((Join-Path $taskOutput $_)) + [IO.Path]::DirectorySeparatorChar
    }
    if (-not ($allowed | Where-Object { $file.StartsWith($_, [StringComparison]::OrdinalIgnoreCase) })) { throw '不属于本轮退役结果范围。' }
    $identity = $group.originalIdentity
    if ($identity.sha256 -ne $group.original.sha256 -or $identity.bytes -ne $group.original.bytes) { throw '归档身份不符。' }
    if (@($taskAudit.results | Where-Object { $_.original -eq $group.original.path -and $_.sha256 -eq $group.original.sha256 -and $_.bytes -eq $group.original.bytes }).Count -ne 1) { throw '独立审核绑定不符。' }
    Assert-TaskRecord $group.original
    $links = @(& fsutil.exe hardlink list $file 2>&1)
    if ($LASTEXITCODE -ne 0 -or @($links | Where-Object { "$_".Trim().StartsWith('\') }).Count -ne 1) { throw '硬链接范围变化。' }
    if (Test-Path -LiteralPath (Join-Path (Split-Path -Path $file -Parent) 'EVIDENCE_ARCHIVED_20261031.json')) { throw '恢复标记已存在。' }
}

# 模块二：全部共享对象、保留依赖、工具及十八份原结果通过后才开始删除。
foreach ($record in @($taskPlan.protectedFiles) + @($taskPlan.sourceFiles) + @($taskPlan.historicalPlans) + @($taskPlan.archiveFiles)) { Assert-TaskRecord $record }
foreach ($group in $taskPlan.groups) {
    Assert-TaskRecord $group.manifest
    Assert-TaskOriginal $group
}

# 模块三：先写恢复入口，再逐文件删除；不递归删除，不改写任何共享数据对象。
$taskFreeBefore = (Get-PSDrive -Name ([IO.Path]::GetPathRoot($taskRoot).Substring(0,1))).Free
$taskRemoved = @()
foreach ($group in $taskPlan.groups) {
    if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) { throw '计算进程已启动，停止维护。' }
    Assert-TaskOriginal $group
    $file = Get-TaskPath $group.original.path
    $marker = Join-Path (Split-Path -Path $file -Parent) 'EVIDENCE_ARCHIVED_20261031.json'
    $entry = [ordered]@{status='losslessly-archived-exact-original-SHA-verified'; original=$group.original.path; originalSHA256=$group.original.sha256; manifest=$group.manifest; store=$taskPlan.store; allScientificDataRetained=$true; directlyLoadableMATAtOldPath=$false; restoreScript='scripts/archiveMat5History.py'; note='全部结果原字节去重保留；旧MAT先恢复到新目录；对象库仅在本地，GitHub没有备份大型数据。'}
    [IO.File]::WriteAllText($marker, (($entry | ConvertTo-Json -Depth 12) + "`n"), [Text.UTF8Encoding]::new($false))
    Remove-Item -LiteralPath $file -Force
    $taskRemoved += $group.original.path
    Add-Content -LiteralPath (Join-Path $taskFolder 'operations.jsonl') -Value (@{path=$group.original.path; sha256=$group.original.sha256; removedUTC=[DateTime]::UtcNow.ToString('o')} | ConvertTo-Json -Compress) -Encoding UTF8
}
$taskFreeAfter = (Get-PSDrive -Name ([IO.Path]::GetPathRoot($taskRoot).Substring(0,1))).Free
$taskReceipt = [ordered]@{status='executed'; originalsRemoved=$taskRemoved; freeBytesBeforeDeletion=$taskFreeBefore; freeBytesAfterDeletion=$taskFreeAfter; deletionWindowFreeIncrease=$taskFreeAfter-$taskFreeBefore; freeBytesBeforeArchive=$taskPlan.freeBytesBeforeArchive; freeIncreaseSinceArchiveStart=$taskFreeAfter-$taskPlan.freeBytesBeforeArchive; netContentBytesReclaimed=$taskPlan.netContentBytes; note='删除窗口含先前另存归档，不能当作净回收；保留全部数值和当前依赖。'}
[IO.File]::WriteAllText((Join-Path $taskFolder 'execution.json'), (($taskReceipt | ConvertTo-Json -Depth 12) + "`n"), [Text.UTF8Encoding]::new($false))
$taskReceipt | ConvertTo-Json -Depth 12
