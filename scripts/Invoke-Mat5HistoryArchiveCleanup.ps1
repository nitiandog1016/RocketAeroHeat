# 模块一：仅处理已通过独立原文件SHA恢复验证的二十份旧结果。
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$taskOutput = [IO.Path]::GetFullPath((Join-Path $taskRoot 'outputs'))
$taskFolder = Join-Path $taskRoot 'docs/validation/storage_20261031/output_review'
$taskPlan = Get-Content -LiteralPath (Join-Path $taskFolder 'plan.json') -Raw | ConvertFrom-Json
$taskAudit = Get-Content -LiteralPath (Join-Path $taskFolder 'predelete_archive_audit.json') -Raw | ConvertFrom-Json
if ($taskAudit.status -ne 'passed-independent-exact-original-MAT-SHA' -or $taskAudit.verifiedOriginals -ne 20) { throw '独立恢复验证未通过。' }
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

# 模块二：所有对象、当前依赖及构造工具先核对，任何失败都不开始删除。
foreach ($record in @($taskPlan.protectedFiles) + @($taskPlan.sourceFiles) + @($taskPlan.archiveFiles)) { Assert-TaskRecord $record }
foreach ($group in $taskPlan.groups) {
    Assert-TaskRecord $group.original
    Assert-TaskRecord $group.manifest
    $file = Get-TaskPath $group.original.path
    $allowed = Join-Path $taskOutput 'retained_line_results_20261031'
    if (-not $file.StartsWith($allowed + [IO.Path]::DirectorySeparatorChar, [StringComparison]::OrdinalIgnoreCase)) { throw '不属于本轮旧结果范围。' }
    if (@($taskAudit.results | Where-Object {$_.original -eq $group.original.path -and $_.sha256 -eq $group.original.sha256}).Count -ne 1) { throw '独立审计绑定不符。' }
    $links = @(& fsutil.exe hardlink list $file 2>&1)
    if ($LASTEXITCODE -ne 0 -or @($links | Where-Object {"$_".Trim().StartsWith('\')}).Count -ne 1) { throw '原对象硬链接范围变化。' }
}

# 模块三：逐个写清恢复入口，再用原生PowerShell删除精确文件；不递归删除目录。
$taskFreeBefore = (Get-PSDrive -Name ([IO.Path]::GetPathRoot($taskRoot).Substring(0,1))).Free
$taskRemoved = @()
foreach ($group in $taskPlan.groups) {
    Assert-TaskRecord $group.original
    $file = Get-TaskPath $group.original.path
    $marker = Join-Path (Split-Path -Path $file -Parent) 'EVIDENCE_ARCHIVED_20261031.json'
    if (Test-Path -LiteralPath $marker) { throw '恢复标记已存在。' }
    $entry = [ordered]@{status='losslessly-archived-exact-original-SHA-verified'; original=$group.original.path; originalSHA256=$group.original.sha256; manifest=$group.manifest; store=$taskPlan.store; allScientificDataRetained=$true; directlyLoadableMATAtOldPath=$false; restoreScript='scripts/archiveMat5History.py'; note='完整数值按原字节去重保留，旧MAT路径须先恢复到新目录；对象库仍仅存本地，未由GitHub备份。'}
    [IO.File]::WriteAllText($marker, (($entry | ConvertTo-Json -Depth 12) + "`n"), [Text.UTF8Encoding]::new($false))
    Remove-Item -LiteralPath $file -Force
    $taskRemoved += $group.original.path
    Add-Content -LiteralPath (Join-Path $taskFolder 'operations.jsonl') -Value (@{path=$group.original.path; sha256=$group.original.sha256; removedUTC=[DateTime]::UtcNow.ToString('o')} | ConvertTo-Json -Compress) -Encoding utf8
}
$taskFreeAfter = (Get-PSDrive -Name ([IO.Path]::GetPathRoot($taskRoot).Substring(0,1))).Free
$taskReceipt = [ordered]@{status='executed'; originalsRemoved=$taskRemoved; freeBytesBeforeDeletion=$taskFreeBefore; freeBytesAfterDeletion=$taskFreeAfter; deletionWindowFreeIncrease=$taskFreeAfter-$taskFreeBefore; freeBytesBeforeArchive=$taskPlan.freeBytesBeforeArchive; freeIncreaseSinceArchiveStart=$taskFreeAfter-$taskPlan.freeBytesBeforeArchive; netContentBytesReclaimed=$taskPlan.netContentBytes; note='删除窗口包含先前另存的归档，不得当作净回收；全部数值保留，当前文件未删。'}
[IO.File]::WriteAllText((Join-Path $taskFolder 'execution.json'), (($taskReceipt | ConvertTo-Json -Depth 12) + "`n"), [Text.UTF8Encoding]::new($false))
$taskReceipt | ConvertTo-Json -Depth 12
