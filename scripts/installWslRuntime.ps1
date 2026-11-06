# 模块一：仅安装已校验的WSL；数据目录位于E盘，不启动科研计算或自动重启。
$ErrorActionPreference = 'Stop'
$rahRoot = 'E:\RocketAeroHeat\tools\runtime'
$rahSetup = Join-Path $rahRoot 'setup-v177'
$rahMsi = Join-Path $rahSetup 'downloads\wsl.3.0.3.0.x64.msi'
$rahInstall = Join-Path $rahRoot 'windows-wsl'
$rahResult = [ordered]@{ phase='starting'; started=(Get-Date).ToString('o'); restartNeeded=$false; installLocation=$rahInstall }
$rahResultPath = Join-Path $rahSetup 'windows-install-result.json'
function Save-RahResult { $rahResult | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $rahResultPath -Encoding UTF8 }
Save-RahResult
Start-Transcript -Path (Join-Path $rahSetup 'windows-install-transcript.txt') -NoClobber
try {
    $rahIdentity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $rahPrincipal = New-Object Security.Principal.WindowsPrincipal($rahIdentity)
    if (-not $rahPrincipal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) { throw 'Administrator consent is required.' }
    if ((Get-FileHash -LiteralPath $rahMsi -Algorithm SHA256).Hash.ToLower() -ne 'ae977286974255bd26bf925f00f40f5692c6b2084155cbc0563c09aea091f211') { throw 'MSI SHA256 mismatch.' }
    $rahSignature = Get-AuthenticodeSignature -LiteralPath $rahMsi
    if ($rahSignature.Status -ne 'Valid' -or $rahSignature.SignerCertificate.Subject -notmatch 'O=Microsoft Corporation') { throw 'Microsoft installer signature validation failed.' }
    $rahResult.signature = $rahSignature.SignerCertificate.Subject
    # 模块二：WSL2需要虚拟机平台；不为WSL1额外启用旧组件。
    $rahFeature = Get-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform
    $rahResult.originalFeatureState = [string]$rahFeature.State
    if ($rahFeature.State -ne 'Enabled') {
        $rahEnabled = Enable-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform -All -NoRestart
        $rahResult.restartNeeded = [bool]$rahEnabled.RestartNeeded
    }
    if (Test-Path -LiteralPath $rahInstall) { throw 'Refusing to overwrite an existing WSL program directory.' }
    New-Item -ItemType Directory -Path $rahInstall | Out-Null
    # 新建服务程序目录不继承项目的可写权限，管理员/系统可写，普通用户只读执行。
    & icacls.exe $rahInstall /inheritance:r /grant:r '*S-1-5-18:(OI)(CI)F' '*S-1-5-32-544:(OI)(CI)F' '*S-1-5-32-545:(OI)(CI)RX'
    if ($LASTEXITCODE -ne 0) { throw 'Unable to set service-program directory permissions.' }
    $rahResult.phase = 'installing-msi'; Save-RahResult
    # 模块三：使用官方安装目录参数，把可移动WSL程序也安装到E盘。
    $rahArgs = @('/i', ('"' + $rahMsi + '"'), '/qn', '/norestart', 'REBOOT=ReallySuppress',
        ('INSTALLDIR="' + $rahInstall + '"'), '/L*v', ('"' + (Join-Path $rahSetup 'windows-msi.log') + '"'))
    $rahProcess = Start-Process -FilePath 'msiexec.exe' -ArgumentList $rahArgs -PassThru -Wait -WindowStyle Hidden
    $rahResult.msiExitCode = $rahProcess.ExitCode
    if ($rahProcess.ExitCode -notin @(0,3010)) { throw "MSI failed: $($rahProcess.ExitCode)" }
    if ($rahProcess.ExitCode -eq 3010) { $rahResult.restartNeeded = $true }
    $rahRegistered = Get-ItemProperty -LiteralPath 'HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Lxss\MSI'
    $rahResult.registeredLocation = $rahRegistered.InstallLocation
    if ($rahRegistered.InstallLocation.TrimEnd('\') -ne $rahInstall) { throw 'WSL did not use the required E drive location.' }
    $rahResult.phase = 'windows-installed'
    $rahResult.finished = (Get-Date).ToString('o')
    Save-RahResult
} catch {
    $rahResult.phase = 'failed'
    $rahResult.error = $_.Exception.Message
    $rahResult.finished = (Get-Date).ToString('o')
    Save-RahResult
    throw
} finally {
    Stop-Transcript
}
