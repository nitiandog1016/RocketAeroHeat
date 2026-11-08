[CmdletBinding()]
param(
    [string]$MatlabPath = 'D:\matlab\bin\matlab.exe',
    [string]$DesktopPath = [Environment]::GetFolderPath('Desktop')
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$iconPath = Join-Path $projectRoot 'assets\icons\rocket-aero-heat.ico'
$entryPath = Join-Path $projectRoot 'launchRocketAeroHeat.m'
foreach ($requiredPath in @($MatlabPath, $iconPath, $entryPath)) {
    if (-not (Test-Path -LiteralPath $requiredPath -PathType Leaf)) {
        throw "Required file not found: $requiredPath"
    }
}
if (-not (Test-Path -LiteralPath $DesktopPath -PathType Container)) {
    throw "Desktop folder not found: $DesktopPath"
}
$MatlabPath = (Resolve-Path -LiteralPath $MatlabPath).Path
$shortcutPath = Join-Path $DesktopPath 'RocketAeroHeat.lnk'
$shell = New-Object -ComObject WScript.Shell
try {
    # Preserve a pre-existing shortcut instead of silently replacing user changes.
    if (Test-Path -LiteralPath $shortcutPath) {
        $backupPath = "$shortcutPath.backup-$(Get-Date -Format 'yyyyMMdd-HHmmss-fff')"
        Copy-Item -LiteralPath $shortcutPath -Destination $backupPath -ErrorAction Stop
    }
    $shortcut = $shell.CreateShortcut($shortcutPath)
    $shortcut.TargetPath = $MatlabPath
    $shortcut.WorkingDirectory = $projectRoot
    # Interactive MATLAB session: the workbench stays usable, MATLAB starts minimized.
    # Only this newly launched session exits when the workbench window closes.
    $startup = "try,launchRocketAeroHeat;catch err,disp(getReport(err,'extended'));uiwait(errordlg(err.message,'RocketAeroHeat'));exit(1);end;exit(0);"
    $shortcut.Arguments = '-nosplash -automation -sd "' + $projectRoot + '" -r "' + $startup + '"'
    $shortcut.IconLocation = "$iconPath,0"
    $shortcut.Description = 'RocketAeroHeat research workbench'
    $shortcut.WindowStyle = 1
    $shortcut.Save()
    Write-Output $shortcutPath
}
finally {
    if ($null -ne $shortcut) { [void][Runtime.InteropServices.Marshal]::ReleaseComObject($shortcut) }
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($shell)
}
