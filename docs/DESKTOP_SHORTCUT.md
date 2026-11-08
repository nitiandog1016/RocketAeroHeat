# 桌面启动

双击桌面的 **RocketAeroHeat** 快捷方式，打开现有五步骤科研工作台。需要本机已安装且可正常使用的 MATLAB；此入口不是独立编译版。

快捷方式启动独立的交互式 MATLAB 会话，MATLAB 主界面最小化，工作台正常显示。关闭工作台后，该会话退出；工作台原有取消计算与关闭回调保持。重复双击会打开多个独立会话。

启动入口 `launchRocketAeroHeat.m` 调用原 `main.m`，没有自动导入模型或启动求解。启动日志和窗口打开记录保存在 `outputs/desktop_launcher/`；记录只说明界面启动状态，不代表求解或物理验证。

原始图标由本项目绘制：深蓝底、白色火箭、青色舷窗与气流、橙色尾焰。`assets/icons/` 同时提供可编辑 SVG、窗口 PNG 及含 16—256 像素尺寸的 Windows ICO。

重新创建快捷方式（PowerShell）：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File E:\RocketAeroHeat\scripts\installDesktopShortcut.ps1
```

MATLAB 安装在其他位置时传入 `-MatlabPath '实际路径\matlab.exe'`。安装脚本从自身位置解析工程目录，并使用 Windows 当前用户桌面位置；若已有同名快捷方式，会先保留备份。移动工程后重新执行安装脚本。

图标可用 `python scripts/createDesktopIcon.py` 重建，需要 matplotlib 和 Pillow；平时启动程序不依赖 Python。

2026-11-08 已从实际桌面 `.lnk` 启动 MATLAB R2023a：工作台可见，5 个流程按钮，当前第 1 步，未加载模型、未执行计算。ICO 的 7 个尺寸及 SVG 格式检查通过。此次只核对启动，不代表工作台全部功能的回归测试；原科学计算源码未修改。
