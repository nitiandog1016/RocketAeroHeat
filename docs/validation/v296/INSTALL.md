# 独立安装与本机缓存

本次复算环境为 64 位 Windows、MATLAB R2023a、Python 3.7.7、NumPy 1.21.6、SciPy 1.7.3、psutil 5.9.8 和官方 SU2 8.6.0。MATLAB 通过 `python` 命令调用这些库。SciPy 参与壁面映射。安装在支持硬链接的磁盘上，体网、缓存和任务目录位于同一项目内。

在项目根目录运行 `python scripts/installSU2Runtime.py`，可安装或核对官方求解器。迁移已有安装时，该命令核验文件并更新本机可执行路径。随后准备下述热流观察器，在 MATLAB 运行 `rocketAeroHeatWorkbench`，选择 **选用 B 型热敏珠 → 6 m/s · 360° 全周 · 56 层 → 准备 56 层全周体网**。确认边界后选择 **从构造初场计算两个半步**，到达指定时刻后结果页自动打开。

**网格步骤**按名义尺寸生成 48 层体网，细分最内八个区间，生成 56 层网和节点映射。结果缓存在项目的 `outputs/thermistor_inputs/B6-56`。此时只产生几何；侧栏实际开始计算后才准备数值初场。

**计算步骤**按实际指纹读取或恢复输入。初场来自程序内部的冷态前驱、同速扇区基流和全周松弛。已安装的观察器有固定的二进制指纹。首次构建使用指定 SU2 构建环境的源码、编译器、两个静态库及相应对象，沿用原观察位置。标准工具目录为 `tools/runtime/su2-thermal-build`。若这套构建工具已装在另一处，可将环境变量 `RAH_THERMISTOR_SU2_BUILD` 设为含 `source`、`toolchain`、`build-v117` 的实际目录。仅安装官方 SU2 可执行文件的设备还需要这套工具，或安装已校验的专用观察器缓存。[观察器源码与许可](../../../cases/thermistor_b/startup6/README.md)

升级旧安装时，程序可一次性迁移指纹相同的本项目旧初场和观察器。换机器时可选择复用同一份内部初场来节省准备时间。如下命令只接受 SHA-256 为 `ea94eb8a3b668a47e1d3e66e59065e704ab05c57ae752fcd6ab75fa112f42445` 的 48 层启动检查点，传入文件保持原样。工作台仍会对该初场延拓，实际推进两个物理半步并复核每个原面的壁压和热流。

```powershell
python scripts/thermistorStartupCache.py --install-seed E:\已有安装\outputs\thermistor_inputs\B6-56\initial\seed_00000.csv
```

输入缓存保留实际源文件指纹。不同内容的旧态或观察器在使用前会被拒绝。准备失败时，可以按异常消息中的工作目录核对终态、日志和输入。

已有另一套经过验收的热敏珠观察器时，可用下面的安装命令省去本机的编译环境。该命令校验固定二进制和随项目发布的对应源码，在新安装中写入独立缓存。

```powershell
python scripts/installThermistorObserver.py E:\已有安装\outputs\thermistor_inputs\B6-56\observer\observe_bead.exe
```
