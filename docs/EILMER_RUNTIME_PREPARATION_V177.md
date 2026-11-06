# Eilmer固定运行时（v177安装完成）

**2026-11-07：WSL3.0.3、Ubuntu24.04.5、LDC1.41.0与固定Eilmer5.0已安装；原生启动、气体文件生成和三温度热物性调用通过。无需重启。** 本次只完成运行环境，不代表箭体、混合网格、热边界或公开基准已验证。详见[安装回执](validation/v177/ENVIRONMENT_SETUP.md)。按用户安排，完成后等待其继续目标指令。

为真实20°完整域准备一个成熟替代后端：GDTK/Eilmer v5.0.0，固定提交`f53f4609a0331d48efee69a4e4f3c3598378cc03`。它原生提供热完美气体、三维混合网格、热壁与Wilcox k–ω；该湍流模型不是现有SST，结果不可同算子接续，也不继承SST资格。SU2 NEMO在准确源码中拒绝任何湍流模型，已排除。

安装前本机Windows的WSL不可用。按[官方安装前提](https://gdtk.uqcloud.net/docs/getting-started/prerequisites/)，Windows部署使用WSL2；固定版本makefile只有Linux/macOS路径，不开展自制Windows移植。现已通过独立WSL环境完成构建与安装健康检查，实际完整域入口仍待后续验证。

## 可审阅的环境安排

1. 2026-11-07用户已明确授权下载、安装及必要重启；此前等待确认/不自动重启的安排由本次授权取代。启用Windows的WSL2组件；需要重启时先保存状态与接续说明。Windows管理员确认仍由用户完成，不自动操作安全确认窗口。
2. 安装独立Ubuntu 24.04环境，发行版数据优先放在E盘的`E:\RocketAeroHeat\tools\runtime\wsl`。启用后按实际WSL版本核对`--location`支持；不隐式退回C盘大安装。
3. 在该Linux环境中安装官方列出的编译依赖，使用新的`/opt/rocket-aeroheat/eilmer-v5.0.0`（或明确的用户目录）作为独立前缀。虚拟磁盘位于E盘，源码/构建/依赖仍需单列空间；预留12GiB为规划值，未实测安装量。
4. 执行项目内`scripts/setupEilmerRuntime.sh`。固定LDC1.41.0 Linux x86_64编译器压缩包89,108,296字节，SHA256`4a439457f0fe59e69d02fd6b57549fc3c87ad0f55ad9fb9e42507b6f8e327c8f`；下载后校验，固定GDTK提交，二线程/单次构建1200秒/8GiB地址空间上限。总安装预算需在实际启用前另登记，下载及系统依赖不假称在构建1200秒内。
5. 完成后先核准确来源、气体模型、原生混合单元与边界/热流小合同，再登记完整域候选；不直接启动长算。MATLAB与工作台继续自动组织调用，日常使用者不需要导入外部CFD。

Linux编译脚本不调用sudo、不修改Windows组件、不重启、不改用户shell配置；拒绝覆盖已有前缀。系统安装和依赖安装由分开的具名脚本完成。构建与健康检查已实际完成；混合网格读取、低温物性适用域、热壁接口与三维RANS热验证仍待完成。

## 2026-11-07环境安装阶段（用户已授权）

本阶段仅完成运行环境和安装健康检查，完成后等待用户指示再继续科研目标。安装包、发行版虚拟盘、WSL交换文件和可移动运行时均放在`E:\RocketAeroHeat\tools\runtime`；Windows固定系统组件保留系统位置。官方安装包先核对发布哈希，WSL安装包另检查Microsoft签名。

下载/系统安装/依赖安装/编译/检查合计预留90分钟（不含重启后用户重新登录的等待），单次编译仍为1200秒、二线程、8GiB地址空间；下载单文件最多1200秒。新增磁盘规划约12GiB，32GiB时停止继续扩展并记录原因。该环境预算独立于原v177小源码/CFD预算，不启动整箭或公开基准CFD，不授予新物理资格。

持久状态与安装日志保存在`tools/runtime/setup-v177`，接续入口记录在`docs/validation/v177/ENVIRONMENT_SETUP.md`。实际版本、路径、退出状态和是否需要重启随执行记录；尚未验证的阶段不能写成已就绪。

来源：[固定发布](https://github.com/gdtk-uq/gdtk/tree/f53f4609a0331d48efee69a4e4f3c3598378cc03)、[官方功能](https://gdtk.uqcloud.net/docs/eilmer/releases/)、[准确makefile](https://github.com/gdtk-uq/gdtk/blob/f53f4609a0331d48efee69a4e4f3c3598378cc03/src/lmr/makefile)、[固定编译器](https://github.com/ldc-developers/ldc/releases/tag/v1.41.0)。许可随对应官方源码保存；公开气动热控制另外记录实验来源与模型差异。
