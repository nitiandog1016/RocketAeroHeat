#!/usr/bin/env bash
# 模块一：仅供新建的RocketAeroHeat Ubuntu环境使用，安装官方构建依赖。
set -euo pipefail
[[ "$(id -u)" == 0 ]] || { echo 'Run as root in the dedicated WSL distribution.' >&2; exit 2; }
[[ "$(. /etc/os-release; printf '%s' "$ID:$VERSION_ID")" == ubuntu:24.04 ]] || exit 2
export DEBIAN_FRONTEND=noninteractive
timeout 600 apt-get update -o APT::Update::Error-Mode=any
timeout 1200 apt-get install --no-install-recommends -y \
  build-essential gfortran git curl xz-utils ca-certificates \
  libreadline-dev libncurses-dev libopenmpi-dev libplot-dev \
  python3 python3-cffi python3-numpy ruby ruby-ffi
# 模块二：建立专用普通用户，编译不以root运行；不设置或索取Windows密码。
if ! id rah >/dev/null 2>&1; then
    useradd --create-home --user-group --shell /bin/bash rah
fi
install -d -o rah -g rah /opt/rocket-aeroheat
cat > /etc/wsl.conf <<'EOF'
[boot]
systemd=false
[user]
default=rah
EOF
# 模块三：保存系统包和编译入口版本，后续仅调用项目隔离运行时。
dpkg-query -W > /opt/rocket-aeroheat/debian-packages.txt
gcc --version | head -n 1
gfortran --version | head -n 1
mpicc --showme:version
python3 --version
printf '%s\n' 'Linux dependencies installed. No CFD was run.'
