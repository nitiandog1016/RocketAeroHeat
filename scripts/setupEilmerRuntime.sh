#!/usr/bin/env bash
# 模块一：在已经可用的Linux/WSL中建立独立固定版本运行时，不启用Windows组件。
# 需事先安装：build-essential gfortran git curl xz-utils libreadline-dev
# libncurses-dev libopenmpi-dev libplot-dev python3 python3-cffi ruby ruby-ffi。
# 先读取docs/EILMER_RUNTIME_PREPARATION_V177.md；实际运行状态以安装回执为准。
set -euo pipefail
if [[ "${1:-}" == "--help" ]]; then
    printf '%s\n' 'Usage: setupEilmerRuntime.sh NEW_ABSOLUTE_PREFIX [COMPILER_ARCHIVE]' \
      'Pins Eilmer v5.0.0 and LDC 1.41.0; no sudo, OS feature changes, or CFD.'
    exit 0
fi
if [[ $# -lt 1 || $# -gt 2 || "$1" != /* || -e "$1" ]]; then
    printf '%s\n' 'Provide a new absolute prefix; existing paths are never overwritten.' >&2
    exit 2
fi
rah_prefix="$1"
[[ "$(uname -s)" == Linux && "$(uname -m)" == x86_64 ]] || { printf '%s\n' 'Linux x86_64 required.' >&2; exit 2; }
for tool in git curl tar xz make gcc gfortran python3 timeout sha256sum mpicc; do
    command -v "$tool" >/dev/null || { printf 'Missing prerequisite: %s\n' "$tool" >&2; exit 2; }
done
mkdir -p "$rah_prefix"
exec > >(tee "$rah_prefix/setup.log") 2>&1
# 模块二：准确来源固定；不执行远程脚本，下载编译器先校验公开SHA256。
rah_commit=f53f4609a0331d48efee69a4e4f3c3598378cc03
rah_compiler=ldc2-1.41.0-linux-x86_64
rah_sha=4a439457f0fe59e69d02fd6b57549fc3c87ad0f55ad9fb9e42507b6f8e327c8f
rah_archive="${2:-$rah_prefix/$rah_compiler.tar.xz}"
if [[ $# == 1 ]]; then
    timeout 300 curl --fail --location --retry 1 --output "$rah_archive" \
      "https://github.com/ldc-developers/ldc/releases/download/v1.41.0/$rah_compiler.tar.xz"
fi
printf '%s  %s\n' "$rah_sha" "$rah_archive" | sha256sum --check --status
tar -xJf "$rah_archive" -C "$rah_prefix"
git init "$rah_prefix/source"
git -C "$rah_prefix/source" remote add origin https://github.com/gdtk-uq/gdtk.git
timeout 300 git -C "$rah_prefix/source" fetch --depth 1 origin "$rah_commit"
git -C "$rah_prefix/source" checkout --detach FETCH_HEAD
[[ "$(git -C "$rah_prefix/source" rev-parse HEAD)" == "$rah_commit" ]]
# 模块三：官方完整热物性/湍流/NK构建，二线程；不修改源码或用户shell配置。
export PATH="$rah_prefix/$rah_compiler/bin:$PATH"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
ulimit -v 8388608
ldc2 --version
timeout 1200 make -C "$rah_prefix/source/src/lmr" -j2 install \
  DMD=ldc2 FLAVOUR=fast WITH_FPE=1 INSTALL_DIR="$rah_prefix/install"
git -C "$rah_prefix/source" diff --exit-code
printf '%s\n' "$rah_commit" > "$rah_prefix/source_commit.txt"
sha256sum "$rah_prefix/install/bin/lmr" "$rah_prefix/install/bin/lmr-run" \
  "$rah_prefix/install/bin/prep-gas" > "$rah_prefix/executables.sha256"
# 模块四：独立激活脚本，只供项目调用，不修改用户shell配置。
cat > "$rah_prefix/environment.sh" <<EOF
export DGD='$rah_prefix/install'
export DGD_REPO='$rah_prefix/source'
export PATH='$rah_prefix/install/bin:$rah_prefix/$rah_compiler/bin':\$PATH
export DGD_LUA_PATH='$rah_prefix/install/lib/?.lua'
export DGD_LUA_CPATH='$rah_prefix/install/lib/?.so'
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
EOF
printf '%s\n' 'Build completed; runtime interface, mesh and gas checks are still required. No CFD was run.'
