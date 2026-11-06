#!/usr/bin/env bash
# 模块一：核对固定源码、可执行程序及动态库，不运行任何流场求解。
set -euo pipefail
rah_prefix=/opt/rocket-aeroheat/eilmer-v5.0.0
source "$rah_prefix/environment.sh"
[[ "$(git -C "$DGD_REPO" rev-parse HEAD)" == f53f4609a0331d48efee69a4e4f3c3598378cc03 ]]
git -C "$DGD_REPO" diff --exit-code
sha256sum --check "$rah_prefix/executables.sha256"
test -x "$DGD/bin/prep-gas"
for rah_bin in lmr lmr-run; do
    test -x "$DGD/bin/$rah_bin"
    ldd "$DGD/bin/$rah_bin" > "$rah_prefix/$rah_bin-libraries.txt"
    if grep -q 'not found' "$rah_prefix/$rah_bin-libraries.txt"; then exit 3; fi
done
lmr --version-long
lmr help > "$rah_prefix/lmr-help.txt"
lmr help prep-gas > "$rah_prefix/prep-gas-help.txt"
# 模块二：仅调用原生热物性库，验证气体文件生成和Lua接口是否可用。
rah_label="${1:-installation-health}"
[[ "$rah_label" =~ ^installation-health(-[0-9]+)?$ ]] || exit 2
rah_check="$rah_prefix/$rah_label"
test ! -e "$rah_check"
mkdir "$rah_check"
cd "$rah_check"
printf '%s\n' 'model = "thermally perfect gas"' 'species = {"air"}' > air.inp
prep-gas air.inp air.lua
cat > gas-health.lua <<'EOF'
-- 安装健康检查：不读取箭体，不生成流场，不授予任何物理验证资格。
local gm = GasModel:new{"air.lua"}
local gs = GasState:new{gm}
gs.massf = {air=1.0}
local file = assert(io.open("gas-health.csv", "w"))
file:write("T_K,p_Pa,rho_kg_m3,a_m_s,cp_J_kgK,mu_Pas,k_W_mK\n")
local cpFirst = nil
for _, T in ipairs({226.50908361133003, 300.0, 1000.0}) do
    gs.p = 1197.0279007383137
    gs.T = T
    gm:updateThermoFromPT(gs)
    gm:updateSoundSpeed(gs)
    gm:updateTransCoeffs(gs)
    local cp = gm:Cp(gs)
    assert(gs.rho > 0 and gs.a > 0 and cp > 0 and gs.mu > 0 and gs.k > 0)
    file:write(string.format("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
        gs.T, gs.p, gs.rho, gs.a, cp, gs.mu, gs.k))
    if cpFirst == nil then cpFirst = cp end
    if T == 1000.0 then assert(math.abs(cp - cpFirst) > 1.0) end
end
file:close()
print("Native gas preparation and thermally perfect property calls passed; no CFD was run.")
EOF
lmr custom-script --job=gas-health.lua
cat gas-health.csv
uname -sr
cat /etc/os-release
nproc
df -h / /mnt/e
printf '%s\n' 'Installation health passed. Mesh/thermal-boundary/CFD validation remains a future task.'
