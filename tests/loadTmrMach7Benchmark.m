function benchmark = loadTmrMach7Benchmark(dataRoot)
% 模块用途：读取近Mach7、20度轴对称外扩面的公开实验合同并核对原始字节。
% 本模块只提供测点、单位和边界身份，不求解流场，不代表计算对实验验证。
% 输入 dataRoot：含contract.json和source_manifest.json的目录，可省略。
% 输出保留原始重复行；压力比较只能使用independentComparisonPoint为真的行。

%% 模块一：定位单一基准合同，拒绝混入已经完成物理验证的虚假状态
if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(fileparts(mfilename('fullpath')), 'data', 'tmr_mach7');
end
if ~(ischar(dataRoot) || (isstring(dataRoot) && isscalar(dataRoot)))
    error('RocketAeroHeat:TmrContract:Input', '数据目录必须为文本路径。');
end
dataRoot = char(dataRoot);
contract = jsondecode(fileread(fullfile(dataRoot, 'contract.json')));
if ~strcmp(contract.schema, 'RocketAeroHeat.tmr-mach7-acc.v1') ...
        || contract.scope.solver_validation_complete ...
        || contract.scope.new_physical_area_m2 ~= 0
    error('RocketAeroHeat:TmrContract:Scope', '合同版本或证据身份不符合本数据切片。');
end
if ~strcmp(contract.case_id, 'NASA-TMR-ACC-SBLI-M7-20deg') ...
        || ~strcmp(contract.status, 'source_data_and_unit_contract_only') ...
        || ~contract.scope.axisymmetric || contract.scope.planar_2d_equivalent ...
        || contract.scope.angle_of_attack_deg ~= 0 || contract.scope.flare_angle_deg ~= 20
    error('RocketAeroHeat:TmrContract:Scope', '本合同仅适用于零攻角、20度轴对称外扩面数据。');
end
manifestFile = fullfile(dataRoot, 'source_manifest.json');
if ~strcmp(sha256File(manifestFile), contract.source_manifest_sha256)
    error('RocketAeroHeat:TmrContract:Hash', '原始来源清单与合同绑定的SHA256不同。');
end
manifest = jsondecode(fileread(manifestFile));

%% 模块二：核对每份来源文件的SHA256、长度及归档提取身份
for j = 1:numel(manifest.files)
    if iscell(manifest.files)
        item = manifest.files{j};
    else
        item = manifest.files(j);
    end
    % 来源清单只允许在本数据目录内的相对路径，避免意外读取其它文件。
    if ~isempty(regexp(item.file, '(^[/\\])|(^[A-Za-z]:)|(^|[/\\])\.\.([/\\]|$)', 'once'))
        error('RocketAeroHeat:TmrContract:Path', '来源文件路径不在合同目录内。');
    end
    sourcePath = fullfile(dataRoot, item.file);
    info = dir(sourcePath);
    if isempty(info) || info.bytes ~= item.bytes || ~strcmp(sha256File(sourcePath), item.sha256)
        error('RocketAeroHeat:TmrContract:Hash', '原始来源已改变或缺失：%s', item.file);
    end
end

%% 模块三：逐行核对TMR公开数字，保留实验重复行，禁止将占位列当测量
upstream = readNumericRows(fullfile(dataRoot, 'raw', 'upstream_exp.dat'), 3);
pressure = readNumericRows(fullfile(dataRoot, 'raw', 'pw_exp.dat'), 2);
heating = readNumericRows(fullfile(dataRoot, 'raw', 'qw_exp.dat'), 2);
if ~isequal(size(upstream), [21 3]) || ~isequal(size(pressure), [27 2]) ...
        || ~isequal(size(heating), [27 2])
    error('RocketAeroHeat:TmrContract:Rows', '原始测点数量不符合本切片。');
end
p = contract.measurements.pressure_rows;
q = contract.measurements.heating_rows;
u = contract.measurements.upstream_profile;
checkSame([p.source_row]', (1:27)');
checkSame([q.source_row]', (1:27)');
checkSame([u.source_row]', (1:21)');
checkSame([p.S_m]', pressure(:, 1)*0.01);
checkSame([p.ratio]', pressure(:, 2));
checkSame([q.S_m]', heating(:, 1)*0.01);
checkSame([q.ratio]', heating(:, 2));
checkSame([u.y_from_wall_m]', upstream(:, 1)*0.01);
checkSame([u.velocity_ratio]', upstream(:, 2));
checkSame([u.total_temperature_ratio]', upstream(:, 3));
originalP = readTecplotZone(fullfile(dataRoot, 'raw', 'KussoyCompCorner.dat'), 'KussoyP20', 4);
originalQ = readTecplotZone(fullfile(dataRoot, 'raw', 'KussoyCompCorner.dat'), 'KussoyQ20', 4);
originalU = readTecplotZone(fullfile(dataRoot, 'raw', 'KussoyBLProfile.dat'), 'KussoyBL20_IncomingBL', 5);
checkSame([u.measured_M]', originalU(:, 2));
checkSame([u.pressure_ratio]', originalU(:, 3));
if ~isequal(size(originalP), [24 4]) || ~isequal(size(originalQ), [27 4]) ...
        || ~isequal(size(originalU), [21 5])
    error('RocketAeroHeat:TmrContract:Rows', '官方数据库中的对应区域尺寸不正确。');
end
uniqueMask = true(27, 1);
for j = 2:27
    uniqueMask(j) = ~any(all(pressure(1:j-1, :) == pressure(j, :), 2));
end
if ~isequal(logical([p.independent_comparison_point]'), uniqueMask) || nnz(uniqueMask) ~= 24
    error('RocketAeroHeat:TmrContract:Duplicates', '压力重复行不得增加独立测点权重。');
end
if any(abs(pressure(uniqueMask, :)-originalP(:, [1 2])) > 1e-6, 'all') ...
        || any(abs(heating-originalQ(:, [1 3])) > 1e-6, 'all') ...
        || any(abs(upstream-originalU(:, [1 4 5])) > 1e-6, 'all')
    error('RocketAeroHeat:TmrContract:Mismatch', 'TMR与原作者数据库超过转录精度差异。');
end

%% 模块四：核对量纲、沿壁坐标与项目坐标转换；输入和测量参考保持独立
ref = contract.experimental_local_reference;
inlet = contract.computational_truncated_inflow;
checkSame(ref.Mach, 7.05);
checkSame(inlet.Mach, 7.11);
checkSame(ref.location_source_x_m, -0.06);
checkSame(inlet.source_x_m, -0.7839);
checkSame(contract.geometry.wall_condition.temperature_K, 311);
checkSame(ref.pressure_Pa, 576);
checkSame(ref.reported_wall_heating_magnitude_W_m2, 9300);
checkSame(contract.geometry.cylinder_radius_m, 0.1015);
checkSame(contract.geometry.cylinder_diameter_m, 0.203);
checkSame(contract.geometry.flare_start_distance_from_nose_m, 1.39);
checkSame(contract.geometry.nose_cone_half_angle_deg, 10);
if ~strcmp(contract.geometry.wall_condition.type, 'prescribed_no_slip_isothermal') ...
        || ~isempty(contract.geometry.nose_arc_radius_m) || ~isempty(contract.geometry.flare_end_S_m)
    error('RocketAeroHeat:TmrContract:Missing', '给定壁面身份或未知完整前体尺寸不得改写。');
end
conditions = fileread(fullfile(dataRoot, 'raw', 'KussoyConditions.txt'));
checkSame(ref.Mach, conditionValue(conditions, 'Minf'));
checkSame(ref.temperature_K, conditionValue(conditions, 'Tinf'));
checkSame(ref.pressure_Pa, conditionValue(conditions, 'Pinf'));
checkSame(ref.velocity_m_s, conditionValue(conditions, 'Uinf'));
checkSame(ref.wall_temperature_K, conditionValue(conditions, 'TWall'));
checkSame(ref.reported_wall_heating_magnitude_W_m2, conditionValue(conditions, 'QWallInf'));
checkSame(ref.density_kg_m3, conditionValue(conditions, 'RhoInf'));
checkSame(ref.wall_shear_Pa, conditionValue(conditions, 'TauWallInf'));
checkSame(ref.unit_Re_m_inv, conditionValue(conditions, 'Re/m'));
checkSame(ref.boundary_layer_thickness_m, conditionValue(conditions, 'Delta0')*0.01);
checkSame(ref.displacement_thickness_m, conditionValue(conditions, 'Delta\*0')*0.01);
checkSame(ref.momentum_thickness_m, conditionValue(conditions, 'Theta0')*0.01);
checkSame(ref.Cf, conditionValue(conditions, 'CfInf'));
windInput = fileread(fullfile(dataRoot, 'raw', 'wind_sa.inp'));
windValues = regexp(windInput, '(?im)^freestream\s+total\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)', 'tokens', 'once');
if numel(windValues) ~= 5
    error('RocketAeroHeat:TmrContract:Conditions', '原始WIND总量输入缺失。');
end
windValues = cellfun(@str2double, windValues);
checkSame([inlet.Mach inlet.total_pressure_psi_WIND_input inlet.total_temperature_rankine_WIND_input ...
    inlet.angle_of_attack_deg inlet.sideslip_deg], windValues);
checkSame([p.pressure_Pa]', pressure(:, 2)*ref.pressure_Pa);
checkSame([q.reported_heating_magnitude_W_m2]', heating(:, 2)*ref.reported_wall_heating_magnitude_W_m2);
checkSame(inlet.total_pressure_Pa_WIND_input, inlet.total_pressure_psi_WIND_input*6894.757293168);
checkSame(inlet.total_temperature_K_WIND_input, inlet.total_temperature_rankine_WIND_input*5/9);
checkSame(inlet.wall_temperature_rankine_WIND_input*5/9, 311);
checkCoordinates(p, contract.geometry.cylinder_radius_m);
checkCoordinates(q, contract.geometry.cylinder_radius_m);
if ~isempty(ref.total_temperature_K) || ~isempty(inlet.incoming_turbulence_profiles)
    error('RocketAeroHeat:TmrContract:Missing', '本切片缺失的总温或湍流剖面不得自行补值。');
end

%% 模块五：输出可比较数组和审计状态，不建立生产面元覆盖
benchmark = struct('contract', contract, 'sourceManifest', manifest);
benchmark.pressure = struct('S_m', [p.S_m]', 'sourceX_m', [p.source_x_m]', ...
    'projectX_m', [p.project_X_m]', 'radius_m', [p.radius_m]', ...
    'ratio', [p.ratio]', 'pressure_Pa', [p.pressure_Pa]', ...
    'independentComparisonPoint', uniqueMask, 'rawSourceRows', (1:27)');
benchmark.heating = struct('S_m', [q.S_m]', 'sourceX_m', [q.source_x_m]', ...
    'projectX_m', [q.project_X_m]', 'radius_m', [q.radius_m]', ...
    'ratio', [q.ratio]', 'reportedMagnitude_W_m2', [q.reported_heating_magnitude_W_m2]');
benchmark.upstream = struct('projectX_m', -ref.location_source_x_m, ...
    'yFromWall_m', [u.y_from_wall_m]', 'velocityRatio', [u.velocity_ratio]', ...
    'totalTemperatureRatio', [u.total_temperature_ratio]', ...
    'measuredMach', [u.measured_M]', 'pressureRatio', [u.pressure_ratio]');
benchmark.audit = struct('sourceHashesVerified', true, 'unitAndCoordinateContractVerified', true, ...
    'sourceFileCount', numel(manifest.files), 'independentPressurePoints', 24, ...
    'heatingPoints', 27, 'upstreamPoints', 21, 'simulationValidationComplete', false, ...
    'newPhysicalArea_m2', 0);
end

function checkCoordinates(rows, radius)
% 模块用途：沿壁弧长不等于外扩面轴向距离，项目X方向与实验方向相反。
s = [rows.S_m]';
x = s;
flare = s > 0;
x(flare) = s(flare)*cosd(20);
r = radius+max(s, 0)*sind(20);
checkSame([rows.source_x_m]', x);
checkSame([rows.project_X_m]', -x);
checkSame([rows.radius_m]', r);
end

function value = conditionValue(text, name)
% 模块用途：直接从官方条件文件核对标量，不从云图或别的工况猜测参考量。
token = regexp(text, ['(?im)^\s*' name '\s*=\s*([-+0-9.Ee]+)'], 'tokens', 'once');
if isempty(token)
    error('RocketAeroHeat:TmrContract:Conditions', '官方条件文件缺少所需量：%s', name);
end
value = str2double(token{1});
end

function rows = readTecplotZone(path, zoneName, ncol)
% 模块用途：只取单一20度数据区，不把其它外扩角或占位列并成测量。
lines = regexp(fileread(path), '\r\n|\n|\r', 'split');
rows = zeros(0, ncol);
active = false;
for j = 1:numel(lines)
    line = strtrim(lines{j});
    token = regexp(line, '(?i)^zone\s+t="([^"]+)"', 'tokens', 'once');
    if ~isempty(token)
        if active, break, end
        active = strcmp(token{1}, zoneName);
    elseif active && ~isempty(regexp(line, '^[-+0-9.]', 'once'))
        values = sscanf(line, '%f')';
        if numel(values) ~= ncol || any(~isfinite(values))
            error('RocketAeroHeat:TmrContract:Numeric', '官方数据区的数字列不正确。');
        end
        rows(end+1, :) = values; %#ok<AGROW>
    end
end
end

function rows = readNumericRows(path, ncol)
% 模块用途：只读取TMR测量文件里的数字行，拒绝列数错误及非有限值。
lines = regexp(fileread(path), '\r\n|\n|\r', 'split');
rows = zeros(0, ncol);
for j = 1:numel(lines)
    line = strtrim(lines{j});
    if isempty(line) || isempty(regexp(line, '^[-+0-9.]', 'once'))
        continue
    end
    values = sscanf(line, '%f')';
    if numel(values) ~= ncol || any(~isfinite(values))
        error('RocketAeroHeat:TmrContract:Numeric', '测量原文件存在列数错误或非有限数。');
    end
    rows(end+1, :) = values; %#ok<AGROW>
end
end

function checkSame(actual, expected)
% 模块用途：核对转录值和单位换算，容差只用于浮点舍入。
if ~isequal(size(actual), size(expected)) || ~isreal(actual) || any(~isfinite(actual(:))) ...
        || any(abs(actual(:)-expected(:)) > 2e-12*max(1, abs(expected(:))))
    error('RocketAeroHeat:TmrContract:Mismatch', '合同与原始数字或单位/坐标转换不一致。');
end
end

function value = sha256File(path)
% 模块用途：核对原始文件字节；SHA通过只证明来源未变，不证明算法正确。
fid = fopen(path, 'rb');
if fid < 0
    error('RocketAeroHeat:TmrContract:Hash', '无法读取来源文件：%s', path);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes = fread(fid, Inf, '*uint8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
value = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2)', 1, []));
end
