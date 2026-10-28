function result = testTmrMach7Benchmark(outputDirectory)
% 模块用途：验证公开基准的来源绑定、转录、单位、坐标及未知量身份。
% 这些断言不执行NS/RANS，也不能作为计算与实验已经吻合的证据。
if nargin < 1, outputDirectory = ''; end
addpath(fileparts(mfilename('fullpath')));
b = loadTmrMach7Benchmark();
checks = 0;

%% 模块一：核对独立测量数量和两份官方数据之间的对应
assert(b.audit.sourceHashesVerified && b.audit.unitAndCoordinateContractVerified);
assert(size(b.pressure.S_m, 1) == 27 && nnz(b.pressure.independentComparisonPoint) == 24);
assert(isequal(find(~b.pressure.independentComparisonPoint), [25;26;27]));
assert(numel(b.heating.S_m) == 27 && numel(b.upstream.yFromWall_m) == 21);
checks = checks+4;

%% 模块二：明确实验当地状态、截短计算入口和沿壁测点的不同身份
c = b.contract;
assert(c.experimental_local_reference.Mach == 7.05 && c.computational_truncated_inflow.Mach == 7.11);
assert(c.experimental_local_reference.temperature_K == 81.2 && c.computational_truncated_inflow.temperature_K_figure == 80);
assert(abs(b.upstream.projectX_m-0.06) < 1e-15);
assert(abs(c.computational_truncated_inflow.total_temperature_K_WIND_input-888.8888888888889) < 1e-10);
assert(c.geometry.wall_condition.temperature_K == 311 && ~c.scope.planar_2d_equivalent);
checks = checks+5;

%% 模块三：检查SI转换、坐标反向及外扩段弧长，避免错配压力/热流测点
assert(abs(b.pressure.pressure_Pa(1)-b.pressure.ratio(1)*576) < 1e-10);
assert(abs(b.heating.reportedMagnitude_W_m2(1)-b.heating.ratio(1)*9300) < 1e-10);
assert(b.heating.S_m(1) < 0 && abs(b.heating.S_m(1)+0.1206) < 1e-8);
assert(abs(b.heating.sourceX_m(end)-b.heating.S_m(end)*cosd(20)) < 1e-14);
assert(abs(b.heating.radius_m(end)-0.1015-b.heating.S_m(end)*sind(20)) < 1e-14);
assert(all(b.pressure.projectX_m == -b.pressure.sourceX_m));
assert(abs(b.upstream.totalTemperatureRatio(1)-0.35) < 1e-7 && b.upstream.velocityRatio(1) == 0);
assert(~isfield(b.pressure, 'wallShear_Pa') && ~isfield(b.heating, 'pressure_Pa'));
checks = checks+8;

%% 模块四：缺失量和证据身份不得被数值零或计算入口参考量替代
assert(isempty(c.experimental_local_reference.total_temperature_K));
assert(isempty(c.computational_truncated_inflow.incoming_turbulence_profiles));
assert(~b.audit.simulationValidationComplete && b.audit.newPhysicalArea_m2 == 0);
checks = checks+3;

%% 模块五：通过篡改副本检查来源绑定和测点合同确实会拒绝失配
sourceRoot = fullfile(fileparts(mfilename('fullpath')), 'data', 'tmr_mach7');
temporaryRoot = [tempname '_RocketAeroHeat_TMR_contract'];
mkdir(temporaryRoot);
cleanup = onCleanup(@() removeTestCopy(temporaryRoot)); %#ok<NASGU>
copyfile(fullfile(sourceRoot, '*'), temporaryRoot);
fid = fopen(fullfile(temporaryRoot, 'raw', 'pw_exp.dat'), 'a');
assert(fid >= 0); fprintf(fid, '\n'); fclose(fid);
expectError(@() loadTmrMach7Benchmark(temporaryRoot), 'RocketAeroHeat:TmrContract:Hash');
copyfile(fullfile(sourceRoot, 'raw', 'pw_exp.dat'), fullfile(temporaryRoot, 'raw', 'pw_exp.dat'));
changed = c;
changed.measurements.pressure_rows(1).pressure_Pa = 0;
fid = fopen(fullfile(temporaryRoot, 'contract.json'), 'w');
assert(fid >= 0); fprintf(fid, '%s', jsonencode(changed)); fclose(fid);
expectError(@() loadTmrMach7Benchmark(temporaryRoot), 'RocketAeroHeat:TmrContract:Mismatch');
changed = c;
changed.measurements.pressure_rows(25).independent_comparison_point = true;
fid = fopen(fullfile(temporaryRoot, 'contract.json'), 'w');
assert(fid >= 0); fprintf(fid, '%s', jsonencode(changed)); fclose(fid);
expectError(@() loadTmrMach7Benchmark(temporaryRoot), 'RocketAeroHeat:TmrContract:Duplicates');
changed = c;
changed.scope.planar_2d_equivalent = true;
writeTestContract(temporaryRoot, changed);
expectError(@() loadTmrMach7Benchmark(temporaryRoot), 'RocketAeroHeat:TmrContract:Scope');
changed = c;
changed.computational_truncated_inflow.total_pressure_psi_WIND_input = 400;
changed.computational_truncated_inflow.total_pressure_Pa_WIND_input = 400*6894.757293168;
writeTestContract(temporaryRoot, changed);
expectError(@() loadTmrMach7Benchmark(temporaryRoot), 'RocketAeroHeat:TmrContract:Mismatch');
changed = c;
changed.experimental_local_reference.total_temperature_K = 888.8888888888889;
writeTestContract(temporaryRoot, changed);
expectError(@() loadTmrMach7Benchmark(temporaryRoot), 'RocketAeroHeat:TmrContract:Missing');
checks = checks+6;

%% 模块六：保存本次数据合同测试终态，与将来的物理验证结果分开
result = b.audit;
result.status = 'passed_data_contract_only';
result.checksPassed = checks;
result.solverRun = false;
if ~isempty(outputDirectory)
    if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
    fid = fopen(fullfile(outputDirectory, 'matlab_contract_test.json'), 'w');
    assert(fid >= 0); fprintf(fid, '%s\n', jsonencode(result)); fclose(fid);
end
fprintf('TMR Mach7 data contract: %d checks passed; no CFD/experimental validation claimed.\n', checks);
end

function writeTestContract(root, contract)
% 模块用途：只修改临时副本，构造应被拒绝的实验身份或边界输入。
fid = fopen(fullfile(root, 'contract.json'), 'w');
assert(fid >= 0); fprintf(fid, '%s', jsonencode(contract)); fclose(fid);
end

function expectError(action, identifier)
% 模块用途：确认合同拒绝真实失配，避免只检测正常路径。
try
    action();
catch ex
    assert(strcmp(ex.identifier, identifier), '出现非预期错误：%s', ex.identifier);
    return
end
error('RocketAeroHeat:TmrContract:Test', '失配合同未被拒绝：%s', identifier);
end

function removeTestCopy(path)
% 模块用途：只删除本测试在系统临时目录下创建的明确副本。
fullPath = char(java.io.File(path).getCanonicalPath());
tempRoot = char(java.io.File(tempdir).getCanonicalPath());
if startsWith(lower(fullPath), [lower(tempRoot) filesep]) ...
        && endsWith(fullPath, '_RocketAeroHeat_TMR_contract') && exist(fullPath, 'dir')
    rmdir(fullPath, 's');
end
end
