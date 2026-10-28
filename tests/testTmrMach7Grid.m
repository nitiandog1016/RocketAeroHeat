function result = testTmrMach7Grid(outputDirectory)
%TESTTMRMACH7GRID 冻结子午面节点、嵌套、接口、坐标及来源拒绝合同测试。
% 模块用途：只测试网格读取，不求解NS/RANS，也不宣称近壁分辨率足够。
if nargin < 1, outputDirectory = ''; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'regionalFlow', 'axisymmetric'));
dataRoot = fullfile(root, 'tests', 'data', 'tmr_mach7');
grids = cell(5, 1);
checks = 0;

%% 模块一：五档原节点、单位、严格接口和边界逐一核对
for level = 0:4
    g = loadTmrMach7Grid(level, dataRoot);
    ni = 320/2^level+1;
    nj = 400/2^level+1;
    assert(g.ni == 2*ni-1 && g.nj == nj && numel(g.blocks) == 2);
    assert(isequal(size(g.X), [2*ni-1 nj]) && isequal(size(g.R), [2*ni-1 nj]));
    for b = 1:2
        assert(isequal(size(g.blocks(b).sourceX), [ni nj]));
        assert(isequal(g.blocks(b).sourceX, .01*g.blocks(b).sourceX_cm));
        assert(isequal(g.blocks(b).R, .01*g.blocks(b).R_cm));
        assert(isequal(g.blocks(b).X, -g.blocks(b).sourceX));
    end
    assert(isequal(g.blocks(1).X(end, :), g.blocks(2).X(1, :)));
    assert(isequal(g.blocks(1).R(end, :), g.blocks(2).R(1, :)));
    assert(g.interface.exactEquality && g.interface.joinedRow == ni && g.interface.nodeCount == nj);
    assert(isequal(g.sourceX, [g.blocks(1).sourceX; g.blocks(2).sourceX(2:end, :)]));
    assert(isequal(g.R, [g.blocks(1).R; g.blocks(2).R(2:end, :)]));
    assert(isequal(g.X, -g.sourceX));
    assert(all(g.R(:) > 0) && all(isfinite(g.X(:))) && all(isfinite(g.R(:))));
    assert(all(diff(g.X, 1, 1) < 0, 'all') && all(diff(g.R, 1, 2) > 0, 'all'));
    assert(abs(g.sourceX(1, 1)+.78389988) < 5e-16);
    assert(abs(g.R(1, 1)-.10150094) < 5e-16 && abs(g.R(1, end)-.50790094) < 5e-16);
    assert(abs(g.sourceX(end, 1)-.23492206) < 5e-16 && abs(g.R(end, 1)-.18700496) < 5e-16);
    verifyBoundary(g.boundaries.iMin, 'inlet', ones(nj, 1), (1:nj)', size(g.X));
    verifyBoundary(g.boundaries.iMax, 'outlet', g.ni*ones(nj, 1), (1:nj)', size(g.X));
    verifyBoundary(g.boundaries.jMin, 'wall', (1:g.ni)', ones(g.ni, 1), size(g.X));
    verifyBoundary(g.boundaries.jMax, 'farfield', (1:g.ni)', nj*ones(g.ni, 1), size(g.X));
    assert(g.axisymmetric && ~g.planarTwoDimensionalFlow && ~g.physicalSolution);
    assert(~g.eligibleForSurfaceCoverage && g.newSolvedAreaM2 == 0);
    assert(strcmp(g.units, 'm') && strcmp(g.sourceUnits, 'cm') && g.source.originalNodesPreserved);
    assert(~g.source.externalCFDSolutionImported && g.source.geometryOnly);
    grids{level+1} = g;
    checks = checks+8;
end

%% 模块二：只抽取原节点构成下一档，禁止重采样或微调20度名义几何
for fine = 1:4
    a = grids{fine}; b = grids{fine+1};
    assert(isequal(b.sourceX, a.sourceX(1:2:end, 1:2:end)));
    assert(isequal(b.X, a.X(1:2:end, 1:2:end)));
    assert(isequal(b.R, a.R(1:2:end, 1:2:end)));
    for block = 1:2
        assert(isequal(b.blocks(block).sourceX_cm, a.blocks(block).sourceX_cm(1:2:end, 1:2:end)));
        assert(isequal(b.blocks(block).R_cm, a.blocks(block).R_cm(1:2:end, 1:2:end)));
    end
    checks = checks+2;
end
default = loadTmrMach7Grid();
assert(default.levelDown == 4 && isequal(default.X, grids{5}.X) && isequal(default.R, grids{5}.R));
checks = checks+1;

%% 模块三：非法档位、路径和来源篡改必须先拒绝，不放开冻结SHA来测试坏网格
invalid = {-1, 5, 1.5, Inf, NaN, 1i, [0 1], true, '0', struct()};
for k = 1:numel(invalid)
    expectError(@() loadTmrMach7Grid(invalid{k}, dataRoot), 'RocketAeroHeat:TmrGrid:Input');
end
expectError(@() loadTmrMach7Grid(4, 42), 'RocketAeroHeat:TmrGrid:Input');
temporaryRoot = [tempname '_RocketAeroHeat_TmrMach7GridTest'];
mkdir(temporaryRoot); mkdir(fullfile(temporaryRoot, 'raw'));
cleanup = onCleanup(@() removeTestCopy(temporaryRoot)); %#ok<NASGU>
copyfile(fullfile(dataRoot, 'source_manifest.json'), temporaryRoot);
copyfile(fullfile(dataRoot, 'raw', 'Axiswblim7_grids.zip'), fullfile(temporaryRoot, 'raw'));
archive = fullfile(temporaryRoot, 'raw', 'Axiswblim7_grids.zip');
fid = fopen(archive, 'ab'); assert(fid >= 0); fprintf(fid, 'tampered'); fclose(fid);
expectError(@() loadTmrMach7Grid(4, temporaryRoot), 'RocketAeroHeat:TmrGrid:Hash');
copyfile(fullfile(dataRoot, 'raw', 'Axiswblim7_grids.zip'), archive);
fid = fopen(fullfile(temporaryRoot, 'source_manifest.json'), 'ab');
assert(fid >= 0); fprintf(fid, '\n'); fclose(fid);
expectError(@() loadTmrMach7Grid(4, temporaryRoot), 'RocketAeroHeat:TmrGrid:Hash');
copyfile(fullfile(dataRoot, 'source_manifest.json'), temporaryRoot);
delete(archive);
expectError(@() loadTmrMach7Grid(4, temporaryRoot), 'RocketAeroHeat:TmrGrid:Source');
expectError(@() loadTmrMach7Grid(4, fullfile(temporaryRoot, 'absent')), 'RocketAeroHeat:TmrGrid:Source');
checks = checks+15;

%% 模块四：保存真实MATLAB节点供独立Python逐节点复核，证据仍限几何合同
result = struct('status', 'passed_frozen_grid_contract_only', 'checkGroupsPassed', checks, ...
    'levelsChecked', 0:4, 'exactNestedGrids', true, 'exactSharedInterface', true, ...
    'sourceNodesPreserved', true, 'axisymmetricRequired', true, 'solverRun', false, ...
    'physicalSolution', false, 'newSolvedAreaM2', 0, 'rejectionChecks', 15);
result.rejectionScope = '来源字节篡改在解析前被SHA拒绝；未通过改写SHA绕过冻结来源合同';
result.sourceBindings = [binding(fullfile(root, 'regionalFlow', 'axisymmetric', 'loadTmrMach7Grid.m'), ...
    'regionalFlow/axisymmetric/loadTmrMach7Grid.m'), binding(mfilename('fullpath'), 'tests/testTmrMach7Grid.m')];
if ~isempty(outputDirectory)
    if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
    save(fullfile(outputDirectory, 'grid_contract_evidence.mat'), 'grids', 'result', '-v7');
    fid = fopen(fullfile(outputDirectory, 'matlab_grid_contract.json'), 'w');
    assert(fid >= 0); fprintf(fid, '%s\n', jsonencode(result)); fclose(fid);
end
fprintf('TMR grid contract passed: five levels, exact nodes/interfaces/nesting, %d refusal cases; no solver.\n', result.rejectionChecks);
end

function verifyBoundary(item, kind, i, j, shape)
% 模块用途：入口/出口/壁面/外域标识与拼接后的节点索引应一一一致。
assert(strcmp(item.kind, kind) && isequal(item.nodeIJ, [i j]));
assert(isequal(item.linearNodeIndices, sub2ind(shape, i, j)));
end

function expectError(action, identifier)
% 模块用途：核对明确拒绝，不把任意异常都当成预期成功。
try
    action();
catch ex
    assert(strcmp(ex.identifier, identifier), '预期%s，实际%s。', identifier, ex.identifier);
    return
end
error('RocketAeroHeat:TmrGrid:Test', '非法网格合同未拒绝：%s', identifier);
end

function record = binding(path, relative)
% 模块用途：保存实际执行源文件SHA，不能只依靠文件名声明版本。
if ~endsWith(path, '.m'), path = [path '.m']; end
fid = fopen(path, 'rb'); assert(fid >= 0);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes = fread(fid, Inf, '*uint8'); digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
record = struct('path', relative, 'sha256', lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2)', 1, [])));
end

function removeTestCopy(path)
% 模块用途：只删除本测试创建、绝对路径经确认位于系统临时根目录内的副本。
absolute = char(java.io.File(path).getCanonicalPath());
temporary = char(java.io.File(tempdir).getCanonicalPath());
if startsWith(lower(absolute), [lower(temporary) filesep]) ...
        && endsWith(absolute, '_RocketAeroHeat_TmrMach7GridTest') && exist(absolute, 'dir')
    rmdir(absolute, 's');
end
end
