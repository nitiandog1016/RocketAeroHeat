function grid = loadTmrMach7Grid(levelDown, dataRoot)
%LOADTMRMACH7GRID 读取冻结NASA/TMR近Mach7轴对称外扩面原始网格。
% 模块用途：只提供几何节点与边界身份，不导入CFD解、不生成物理解。
% 输入 levelDown：0为最细321×401，4为最粗21×26，缺省为4。
% 输入 dataRoot：tests/data/tmr_mach7，可省略；必须保留冻结原始归档。
% 输出 blocks(1:2).sourceX/X/R：两块原始节点，单位m，X=-sourceX。
% 输出 sourceX/X/R：沿流向i拼接，并只保留一次完全一致的共享接口。
% i随实验sourceX增加、随项目X减小；j从壁面指向外域，R始终为正。
% 本二维网格是轴对称子午面，不能作为平面二维流动使用。

%% 模块一：固定入口与原始来源身份，拒绝任意改写的归档/清单
if nargin < 1 || isempty(levelDown), levelDown = 4; end
if ~isnumeric(levelDown) || islogical(levelDown) || ~isscalar(levelDown) ...
        || ~isreal(levelDown) || ~isfinite(levelDown) ...
        || levelDown ~= fix(levelDown) || levelDown < 0 || levelDown > 4
    error('RocketAeroHeat:TmrGrid:Input', '网格档位必须是0至4的有限实整数。');
end
levelDown = double(levelDown);
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 2 || isempty(dataRoot)
    dataRoot = fullfile(projectRoot, 'tests', 'data', 'tmr_mach7');
end
if ~(ischar(dataRoot) && isrow(dataRoot)) ...
        && ~(isstring(dataRoot) && isscalar(dataRoot) && ~ismissing(dataRoot))
    error('RocketAeroHeat:TmrGrid:Input', '数据根目录必须是单一文本路径。');
end
dataRoot = char(dataRoot);
if ~exist(dataRoot, 'dir')
    error('RocketAeroHeat:TmrGrid:Source', '未找到冻结TMR数据目录。');
end
archivePath = fullfile(dataRoot, 'raw', 'Axiswblim7_grids.zip');
manifestPath = fullfile(dataRoot, 'source_manifest.json');
if ~exist(archivePath, 'file') || ~exist(manifestPath, 'file')
    error('RocketAeroHeat:TmrGrid:Source', '缺少原始网格ZIP或其冻结来源清单。');
end
archiveSHA = '6d12b620b0e1fff13c0567a9ad4f18c88133ed9d87ee332be8c9d775c7cdc4aa';
manifestSHA = '55606c2eb77eaf1e83f051bed47a639aba3d0f18abb171e3bc2705d49970352d';
if ~strcmp(sha256File(manifestPath), manifestSHA) || ~strcmp(sha256File(archivePath), archiveSHA)
    error('RocketAeroHeat:TmrGrid:Hash', '网格归档或来源清单SHA256与冻结版本不同。');
end
manifest = jsondecode(fileread(manifestPath));
matches = 0;
record = struct();
for k = 1:numel(manifest.files)
    if iscell(manifest.files), item = manifest.files{k}; else, item = manifest.files(k); end
    if strcmp(item.file, 'raw/Axiswblim7_grids.zip'), matches = matches+1; record = item; end
end
if matches ~= 1 || ~strcmp(record.sha256, archiveSHA) || record.bytes ~= 8826962
    error('RocketAeroHeat:TmrGrid:Hash', '原始ZIP在来源清单中的身份不一致。');
end

%% 模块二：只选择明确的子午面成员，检查名称、唯一性及压缩成员SHA
memberRoot = 'u/piyer/nasa_tmr/gitlab/turbmodels/Axiswblim7_grids/';
memberName = sprintf('%skussoy_2zn_%dlevdn.p2dfmt.gz', memberRoot, levelDown);
memberSHA = { ...
    '4466f018d1c3a00f35deb75ea473445f01f517bec16b475191d9cb1c3892b323', ...
    '8ae30fefc15d7c58d55305d9820d13ff7de603a670c15ffd670738befbcac522', ...
    '8980c3ba3c2ebdfdd661617dc27380df2af9a3ed45822fa8b4a1dba1596f11d8', ...
    '82c8c16087cd77993b48211e816c0d53501a26f4c1569d40c3baf62969a24f70', ...
    '29486068c3a8f52a06dd95bec9e4b77805036502cc7c654ac92ff163da1ead4f'};
plainSHA = { ...
    '57809884702937a96feb1946c093d3f5f6dc3ed02747f0904ebb4455a6632408', ...
    'da1bbf898f1452094e11ba0e78f9a4efea8faaa5a5a603017b74ef90d7429655', ...
    '6caf17da534593f853dac58ad6a5b8db8c6bcd82a03077818a85cbc2cb54fcde', ...
    '64249892a11ff710c8f9cb9f08a738ff4ac6591b523d3262c0df485ffa1d2dc4', ...
    'd1cde541fa49afbdf385efdba9d61e55f23da6f5adb7ee3dd38e716cae8fe59c'};
inspectArchive(archivePath, memberName, memberRoot);
temporaryRoot = [tempname '_RocketAeroHeat_TmrMach7Grid'];
mkdir(temporaryRoot);
cleanup = onCleanup(@() removeTemporary(temporaryRoot)); %#ok<NASGU>
try
    unzip(archivePath, temporaryRoot);
    compressed = fullfile(temporaryRoot, strrep(memberName, '/', filesep));
    if ~exist(compressed, 'file') || ~strcmp(sha256File(compressed), memberSHA{levelDown+1})
        error('RocketAeroHeat:TmrGrid:Hash', '选中网格压缩成员缺失或SHA不一致。');
    end
    uncompressed = gunzip(compressed, temporaryRoot);
catch ex
    if startsWith(ex.identifier, 'RocketAeroHeat:TmrGrid:'), rethrow(ex); end
    error('RocketAeroHeat:TmrGrid:Archive', '原始ZIP/GZIP读取失败：%s', ex.message);
end
if numel(uncompressed) ~= 1 || ~strcmp(sha256File(uncompressed{1}), plainSHA{levelDown+1})
    error('RocketAeroHeat:TmrGrid:Hash', '解压后的原始PLOT3D字节不一致。');
end

%% 模块三：按PLOT3D列主序读原节点，保留cm原值及SI/项目坐标
text = fileread(uncompressed{1});
[values, ~, message, nextIndex] = sscanf(text, '%f');
if ~isempty(message) || (nextIndex <= numel(text) && ~isempty(strtrim(text(nextIndex:end)))) ...
        || ~isreal(values) || any(~isfinite(values)) || numel(values) < 5
    error('RocketAeroHeat:TmrGrid:Format', 'PLOT3D数据存在非法、非有限或未读取字符。');
end
ni = 320/2^levelDown+1;
nj = 400/2^levelDown+1;
if values(1) ~= 2 || ~isequal(values(2:5)', [ni nj ni nj]) ...
        || numel(values) ~= 5+4*ni*nj
    error('RocketAeroHeat:TmrGrid:Format', 'PLOT3D块数、尺寸或数字总数不符合冻结档位。');
end
blocks = repmat(struct('sourceX_cm', [], 'R_cm', [], 'sourceX', [], 'X', [], 'R', [], ...
    'ni', ni, 'nj', nj), 1, 2);
cursor = 6;
for b = 1:2
    number = ni*nj;
    blocks(b).sourceX_cm = reshape(values(cursor:cursor+number-1), ni, nj);
    blocks(b).R_cm = reshape(values(cursor+number:cursor+2*number-1), ni, nj);
    blocks(b).sourceX = blocks(b).sourceX_cm*0.01;
    blocks(b).X = -blocks(b).sourceX;
    blocks(b).R = blocks(b).R_cm*0.01;
    cursor = cursor+2*number;
    if any(blocks(b).R(:) <= 0) || any(diff(blocks(b).sourceX, 1, 1) <= 0, 'all') ...
            || any(diff(blocks(b).R, 1, 2) <= 0, 'all')
        error('RocketAeroHeat:TmrGrid:Geometry', '原节点半径或流向/壁外方向不符合本子午面身份。');
    end
end

%% 模块四：共享接口必须逐值相同；仅删除重复接口行，不重插值或修几何
if ~isequal(blocks(1).sourceX_cm(end, :), blocks(2).sourceX_cm(1, :)) ...
        || ~isequal(blocks(1).R_cm(end, :), blocks(2).R_cm(1, :))
    error('RocketAeroHeat:TmrGrid:Interface', '两块原网格共享接口不完全一致。');
end
grid = struct();
grid.schema = 'RocketAeroHeat.axisymmetric-tmr-grid.v1';
grid.levelDown = levelDown;
grid.blocks = blocks;
grid.sourceX = [blocks(1).sourceX; blocks(2).sourceX(2:end, :)];
grid.X = -grid.sourceX;
grid.R = [blocks(1).R; blocks(2).R(2:end, :)];
grid.ni = 2*ni-1;
grid.nj = nj;
grid.niPerBlock = ni;
grid.interface = struct('joinedRow', ni, 'block1Row', ni, 'block2Row', 1, ...
    'nodeCount', nj, 'exactEquality', true);
grid.boundaries.iMin = boundary('inlet', ones(nj, 1), (1:nj)', [grid.ni nj]);
grid.boundaries.iMax = boundary('outlet', grid.ni*ones(nj, 1), (1:nj)', [grid.ni nj]);
grid.boundaries.jMin = boundary('wall', (1:grid.ni)', ones(grid.ni, 1), [grid.ni nj]);
grid.boundaries.jMax = boundary('farfield', (1:grid.ni)', nj*ones(grid.ni, 1), [grid.ni nj]);
grid.units = 'm';
grid.sourceUnits = 'cm';
grid.coordinateConvention = 'sourceX随流向增加；项目X=-sourceX，鼻尖+X、零攻角来流-X；R为轴距';
grid.indexConvention = 'i=inlet至outlet；j=wall至farfield；节点不是单元中心';
grid.axisymmetric = true;
grid.planarTwoDimensionalFlow = false;
grid.physicalSolution = false;
grid.eligibleForSurfaceCoverage = false;
grid.newSolvedAreaM2 = 0;
grid.source = struct('archive', 'raw/Axiswblim7_grids.zip', 'archiveSHA256', archiveSHA, ...
    'manifestSHA256', manifestSHA, 'member', memberName, 'memberSHA256', memberSHA{levelDown+1}, ...
    'uncompressedSHA256', plainSHA{levelDown+1}, 'originalNodesPreserved', true, ...
    'geometryOnly', true, 'externalCFDSolutionImported', false);
end

function item = boundary(kind, i, j, shape)
% 模块用途：边界只标明几何身份与节点索引，不自行给定压力/湍流或出口条件。
item = struct('kind', kind, 'nodeIJ', [i j], 'linearNodeIndices', sub2ind(shape, i, j));
end

function inspectArchive(path, selectedMember, memberRoot)
% 模块用途：在解压前检查成员安全与唯一性，防止读取相似命名的其它网格。
archive = java.util.zip.ZipFile(path);
cleanup = onCleanup(@() archive.close()); %#ok<NASGU>
entries = archive.entries();
names = cell(0, 1);
while entries.hasMoreElements()
    entry = entries.nextElement();
    name = char(entry.getName());
    if ~startsWith(name, memberRoot) || ~isempty(regexp(name, '(^[/\\])|:|\\|(^|/)\.\.(/|$)', 'once'))
        error('RocketAeroHeat:TmrGrid:Archive', 'ZIP成员路径不在冻结来源目录中。');
    end
    names{end+1, 1} = name; %#ok<AGROW>
end
if numel(names) ~= 17 || numel(unique(names)) ~= 17 || nnz(strcmp(names, selectedMember)) ~= 1
    error('RocketAeroHeat:TmrGrid:Archive', 'ZIP成员数、唯一性或目标网格身份不正确。');
end
end

function hash = sha256File(path)
% 模块用途：对原始字节核验SHA，不将文件名或ZIP内部时间戳当成来源证明。
fid = fopen(path, 'rb');
if fid < 0, error('RocketAeroHeat:TmrGrid:Source', '不能读取来源文件：%s', path); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
bytes = fread(fid, Inf, '*uint8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
hash = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2)', 1, []));
end

function removeTemporary(path)
% 模块用途：只清理本函数在系统临时目录创建、且绝对路径经检查的工作副本。
absolute = char(java.io.File(path).getCanonicalPath());
temporary = char(java.io.File(tempdir).getCanonicalPath());
if startsWith(lower(absolute), [lower(temporary) filesep]) ...
        && endsWith(absolute, '_RocketAeroHeat_TmrMach7Grid') && exist(absolute, 'dir')
    rmdir(absolute, 's');
end
end
