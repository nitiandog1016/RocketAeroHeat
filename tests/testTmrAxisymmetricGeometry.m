function report = testTmrAxisymmetricGeometry(outputDirectory, levels)
%TESTTMRAXISYMMETRICGEOMETRY 在冻结TMR原网格核对轴对称度量及静止残差。
% 模块一：缺省检查4/3/2三档，保留原节点、两块共享接口及四条边界。
% 模块二：用独立三角形积分核对每弧度真体积、平面面积和体积质心。
% 模块三：给定常压/常温静止解析迹，按局部声学与扩散时间归一残差。
% 模块四：只保存最粗被测网格及逐档摘要，不保存细网格大MAT文件。
% 这是几何与空间算子一致性测试；没有推进Mach7流场，也没有实验验证。
% 调用示例：testTmrAxisymmetricGeometry(folder,4:-1:0)可另查全部五档。

%% 一、固定来源、接受门槛及等温静止边界；这些数值不是Mach7求解入口
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root, 'regionalFlow', 'axisymmetric'));
if nargin < 1 || isempty(outputDirectory)
    outputDirectory = fullfile(root, 'outputs', 'v62_validation', 'tmr_geometry');
end
if nargin < 2 || isempty(levels), levels = 4:-1:2; end
assert((ischar(outputDirectory) && isrow(outputDirectory)) || ...
    (isstring(outputDirectory) && isscalar(outputDirectory) && ~ismissing(outputDirectory)), ...
    'TmrAxisFVTest:Input', '输出目录须为单一路径。');
validateattributes(levels, {'numeric'}, {'real','finite','vector','integer','>=',0,'<=',4});
levels = reshape(double(levels), 1, []);
assert(numel(unique(levels)) == numel(levels), 'TmrAxisFVTest:Input', '网格档位不能重复。');
if ~exist(outputDirectory, 'dir'), mkdir(outputDirectory); end
gas = struct('gamma',1.4,'R',287.05,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]);
gas = axisymmetricFVGas(gas);
temperature = 311;
pressure = 576;
rho = pressure/(gas.R*temperature);
primitive = [rho 0 0 temperature];
trace = struct('type','prescribedTrace','primitive',primitive);
boundary = struct('iMin',trace,'iMax',trace, ...
    'jMin',struct('type','wall','temperature',temperature,'velocity',[0 0]),'jMax',trace);
acceptance = struct('metricClosureRelativeByPlanarArea',1e-9, ...
    'triangleAreaAndVolumeRelative',1e-10,'triangleCentroidRelativeToCellDiameter',1e-8, ...
    'localTimeScaledResidual',1e-11,'wallHeatFluxRelativeToConductiveScale',1e-11, ...
    'conservationDefectRelativeToBoundaryFluxScale',1e-12);
sourcePaths = {'regionalFlow/axisymmetric/loadTmrMach7Grid.m', ...
    'regionalFlow/axisymmetric/buildAxisymmetricFVGeometry.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVGas.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVState.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVConservative.m', ...
    'regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m', ...
    'tests/testTmrAxisymmetricGeometry.m'};
sourceBindings = repmat(struct('path','','sha256',''), 1, numel(sourcePaths));
for k = 1:numel(sourcePaths), sourceBindings(k) = binding(root,sourcePaths{k}); end
report = struct('status','running-not-yet-tested','levelsRequested',levels,'levelsCompleted',[], ...
    'acceptance',acceptance,'sourceBindings',sourceBindings,'levels',struct([]), ...
    'staticPrimitiveInput',primitive,'staticPressureInputPa',pressure,'gasInput',gas, ...
    'spatialResidualEvaluated',false,'timeIntegrationPerformed',false,'Mach7FlowComputed',false, ...
    'experimentValidation',false,'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'newSolvedAreaM2',0,'wallTemperatureValueKind','prescribed', ...
    'boundaryIdentity',struct('iMin','prescribed static analytic trace', ...
        'iMax','prescribed static analytic trace','jMin','prescribed no-slip isothermal wall', ...
        'jMax','prescribed static analytic trace'));
writeJSON(fullfile(outputDirectory,'tmr_geometry_summary.json'),report);

try
    for level = levels
        started = tic;
        grid = loadTmrMach7Grid(level);
        mesh = buildAxisymmetricFVGeometry(grid.X,grid.R);
        f = mesh.faces;
        inside = f.neighbour > 0;
        nc = mesh.cellCount;

        %% 二、原节点、边界及两块接口只拼接一次；共享通量必须只存一列
        assert(isequal(mesh.X,grid.X) && isequal(mesh.R,grid.R));
        assert(isequal(grid.blocks(1).X(end,:),grid.blocks(2).X(1,:)) && ...
            isequal(grid.blocks(1).R(end,:),grid.blocks(2).R(1,:)));
        assert(isequal(mesh.shape,[grid.ni-1 grid.nj-1]) && nc == prod(mesh.shape));
        assert(~mesh.periodicI && ~any(f.periodic) && ~any(f.neighbourShift(:)));
        assert(all(f.boundarySide(inside)==0) && all(f.boundarySide(~inside)>0));
        counts = accumarray(f.boundarySide(~inside),1,[4 1]);
        assert(isequal(counts,[grid.nj-1;grid.nj-1;grid.ni-1;grid.ni-1]));
        columnCount = full(sum(abs(mesh.incidence),1))';
        assert(all(columnCount(inside)==2) && all(columnCount(~inside)==1));
        assert(all(full(sum(mesh.incidence(:,inside),1))==0));
        assert(all(full(sum(mesh.incidence(:,~inside),1))==1));
        [ownerI,ownerJ] = ind2sub(mesh.shape,f.owner);
        neighbourI = zeros(size(f.owner)); neighbourJ = neighbourI;
        [neighbourI(inside),neighbourJ(inside)] = ind2sub(mesh.shape,f.neighbour(inside));
        interface = find(inside & ((ownerI==grid.niPerBlock-1 & neighbourI==grid.niPerBlock) | ...
            (ownerI==grid.niPerBlock & neighbourI==grid.niPerBlock-1)));
        assert(numel(interface)==grid.nj-1 && all(ownerJ(interface)==neighbourJ(interface)));
        assert(isequal(sort(ownerJ(interface)),(1:grid.nj-1)'));
        assert(all(f.boundarySide(interface)==0));
        [~,order] = sort(ownerJ(interface)); interface = interface(order);
        interfaceNodes = [grid.X(grid.interface.joinedRow,:)',grid.R(grid.interface.joinedRow,:)'];
        interfaceDelta = diff(interfaceNodes,1,1);
        interfaceArea = mean([interfaceNodes(1:end-1,2),interfaceNodes(2:end,2)],2) .* ...
            vecnorm(interfaceDelta,2,2);
        interfaceAreaError = max(abs(f.areaPerRadian(interface)-interfaceArea)./interfaceArea);
        assert(interfaceAreaError<1e-12);

        %% 三、独立拆分三角形积分真矩；不以几何模块自身返回值作为参考
        geometryReference = triangleReference(mesh);
        areaError = max(abs(mesh.planarArea-geometryReference.area)./geometryReference.area);
        volumeError = max(abs(mesh.volumePerRadian-geometryReference.volume)./geometryReference.volume);
        centroidError = max(vecnorm(mesh.volumeCentroid-geometryReference.volumeCentroid,2,2) ./ ...
            geometryReference.cellDiameter);
        metric = mesh.incidence*f.areaVectorPerRadian-[zeros(nc,1),mesh.planarArea];
        metricRelative = max(abs(metric)./mesh.planarArea,[],'all');
        assert(all(mesh.planarArea>0) && all(mesh.volumePerRadian>0) && all(f.normalDistance>0));
        assert(areaError<acceptance.triangleAreaAndVolumeRelative && ...
            volumeError<acceptance.triangleAreaAndVolumeRelative);
        assert(centroidError<acceptance.triangleCentroidRelativeToCellDiameter);
        assert(metricRelative<acceptance.metricClosureRelativeByPlanarArea);

        %% 四、常温静止空间残差按各单元真实稳定频率归一，避免薄层原始率误读
        Q = axisymmetricFVConservative(repmat(primitive,nc,1),gas);
        [rate,info] = axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);
        localCharacteristicTime = 1./(info.convectiveFrequency+info.diffusiveFrequency);
        assert(all(isfinite(localCharacteristicTime) & localCharacteristicTime>0));
        acousticSpeed = sqrt(gas.gamma*gas.R*temperature);
        residualScale = [rho,rho*acousticSpeed,rho*acousticSpeed,rho*gas.cv*temperature];
        scaledResidual = abs(rate).*localCharacteristicTime./residualScale;
        maximumScaledResidual = max(scaledResidual,[],'all');
        assert(maximumScaledResidual<acceptance.localTimeScaledResidual, ...
            'TmrAxisFVTest:StaticResidual','静止状态的局部时间归一残差超限。');
        w = info.wall.faceID;
        assert(numel(w)==grid.ni-1 && all(f.boundarySide(w)==3));
        assert(all(info.wall.temperatureInput==temperature) && ~any(info.wall.velocityInput(:)));
        assert(all(info.faceFluxIntegralPerRadian(w,1)==0));
        conductiveScale = gas.k*temperature./f.normalDistance(w);
        heatRelative = max(abs(info.wall.heatFluxIntoWall)./conductiveScale);
        assert(heatRelative<acceptance.wallHeatFluxRelativeToConductiveScale, ...
            'TmrAxisFVTest:WallHeat','等温静止墙热流超过局部导热尺度的舍入级阈值。');
        assert(~any(info.wall.mechanicalPowerIntoFluidPerRadian));
        boundaryArea = sum(f.areaPerRadian(mesh.boundaryFaces));
        fluxScale = residualScale*acousticSpeed*boundaryArea;
        ledger = sum(rate.*mesh.volumePerRadian,1)+info.netBoundaryFlux-info.sourceIntegral;
        ledgerRelative = max(abs(ledger)./fluxScale);
        assert(ledgerRelative<acceptance.conservationDefectRelativeToBoundaryFluxScale);
        assert(~info.physicalSolution && ~info.eligibleForSurfaceCoverage && ...
            strcmp(info.wallTemperatureValueKind,'prescribed'));

        summary = struct('levelDown',level,'status','passed-geometry-and-static-operator-only', ...
            'nodeShape',size(grid.X),'cellShape',mesh.shape,'cellCount',nc, ...
            'faceCount',numel(f.owner),'interiorFaceCount',nnz(inside), ...
            'boundaryFaceCounts',counts','sharedBlockInterfaceFaceCount',numel(interface), ...
            'originalNodesPreserved',true,'oneStoredFacePerSharedInterface',true, ...
            'sharedBlockInterfaceAreaRelativeError',interfaceAreaError, ...
            'minimumPlanarAreaM2',min(mesh.planarArea),'minimumVolumePerRadianM3',min(mesh.volumePerRadian), ...
            'volumePerRadianM3',sum(mesh.volumePerRadian), ...
            'minimumNormalDistanceM',min(f.normalDistance), ...
            'metricClosureRelativeByPlanarArea',metricRelative, ...
            'triangleAreaRelativeError',areaError,'triangleVolumeRelativeError',volumeError, ...
            'triangleCentroidRelativeToCellDiameter',centroidError, ...
            'rawStaticResidualMaximumByEquation',max(abs(rate),[],1), ...
            'localTimeScaledResidualMaximumByEquation',max(scaledResidual,[],1), ...
            'localTimeScaledResidualMaximum',maximumScaledResidual, ...
            'localCharacteristicTimeRangeS',[min(localCharacteristicTime),max(localCharacteristicTime)], ...
            'residualScale',residualScale,'wallHeatFluxMaximumAbsWm2',max(abs(info.wall.heatFluxIntoWall)), ...
            'wallHeatFluxRelativeToConductiveScale',heatRelative, ...
            'wallIntegratedHeatPerRadianW',sum(info.wall.heatFluxIntoWall.*info.wall.areaPerRadian), ...
            'conservationDefect',ledger,'conservationDefectRelativeToBoundaryFluxScale',ledgerRelative, ...
            'elapsedSeconds',toc(started),'source',grid.source, ...
            'wallTemperatureValueKind','prescribed','physicalSolution',false, ...
            'eligibleForSurfaceCoverage',false,'newSolvedAreaM2',0);
        if isempty(report.levels), report.levels=summary; else, report.levels(end+1)=summary; end
        report.levelsCompleted(end+1) = level;
        report.spatialResidualEvaluated = true;
        writeJSON(fullfile(outputDirectory,sprintf('level_%d_summary.json',level)),summary);
        if level == max(levels)
            save(fullfile(outputDirectory,'coarsest_mesh_evidence.mat'),'grid','mesh','Q','rate','info', ...
                'gas','boundary','localCharacteristicTime','residualScale','geometryReference','summary','-v7');
        end
        writeJSON(fullfile(outputDirectory,'tmr_geometry_summary.json'),report);
        fprintf('TMR level %d: %d cells; metric %.3g; local-time residual %.3g; wall |q| %.3g W/m2.\n', ...
            level,nc,metricRelative,maximumScaledResidual,summary.wallHeatFluxMaximumAbsWm2);
    end
    for k = 1:numel(sourceBindings)
        currentBinding = binding(root,sourcePaths{k});
        assert(strcmp(sourceBindings(k).sha256,currentBinding.sha256), ...
            'TmrAxisFVTest:SourceChanged','测试期间源文件改变，证据不予验收。');
    end
    report.status = 'passed-frozen-TMR-geometry-and-static-operator-only';
    report.normalization = ['abs(dQ/dt) * [convectiveFrequency+diffusiveFrequency]^-1 / ' ...
        '[rho,rho*a,rho*a,rho*cv*T]; local to each cell, no time integration'];
    report.limitations = {'static analytic traces are numerical test inputs, not Mach7 experimental boundary data'; ...
        'raw rate grows in thin cells; only the explicitly reported local-time normalization is accepted'; ...
        'no shock, turbulence, separation, reattachment, y-plus or Mach7 heat-flux accuracy established'; ...
        'all volumes and integrated fluxes are per radian; multiply by 2*pi for a full revolution'};
    save(fullfile(outputDirectory,'per_level_summaries.mat'),'report','-v7');
    writeJSON(fullfile(outputDirectory,'tmr_geometry_summary.json'),report);
catch exception
    report.status = 'failed-not-qualified';
    report.failure = struct('identifier',exception.identifier,'message',exception.message);
    writeJSON(fullfile(outputDirectory,'tmr_geometry_summary.json'),report);
    rethrow(exception);
end
end

function reference = triangleReference(mesh)
% 模块用途：每个四边形拆成两三角形，采用二次精确三点积分计算r权真矩。
nc = mesh.cellCount;
area = zeros(nc,1); volume = area; moment = zeros(nc,2); diameter = area;
barycentric = [2/3 1/6 1/6;1/6 2/3 1/6;1/6 1/6 2/3];
for c = 1:nc
    p = mesh.vertices(mesh.cells(c,:),:);
    origin = p(1,:); local = p-origin;
    diameter(c) = max([vecnorm(local(3,:)),vecnorm(local(4,:)-local(2,:))]);
    for index = {[1 2 3],[1 3 4]}
        triangle = local(index{1},:);
        a = triangle(2,:)-triangle(1,:); b = triangle(3,:)-triangle(1,:);
        triangleArea = (a(1)*b(2)-a(2)*b(1))/2;
        assert(triangleArea>0,'TmrAxisFVTest:Triangle','独立三角形积分发现非正面积。');
        points = barycentric*triangle;
        radius = points(:,2)+origin(2);
        area(c) = area(c)+triangleArea;
        volume(c) = volume(c)+triangleArea*mean(radius);
        moment(c,:) = moment(c,:)+triangleArea*mean(radius.*points,1);
    end
    moment(c,:) = origin+moment(c,:)/volume(c);
end
reference = struct('area',area,'volume',volume,'volumeCentroid',moment,'cellDiameter',diameter, ...
    'method','two triangles; degree-two three-point quadrature; SI per radian');
end

function record = binding(root,relative)
% 模块用途：绑定实际执行源码的完整字节，不只记录文件名。
fid = fopen(fullfile(root,relative),'rb'); assert(fid>=0);
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes = fread(fid,Inf,'*uint8'); digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
record = struct('path',relative,'sha256',lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[])));
end

function writeJSON(path,value)
% 模块用途：逐档保存摘要和终态，失败运行不得继续保留passed总状态。
fid = fopen(path,'w','n','UTF-8'); assert(fid>=0);
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end
