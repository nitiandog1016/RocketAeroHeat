function [current,report,page]=runThermistorBeadCase(folder,mode,windowOpen,inputMesh,thermalInput)
%RUNTHERMISTORBEADCASE 读取 B 型热敏珠的冻结结果，或从冷初场计算同一工况。
% [mesh,report]=runThermistorBeadCase() 只读随项目发布的原面表。
% [mesh,report,page]=runThermistorBeadCase([], 'solve') 运行内嵌求解器并导出离线页。
% 'archive7p5' 和 'solve7p5' 使用 7.5° 周向单元，并与 15° 归档按面配对。
% 'startup6' 读取 6 m/s 完整周向域的首个已接受时刻，采用该时刻的固定热边界。
% 'timeSteps6' 在相同物理时刻加入两个半步的全表面结果及逐面差值。
% 'nearWall6' 读取两个半步后的 56 层全表面结果，与 48 层同刻结果对照。
% runThermistorBeadCase([], 'archive', [], [], 299.15) 使用参照场计算给定壁温的稳态结果。
% 第五参数也可为 struct('netHeatingPowerW',0.004)，按整珠净功率反求稳态等温壁温。
% 显式指定 outputs 下的新目录，或请求第三个输出，都会在那里导出离线页。
root=fileparts(mfilename('fullpath'));outputs=fullfile(root,'outputs');
if nargin<1,folder='';end
if nargin<2,mode='archive';end
if nargin<3,windowOpen=[];end
mode=validatestring(mode,{'archive','solve','archive7p5','solve7p5','startup6','timeSteps6','nearWall6'});
fresh=startsWith(mode,'solve');refined=endsWith(mode,'7p5');
angularCells=2+2*refined;
if nargin<4,inputMesh=[];end
startup=any(strcmp(mode,{'startup6','timeSteps6','nearWall6'}));
if startup
    assert(isempty(inputMesh) && nargin<5,'RAH:Bead:Startup','6 m/s 快照使用绑定的完整体网和 303.15 K 给定壁温。');
else
    boundary=thermistorWorkbenchBoundary();
    if nargin>=5
        if isstruct(thermalInput)
            if isscalar(thermalInput) && isequal(fieldnames(thermalInput),{'netHeatingPowerW'})
                boundary=rmfield(boundary,'wallTemperatureK');
                boundary.netHeatingPowerW=thermalInput.netHeatingPowerW;
            elseif isscalar(thermalInput) && isequal(fieldnames(thermalInput),{'wallTemperatureK'})
                boundary.wallTemperatureK=thermalInput.wallTemperatureK;
            else
                boundary=thermalInput;
            end
        else
            boundary.wallTemperatureK=thermalInput;
        end
    end
    boundary=thermistorWorkbenchBoundary(boundary,angularCells);
    powerMode=isfield(boundary,'netHeatingPowerW');
    if ~isempty(inputMesh)
        % 将显式计算或归档绑定到主侧栏已确认的米制原面与实际体网。
        verifyThermistorWorkbenchMesh(inputMesh,angularCells);
    end
end
writePage=fresh || ~isempty(folder) || nargout>=3;page='';
if writePage
    if isempty(folder),folder=fullfile(outputs,['thermistor_' char(java.util.UUID.randomUUID())]);end
    targetDir=java.io.File(char(folder));
    if ~targetDir.isAbsolute(),targetDir=java.io.File(pwd,char(folder));end
    folder=char(targetDir.getCanonicalPath());
    assert(startsWith(folder,[char(java.io.File(outputs).getCanonicalPath()) filesep],'IgnoreCase',ispc) && ...
        java.io.File(folder).mkdirs(),'RAH:Bead:Output','请选择项目 outputs 下尚未使用的目录。');
end

% 已接受的全周启动步按自身的边界与能量账读取，并从同一入口导出三维工作台。
if startup
    if strcmp(mode,'timeSteps6')
        [current,report]=readThermistorStartup();
        [current,report]=attachThermistorStepComparison(current,report);
    else
        [current,report]=readThermistorStartup(mode);
    end
    if writePage,page=exportAerothermalWeb(current,fullfile(folder,'workbench'));end
    return;
end
% 模块一：新计算经过同网格三维动量、强温差和反力核对；浏览时校验归档后再解包。
if ~fresh || refined
    [baseline,baseCleanup]=openThermistorArchive(2); %#ok<NASGU>
end
if fresh
    native=fullfile(folder,'native');
    pulse=[];
    if ~isempty(windowOpen),pulse=@checkWindow;end
    args={fullfile(root,'scripts','runThermistorCase.py'),native};
    if refined,args=[args,{'--angular-cells','4'}];end
    if ~isempty(inputMesh)
        p=inputMesh.thermistorPrepared;
        args=[args,{'--prepared',p.directory,'--prepared-sha',p.manifestSHA256}];
    end
    runCM011Python(args,fullfile(folder,'solve.log'),pulse);
    geo=fullfile(native,'strong');result=fullfile(native,'result');
else
    geo=baseline;result=baseline;
    if refined
        [result,refinedCleanup]=openThermistorArchive(4); %#ok<NASGU>
    end
end
spec=jsondecode(fileread(fullfile(geo,'case.json')));
stats=jsondecode(fileread(fullfile(result,'summary.json')));
rowsFile=fullfile(result,'wall.csv');bodyFile=fullfile(geo,'body.json');
assert(strcmp(spec.grid,'outer') && strcmp(stats.grid,'outer') && ...
    strcmp(cm011FileSHA256(fullfile(geo,'case.json')),'46d981fce8d8b3f456f0a36cfde7ccd84cd47cb721d2d71a31ceea87fded038c') && ...
    strcmp(cm011FileSHA256(bodyFile),stats.bodySHA256) && strcmp(cm011FileSHA256(rowsFile),stats.wallTableSHA256), ...
    'RAH:Bead:Binding','工况、原面表或三角面与数值回执不符。');
body=jsondecode(fileread(bodyFile));mesh=struct('V',body.V,'F',body.F);
rows=readtable(rowsFile);n=size(mesh.F,1);
assert(stats.fullSurfaceComputed && stats.originalFaces==n && n==5216 && ...
    isequal(rows.faceID,(0:n-1)') && stats.boundaryRelativeDifference<.005 && stats.heatChange20<.002 && ...
    all(structfun(@(v)v < -8,stats.last20MaxLogRMS)) && all(structfun(@(v)v < -8,stats.continued20MaxLogRMS)), ...
    'RAH:Bead:Receipt','原面编号、热量账或收敛窗未通过。');
assert(all(rows.pressurePa>0 & rows.heatFluxWm2<0) && all(rows.nearWallCellTemperatureK>spec.ambientK & ...
    rows.nearWallCellTemperatureK<spec.wallK) && max(abs(rows.pressurePa-spec.pressurePa-rows.pressureGaugePa))<2e-10, ...
    'RAH:Bead:Fields','压力、热流或首层温度读数不符。');
f=mesh.F;v=mesh.V;area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(all(abs(rows.areaM2-area)<1e-10*area),'RAH:Bead:Area','逐面面积与原几何不符。');

% 恒物性、无源被动温差方程在固定速度场中是线性的，压力和面积直接沿用。
% 原归档的温度界限与热反力先按 12 K 参照工况审核，再换算到给定边界。
referenceDelta=spec.wallK-spec.ambientK;
referenceRelease=-stats.powerW;
conductance=referenceRelease/referenceDelta;
if powerMode
    % 只把扣除引线与辐射等损失后的净功率交给空气侧；采用本次所选场的实际热导。
    assert(boundary.netHeatingPowerW<=referenceRelease,'RAH:Bead:Power','净功率超过本次参照场的范围。');
    scale=boundary.netHeatingPowerW/referenceRelease;
    deltaT=referenceDelta*scale;
    wallTemperatureK=spec.ambientK+deltaT;
else
    wallTemperatureK=boundary.wallTemperatureK;
    deltaT=wallTemperatureK-spec.ambientK;
    scale=deltaT/referenceDelta;
end
if scale==0
    rows.heatFluxWm2(:)=0;
    rows.nearWallCellTemperatureK(:)=spec.ambientK;
elseif scale~=1
    rows.heatFluxWm2=scale*rows.heatFluxWm2;
    rows.nearWallCellTemperatureK=spec.ambientK+scale*(rows.nearWallCellTemperatureK-spec.ambientK);
end

% 模块二：同一名义等温珠体的全部原面交给区域求解入口，保留每个变量的算法来源。
condition=struct('caseLabel','B 型热敏珠 · 名义长椭球 · 沿长轴来流','velocityMps',[-spec.freestreamMps 0 0], ...
    'ambientPressurePa',spec.pressurePa,'ambientTemperatureK',spec.ambientK,'wallTemperatureK',wallTemperatureK, ...
    'alphaDeg',0,'betaDeg',0,'diameterM',spec.diameterM,'geometryScope','光滑等温珠体，未建引线');
setup=struct('condition',condition,'region',repmat("热敏珠完整外表面",n,1), ...
    'wallTemperatureK',repmat(wallTemperatureK,n,1),'engineeringSolver',@emptyWall,'localSolver',@localWall);
if powerMode
    setup.condition=rmfield(condition,'wallTemperatureK');
    setup.condition.netHeatingPowerW=boundary.netHeatingPowerW;
end
[current,~,report]=runAeroHeatSingleCase(mesh,struct('analysisMode','hybridSurface','hybridSurface',setup));
targetPower=scale*stats.powerW;
if powerMode
    targetPower=-boundary.netHeatingPowerW;
    current.hybridSolvedWallTemperatureK=current.hybridWallTemperatureK;
    current=rmfield(current,'hybridWallTemperatureK');
end
assert(report.fullSurfaceComputed && isequal(current.hybridWall.faceID,rows.faceID) && ...
    abs(report.wholeSurfaceHeatPowerW-targetPower)<=1e-11*abs(targetPower),'RAH:Bead:Power','完整原面功率与给定壁温的反力不符。');
current.hybridWall.pressureGaugePa=rows.pressureGaugePa;
identity=struct('kind','SU2-8.6.0-thermistor-B','sha256',stats.sourceFrameSHA256, ...
    'reactionSHA256',stats.sourceReactionSHA256,'wallSHA256',stats.wallTableSHA256,'bodySHA256',stats.bodySHA256);
current.evidenceIdentity.sourceFrame=identity;current.evidenceIdentity.sourceTraceable=true;
current.evidenceIdentity.productionQualified=false;
current.evidenceIdentity.scope='名义等温珠体的数值解。实验只给整珠平均换热系数。';
thermal=struct('referenceWallTemperatureK',spec.wallK,'wallTemperatureK',wallTemperatureK, ...
    'ambientTemperatureK',spec.ambientK,'deltaTemperatureK',deltaT,'scale',scale, ...
    'densityKgM3',spec.rho,'viscosityPaS',spec.mu,'specificHeatJkgK',spec.cp,'conductivityWmK',spec.k, ...
    'nativeFrameSHA256',stats.sourceFrameSHA256,'nativeReactionSHA256',stats.sourceReactionSHA256);
if powerMode
    thermal.netHeatingPowerW=boundary.netHeatingPowerW;
    thermal.conductanceWK=conductance;
    thermal.referenceReleaseW=referenceRelease;
    thermal.surfaceReleaseW=-sum(area.*current.hybridWall.heatFluxWm2);
    thermal.balanceErrorW=thermal.surfaceReleaseW-boundary.netHeatingPowerW;
    thermal.wallMethod='稳态等温集总热平衡：环境温度加净加热功率除以同网格空气侧热导';
    thermal.scope='名义等温珠体，恒物性强迫对流。输入功率已扣除引线、辐射等其他损失。';
    assert(abs(thermal.balanceErrorW)<=1e-11*boundary.netHeatingPowerW,'RAH:Bead:Power','逐面放热积分与输入净功率不符。');
end
if scale~=1 || powerMode
    current.evidenceIdentity.sourceFrame.thermalConversion=thermal;
end

% 模块三：按论文的名义面积比较总功率；逐面验证面积、给定温度与首层温度独立保存。
meanH=NaN;
if deltaT>0
    meanH=-report.wholeSurfaceHeatPowerW/(deltaT*spec.areaSourceM2);
    assert(abs(meanH/stats.nominalExperimentAreaCoefficientWm2K-1)<1e-11,'RAH:Bead:Mean','名义面积换热积分不符。');
end
current.thermistorBead=struct('meanH',meanH,'verifiedJointAreaM2',0,'thermal',thermal);
if scale==1 && ~powerMode
    comparison=jsondecode(fileread(fullfile(root,'docs','validation','v282','evidence','comparison.json')));
    current.thermistorBead.experiment=comparison.experiment;
    current.thermistorBead.relativeDifference=meanH/comparison.experiment.hWm2K-1;
end
current.thermistorBead.cellAngleDeg=15;
if refined
    % 模块四：根据两份各自绑定的场计算同面差值，主热流和主压力保留细网原值。
    old=jsondecode(fileread(fullfile(baseline,'summary.json')));
    oldWall=readtable(fullfile(baseline,'wall.csv'));
    assert(stats.angularCells==4 && stats.cellAngleDeg==7.5 && ...
        strcmp(stats.meshSHA256,'f59482ab24bcdc3d886fac3e442ff474a13c2d34c0596c822ebfa2e7cb81e723') && ...
        strcmp(stats.bodySHA256,old.bodySHA256) && strcmp(stats.configSHA256,old.configSHA256) && ...
        isequal(oldWall.faceID,rows.faceID) && isequal(oldWall.areaM2,rows.areaM2), ...
        'RAH:Bead:Angle','周向网格、配置或原面身份与对照不符。');
    oldScale=scale;
    if powerMode,oldScale=boundary.netHeatingPowerW/(-old.powerW);end
    if oldScale<=1
        current.beadAngularDelta=struct('pressurePa',rows.pressureGaugePa-oldWall.pressureGaugePa, ...
            'heatFluxWm2',rows.heatFluxWm2-oldScale*oldWall.heatFluxWm2);
        if powerMode
            current.thermistorBead.baselineSolvedWallTemperatureK=spec.ambientK+referenceDelta*oldScale;
        end
        if deltaT>0,current.thermistorBead.baselineH=old.nominalExperimentAreaCoefficientWm2K;end
    else
        current.thermistorBead.comparisonUnavailable='输入功率超过 15° 网的支持上限，同功率两网对照不可用。';
    end
    current.thermistorBead.cellAngleDeg=stats.cellAngleDeg;
    current.thermistorBead.baselineFrameSHA256=old.sourceFrameSHA256;
    current.evidenceIdentity.sourceFrame.comparisonFrameSHA256=old.sourceFrameSHA256;
end
if scale==1 && ~powerMode,report.native=stats;else,report.referenceNative=stats;end
report.thermal=thermal;report.reusedArchive=~fresh;report.comparison=current.thermistorBead;
if ~isempty(inputMesh)
    p=inputMesh.thermistorPrepared;
    assert(isequal(mesh.V,inputMesh.V) && isequal(mesh.F,inputMesh.F),'RAH:Bead:Geometry','返回原面与当前输入不符。');
    if fresh
        assert(strcmp(cm011FileSHA256(fullfile(geo,'mesh.su2')),p.meshSHA256),'RAH:Bead:Mesh','求解终场与确认体网不符。');
        report.outputDirectory=folder;
    end
    current.evidenceIdentity.meshSHA256=p.meshSHA256;
    current.evidenceIdentity.preparedSHA256=p.manifestSHA256;
    current.evidenceIdentity.reusedArchive=~fresh;
end
if writePage,page=exportAerothermalWeb(current,fullfile(folder,'workbench'));end

    function checkWindow()
        % 关闭主窗口后通知当前原生监督器收回本次计算，其他工况继续独立运行。
        drawnow limitrate;
        try
            if ~windowOpen(),error('RAH:Bead:Closed','热敏珠窗口已关闭。');end
        catch ex
            fid=fopen(fullfile(folder,'thermistor_cancel.request'),'w');
            if fid>=0,fclose(fid);end
            rethrow(ex);
        end
    end
    function empty=emptyWall(~,~,~)
        % 对这个热敏珠，采用同一连通外域的数值解覆盖全部表面。
        empty=rows([],{'faceID','pressurePa','heatFluxWm2'});empty.method=strings(0,1);
    end
    function solution=localWall(~,~,~,~,target,~)
        wall=rows(target+1,{'faceID','pressurePa','heatFluxWm2','nearWallCellTemperatureK'});
        count=height(wall);
        wall.method=repmat("SU2 8.6.0 三维层流，30° 对称扇区",count,1);
        wall.pressureMethod=repmat("原生壁节点相对压力按面内权重平均，加环境静压",count,1);
        wall.heatFluxMethod=repmat("强等温对流扩散方程的壁反力，按原三角面积守恒分配",count,1);
        wall.nearWallTemperatureMethod=repmat("首层流体节点的温差标量换算为 K，再按面内权重平均",count,1);
        if scale~=1
            wall.heatFluxMethod(:)="恒物性线性换算：已收敛的强等温壁反力乘以当前温差与参照温差之比，原面守恒分配";
            wall.nearWallTemperatureMethod(:)="环境温度加上当前壁温差乘以原面对应的首层无量纲温差";
        end
        solution=struct('accepted',true,'reason','六个分量的末窗残差、温界与边界热量账通过', ...
            'wall',wall,'diagnostics',struct('sourceFrameSHA256',stats.sourceFrameSHA256));
    end
end
