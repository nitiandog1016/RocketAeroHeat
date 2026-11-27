function [result,c,report]=runSimpleLauncherWindTunnel(mesh,cfg,progressFcn)
% 模块一：将当前三角形及来流与固定风洞试件逐面绑定，保留实际原面顺序。
if nargin<3,progressFcn=[];end
root=fileparts(mfilename('fullpath'));folder=fullfile(root,'cases','simple_launcher');
wake=isfield(cfg.launcherBoundary,'elapsedPhysicalTimeS');
mode='steady';if wake,mode='wake';end
c=simpleLauncherBoundary(mode);assert(isequaln(cfg.launcherBoundary,c),'RAH:Launcher:Boundary','风洞工况不一致。');
def=jsondecode(fileread(fullfile(folder,'geometry.json')));
input=struct('geometry',struct('stlFile',fullfile(folder,def.file),'scale',1));
assert(strcmp(cm011FileSHA256(input.geometry.stlFile),def.fileSHA256) && strcmp(def.fileSHA256,c.surfaceSHA256));
known=loadSTLMesh(input,true);n=size(known.F,1);
assert(size(mesh.F,1)==n && n==def.faces,'RAH:Launcher:Mesh','应使用当前完整原面网格。');
for j=1:3
    assert(max(abs(mesh.V(mesh.F(:,j),:)-known.V(known.F(:,j),:)),[],'all')<1e-10, ...
        'RAH:Launcher:Mesh','风洞模型的面序、单位或朝向有变。');
end
% 模块二：默认从当前模型独立新算。明确读取归档时，核对发布包、工况和壁面表指纹。
job=fullfile(root,'outputs',['launcher_' char(java.util.UUID.randomUUID())]);
action='solve';if isfield(cfg,'launcherAction'),action=validatestring(cfg.launcherAction,{'solve','archive'});end
seal=struct;
caseFolder=folder;script='runSimpleLauncherCase.py';
if wake,caseFolder=fullfile(folder,'wake');script='runSimpleLauncherWakeCase.py';end
receiptFile=fullfile(caseFolder,'receipt.json');
if isfile(receiptFile),seal=jsondecode(fileread(receiptFile));end
if strcmp(action,'archive')
    zip=fullfile(caseFolder,'wind_tunnel.zip');
    assert(isfield(seal,'sha256') && strcmp(cm011FileSHA256(zip),seal.sha256),'RAH:Launcher:Archive','归档指纹不符。');
    mkdir(job);
else
    logPath=[job '.log'];
    runCM011Python({fullfile(root,'scripts',script),job},logPath,@heartbeat);
    zip=fullfile(job,'result','wind_tunnel.zip');
end
data=fullfile(job,'wall');unzip(zip,data);
s=jsondecode(fileread(fullfile(data,'summary.json')));s.condition.velocityMps=reshape(s.condition.velocityMps,1,3);
assert(s.accepted && isequaln(s.condition,c) && strcmp(s.surfaceSHA256,def.fileSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(data,'wall.csv')),s.wallSHA256),'RAH:Launcher:Snapshot','工况或面表指纹不符。');
% 模块三：从原面面积独立复算整件功率，科研工作台按变量显示各自的来源。
w=readtable(fullfile(data,'wall.csv'));assert(isequal(w.faceID,(1:n)') && height(w)==n);
f=mesh.F;v=mesh.V;area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(max(abs(area-w.areaM2))<1e-10 && all(w.pressurePa>0) && all(isfinite(w.nearWallK)));
report=summarizeComputedWall(area,isfinite(w.pressurePa),isfinite(w.heatFluxWm2),w.heatFluxWm2);
assert(report.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW+s.outwardW)<1e-10 && ...
    abs(report.totalArea-s.areaM2)<1e-10,'RAH:Launcher:Coverage','当前原面覆盖或功率不符。');
regions=strings(n,1);
names=struct('sphericalNose','球形鼻端','cone','锥形头罩','cylinder','圆柱箭身', ...
    'contraction','收缩肩','afterbody','尾段','base','平底');
for g=def.regions',regions(g.firstFace:g.lastFace)=names.(g.name);end
assert(all(strlength(regions)>0));
% 模块四：按每个原面实际覆盖的范围标注读数，另列同一终帧的平底径向静压。
% 壁面表始终保存面积平均值。原生流场中的测点值只用于探针说明。
w.sampleScope=repmat("原三角面面积平均。",n,1);
nose=regions==names.sphericalNose;
axialDistanceMm=(max(v(:,1))-reshape(v(f(nose,:),1),[],3))*1000;
w.sampleScope(nose)=compose("鼻端面平均；距鼻尖的轴向距离 %.3f 至 %.3f mm。", ...
    min(axialDistanceMm,[],2),max(axialDistanceMm,[],2));
base=regions==names.base;
baseVertices=v(unique(f(base,:)),:);
radiusMm=1000*max(hypot(baseVertices(:,2),baseVertices(:,3)));
w.sampleScope(base)=sprintf('平底扇片面平均；径向范围为底心至外缘（0 至 %.1f mm）。',radiusMm);
w.pressureSampleScope=w.sampleScope;
studyFile=fullfile(folder,'wind_tunnel_stations.json');
if isfile(studyFile) && isfield(seal,'stationsSHA256') && strcmp(cm011FileSHA256(studyFile),seal.stationsSHA256)
    study=jsondecode(fileread(studyFile));grid=study.grids.outer;
    if strcmp(study.surfaceSHA256,s.surfaceSHA256) && strcmp(grid.wallSHA256,s.wallSHA256) && ...
            strcmp(grid.restartSHA256,s.runs.outer.restartSHA256)
        stations=compose('%.3f R 处 %+.5f',study.tapRadiusOverBaseRadius(:),grid.tapGaugeMilliPa(:));
        w.pressureSampleScope(base)=w.sampleScope(base)+" 同一流场的径向点值（相对环境，mPa）为 "+ ...
            strjoin(stations,'，')+"。R 为底面半径。";
    end
end
w.faceID=w.faceID-1;w.pressureGaugePa=w.pressurePa-c.ambientPressurePa;
w.method=repmat("固定低速风洞数值解",n,1);
w.pressureMethod=repmat("SU2 8.6.0 不可压恒密度层流壁压",n,1);
w.heatFluxMethod=repmat("轴对称稳态显热有限元、守恒上风通量及等温壁反力",n,1);
w.nearWallCellTemperatureK=w.nearWallK;
w.nearWallTemperatureMethod=repmat("显热场插值至第一层流体节点后按原面面积平均",n,1);
% 瞬态快照的流动、热场、时间步及储热量有独立身份；两半步末值作为本次主值。
kind='fixed-low-speed-wind-tunnel';
if wake
    flowFrame=s.sourceFrameSHA256;kind='fixed-low-speed-wake-bdf1';
    w.pressureMethod(:)="SU2 8.6.0 三维恒密度无滑移动量，BDF1 两半步末壁压";
    w.heatFluxMethod(:)="三维守恒质量通量、周期周向输运、BDF1 流体储热和强等温壁反力";
    w.nearWallTemperatureMethod(:)="同一末时刻的三维显热场插值至首层流体节点后映射到原面";
else
    flowFrame=s.runs.outer.restartSHA256;
end
result=struct('V',v,'F',f,'localGeometry',struct('region',{cellstr(regions)}), ...
    'hybridWall',w,'hybridCondition',c,'hybridWallTemperatureK',repmat(c.wallTemperatureK,n,1), ...
    'launcherSummary',s,'launcherAction',action,'evidenceIdentity',struct('origin','research','physicalCondition',c, ...
    'sourceTraceable',true,'sourceFrame',struct('kind',kind, ...
    'wallSHA256',s.wallSHA256,'flowSHA256',flowFrame,'stlSHA256',def.fileSHA256), ...
    'scope',c.scope,'productionQualified',false));
if wake
    result.evidenceIdentity.sourceFrame.thermalSHA256=s.thermalFrameSHA256;
    result.evidenceIdentity.sourceFrame.initialFlowSHA256=s.initialFlowSHA256;
    result.evidenceIdentity.sourceFrame.initialThermalSHA256=s.initialThermalSHA256;
end
% 模块五：按本次速度终帧和面表绑定全周热网复核。旧归档仍然可读。
comparisonFolder=fileparts(zip);
comparisonSeal=jsondecode(fileread(fullfile(comparisonFolder,'receipt.json')));
if isfield(comparisonSeal,'angularHeatSHA256')
    pack=fullfile(comparisonFolder,'angular_heat.zip');
    assert(isfile(pack) && strcmp(cm011FileSHA256(pack),comparisonSeal.angularHeatSHA256), ...
        'RAH:Launcher:AngularArchive','全周复核包与本次工况回执不匹配。');
    sub=fullfile(job,'angular_wall');unzip(pack,sub);
    a=jsondecode(fileread(fullfile(sub,'comparison.json')));
    a.condition.velocityMps=reshape(a.condition.velocityMps,1,3);
    path=fullfile(sub,'angular_wall.csv');
    assert(a.accepted && a.sectors==48 && a.faces==n && isequaln(a.condition,c) && ...
        strcmp(a.surfaceSHA256,def.fileSHA256) && strcmp(a.flowSHA256,s.runs.outer.restartSHA256) && ...
        strcmp(a.primaryWallSHA256,s.wallSHA256) && strcmp(a.csvSHA256,cm011FileSHA256(path)) && ...
        a.energyImbalance<0.005 && a.relativeFreeRowResidual<1e-8,'RAH:Launcher:AngularSource', ...
        '全周复核的几何、工况、流场或热量收支不一致。');
    tab=readtable(path);
    assert(isequal(tab.faceID,(1:n)') && max(abs(tab.areaM2-area))<1e-10 && ...
        all(isfinite(tab.heatFluxWm2)) && abs(area'*tab.heatFluxWm2+a.outwardW)<1e-12, ...
        'RAH:Launcher:AngularFaces','全周复核面表或积分功率不一致。');
    a.heatFluxWm2=tab.heatFluxWm2;
    result.launcherAngularHeat=a;
end
report.status='completed';report.runDirectory=job;
    function heartbeat()
        if ~isempty(progressFcn),progressFcn(struct('status','running','message','低速风洞计算中'));end
        drawnow limitrate
    end
end
