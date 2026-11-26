function [result,c,report]=runSimpleLauncherWindTunnel(mesh,cfg,progressFcn)
% 模块一：将当前三角形及来流与固定风洞试件逐面绑定，保留实际原面顺序。
if nargin<3,progressFcn=[];end
root=fileparts(mfilename('fullpath'));folder=fullfile(root,'cases','simple_launcher');
c=simpleLauncherBoundary();assert(isequaln(cfg.launcherBoundary,c),'RAH:Launcher:Boundary','风洞工况不一致。');
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
if strcmp(action,'archive')
    zip=fullfile(folder,'wind_tunnel.zip');seal=jsondecode(fileread(fullfile(folder,'receipt.json')));
    assert(strcmp(cm011FileSHA256(zip),seal.sha256),'RAH:Launcher:Archive','归档指纹不符。');
    mkdir(job);
else
    logPath=[job '.log'];
    runCM011Python({fullfile(root,'scripts','runSimpleLauncherCase.py'),job},logPath,@heartbeat);
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
w.faceID=w.faceID-1;w.pressureGaugePa=w.pressurePa-c.ambientPressurePa;
w.method=repmat("固定低速风洞数值解",n,1);
w.pressureMethod=repmat("SU2 8.6.0 不可压恒密度层流壁压",n,1);
w.heatFluxMethod=repmat("轴对称稳态显热有限元、守恒上风通量及等温壁反力",n,1);
w.nearWallCellTemperatureK=w.nearWallK;
w.nearWallTemperatureMethod=repmat("显热场插值至第一层流体节点后按原面面积平均",n,1);
result=struct('V',v,'F',f,'localGeometry',struct('region',{cellstr(regions)}), ...
    'hybridWall',w,'hybridCondition',c,'hybridWallTemperatureK',repmat(c.wallTemperatureK,n,1), ...
    'launcherSummary',s,'launcherAction',action,'evidenceIdentity',struct('origin','research','physicalCondition',c, ...
    'sourceTraceable',true,'sourceFrame',struct('kind','fixed-low-speed-wind-tunnel', ...
    'wallSHA256',s.wallSHA256,'flowSHA256',s.runs.outer.restartSHA256,'stlSHA256',def.fileSHA256), ...
    'scope',c.scope,'productionQualified',false));
report.status='completed';report.runDirectory=job;
    function heartbeat()
        if ~isempty(progressFcn),progressFcn(struct('status','running','message','低速风洞计算中'));end
        drawnow limitrate
    end
end
