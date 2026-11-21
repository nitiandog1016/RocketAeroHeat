function [current,report]=runThermistorSidewindWorkbench(input,cfg,progress)
% 模块一：绑定实际全周体网、试件原面及确认过的 3 m/s 侧吹边界。
root=fileparts(mfilename('fullpath'));p=input.thermistorPrepared;
b=thermistorWorkbenchBoundary(cfg.thermistorBoundary,cfg.thermistorAngularCells,cfg);
key=jsondecode(fileread(fullfile(root,'cases','thermistor_b','sidewind3.json')));
assert(strcmp(cm011FileSHA256(p.manifest),p.manifestSHA256) && isequal(jsondecode(fileread(p.manifest)),rmfield(p,{'manifest','manifestSHA256'})));
assert(p.angularCells==48 && strcmp(cm011FileSHA256(fullfile(p.directory,'mesh.su2')),key.meshSHA256) && ...
    strcmp(p.meshSHA256,key.meshSHA256) && strcmp(cm011FileSHA256(fullfile(p.directory,'nodes.npz')),key.nodeMapSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(p.directory,'mesh_audit.json')),p.meshAuditSHA256),'RAH:Lead:Mesh','侧吹体网、原面映射或验收回执已改变。');
source='solve';if isfield(cfg,'thermistorResultSource'),source=validatestring(cfg.thermistorResultSource,{'archive','solve'});end
runFolder='';
if strcmp(source,'archive')
    [folder,cleanup]=openThermistorArchive('sidewind3'); %#ok<NASGU>
    geometry=folder;
else
    if isfield(cfg,'thermistorRunDirectory'),runFolder=cfg.thermistorRunDirectory;
    else,runFolder=fullfile(root,'outputs',['thermistor_sidewind_run_' char(java.util.UUID.randomUUID())]);end
    timer=tic;
    runCM011Python({fullfile(root,'scripts','thermistorSidewindWorkbench.py'),'calculate',runFolder,'--prepared',p.directory}, ...
        [runFolder '.log'],@pulse,true);
    done=jsondecode(fileread(fullfile(runFolder,'completed.json')));
    assert(ismember(string(done.stage),["steady01","steady02"]));
    geometry=fullfile(runFolder,done.stage);folder=fullfile(geometry,'surface');
end
% 模块二：复核边界能量账、源帧和全部七区的面元面积分，再写入混合计算结果。
c=jsondecode(fileread(fullfile(geometry,'case.json')));body=jsondecode(fileread(fullfile(geometry,'body.json')));
s=jsondecode(fileread(fullfile(folder,'summary.json')));w=readtable(fullfile(folder,'wall.csv'));
assert(strcmp(c.caseKind,'thermistor-B-sidewind-two-leads') && isequal(c.velocityMps(:)',[0 -3 0]) && ...
    c.airK==b.airTemperatureK && c.airPa==b.airPressurePa && c.wallK==b.wallTemperatureK && c.sectorDeg==360);
assert(strcmp(cm011FileSHA256(fullfile(geometry,'case.json')),key.caseSHA256) && strcmp(s.caseSHA256,key.caseSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(geometry,'body.json')),p.bodySHA256) && strcmp(s.bodySHA256,p.bodySHA256) && ...
    strcmp(s.meshSHA256,p.meshSHA256) && strcmp(cm011FileSHA256(fullfile(folder,'wall.csv')),s.wallSHA256));
if ~isempty(runFolder)
    assert(strcmp(s.sourceFrameSHA256,cm011FileSHA256(fullfile(geometry,'restart.csv'))) && ...
        strcmp(s.sourceReactionSHA256,cm011FileSHA256(fullfile(geometry,'reaction0','restart.csv'))) && ...
        strcmp(done.wallSHA256,s.wallSHA256) && strcmp(done.frameSHA256,s.sourceFrameSHA256));
end
v=input.V;f=input.F;n=size(f,1);
a=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(isequal(body.F,f) && isequal(body.V,v) && isequal(w.faceID,(0:n-1)') && isequal(w.regionID,body.regionID) && ...
    n==10272 && s.faces==n && all(abs(a-w.areaM2)<1e-10*a) && all(w.pressurePa>0) && ...
    all(w.nearWallCellTemperatureK>=c.airK & w.nearWallCellTemperatureK<=c.wallK));
assert(all(structfun(@(x)x<-8,s.residualWindow)) && all(structfun(@(x)x<-8,s.plus20Window)) && ...
    s.relativeChange20<.002 && abs(s.energy.relativeDifference)<.005 && s.energy.domainCopies==1 && ...
    strcmp(s.energy.inputFrameSHA256,s.sourceFrameSHA256) && abs(a'*w.heatFluxWm2+s.energy.wallReleaseW)<1e-10);
for j=1:7
    use=w.regionID==j;assert(s.regions(j).faces==nnz(use) && abs(a(use)'*w.heatFluxWm2(use)+s.regions(j).outwardW)<1e-12);
end
condition=struct('caseLabel','带线 B 型热敏珠 · 全周侧吹','velocityMps',[0 -3 0],'ambientPressurePa',c.airPa, ...
    'ambientTemperatureK',c.airK,'wallTemperatureK',c.wallK,'alphaDeg',0,'betaDeg',90,'geometryScope','完整 360° 名义带线试件');
names=string(body.regionNames);
setup=struct('condition',condition,'region',names(body.regionID),'wallTemperatureK',c.wallK*ones(n,1), ...
    'engineeringSolver',@emptyWall,'localSolver',@wallData);
[current,~,report]=runAeroHeatSingleCase(struct('F',f,'V',v),struct('analysisMode','hybridSurface','hybridSurface',setup));
assert(report.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW+s.totalOutwardW)<1e-12);
current.thermistorSidewind=s;current.hybridWall.pressureGaugePa=w.pressureGaugePa;
current.hybridWall.pressureGaugeMethod=repmat("侧吹全周壁节点相对环境的表压，按原三角形平均并换算为 Pa",n,1);
current.evidenceIdentity.sourceFrame=struct('kind','SU2-8.6.0-thermistor-B-full-sidewind-3ms', ...
    'sha256',s.sourceFrameSHA256,'bodySHA256',s.bodySHA256,'reactionSHA256',s.sourceReactionSHA256,'wallSHA256',s.wallSHA256);
current.evidenceIdentity.meshSHA256=p.meshSHA256;current.evidenceIdentity.preparedSHA256=p.manifestSHA256;
current.evidenceIdentity.runDirectory=runFolder;current.evidenceIdentity.reusedArchive=strcmp(source,'archive');
current.evidenceIdentity.sourceTraceable=true;current.evidenceIdentity.productionQualified=false;
current.evidenceIdentity.scope='3 m/s、−Y 侧向来流，20D 外域的三维不可压层流稳态解，侧吹实验精度未知。';
report.outputDirectory=runFolder;report.regions=s.regions;
    function pulse()
        if ~isempty(progress),progress(struct('phase','thermistorSidewind','status','running','fraction',NaN, ...
                'message','','elapsedSeconds',toc(timer),'errorIdentifier','','errorMessage',''));end
    end
    function t=emptyWall(~,~,~)
        t=w([],{'faceID','pressurePa','heatFluxWm2'});t.method=strings(0,1);
    end
    function result=wallData(~,~,~,~,target,~)
        t=w(target+1,{'faceID','pressurePa','heatFluxWm2','nearWallCellTemperatureK'});count=height(t);
        t.method=repmat("SU2 8.6.0 全周侧吹不可压层流",count,1);
        t.pressureMethod=repmat("侧吹全周壁节点表压按面平均，加来流静压",count,1);
        t.heatFluxMethod=repmat("6 K 强温差壁节点反力按本方位邻接面积守恒分配，正热流入壁",count,1);
        t.nearWallTemperatureMethod=repmat("本方位三个首层流体节点的温差标量按 6 K 换算并取平均",count,1);
        result=struct('accepted',true,'reason','侧吹末窗、温界、独立反力、热稳定和远场能量账通过', ...
            'wall',t,'diagnostics',struct('sourceFrameSHA256',s.sourceFrameSHA256));
    end
end
