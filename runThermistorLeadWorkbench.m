function [current,report]=runThermistorLeadWorkbench(input,cfg,progress)
% 模块一：绑定已确认的米制原面、求解体网、冻结热边界与所选数据来源。
root=fileparts(mfilename('fullpath'));p=input.thermistorPrepared;
contract=jsondecode(fileread(fullfile(root,'cases','thermistor_b','leads3.json')));
assert(isequal(cfg.thermistorBoundary,thermistorWorkbenchBoundary()) && cfg.thermistorAngularCells==4, ...
    'RAH:Lead:Boundary','带线工况使用 3 m/s 空气和 305.15 K 给定壁温。');
saved=jsondecode(fileread(p.manifest));
assert(isequal(saved,rmfield(p,{'manifest','manifestSHA256'})) && strcmp(cm011FileSHA256(p.manifest),p.manifestSHA256) && ...
    strcmp(p.meshSHA256,contract.meshSHA256) && strcmp(p.bodySHA256,contract.bodySHA256) && ...
    strcmp(cm011FileSHA256(fullfile(p.directory,'mesh.su2')),contract.meshSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(p.directory,'nodes.npz')),contract.nodeMapSHA256), ...
    'RAH:Lead:Mesh','已确认的体网或真实节点映射发生变化。');
mode='solve';if isfield(cfg,'thermistorResultSource'),mode=validatestring(cfg.thermistorResultSource,{'archive','solve'});end
if strcmp(mode,'archive')
    [folder,guard]=openThermistorArchive('leads3'); %#ok<NASGU>
    geometry=folder;runFolder='';
else
    if isfield(cfg,'thermistorRunDirectory'),runFolder=cfg.thermistorRunDirectory;
    else,runFolder=fullfile(root,'outputs',['thermistor_leads_run_' char(java.util.UUID.randomUUID())]);end
    logFile=[runFolder '.log'];timer=tic;
    runCM011Python({fullfile(root,'scripts','thermistorLeadWorkbench.py'),'calculate',runFolder, ...
        '--prepared',p.directory},logFile,@pulse,true);
    folder=fullfile(runFolder,'surface');geometry=runFolder;
end

% 模块二：按七个暴露区域接入同一逐面台账，独立积分面元功率并核对温度身份。
s=jsondecode(fileread(fullfile(folder,'summary.json')));bodyFile=fullfile(geometry,'body.json');
c=jsondecode(fileread(fullfile(geometry,'case.json')));body=jsondecode(fileread(bodyFile));
wallFile=fullfile(folder,'wall.csv');rows=readtable(wallFile);
assert(strcmp(s.wallSHA256,cm011FileSHA256(wallFile)) && strcmp(s.meshSHA256,p.meshSHA256) && ...
    strcmp(s.bodySHA256,p.bodySHA256) && strcmp(s.bodySHA256,cm011FileSHA256(bodyFile)) && ...
    strcmp(s.caseSHA256,cm011FileSHA256(fullfile(geometry,'case.json'))) && ...
    isequal(input.V,body.V) && isequal(input.F,body.F),'RAH:Lead:Source','原面表与本次模型、体网或工况不符。');
if ~isempty(runFolder)
    assert(strcmp(s.sourceFrameSHA256,cm011FileSHA256(fullfile(runFolder,'restart.csv'))) && ...
        strcmp(s.sourceReactionSHA256,cm011FileSHA256(fullfile(runFolder,'reaction0','restart.csv'))), ...
        'RAH:Lead:Source','壁面结果与本次保存帧不符。');
end
n=size(body.F,1);f=body.F;v=body.V;
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(c.speedMps==3 && c.airPa==101325 && c.airK==293.15 && c.wallK==305.15 && ...
    isequal(rows.faceID,(0:n-1)') && isequal(rows.regionID,body.regionID) && ...
    all(abs(rows.areaM2-area)<1e-10*area) && all(rows.pressurePa>0) && ...
    all(rows.nearWallCellTemperatureK>=c.airK & rows.nearWallCellTemperatureK<=c.wallK) && ...
    abs(sum(area.*rows.heatFluxWm2)+s.totalOutwardW)<1e-12 && abs(s.energy.relativeDifference)<.005 && ...
    s.relativeChange20<.002,'RAH:Lead:Wall','边界、面积、功率或热稳定性核对失败。');
for j=1:numel(s.regions)
    ix=rows.regionID==j;
    assert(nnz(ix)==s.regions(j).faces && abs(sum(area(ix).*rows.heatFluxWm2(ix))+s.regions(j).outwardW)<1e-12);
end
labels=string(body.regionNames);
condition=struct('caseLabel','带引线 B 型热敏珠','velocityMps',[-3 0 0],'ambientPressurePa',c.airPa, ...
    'ambientTemperatureK',c.airK,'wallTemperatureK',c.wallK,'alphaDeg',0,'betaDeg',0, ...
    'geometryScope','两端带线名义试件，30° 对称稳态，完整暴露面含接头及端面');
setup=struct('condition',condition,'region',labels(rows.regionID),'wallTemperatureK',repmat(c.wallK,n,1), ...
    'engineeringSolver',@emptyWall,'localSolver',@wallData);
[current,~,report]=runAeroHeatSingleCase(struct('F',f,'V',v),struct('analysisMode','hybridSurface','hybridSurface',setup));
assert(report.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW+s.totalOutwardW)<1e-12);
current.hybridWall.pressureGaugePa=rows.pressureGaugePa;current.thermistorLead=s;
current.evidenceIdentity.sourceFrame=struct('kind','SU2-8.6.0-thermistor-B-leads-3ms','sha256',s.sourceFrameSHA256, ...
    'reactionSHA256',s.sourceReactionSHA256,'bodySHA256',s.bodySHA256,'wallSHA256',s.wallSHA256);
current.evidenceIdentity.meshSHA256=p.meshSHA256;current.evidenceIdentity.preparedSHA256=p.manifestSHA256;
current.evidenceIdentity.runDirectory=runFolder;current.evidenceIdentity.reusedArchive=strcmp(mode,'archive');
current.evidenceIdentity.sourceTraceable=true;current.evidenceIdentity.productionQualified=false;
current.evidenceIdentity.scope='轴向来流下两端带线等温试件的对称稳态分支，实验对照为珠体平均换热。';
current=attachThermistorLeadNearWall(current,runFolder);
report.outputDirectory=runFolder;report.regions=s.regions;

    function pulse()
        if ~isempty(progress),progress(struct('phase','thermistorLead','status','running','fraction',NaN, ...
                'message','','elapsedSeconds',toc(timer),'errorIdentifier','','errorMessage',''));end
    end
    function w=emptyWall(~,~,~)
        w=rows([],{'faceID','pressurePa','heatFluxWm2'});w.method=strings(0,1);
    end
    function result=wallData(~,~,~,~,target,~)
        w=rows(target+1,{'faceID','pressurePa','heatFluxWm2','nearWallCellTemperatureK'});m=height(w);
        w.method=repmat("SU2 8.6.0 三维对称稳态不可压层流",m,1);
        w.pressureMethod=repmat("壁节点表压按三角形平均，加来流静压",m,1);
        w.heatFluxMethod=repmat("强温差标量的壁节点反力按邻接面积守恒分配，正热流入壁",m,1);
        w.nearWallTemperatureMethod=repmat("对应三个壁节点的首层流体温差标量换算为 K 后平均",m,1);
        result=struct('accepted',true,'reason','六项残差、温界、独立壁反力、热稳定性与远场收支通过', ...
            'wall',w,'diagnostics',struct('sourceFrameSHA256',s.sourceFrameSHA256));
    end
end
