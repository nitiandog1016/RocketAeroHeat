function [current,report]=runThermistorLeadWorkbench(input,cfg,progress)
% 模块一：绑定已确认的米制原面、求解体网、冻结热边界与所选数据来源。
root=fileparts(mfilename('fullpath'));p=input.thermistorPrepared;
midbody=strcmp(cfg.thermistorVariant,'midbody');options={};
refineLeads=isfield(cfg,'thermistorRefineLeads') && cfg.thermistorRefineLeads;
if midbody
    archive='midbody3';options={'--midbody'};caseLabel='光滑 B 型热敏珠 · 中段加密';
    scope='光滑名义珠体，30° 对称稳态，中段近壁过渡区局部加密';
else
    archive='leads3';caseLabel='带引线 B 型热敏珠';scope='两端带线名义试件，30° 对称稳态，完整暴露面含接头及端面';
    if refineLeads
        archive='leadsMidbody3';options={'--refine-leads'};caseLabel='带引线 B 型热敏珠 · 中段加密';
    end
end
contract=jsondecode(fileread(fullfile(root,'cases','thermistor_b',[archive '.json'])));
b=thermistorWorkbenchBoundary(cfg.thermistorBoundary,cfg.thermistorAngularCells,cfg);
powerMode=isfield(b,'netHeatingPowerW');
saved=jsondecode(fileread(p.manifest));
assert(isequal(saved,rmfield(p,{'manifest','manifestSHA256'})) && strcmp(cm011FileSHA256(p.manifest),p.manifestSHA256) && ...
    strcmp(p.meshSHA256,contract.meshSHA256) && strcmp(p.bodySHA256,contract.bodySHA256) && ...
    strcmp(cm011FileSHA256(fullfile(p.directory,'mesh.su2')),contract.meshSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(p.directory,'nodes.npz')),contract.nodeMapSHA256), ...
    'RAH:Lead:Mesh','已确认的体网或真实节点映射发生变化。');
mode='solve';if isfield(cfg,'thermistorResultSource'),mode=validatestring(cfg.thermistorResultSource,{'archive','solve'});end
if strcmp(mode,'archive')
    [folder,guard]=openThermistorArchive(archive); %#ok<NASGU>
    geometry=folder;runFolder='';
else
    if isfield(cfg,'thermistorRunDirectory'),runFolder=cfg.thermistorRunDirectory;
    else,runFolder=fullfile(root,'outputs',['thermistor_leads_run_' char(java.util.UUID.randomUUID())]);end
    logFile=[runFolder '.log'];timer=tic;
    runCM011Python({fullfile(root,'scripts','thermistorLeadWorkbench.py'),'calculate',runFolder, ...
        '--prepared',p.directory,options{:}},logFile,@pulse,true);
    folder=fullfile(runFolder,'surface');geometry=runFolder;
end

% 模块二：按实际暴露区域接入逐面台账，独立积分面元功率并核对温度身份。
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
% 模块三：净功率由整件七区共同承担，逐面换算后独立积分各区与总功率。
wallK=c.wallK;summary=s;
if powerMode
    [rows,scale,deltaT,wallK]=scaleThermistorWallField(rows,b,c.airK,c.wallK,s.totalOutwardW);
    regions=s.regions;
    for j=1:numel(regions)
        ix=rows.regionID==j;
        regions(j).outwardW=-sum(area(ix).*rows.heatFluxWm2(ix));
        regions(j).intoWallWm2=[min(rows.heatFluxWm2(ix));max(rows.heatFluxWm2(ix))];
        if scale==0,regions(j).peakFaceID=[];end
    end
    release=[regions.outwardW];
    thermal=struct('netHeatingPowerW',b.netHeatingPowerW,'referenceReleaseW',s.totalOutwardW, ...
        'conductanceWK',s.totalOutwardW/(c.wallK-c.airK),'ambientTemperatureK',c.airK,'referenceWallTemperatureK',c.wallK, ...
        'deltaTemperatureK',deltaT,'wallTemperatureK',wallK,'scale',scale,'surfaceReleaseW',sum(release), ...
        'balanceErrorW',sum(release)-b.netHeatingPowerW,'nativeFrameSHA256',s.sourceFrameSHA256,'nativeReactionSHA256',s.sourceReactionSHA256, ...
        'wallMethod','整件稳态等温集总热平衡，环境温度加净对流功率除以七区总热导', ...
        'scope','名义珠体和两根引线共用一个壁温。输入为全部暴露面传给空气的净对流功率。');
    assert(abs(thermal.balanceErrorW)<=1e-11*b.netHeatingPowerW,'RAH:Lead:Power','七区对流放热与整件输入功率不符。');
    summary=struct('thermal',thermal,'referenceNative',s,'regions',regions,'faces',n,'areaM2',sum(area), ...
        'totalOutwardW',sum(release),'beadAndJointsW',sum(release(3:5)), ...
        'wireSidesW',sum(release([2 6])),'wireEndsW',sum(release([1 7])),'jointExperimentVerifiedM2',0);
end
labels=string(body.regionNames);
condition=struct('caseLabel',caseLabel,'velocityMps',[-3 0 0],'ambientPressurePa',c.airPa, ...
    'ambientTemperatureK',c.airK,'wallTemperatureK',c.wallK,'alphaDeg',0,'betaDeg',0, ...
    'geometryScope',scope);
if powerMode
    condition=rmfield(condition,'wallTemperatureK');condition.netHeatingPowerW=b.netHeatingPowerW;
end
setup=struct('condition',condition,'region',labels(rows.regionID),'wallTemperatureK',repmat(wallK,n,1), ...
    'engineeringSolver',@emptyWall,'localSolver',@wallData);
[current,~,report]=runAeroHeatSingleCase(struct('F',f,'V',v),struct('analysisMode','hybridSurface','hybridSurface',setup));
assert(report.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW+summary.totalOutwardW)<1e-12);
if powerMode
    current.hybridSolvedWallTemperatureK=current.hybridWallTemperatureK;
    current=rmfield(current,'hybridWallTemperatureK');
    assert(abs(report.wholeSurfaceHeatPowerW+b.netHeatingPowerW)<=1e-11*b.netHeatingPowerW);
    report.thermal=thermal;
end
current.hybridWall.pressureGaugePa=rows.pressureGaugePa;
if midbody,current.thermistorMidbody=s;else,current.thermistorLead=summary;end
current.evidenceIdentity.sourceFrame=struct('kind',['SU2-8.6.0-thermistor-B-' cfg.thermistorVariant '-3ms'],'sha256',s.sourceFrameSHA256, ...
    'reactionSHA256',s.sourceReactionSHA256,'bodySHA256',s.bodySHA256,'wallSHA256',s.wallSHA256);
current.evidenceIdentity.meshSHA256=p.meshSHA256;current.evidenceIdentity.preparedSHA256=p.manifestSHA256;
current.evidenceIdentity.runDirectory=runFolder;current.evidenceIdentity.reusedArchive=strcmp(mode,'archive');
current.evidenceIdentity.sourceTraceable=true;current.evidenceIdentity.productionQualified=false;
current.evidenceIdentity.scope=[scope '。实验对照为珠体平均换热。'];
current=attachThermistorLeadNearWall(current,runFolder,folder);
if powerMode
    current.evidenceIdentity.sourceFrame.thermalConversion=thermal;
    rev=logical(current.thermistorLeadFlow.reverse);
    current.thermistorLeadFlow.summary.reverseReleaseW=-sum(area(rev).*rows.heatFluxWm2(rev));
elseif midbody,current=attachThermistorPairedControl(current,'midbody');
elseif refineLeads
    current=attachThermistorPairedControl(current,'wiredRefined');
    current=attachThermistorLeadGridComparison(current);
else,current=attachThermistorPairedControl(current,'wired');end
if refineLeads,current=attachThermistorOuterDomain(current);end
report.outputDirectory=runFolder;report.regions=summary.regions;

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
        if powerMode
            w.heatFluxMethod(:)="整件集总热平衡求壁温，以温差比缩放同网格强等温壁反力，按原面积分";
            w.nearWallTemperatureMethod(:)="同网格首层流体超温按整件热平衡温差比换算，再加环境温度";
        end
        result=struct('accepted',true,'reason','六项残差、温界、独立壁反力、热稳定性与远场收支通过', ...
            'wall',w,'diagnostics',struct('sourceFrameSHA256',s.sourceFrameSHA256));
    end
end
