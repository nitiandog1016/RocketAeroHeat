function [current,report]=runThermistorWireWorkbench(input,cfg,progress)
% 模块一：核对冻结的试件、真实体网及长引线边界，选择原生新算、接续或已算全场。
if nargin<3,progress=[];end
root=fileparts(mfilename('fullpath'));p=input.thermistorPrepared;
b=thermistorWorkbenchBoundary(cfg.thermistorBoundary,cfg.thermistorAngularCells,cfg);
assert(strcmp(b.leadThermalModel,'longWire') && p.angularCells==96 && ...
    strcmp(cm011FileSHA256(p.manifest),p.manifestSHA256) && isequal(jsondecode(fileread(p.manifest)),rmfield(p,{'manifest','manifestSHA256'})));
action='solve';if isfield(cfg,'thermistorResultSource'),action=validatestring(cfg.thermistorResultSource,{'solve','archive','resumeWire'});end
job='';reading=isfield(b,'resistanceTemperatureK');
if strcmp(action,'archive') || reading
    tag='longWire12';
    if reading
        assert(any(strcmp(action,{'solve','archive'})));
        if strcmp(action,'archive'),tag='resistance12';end
    end
    [folder,guard]=openThermistorArchive(tag); %#ok<NASGU>
    key=jsondecode(fileread(fullfile(root,'cases','thermistor_b',[tag '.json'])));
    assert(strcmp(cm011FileSHA256(fullfile(folder,'summary.json')),key.summarySHA256));
else
    if isfield(cfg,'thermistorRunDirectory'),job=cfg.thermistorRunDirectory;
    elseif strcmp(action,'solve'),job=fullfile(root,'outputs',['thermistor_wire_' char(java.util.UUID.randomUUID())]);
    else,error('RAH:Wire:Job','请选择引线计算任务目录。');end
    task=fullfile(job,'run');
    if strcmp(action,'resumeWire') && isfile(fullfile(job,'prepared','profile.json')),task=job;end
    if ~strcmp(action,'resumeWire') || ~isfile(fullfile(task,'completed.json'))
        token=char(java.util.UUID.randomUUID());if isfield(cfg,'thermistorInvocation'),token=cfg.thermistorInvocation;end
        args={fullfile(root,'scripts','thermistorWireWorkbench.py'),action,job,p.directory,'--invocation',token};
        if isfield(cfg,'wireSeedDirectory') && ~isempty(cfg.wireSeedDirectory),args=[args {'--base',cfg.wireSeedDirectory}];end
        started=tic;runCM011Python(args,[job '.' char(java.util.UUID.randomUUID()) '.log'],@pulse,true);
        ack=jsondecode(fileread(fullfile(job,['receipt-' token '.json'])));assert(strcmp(ack.invocation,token));
        if strcmp(ack.status,'paused'),error('RAH:Bead:Paused','引线任务已在保存点暂停。%s',job);end
        assert(strcmp(ack.status,'complete'));
    end
    folder=fullfile(task,'surface');
    done=jsondecode(fileread(fullfile(task,'completed.json')));
    assert(~isempty(regexp(done.stage,'^steady0[1-5]$','once')) && strcmp(done.summarySHA256,cm011FileSHA256(fullfile(folder,'summary.json'))) && ...
        strcmp(done.frame,cm011FileSHA256(fullfile(task,done.stage,'restart.csv'))));
end
s=jsondecode(fileread(fullfile(folder,'summary.json')));rows=readtable(fullfile(folder,'wall.csv'));
if ~strcmp(action,'archive') && ~reading
    assert(strcmp(done.frame,s.sourceFrameSHA256) && strcmp(done.wallSHA256,s.wallSHA256) && ...
        strcmp(cm011FileSHA256(fullfile(task,'prepared','profile.dat')),s.boundarySHA256) && ...
        strcmp(cm011FileSHA256(fullfile(task,done.stage,'reaction0','restart.csv')),s.reactionSHA256));
end
body=jsondecode(fileread(fullfile(p.directory,'body.json')));c=jsondecode(fileread(fullfile(p.directory,'case.json')));
assert(strcmp(s.meshSHA256,p.meshSHA256) && strcmp(s.nodeMapSHA256,p.nodeMapSHA256) && strcmp(s.bodySHA256,p.bodySHA256));
names={'mesh.su2','nodes.npz','body.json'};digests={s.meshSHA256,s.nodeMapSHA256,s.bodySHA256};
for k=1:3,assert(strcmp(cm011FileSHA256(fullfile(p.directory,names{k})),digests{k}));end
assert(strcmp(cm011FileSHA256(fullfile(folder,'wall.csv')),s.wallSHA256) && strcmp(s.energy.sourceFrameSHA256,s.sourceFrameSHA256));
assert(all(structfun(@(x)x<-8,s.window)) && all(structfun(@(x)x<-8,s.plus20Window)) && abs(s.energy.relativeDifference)<.005 && s.relativeChange20<.002);
if reading
    [rows,s,measurement,change]=applyThermistorThermometry(rows,s,c,b.resistanceTemperatureK,action);beadK=measurement.surfaceK;
else,beadK=b.wallTemperatureK;end
% 模块二：将全部原面的输入壁温、壁反力、表压和共同首层温度按实际来源接入统一求解出口。
f=input.F;v=input.V;n=size(f,1);a=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(isequal(f,body.F) && isequal(v,body.V) && n==s.faces && isequal(rows.faceID,(0:n-1)') && isequal(rows.regionID,body.regionID) && ...
    all(abs(rows.areaM2-a)<1e-18) && all(rows.pressurePa>0) && all(abs(rows.pressurePa-(b.airPressurePa+rows.pressureGaugePa))<1e-9));
assert(all(isfinite(rows.heatFluxWm2)) && all(rows.nearWallCellTemperatureK>=b.airTemperatureK & rows.nearWallCellTemperatureK<=beadK));
assert(abs(min(rows.wallInputK)-s.profile.endK)<1e-9 && abs(max(rows.wallInputK)-beadK)<1e-9);
condition=struct('caseLabel','B 珠 · 长引线温度边界','velocityMps',[0 -3 0],'ambientPressurePa',b.airPressurePa, ...
    'ambientTemperatureK',b.airTemperatureK,'beadWallK',beadK,'leadThermalModel','longWire', ...
    'boundarySHA256',s.boundarySHA256,'alphaDeg',0,'betaDeg',90,'geometryScope','完整 360° 名义带线试件');
if reading,condition.resistanceTemperatureK=measurement.internalK;condition.temperatureModel=measurement.method;end
regions=string(body.regionNames);setup=struct('condition',condition,'region',regions(body.regionID), ...
    'wallTemperatureK',rows.wallInputK,'engineeringSolver',@emptyWall,'localSolver',@wallData);
[current,~,report]=runAeroHeatSingleCase(struct('F',f,'V',v),struct('analysisMode','hybridSurface','hybridSurface',setup));
assert(report.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW+s.totalOutwardW)<1e-12 && abs(a'*rows.heatFluxWm2+s.energy.wallReleaseW)<1e-10);
for k=1:7
    use=body.regionID==k;
    assert(s.regions(k).faces==nnz(use) && abs(a(use)'*rows.heatFluxWm2(use)+s.regions(k).outwardW)<1e-12);
end
% 模块三：按珠部实验定义比较两种边界，保存定向导热、分区空气侧对流及方法身份。
[ref,refGuard]=openThermistorArchive('windPoint12'); %#ok<NASGU>
paper=jsondecode(fileread(fullfile(ref,'comparison.json')));s.experiment=paper.experiment;
s.action=action;s.controlHWm2K=s.controlBeadW/(c.sourceAreaM2*(beadK-c.airK));
s.bands=thermistorWireBands(v,f,rows,s.profile,b.airTemperatureK,beadK);
assert(abs(sum([s.bands.backSideW])-s.regions(2).outwardW)<1e-12 && abs(sum([s.bands.frontSideW])-s.regions(6).outwardW)<1e-12);
s.relativeDifference=s.beadHWm2K/s.experiment.hWm2K-1;s.experimentEquivalentW=s.experiment.hWm2K*c.sourceAreaM2*12;
assert(abs(s.beadHWm2K-s.beadAndJointsW/(c.sourceAreaM2*(beadK-c.airK)))<1e-9);
if reading
    measurement.originalBias=measurement.referenceBeadW/s.experimentEquivalentW-1;measurement.remainingBias=measurement.apparentHWm2K/s.experiment.hWm2K-1;
    measurement.powerReductionW=measurement.referenceBeadW-s.beadAndJointsW;measurement.remainingGapW=s.beadAndJointsW-s.experimentEquivalentW;
    measurement.fractionOfPowerGap=measurement.powerReductionW/(measurement.referenceBeadW-s.experimentEquivalentW);
    measurement.paperUpperScenarioRatio=.03;
    measurement.upperScenarioGapW=measurement.referenceBeadW/(1+measurement.paperUpperScenarioRatio)-s.experimentEquivalentW;
    measurement.upperScenarioBias=measurement.upperScenarioGapW/s.experimentEquivalentW;
    s.relativeDifference=measurement.remainingBias;s.thermometry=measurement;current.thermometryFluxChange=change;
end
current.thermistorWire=s;current.thermistorSidewind=struct('beadAndJointsW',s.beadAndJointsW,'totalOutwardW',s.totalOutwardW, ...
    'regions',s.regions,'peakFaceID',s.peakFaceID);
current.hybridWall.pressureGaugePa=rows.pressureGaugePa;
current.hybridWall.pressureGaugeMethod=repmat("细面壁节点表压交叠积分至原面",n,1);
current.thermistorWireFluxChange=rows.heatFluxWm2-rows.uniformHeatFluxWm2;
current.evidenceIdentity.sourceFrame=struct('kind','SU2-8.6.0-thermistor-long-wire','sha256',s.sourceFrameSHA256, ...
    'bodySHA256',s.bodySHA256,'reactionSHA256',s.reactionSHA256,'wallSHA256',s.wallSHA256);
current.evidenceIdentity.meshSHA256=p.meshSHA256;current.evidenceIdentity.preparedSHA256=p.manifestSHA256;
current.evidenceIdentity.runDirectory=job;current.evidenceIdentity.reusedArchive=reading || strcmp(action,'archive');
if reading
    current.evidenceIdentity.sourceFrame.kind='SU2-8.6.0-thermistor-resistance-temperature';
    current.evidenceIdentity.sourceFrame.fieldFactor=measurement.fieldFactor;current.evidenceIdentity.sourceFrame.temperatureMethod=measurement.fieldMethod;
end
current.evidenceIdentity.sourceTraceable=true;current.evidenceIdentity.productionQualified=false;
current.evidenceIdentity.scope='3 m/s 恒物性空气，线根到截面的给定壁温来自长线导热模型；局部实验精度未知。';
report.outputDirectory=job;report.regions=s.regions;report.wire=s;
    function pulse()
        if ~isempty(progress),progress(struct('phase','wire','status','running','message','','fraction',NaN,'elapsedSeconds',toc(started)));end
    end
    function t=emptyWall(~,~,~)
        t=rows([],{'faceID','pressurePa','heatFluxWm2'});t.method=strings(0,1);
    end
    function r=wallData(~,~,~,~,ids,~)
        t=rows(ids+1,{'faceID','pressurePa','heatFluxWm2','nearWallCellTemperatureK'});count=height(t);
        t.method=repmat("SU2 8.6.0 不可压层流及全壁逐节点强标量",count,1);
        t.pressureMethod=repmat("壁节点表压交叠积分，加 101325 Pa 环境压力",count,1);
        t.heatFluxMethod=repmat("非均匀强等温壁行反力，按原面交叠面积守恒积分，正值入壁",count,1);
        t.nearWallTemperatureMethod=repmat("同一三个实际首层流体节点的标量按 24 K 单位换算后平均",count,1);
        if reading && strcmp(action,'solve')
            t.heatFluxMethod=repmat("v317 长引线节点壁反力的原面热流乘全壁统一超温比",count,1);
            t.nearWallTemperatureMethod=repmat("固定首层三节点的参照超温乘全壁统一超温比，加环境温度",count,1);
        end
        r=struct('accepted',true,'reason','原生温界、六项末窗、双次反力和同态能量账通过','wall',t);
    end
end
