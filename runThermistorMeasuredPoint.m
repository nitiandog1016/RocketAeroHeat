function [current,report]=runThermistorMeasuredPoint(input,cfg)
% 模块一：核对工作台实际准备的 96 档体网、原面和图 7 工况。
root=fileparts(mfilename('fullpath'));p=input.thermistorPrepared;cases=fullfile(root,'cases','thermistor_b');
b=thermistorWorkbenchBoundary(cfg.thermistorBoundary,cfg.thermistorAngularCells,cfg);
key=jsondecode(fileread(fullfile(cases,'windPoint12.json')));
assert(p.angularCells==96 && strcmp(p.meshSHA256,key.meshSHA256) && strcmp(p.nodeMapSHA256,key.nodeMapSHA256) && ...
    strcmp(cm011FileSHA256(p.manifest),p.manifestSHA256) && isequal(jsondecode(fileread(p.manifest)),rmfield(p,{'manifest','manifestSHA256'})));
assert(strcmp(cm011FileSHA256(fullfile(p.directory,'mesh.su2')),key.meshSHA256) && ...
    strcmp(cm011FileSHA256(fullfile(p.directory,'nodes.npz')),key.nodeMapSHA256) && strcmp(p.bodySHA256,key.bodySHA256));
[native,guard]=openThermistorArchive('windPoint12'); %#ok<NASGU>
assert(strcmp(cm011FileSHA256(fullfile(native,'comparison.json')),key.comparisonSHA256));
record=jsondecode(fileread(fullfile(native,'comparison.json')));exp=record.experiment;
exp.measuredAxisAngleDeg=NaN;exp.reportedPointUncertaintyWm2K=NaN;
assert(exp.nominalDeltaTK==12 && exp.speedMps==b.speed && all(exp.deltaTReadRangesK(:,1)<=12 & exp.deltaTReadRangesK(:,2)>=12));
nativeCase=jsondecode(fileread(fullfile(native,'case.json')));
assert(nativeCase.scalarReferenceDeltaK==24 && nativeCase.airK==b.airTemperatureK && nativeCase.wallK==b.wallTemperatureK && ...
    nativeCase.airPa==b.airPressurePa && isequal(nativeCase.velocityMps(:)',[0 -3 0]));
action='solve';if isfield(cfg,'thermistorResultSource'),action=validatestring(cfg.thermistorResultSource,{'solve','archive'});end
folder=native;fieldKey=key;prefix='';
if strcmp(action,'solve')
    prefix='six_';fieldKey=jsondecode(fileread(fullfile(cases,'sidewindAngular96.json')));
    fieldKey.wallSHA256=key.sixKWallSHA256;fieldKey.sourceFrameSHA256=key.sixKFrame;fieldKey.sourceReactionSHA256=key.sixKReaction;
end
for pair={'body','case','wall'}
    name=pair{1};ext='.json';if strcmp(name,'wall'),ext='.csv';end
    file=[prefix name ext];if strcmp(name,'body'),file='body.json';end
    assert(strcmp(cm011FileSHA256(fullfile(folder,file)),fieldKey.([name 'SHA256'])));
end
spec=jsondecode(fileread(fullfile(folder,[prefix 'case.json'])));body=jsondecode(fileread(fullfile(folder,'body.json')));
s=jsondecode(fileread(fullfile(folder,[prefix 'summary.json'])));rows=readtable(fullfile(folder,[prefix 'wall.csv']));
assert(isequal(body.F,input.F) && isequal(body.V,input.V) && strcmp(s.meshSHA256,p.meshSHA256) && ...
    strcmp(s.sourceFrameSHA256,fieldKey.sourceFrameSHA256) && strcmp(s.sourceReactionSHA256,fieldKey.sourceReactionSHA256) && ...
    strcmp(s.wallSHA256,fieldKey.wallSHA256) && strcmp(s.nodeMapSHA256,p.nodeMapSHA256) && ...
    strcmp(s.energy.inputFrameSHA256,s.sourceFrameSHA256) && abs(s.energy.relativeDifference)<.005 && s.relativeChange20<.002 && ...
    all(structfun(@(v)v<-8,s.residualWindow)) && all(structfun(@(v)v<-8,s.plus20Window)));
assert(spec.angularCells==96 && spec.outerRadiusOverD==20 && spec.sourceAreaM2==exp.areaM2 && ...
    isequal(spec.velocityMps(:)',[0 -3 0]) && spec.airK==b.airTemperatureK && spec.airPa==b.airPressurePa);
% 模块二：恒物性条件下从 6 K 基底求当前场，或显示已独立施加 12 K 热边界的原生结果。
[rows,factor,delta,wallK]=scaleThermistorWallField(rows,b,spec.airK,spec.wallK,s.totalOutwardW);
if strcmp(action,'solve')
    assert(factor==2 && strcmp(record.sixKFrame,s.sourceFrameSHA256));
    heatMethod='96 档 6 K 强等温壁反力按交叠面积积分至原面，再按 12/6 温差比换算';
    firstMethod='保留的三个首层节点的流体超温乘 12/6，再加回环境温度';
else
    assert(factor==1 && strcmp(record.twelveKFrame,s.sourceFrameSHA256));
    heatMethod='96 档原生强标量壁值 0.5，按 24 K 单位换算壁反力并在细面与原面交叠区积分';
    firstMethod='同一三个首层节点的原生标量按 24 K 单位换算为温度后平均';
end
pressureMethod='96 档不可压层流壁节点表压按细面线性形函数积分到原面后加环境压力';
f=input.F;v=input.V;n=size(f,1);
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(isequal(rows.faceID,(0:n-1)') && isequal(rows.regionID,body.regionID) && all(abs(area-rows.areaM2)<1e-18) && ...
    all(rows.pressurePa>0) && all(abs(rows.pressureGaugePa+spec.airPa-rows.pressurePa)<1e-9) && ...
    all(rows.nearWallCellTemperatureK>=spec.airK & rows.nearWallCellTemperatureK<=wallK));
condition=struct('caseLabel','B 型热敏珠 · 图 7 测点 · 名义侧吹','velocityMps',[0 -3 0],'ambientPressurePa',spec.airPa, ...
    'ambientTemperatureK',spec.airK,'wallTemperatureK',wallK,'alphaDeg',0,'betaDeg',90,'geometryScope','完整 360° 名义带线试件');
names=string(body.regionNames);
setup=struct('condition',condition,'region',names(body.regionID),'wallTemperatureK',wallK*ones(n,1), ...
    'engineeringSolver',@emptyWall,'localSolver',@wallData);
[current,~,report]=runAeroHeatSingleCase(struct('F',f,'V',v),struct('analysisMode','hybridSurface','hybridSurface',setup));
watts=-area.*rows.heatFluxWm2;
assert(report.fullSurfaceComputed && abs(report.wholeSurfaceHeatPowerW+sum(watts))<1e-12 && abs(sum(watts)-factor*s.totalOutwardW)<1e-12);
% 模块三：按实验定义比较珠部平均量，保持支线、端面、测点范围与温度身份可追溯。
m=struct('action',action,'wallK',wallK,'deltaTK',delta,'areaM2',sum(area),'meshSHA256',p.meshSHA256, ...
    'experiment',exp,'wholeW',sum(watts),'beadW',sum(watts(ismember(rows.regionID,3:5))), ...
    'jointVerifiedM2',0,'sourceCheck',record.nativeCheck,'filmTemperatureK',(spec.airK+wallK)/2);
m.beadHWm2K=m.beadW/(delta*exp.areaM2);m.relativeDifference=m.beadHWm2K/exp.hWm2K-1;
m.relativeAtGraphLimits=m.beadHWm2K./fliplr(exp.graphReadRangeWm2K(:)')-1;
m.equivalentExperimentPowerW=exp.hWm2K*exp.areaM2*delta;
m.scaledFrom=struct('frame',key.sixKFrame,'wallSHA256',key.sixKWallSHA256,'wallK',299.15,'temperatureRatio',2);
m.originalSixKFrame=record.checkSeedFrame;m.unmatched=record.unmatched;m.sixKPressureDriftPa=record.sixKPressureDriftPa;
m.nativeCheck=struct('frame',key.sourceFrameSHA256,'wallSHA256',key.wallSHA256,'wallK',305.15,'scalarUnitK',24,'wallScalar',.5);
for j=1:numel(names)
    use=rows.regionID==j;
    m.regions(j)=struct('name',char(names(j)),'faces',nnz(use),'areaM2',sum(area(use)),'outwardW',sum(watts(use)));
    assert(abs(m.regions(j).outwardW-factor*s.regions(j).outwardW)<1e-12);
end
[~,m.peakFaceID]=min(rows.heatFluxWm2);
[ax,axGuard]=openThermistorArchive('leadsThird3'); %#ok<NASGU>
axCase=jsondecode(fileread(fullfile(ax,'case.json')));axS=jsondecode(fileread(fullfile(ax,'summary.json')));
axRows=readtable(fullfile(ax,'wall.csv'));
assert(strcmp(axCase.caseKind,'thermistor-B-axial-two-leads') && axCase.wallK==wallK && axCase.airK==spec.airK && ...
    axCase.airPa==spec.airPa && axCase.speedMps==spec.speedMps && strcmp(axS.bodySHA256,p.bodySHA256) && ...
    isequal(axRows.faceID,rows.faceID) && isequal(axRows.regionID,rows.regionID) && strcmp(cm011FileSHA256(fullfile(ax,'wall.csv')),axS.wallSHA256));
bead=ismember(rows.regionID,3:5);axPower=-area(bead)'*axRows.heatFluxWm2(bead);axH=axPower/(exp.areaM2*delta);
assert(abs(axPower-axS.beadAndJointsW)<1e-12);
m.axial=struct('method','轴向第三档 30° 对称稳态','hWm2K',axH,'beadW',axPower,'relativeDifference',axH/exp.hWm2K-1, ...
    'frame',axS.sourceFrameSHA256,'meshSHA256',axS.meshSHA256,'wallSHA256',axS.wallSHA256);
m.correlation=thermistorBeadMeanConvection(spec,wallK);assert(m.correlation.available);
m.correlation.relativeDifference=m.correlation.meanHWm2K/exp.hWm2K-1;
current.thermistorMeasuredPoint=m;
current.thermistorSidewind=struct('beadAndJointsW',m.beadW,'totalOutwardW',m.wholeW,'regions',m.regions,'peakFaceID',m.peakFaceID);
current.hybridWall.pressureGaugePa=rows.pressureGaugePa;current.hybridWall.pressureGaugeMethod=repmat("细面壁节点表压积分回原面",n,1);
current.evidenceIdentity.sourceFrame=struct('kind','SU2-8.6.0-thermistor-B-96-wind-point','sha256',s.sourceFrameSHA256, ...
    'bodySHA256',p.bodySHA256,'reactionSHA256',s.sourceReactionSHA256,'wallSHA256',s.wallSHA256,'referenceWallK',spec.wallK,'ratio',factor);
current.evidenceIdentity.meshSHA256=p.meshSHA256;current.evidenceIdentity.preparedSHA256=p.manifestSHA256;
current.evidenceIdentity.runDirectory='';current.evidenceIdentity.reusedArchive=true;current.evidenceIdentity.sourceTraceable=true;
current.evidenceIdentity.productionQualified=false;current.evidenceIdentity.scope=['B 珠 3 m/s，12 K；恒物性等温名义试件。' exp.mounting];
report.outputDirectory='';report.regions=m.regions;report.measuredPoint=m;
    function t=emptyWall(~,~,~)
        t=rows([],{'faceID','pressurePa','heatFluxWm2'});t.method=strings(0,1);
    end
    function r=wallData(~,~,~,~,ids,~)
        t=rows(ids+1,{'faceID','pressurePa','heatFluxWm2','nearWallCellTemperatureK'});count=height(t);
        t.method=repmat("带线 B 珠 96 档侧吹等温被动标量",count,1);t.pressureMethod=repmat(string(pressureMethod),count,1);
        t.heatFluxMethod=repmat(string(heatMethod),count,1);t.nearWallTemperatureMethod=repmat(string(firstMethod),count,1);
        r=struct('accepted',true,'reason','同网原面、原生强边界与能量账通过','wall',t,'diagnostics',struct('sourceFrameSHA256',s.sourceFrameSHA256));
    end
end
