function report=testPhaseNoseHeatingV227(folder)
%TESTPHASENOSEHEATINGV227 原面保存外缘直接接主工况热剖面/工作台。
d=load('outputs/v104_validation/baseline02/current_case.mat','current','cfg');
src='outputs/v226_validation/matlab02';e=jsondecode(fileread(fullfile(src,'rocket_metadata.json')));
a=readmatrix(fullfile(src,'rocket_faces.csv'));e.pressurePa=a(:,3);e.temperatureK=a(:,4);
e.speedMps=a(:,5);e.densityKgm3=a(:,6);e.velocityMps=a(:,7:9);e.Mach=a(:,10);
e.recoveryTemperatureK=a(:,11);e.polarAngleRad=a(:,12);e.arcLengthM=a(:,13);
e.validFace=logical(a(:,14));e.reasonCode=uint8(a(:,15));e.faceAreaM2=a(:,2);
% CSV不能保存NaN位形；明确恢复原MATLAB nan()位形，再严格核对原SHA。
for key={'pressurePa','temperatureK','speedMps','velocityMps','Mach','recoveryTemperatureK','arcLengthM'}
    x=e.(key{1});x(isnan(x))=nan;e.(key{1})=x;
end
m=d.current;m.phaseAwareNoseEdge=e;cfg=d.cfg;cfg.analysisMode='phaseAwareNoseHeating';cfg.reusePhaseAwareNoseEdge=true;
[m,fs,caseReport]=runAeroHeatSingleCase(m,cfg);h=m.phaseAwareNoseHeating;
write('heating.json',h);write('case_report.json',caseReport);
assert(isequal(m.V,d.current.V)&&isequal(m.F,d.current.F)&&~isfield(m,'heating')&&~isfield(m,'pressure'));
assert(size(m.F,1)==19752&&abs(h.summary.totalAreaM2-1.9646292565715553)<2e-13);
assert(nnz(h.validFace)==1306&&all(~h.validFace|e.validFace),'原1306面未全部取得热候选，查看已保存剖面拒绝。');
payload=buildAerothermalViewerData(m);write('rocket_workbench.json',payload);
assert(numel(payload.fields)==15&&strcmp(payload.defaultField,'heatFlux'));
assert(isequaln(payload.fields(7).values,h.heatFlux)&&all(payload.fields(8).status==2)&&all(isnan(payload.fields(9).values)));
% 实际热/零热驻点控制；只新增本轮剖面，不重做旧外缘或物性表。
hot=solvePhaseAwareLaminarProfile(h.edges(1),h.stations(1),1500);
neutral=solvePhaseAwareLaminarProfile(h.edges(1),h.stations(1),e.stagnation.temperatureK);
write('hot_stagnation.json',hot);write('neutral_stagnation.json',neutral);
assert(hot.available&&neutral.available&&hot.heatFlux<0&&h.profiles{1}.heatFlux>0);
assert(abs(neutral.heatFlux)<1e-4&&~neutral.heatFluxNumericallyResolved);
% 来源、弧线、工况、壁温篡改在昂贵BVP之前拒绝。
rejections=0;id=find(e.validFace,1);
for k=1:5
    bad=m;cc=cfg;
    if k==1,bad.phaseAwareNoseEdge.pressurePa(id)=e.pressurePa(id)+1;
    elseif k==2,bad.phaseAwareNoseEdge.stagnation.enthalpyJkg=e.stagnation.enthalpyJkg+1;
    elseif k==3,bad.V(1,1)=bad.V(1,1)+1e-4;
    elseif k==4,cc.flight.alpha=cfg.flight.alpha+1;
    else,bad.localGeometry.featureFace(id)=true;end
    try,runAeroHeatSingleCase(bad,cc);catch err
        assert(startsWith(err.identifier,'RAH:PhaseNoseHeat:'));rejections=rejections+1;
    end
end
bad=m;bad.phaseAwareNoseHeating.heatFlux(id)=h.heatFlux(id)+1;
try,buildAerothermalViewerData(bad);catch err
    assert(strcmp(err.identifier,'RAH:PhaseNoseHeatViewer:Payload'));rejections=rejections+1;
end
assert(rejections==6);
report=struct('passed',true,'summary',h.summary,'mainEntryUsed',true,'mainViewerUsed',true, ...
    'edgeSnapshotReused',true,'fields',15,'rejections',rejections,'hotHeatFluxWm2',hot.heatFlux, ...
    'neutralHeatFluxWm2',neutral.heatFlux,'physicalQualification',false,'newQualifiedAreaM2',0, ...
    'meshSHA256',e.meshSHA256,'edgeSHA256',e.source.sha256,'heatSHA256',h.sourceSHA256,'freestream',fs);
write('report.json',report);disp(report);
    function write(name,value)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');c=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value));
    end
end
