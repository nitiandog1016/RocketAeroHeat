function report=testPhaseNoseRecoveryV228(folder)
%TESTPHASENOSERECOVERYV228 唯一非零角站壁列实际回归及原面主工作台导出。
% 模块一：复用已有独立恢复数据，不重新计算31站或外缘。
h=jsondecode(fileread('docs/validation/v227/evidence/recovered/heating.json'));
reference=jsondecode(fileread('docs/validation/v227/evidence/recovered/report.json'));
e=jsondecode(fileread('outputs/v226_validation/matlab02/rocket_metadata.json'));
a=readmatrix('outputs/v226_validation/matlab02/rocket_faces.csv');
e.pressurePa=a(:,3);e.temperatureK=a(:,4);e.speedMps=a(:,5);e.densityKgm3=a(:,6);
e.velocityMps=a(:,7:9);e.Mach=a(:,10);e.recoveryTemperatureK=a(:,11);e.arcLengthM=a(:,13);
e.validFace=logical(a(:,14));e.reasonCode=uint8(a(:,15));e.faceAreaM2=a(:,2);
for key={'pressurePa','temperatureK','speedMps','velocityMps','Mach','recoveryTemperatureK','arcLengthM'}
    z=e.(key{1});z(isnan(z))=nan;e.(key{1})=z;
end
fields={'heatFlux','wallShearStress','delta','deltaStar','thetaPhysical','massThickness','skinFriction', ...
    'upperWeight','wallTemperature','lowerStation','upperStation'};
for key=fields,z=h.(key{1});z(isnan(z))=nan;h.(key{1})=z;end
% 模块二：直接执行修复后的30°剖面；参考来自独立物理温度坐标BVP。
j=find(abs(h.anglesRad-pi/6)<1e-12);assert(isscalar(j));
r=solvePhaseAwareLaminarProfile(h.edges(j),h.stations(j),300);write('station.json',r);
assert(r.available&&r.wallShearStress>100);
exact=reference.independentProfiles(2).values;errors=struct;
for key={'heatFlux','wallShearStress','delta','deltaStar','thetaPhysical','massThickness'}
    name=key{1};errors.(name)=abs(r.(name)/exact.(name)-1);
    limit=2e-4;if ismember(name,{'delta','deltaStar','massThickness'}),limit=2e-3;end
    assert(errors.(name)<limit);
end
s=phaseAwareAirState('PT',r.edge.pressurePa,r.edge.temperatureK,'linearTable');
tau=s.dynamicViscosityPas*r.edge.speedMps/r.thicknessScale*r.profile.states(3,1);
assert(abs(tau/r.wallShearStress-1)<1e-12);
% 模块三：实际主导出器恢复15字段/区域面积；保留原面完整分母及未知值。
d=load('outputs/v104_validation/baseline02/current_case.mat','current');m=d.current;
m.phaseAwareNoseEdge=e;m.phaseAwareNoseHeating=h;m.phaseAwareWallTemperature=h.wallTemperature;
m.aerothermalResultMode='phaseAwareNoseHeating';actual=buildAerothermalViewerData(m);
prior=jsondecode(fileread('docs/validation/v227/evidence/recovered/rocket_workbench.json'));
assert(numel(actual.fields)==15&&size(actual.faces,1)==19752);
for k=1:15
    assert(strcmp(actual.fields(k).key,prior.fields(k).key));
    assert(isequaln(actual.fields(k).values,prior.fields(k).values));
    assert(isequal(actual.fields(k).status,prior.fields(k).status));
end
assert(isequal(actual.vertices,prior.vertices)&&isequal(actual.faces,prior.faces));
assert(isequal(actual.region,prior.region)&&numel(actual.regions)==2);
assert(sum([actual.regions.faces])==19752&&abs(sum([actual.regions.area])-actual.totalArea)<1e-12);
for k=1:2
    mask=strcmp(actual.region,actual.regions(k).name);
    assert(nnz(mask)==actual.regions(k).faces&&abs(sum(actual.faceArea(mask))-actual.regions(k).area)<1e-12);
end
assert(nnz(actual.fields(7).status==1)==1306&&all(actual.fields(8).status==2)&&all(actual.fields(9).status==0));
write('workbench.json',actual);
report=struct('passed',true,'actualCorrectedStationExecuted',true,'allStationsRerun',false, ...
    'mainExporterExecuted',true,'fields',15,'faceCount',19752,'fieldValuesEqual',true,'regions',actual.regions, ...
    'relativeErrors',errors,'wallShearStressPa',r.wallShearStress,'heatFluxWm2',r.heatFlux, ...
    'heatCandidateAreaM2',h.summary.heatCandidateAreaM2,'signedPowerW',h.summary.signedPowerW, ...
    'newHeatCandidateAreaM2',0,'newQualifiedAreaM2',0,'physicalQualification',false,'legacyFailuresRetained',true);
write('report.json',report);disp(report);
    function write(name,value)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');c=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value));
    end
end
