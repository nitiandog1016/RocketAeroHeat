function report=testPhaseAwareWaveV225(folder)
%TESTPHASEAWAREWAVEV225 新波系守恒/拒绝/状态继承以及真实箭体入口定位。
% 不重做v224平板或物性全表；原箭体仅用原几何/飞行工况，不读取旧外缘数值。
rows=struct('upstreamPressurePa',{},'upstreamTemperatureK',{},'upstreamSpeedMps',{},'result',{});
p=1197.0279007383137;T=226.50908361133003;
a=phaseAwareAirState('PT',p,T,'linearTable');U=5*a.soundSpeedMps;
cases={'normal',0;'oblique',0;'oblique',1e-5;'oblique',10;'oblique',30;'expansion',0;'expansion',15};
for k=1:size(cases,1)
    w=solvePhaseAwareWave(cases{k,1},p,T,U,cases{k,2});assert(w.available,w.status);
    add(w,p,T,U);
    if ~strcmp(cases{k,1},'normal')&&cases{k,2}==0
        assert(w.pressurePa==p&&w.temperatureK==T&&w.speedMps==U);
    end
end
near=rows(3).result;assert(abs(near.pressurePa/p-1)<2e-6);
refused=struct;
w=solvePhaseAwareWave('oblique',p,T,U,60);assert(~w.available);refused.detached=w.status;
w=solvePhaseAwareWave('expansion',p,T,U,80);assert(~w.available);refused.coldExpansion=w.status;
b=phaseAwareAirState('PT',900000,500,'linearTable');
w=solvePhaseAwareWave('normal',900000,500,3*b.soundSpeedMps);assert(~w.available);refused.jointDomain=w.status;
b=phaseAwareAirState('PT',1000,1500,'linearTable');
w=solvePhaseAwareWave('normal',1000,1500,3*b.soundSpeedMps);assert(~w.available);refused.hotShock=w.status;
w=solvePhaseAwareWave('normal',p,T,.9*a.soundSpeedMps);assert(~w.available);refused.subsonic=w.status;

%% 多转角同一路径：压力/总焓/熵与速度方向继承，不从自由流重置
angles=[10;5;15;15];tangent=[-cosd(angles) zeros(4,1) -sind(angles)];
path=marchPhaseAwareWavePath(p,T,U,[-1 0 0],tangent,[0 -1 0]);assert(all(path.valid));
assert(max(abs(path.turnDeg-[10;-5;10;0]))<1e-12);
in=[p T U];for k=1:4,add(path.waves{k},in(1),in(2),in(3));in=[path.pressurePa(k),path.temperatureK(k),path.speedMps(k)];end
fresh=solvePhaseAwareWave('oblique',p,T,U,5);
assert(abs(path.pressurePa(2)/fresh.pressurePa-1)>1e-4);
bad=[-cosd([10;80;5]) zeros(3,1) -sind([10;80;5])];
blocked=marchPhaseAwareWavePath(p,T,U,[-1 0 0],bad,[0 -1 0]);
assert(blocked.valid(1)&&~any(blocked.valid(2:end))&&all(isnan(blocked.pressurePa(2:end))));
rejected=false;try,marchPhaseAwareWavePath(p,T,U,[-1 0 0],[-sqrt(.5) sqrt(.5) 0],[0 -1 0]);
catch err,rejected=strcmp(err.identifier,'RAH:PhasePath:Planar');end
assert(rejected);
write('waves.json',rows);write('path.json',path);write('blocked_path.json',blocked);

%% 真实原面：先定位鼻部入口障碍，禁止用旧有效面带绕开脱体段
d=load('outputs/v104_validation/baseline02/current_case.mat','current','fs');m=d.current;fs=d.fs;
identity=phaseAwareMeshIdentity(m);area=sqrt(sum(cross(m.V(m.F(:,2),:)-m.V(m.F(:,1),:), ...
    m.V(m.F(:,3),:)-m.V(m.F(:,1),:),2).^2,2))/2;
assert(size(m.F,1)==19752&&abs(sum(area)-1.9646292565715553)<2e-13);
sc=m.surfaceCoord;xx=sc.sectionX(:);ns=numel(xx);radius=zeros(ns,1);deviation=radius;
for k=1:ns
    ids=sc.sectionNodeIDs{k};yz=m.V(ids,2:3);center=sc.sectionCenterYZ(k,:);
    rr=vecnorm(yz-center,2,2);radius(k)=mean(rr);deviation(k)=max(abs(rr-radius(k)));
end
% 这是实际环的轴对称剖面诊断，不声称已求解实际20度横流或映射面元状态。
assert(radius(end)<1e-8,'RAH:PhaseWave:Geometry','该诊断要求单鼻尖截面。');
xx=flipud(xx);radius=flipud(radius);dx=diff(xx);dr=diff(radius);lengths=hypot(dx,dr);
wind=[dx zeros(ns-1,1) -dr]./lengths;lee=[dx zeros(ns-1,1) dr]./lengths;
flow=fs.unitFlowDirection(:)';gas=phaseAwareAirState('PT',fs.pressure,fs.temperature,'linearTable');
speed=fs.Mach*gas.soundSpeedMps;
windPath=marchPhaseAwareWavePath(fs.pressure,fs.temperature,speed,flow,wind,[0 -1 0]);
leePath=marchPhaseAwareWavePath(fs.pressure,fs.temperature,speed,flow,lee,[0 1 0]);
assert(~windPath.valid(1)&&~leePath.valid(1));
assert(~any(windPath.valid)&&~any(leePath.valid));
rocket=struct('meshSHA256',identity,'faces',19752,'areaM2',sum(area),'Mach',fs.Mach, ...
    'alphaDeg',fs.alpha,'betaDeg',fs.beta,'pressurePa',fs.pressure,'temperatureK',fs.temperature, ...
    'newGasSpeedMps',speed,'legacySpeedNotReused',true,'firstWindwardTurnDeg',windPath.turnDeg(1), ...
    'firstLeewardTurnDeg',leePath.turnDeg(1),'windwardReason',windPath.reason{1},'leewardReason',leePath.reason{1}, ...
    'maximumSectionRadialDeviationM',max(deviation),'geometryOnlyMeridianDiagnostic',true, ...
    'originalFaceStatesProduced',false,'newCandidateHeatAreaM2',0,'newQualifiedAreaM2',0, ...
    'boundaryLayerOriginAvailable',false,'existingFallbackUsed',false, ...
    'nextPhysicalRequirement','真实钝鼻脱体激波/亚声速至声速线的同模型外缘与入口历史；二维弱波不能越过缺失入口');
write('rocket_inlet.json',rocket);
fid=fopen(fullfile(folder,'rocket_profile.csv'),'w');clean=onCleanup(@()fclose(fid));
fprintf(fid,'xM,radiusM\n');fprintf(fid,'%.17g,%.17g\n',[xx radius]');clear clean;
report=struct('passed',true,'waveCases',numel(rows),'refusals',refused,'downstreamNoReset',true, ...
    'nonPlanarRejected',true,'rocket',rocket,'CFDRuns',0,'newQualifiedAreaM2',0,'physicalQualification',false);
write('report.json',report);disp(report);
    function add(w,pa,te,sp)
        rows(end+1)=struct('upstreamPressurePa',pa,'upstreamTemperatureK',te,'upstreamSpeedMps',sp,'result',w);
    end
    function write(name,value)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');clean=onCleanup(@()fclose(fid));
        fprintf(fid,'%s\n',jsonencode(value,'ConvertInfAndNaN',true));
    end
end
