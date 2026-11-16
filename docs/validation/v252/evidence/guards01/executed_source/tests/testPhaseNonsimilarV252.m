function report=testPhaseNonsimilarV252(folder)
%TESTPHASENONSIMILARV252 原真实外缘/驻点复用，实际非相似链及壁温历史验证。
% 模块一：读取已绑定数据，不重做原波系/局部相似31站。
d=load('outputs/v104_validation/baseline02/current_case.mat','current');m=d.current;
e=jsondecode(fileread('outputs/v226_validation/matlab02/rocket_metadata.json'));
a=readmatrix('outputs/v226_validation/matlab02/rocket_faces.csv');
e.pressurePa=a(:,3);e.temperatureK=a(:,4);e.speedMps=a(:,5);e.densityKgm3=a(:,6);e.velocityMps=a(:,7:9);
e.Mach=a(:,10);e.recoveryTemperatureK=a(:,11);e.polarAngleRad=a(:,12);e.arcLengthM=a(:,13);
e.validFace=logical(a(:,14));e.reasonCode=uint8(a(:,15));e.faceAreaM2=a(:,2);
for key={'pressurePa','temperatureK','speedMps','velocityMps','Mach','recoveryTemperatureK','arcLengthM'}
    z=e.(key{1});z(isnan(z))=nan;e.(key{1})=z;
end
old=jsondecode(fileread('outputs/v227_validation/audit02/heating.json'));
if iscell(old.profiles),seed=old.profiles{1};else,seed=old.profiles(1);end
angles=old.anglesRad(:)';R=max(m.V(:,1))-e.sphereCenter(1);
[edges,stations]=phaseAwareNoseMarchInputs(e,angles,R);
assert(max(abs([stations.xi]'-[old.stations.xi]'))<1e-15);
domain=struct('geometryKind','rotated-spherical-research','crossflowPresent',false,'separationPresent',false, ...
    'layerState','laminar','sourceRole','engineering-candidate','sourceSHA256',e.source.sha256, ...
    'angleOfAttackDeg',acosd(-e.condition.flowDirection(1)),'radiusM',R);
% 模块二：真实无源沿程解、一次预定细化及同终点壁温的不同上游历史。
base=marchPhaseAwareLaminarProfiles(edges,stations,300*ones(size(angles)),seed,domain);write('base.json',base);accepted(base);
fineAngles=sort(unique([angles .5*(angles(1:end-1)+angles(2:end))]));[fe,fs]=phaseAwareNoseMarchInputs(e,fineAngles,R);
fine=marchPhaseAwareLaminarProfiles(fe,fs,300*ones(size(fineAngles)),seed,domain);write('fine.json',fine);accepted(fine);
Tw=300+zeros(size(angles));deg=rad2deg(angles);mask=deg>15&deg<35;Tw(mask)=Tw(mask)+100*sin(pi*(deg(mask)-15)/20).^2;
pulse=marchPhaseAwareLaminarProfiles(edges,stations,Tw,seed,domain);write('history.json',pulse);accepted(pulse);
% 模块三：同一原面映射和完整原分母，不将球帽研究能力扩大到后体。
hb=mapPhaseAwareMarchToSurface(m,e,base);hf=mapPhaseAwareMarchToSurface(m,e,fine);hp=mapPhaseAwareMarchToSurface(m,e,pulse);
write('surface_base.json',hb);write('surface_fine.json',hf);write('surface_history.json',hp);
assert(nnz(hb.validFace)==1306&&isequal(hb.validFace,hf.validFace)&&isequal(hb.validFace,hp.validFace));
assert(size(m.F,1)==19752&&abs(hb.summary.totalAreaM2-1.9646292565715553)<2e-13);
m.phaseAwareNoseEdge=e;m.phaseAwareWallTemperature=300+zeros(size(m.F,1),1);
viewer=buildPhaseAwareNoseViewerData(m);write('geometry_viewer.json',viewer);
q=cellfun(@(z)z.heatFlux,base.profiles);qf=cellfun(@(z)z.heatFlux,fine.profiles);qh=cellfun(@(z)z.heatFlux,pulse.profiles);
same=ismember(fineAngles,angles);difference=max(abs(qf(same)-q)./max(1,abs(qf(same))));
j=find(deg>=40,1);memory=abs(qh(j)-q(j))/max(1,abs(q(j)));
assert(Tw(j)==300&&memory>1e-5,'RAH:PhaseMarch:Memory','同终点壁温却未分辨上游历史作用。');
% 准入拒绝只检查真实新接口边界，不重做旧物性/平板控制。
rejected=0;bad=domain;bad.crossflowPresent=true;check(@()marchPhaseAwareLaminarProfiles(edges,stations,Tw,seed,bad));
bad=domain;bad.geometryKind='axisymmetric';check(@()marchPhaseAwareLaminarProfiles(edges,stations,Tw,seed,bad));
bad=domain;bad.sourceRole='accepted-edge';check(@()marchPhaseAwareLaminarProfiles(edges,stations,Tw,seed,bad));
badseed=seed;badseed.wallTemperature=400;check(@()marchPhaseAwareLaminarProfiles(edges,stations,Tw,badseed,domain));
badst=stations;badst(2).xi=badst(1).xi;check(@()marchPhaseAwareLaminarProfiles(edges,badst,Tw,seed,domain));
bade=edges;bade(2).temperatureK=bade(2).temperatureK+2;check(@()marchPhaseAwareLaminarProfiles(bade,stations,Tw,seed,domain));
assert(rejected==6);
report=struct('passed',true,'branchCounts',[numel(base.profiles) numel(fine.profiles) numel(pulse.profiles)], ...
    'originalFaces',19752,'candidateFaces',1306,'candidateAreaM2',hb.summary.candidateAreaM2, ...
    'physicalQualification',false,'newQualifiedAreaM2',0,'streamwiseRefinementMaximumQRelativeDifference',difference, ...
    'streamwiseRefinementGate',difference<=.02,'historyEffectAt40Deg',memory,'sameFinalWallTemperature',true, ...
    'rejections',rejected,'MatlabVersion',version,'sourceEdgeSHA256',e.source.sha256,'basePowerW',hb.summary.signedPowerW);
write('report.json',report);disp(report);
    function write(name,obj)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(obj));
    end
    function accepted(c)
        if ~all(c.available),i=find(~c.available,1);disp(c.profiles{i});end
        assert(all(c.available),'RAH:PhaseMarch:Branch','有非相似站拒绝，已保存准确失败。');
    end
    function check(fn)
        try,fn();catch err,assert(startsWith(err.identifier,'RAH:PhaseMarch:'));rejected=rejected+1;end
    end
end
