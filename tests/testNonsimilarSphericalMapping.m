function report=testNonsimilarSphericalMapping(folder)
%TESTNONSIMILARSPHERICALMAPPING 核对冻结条件段重放、真实面元映射及拒绝合同。
% 不重跑球鼻长段；只验证已有数值证据及几何/工况一致性，不是实验验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','nonsimilar'));
if nargin<1,folder=fullfile(root,'outputs','v61_validation','mapping');end
if ~exist(folder,'dir'),mkdir(folder);end
study=loadNonsimilarSphericalStudy;fs=study.freestream;
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=fs.Mach;cfg.flight.altitude=fs.altitude;
cfg.flight.alpha=0;cfg.flight.beta=0;cfg.flight.Twall=study.uniformWallTemperature;
cfg.heating.airThermalModel='thermallyPerfect';
[mesh,meta]=buildCanonicalBody('hemisphereCylinder',struct('resolution',48));mesh=prepare(mesh,cfg);before=mesh;
mapped=mapNonsimilarSphericalStudy(mesh,study,fs,cfg);
assert(isequaln(mesh,before)&&mapped.available&&~mapped.physicalSolution&&~mapped.eligibleForSurfaceCoverage);
assert(~study.historicalPreRunSourceBinding&&study.currentSourceReplayPassed&&numel(study.station)==80);

%% 一、独立解析球面真值：完整顶点均在接受段内，区域外没有补值
origin=[max(mesh.V(:,1))-meta.radius,0,0];r=mesh.V-origin;
angle=acosd(max(-1,min(1,r(:,1)./vecnorm(r,2,2))));
expected=strcmp(meta.faceRegion,'sphericalNose') & ...
    min(angle(mesh.F),[],2)>=study.acceptedIntervalDeg(1)-1e-9 & ...
    max(angle(mesh.F),[],2)<=study.acceptedIntervalDeg(2)+1e-9;
assert(isequal(mapped.faceMask,expected)&&nnz(expected)>20);
assert(all(isnan(mapped.heatFlux(~expected)))&&all(mapped.heatFlux(expected)>0));
assert(all(isnan(mapped.recoveryTemperature))&&all(isnan(mapped.wallTemperatureSolved)));
ids=find(expected);node=mapped.profileNodeIDs(expected,:);w=mapped.profileNodeWeight(expected);
table=[study.station.heatFlux]';truth=table(node(:,1)).*(1-w)+table(node(:,2)).*w;
assert(max(abs(truth-mapped.heatFlux(expected)))<1e-8);
assert(abs(sum(mesh.area(expected).*mapped.heatFlux(expected))-mapped.summary.conditionalHeatPowerW)<1e-9);
assert(mapped.summary.newProductionSolvedAreaM2==0&&mapped.summary.mappedFaces==nnz(expected));
% 顶点全部达标并不能保证边内部达标；用独立闭式反例核验完整三角形界。
lo=40.25;theta=[lo;lo;41];az=[-5;5;0];
v=[cosd(theta),sind(theta).*cosd(az),sind(theta).*sind(az)];
bounds=sphericalTrianglePolarBounds(v,[1 2 3],[0 0 0]);
assert(abs(bounds(1)-atand(tand(lo)*cosd(5)))<1e-12&&bounds(1)<lo);
assert(abs(bounds(2)-41)<1e-12);
az=[0;120;240];v=[cosd(30)*ones(3,1),sind(30)*cosd(az),sind(30)*sind(az)];
bounds=sphericalTrianglePolarBounds(v,[1 2 3],[0 0 0]);
assert(bounds(1)==0&&abs(bounds(2)-30)<1e-12);
assert(mapped.summary.minimumFullFaceAngleDeg>study.acceptedIntervalDeg(1));

%% 二、空间平移不改变同一几何的面元身份、插值及真实热功率
translated=mesh;translated.V=translated.V+[.7 .31 -.22];translated=prepare(translated,cfg);
shifted=mapNonsimilarSphericalStudy(translated,study,fs,cfg);
assert(isequal(mapped.faceMask,shifted.faceMask));
assert(max(abs(mapped.heatFlux(expected)-shifted.heatFlux(expected)))<1e-6);
assert(abs(mapped.summary.conditionalHeatPowerW-shifted.summary.conditionalHeatPowerW)<1e-6);

%% 三、真实条件变化、旧几何、损坏站表和身份升级均须拒绝
refusals=0;changed=fs;changed.alpha=20;changed.unitFlowDirection=[-cosd(20) 0 sind(20)];
reject(@()mapNonsimilarSphericalStudy(mesh,study,changed,cfg));
changed=fs;changed.Mach=fs.Mach+1;reject(@()mapNonsimilarSphericalStudy(mesh,study,changed,cfg));
changed=cfg;changed.flight.Twall=310;reject(@()mapNonsimilarSphericalStudy(mesh,study,fs,changed));
changed=cfg;changed.heating.airThermalModel='caloricallyPerfect';reject(@()mapNonsimilarSphericalStudy(mesh,study,fs,changed));
changed=mesh;changed.V(1,2)=changed.V(1,2)+.001;reject(@()mapNonsimilarSphericalStudy(changed,study,fs,cfg));
changed=study;changed.station(4)=[];reject(@()mapNonsimilarSphericalStudy(mesh,changed,fs,cfg));
changed=study;changed.station(4).heatFlux=complex(changed.station(4).heatFlux,1);reject(@()mapNonsimilarSphericalStudy(mesh,changed,fs,cfg));
changed=study;changed.station(4).thinLayerChecksPassed=NaN;reject(@()mapNonsimilarSphericalStudy(mesh,changed,fs,cfg));
changed=study;changed.freestream.Mach=NaN;reject(@()mapNonsimilarSphericalStudy(mesh,changed,fs,cfg));
changed=study;changed.eligibleForSurfaceCoverage=true;reject(@()mapNonsimilarSphericalStudy(mesh,changed,fs,cfg));
reject(@()loadNonsimilarSphericalStudy(study.source.inletFile,study.source.inletFile));

%% 四、同一WebGL合同显示条件字段，域外断开，不能伪造恢复温度和壁温解
path=exportNonsimilarSphericalStudy(mesh,mapped,fullfile(folder,'web'));
data=jsondecode(fileread(fullfile(folder,'web','results.json')));
heat=data.fields(strcmp({data.fields.key},'conditionalHeatFlux'));
assert(isequal(heat.status,double(expected))&&contains(heat.label,'不计生产覆盖'));
assert(contains(data.caseNote,'条件研究')&&contains(data.caseNote,'新增生产面积0'));
assert(~data.researchIdentity.eligibleForSurfaceCoverage && ...
    ~data.researchIdentity.historicalPreRunSourceBinding && ...
    strcmp(data.researchIdentity.source.studySHA256,study.source.studySHA256));
for j=1:numel(data.regions)
    assert(isequal(sort(fieldnames(data.regions(j).fields)),sort({data.fields.key}')));
    for k=1:numel(data.fields)
        amount=data.regions(j).fields.(data.fields(k).key);
        assert(abs(amount.computedArea+amount.prescribedArea+amount.missingArea-data.regions(j).area)<1e-12);
    end
end
for key={'wallSolved','conditionalRecoveryTemperature'}
    value=data.fields(strcmp({data.fields.key},key{1}));
    assert(all(value.status==0)&&all(isnan(value.values)));
end
bad=mapped;bad.eligibleForSurfaceCoverage=true;reject(@()exportNonsimilarSphericalStudy(mesh,bad,folder));
bad=mapped;bad.heatFlux(ids(1))=NaN;reject(@()exportNonsimilarSphericalStudy(mesh,bad,folder));
bad=mapped;bad.heatFlux(find(~expected,1))=0;reject(@()exportNonsimilarSphericalStudy(mesh,bad,folder));
bad=mapped;bad.faceMask=double(bad.faceMask);reject(@()exportNonsimilarSphericalStudy(mesh,bad,folder));
bad=mapped;bad.wallTemperatureInput(1)=NaN;reject(@()exportNonsimilarSphericalStudy(mesh,bad,folder));
bad=mapped;bad.recoveryTemperature(ids(1))=1000;reject(@()exportNonsimilarSphericalStudy(mesh,bad,folder));
bad=mesh;bad.V(1,1)=bad.V(1,1)+1e-5;reject(@()exportNonsimilarSphericalStudy(bad,mapped,folder));
report=struct('status','passed','acceptedStations',numel(study.station),'mapping',mapped.summary, ...
    'refusals',refusals,'physicalValidation',false,'newProductionSolvedAreaM2',0, ...
    'interpolation',mapped.interpolation,'viewerPath',path,'source',study.source);
save(fullfile(folder,'spherical_mapping.mat'),'study','mapped','mesh','meta','fs','cfg','report');
fid=fopen(fullfile(folder,'spherical_mapping.json'),'w','n','UTF-8');assert(fid>=0);
clean=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('Conditional spherical mapping passed: %d faces, %.9g m2, %.9g W, %d refusals; no production coverage.\n', ...
    nnz(expected),mapped.summary.mappedConditionalAreaM2,mapped.summary.conditionalHeatPowerW,refusals);

    function reject(fun)
        try,fun();catch ex
            if startsWith(ex.identifier,{'RocketAeroHeat:','MATLAB:'}),refusals=refusals+1;return;end
            rethrow(ex);
        end
        error('RocketAeroHeat:Test:NoRefusal','损坏或超域输入未拒绝。');
    end
end

function mesh=prepare(mesh,cfg)
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
end
