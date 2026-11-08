function report=testPhaseAwareLaminarHeating(folder)
%TESTPHASEAWARELAMINARHEATING 验证新相态热链、真实主分派和来源/沿程拒绝。
% 制造平板只作解析/软件控制，不是新箭体外缘、实验数据或物理资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
n=81;s=linspace(0,1,n)';V=zeros(3*n,3);F=reshape(1:3*n,3,n)';
for i=1:n,V(F(i,:),:)=[-s(i)-.004 0 0;-s(i)+.004 .01 0;-s(i)+.004 -.01 0];end
m=struct('V',V,'F',F);
source=struct('id','v224-prescribed-flat-plate','sha256',repmat('a',1,64), ...
    'role','prescribed-control','method','显式均匀外缘制造控制；非实测/非箭体解');
path=struct('predecessorFace',[0;(1:n-1)'],'stepLength',[0;diff(s)], ...
    'distance',s,'pathCost',s,'connectedFace',true(n,1), ...
    'meshSHA256',phaseAwareMeshIdentity(m),'sourceSHA256',source.sha256);
edge=struct('schema','rah-phase-engineering-edge-v1','model','rah-phase-aware-dry-air-of13-v1', ...
    'transportModel','rah-lemmon-dilute-icoTabulated-of13-v1', ...
    'energyReference','hf=0; hs(298.15K)=0; es=hs-RT','meshSHA256',path.meshSHA256, ...
    'originalFaceID',(1:n)','pressurePa',10000*ones(n,1),'temperatureK',220*ones(n,1), ...
    'velocityMps',[-1000*ones(n,1) zeros(n,2)],'source',source,'validFace',true(n,1), ...
    'attachedFace',true(n,1),'laminarFace',true(n,1),'quasiTwoDimensionalFace',true(n,1), ...
    'continuumLengthM',ones(n,1),'freestreamMeanFreePathM',1e-6*ones(n,1),'path',path);
m.phaseAwareEdge=edge;
cfg=struct('flight',struct('Twall',300),'heating',struct('sidewallModel','phaseAwareLaminarIntegral'));
result=calcSurfaceHeating(m,struct,cfg);b=result.phaseAwareLaminar;
assert(nnz(b.validFace)==n-1&&b.reasonCode(1)==8);
nu=b.referenceKinematicViscosityM2s(2);reference=.664./sqrt(1000*s(2:end)/nu);
flatError=max(abs(b.skinFrictionCoefficient(2:end)./reference-1));assert(flatError<2e-4);
assert(all(b.heatFluxWm2(2:end)>0)&&isnan(b.heatFluxWm2(1)));
assert(~b.physicalQualification&&~b.profileReady&&all(isnan(b.energyThicknessM)));
saveRows(result,'uniform.csv');
% 非均匀壁温/外缘路径仍由主入口整链积分，不逐点乘固定比例。
variable=m;variable.phaseAwareEdge.velocityMps(:,1)=-1000*(1+.02*s);
variable.phaseAwareEdge.temperatureK=220+10*s;
vcfg=cfg;vcfg.heating.wallTemperatureByFace=300+20*s;
varResult=calcSurfaceHeating(variable,struct,vcfg);assert(nnz(varResult.phaseAwareLaminar.validFace)==n-1);
saveRows(varResult,'variable.csv');
hot=cfg;hot.flight.Twall=1000;negative=calcSurfaceHeating(m,struct,hot);
assert(all(negative.heating.heatFlux(2:end)<0));saveRows(negative,'negative.csv');
neutral=cfg;neutral.flight.Twall=b.recoveryTemperatureK(2);
zero=calcSurfaceHeating(m,struct,neutral);assert(max(abs(zero.heating.heatFlux(2:end)))<1e-6);
% 下游不能绕过分离、域外/未知外缘、非层流或稀薄状态重新起算。
checks=struct;
cases={'separated','temperature','unknown','notLaminar','rare','research'};
for k=1:numel(cases)
    bad=m;
    switch cases{k}
        case 'separated',bad.phaseAwareEdge.velocityMps(40:end,1)=-200;expected=9;start=40;
        case 'temperature',bad.phaseAwareEdge.temperatureK(40)=50;expected=5;start=40;
        case 'unknown',bad.phaseAwareEdge.temperatureK(40)=NaN;expected=2;start=40;
        case 'notLaminar',bad.phaseAwareEdge.laminarFace(40)=false;expected=3;start=40;
        case 'rare',bad.phaseAwareEdge.continuumLengthM(40)=1e-9;expected=6;start=40;
        case 'research',bad.phaseAwareEdge.source.role='unqualified-research';expected=1;start=1;
    end
    rbad=calcSurfaceHeating(bad,struct,cfg);br=rbad.phaseAwareLaminar;
    assert(br.reasonCode(start)==expected&&all(isnan(br.heatFluxWm2(start:end))));
    if start>1,assert(all(br.reasonCode(start+1:end)==7));end
    checks.(cases{k})=true;
end
% 畸形合同必须直接拒绝，不能通过索引顺序或标量广播“修复”。
rejections=0;
for k=1:9
    bad=m;
    switch k
        case 1,bad.V(1,1)=bad.V(1,1)+1e-4;
        case 2,bad.phaseAwareEdge.originalFaceID([2 3])=[3 2];
        case 3,bad.phaseAwareEdge.energyReference='internal e/h + 300000';
        case 4,bad.phaseAwareEdge.path.pathCost(3)=bad.phaseAwareEdge.path.pathCost(2);
        case 5,bad.phaseAwareEdge.validFace=true;
        case 6,bad.phaseAwareEdge.path.connectedFace(2)=false;
        case 7,bad.phaseAwareEdge.source.role='accepted-edge';
        case 8,bad.phaseAwareEdge.path.distance(1)=.01;
        case 9,bad.phaseAwareEdge.velocityMps=-bad.phaseAwareEdge.velocityMps;
    end
    rejected=false;try,calcSurfaceHeating(bad,struct,cfg);catch,rejected=true;end
    assert(rejected);rejections=rejections+1;
end
% 主工作台数据分派必须使用新链，旧压力对象不能泄漏入本结果。
result.pressure=struct('obsolete','must not be read');
payload=buildAerothermalViewerData(result);assert(numel(payload.fields)==8);
assert(~payload.researchIdentity.productionQualified);
assert(isequaln(payload.fields(1).values,b.heatFluxWm2));
assert(all(payload.fields(6).status==2)&&all(isnan(payload.fields(7).values)));
writeJSON('workbench.json',payload);
rejected=false;try,runAeroHeatSingleCase(m,cfg);catch err,rejected=strcmp(err.identifier,'RAH:PhaseLayer:ProviderNotConnected');end
assert(rejected);
% 新线性选项不改变旧默认输出，兼容回归只覆盖实际修改边界。
old=phaseAwareAirState('PT',[1000;10000],[180;700]);same=phaseAwareAirState('PT',[1000;10000],[180;700],'lemmonAnalytic');
assert(isequaln(old,same));
testBoundaryLayerState();
report=struct('passed',true,'faces',n,'flatPlateCfRelativeDifference',flatError, ...
    'blockerCases',checks,'contractRejections',rejections,'mainEntryUsed',true,'viewerMainEntryUsed',true, ...
    'legacyDefaultCompatible',true,'singleCaseProviderStillRejected',true,'profileReady',false, ...
    'physicalQualification',false,'newQualifiedAreaM2',0);
writeJSON('report.json',report);disp(report);
    function saveRows(mout,name)
        bb=mout.phaseAwareLaminar;ee=mout.phaseAwareEdge;
        rows=[s ee.pressurePa ee.temperatureK sqrt(sum(ee.velocityMps.^2,2)) bb.wallTemperatureK ...
            bb.referenceTemperatureK bb.referenceEnthalpyJkg bb.recoveryTemperatureK bb.thermalDrivingEnthalpyJkg ...
            bb.referencePrandtl bb.referenceKinematicViscosityM2s bb.referenceMomentumThicknessM ...
            bb.skinFrictionCoefficient bb.heatFluxWm2 double(bb.reasonCode)];
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');clean=onCleanup(@()fclose(fid));
        fprintf(fid,'s,p,Te,U,Tw,Tref,href,Taw,dh,Prref,nuref,thetaRef,Cf,q,reason\n');
        fprintf(fid,[repmat('%.17g,',1,14) '%.17g\n'],rows');
    end
    function writeJSON(name,value)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');clean=onCleanup(@()fclose(fid));
        fprintf(fid,'%s\n',jsonencode(value,'ConvertInfAndNaN',true));
    end
end
