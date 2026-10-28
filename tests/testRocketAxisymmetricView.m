function report=testRocketAxisymmetricView(folder,sourceFolder)
%TESTROCKETAXISYMMETRICVIEW 核对真实完整箭体端点的离线显示及原面投影。
% 模块一：只读取真实不可变检查点，完整逐步BE链在导出前重放。
% 模块二：核对共享鼻尖/底心闭合、全部真实墙段、signed热流和真实曲面功率。
% 模块三：原生V/F投影保持、加权p/q及给定Tw分开，未知温度继续NaN。
% 模块四：伪物理解、初猜、假历史、改Q及漏墙拒绝，不创建假端点。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v68_validation','rocket_view');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v68_validation','actual_rocket');end
path=fullfile(sourceFolder,'rocket_axisymmetric_evidence.mat');sourceSHA=fileSHA(path);
source=load(path,'caseDefinition','implicit');item=source.caseDefinition;result=source.implicit;
opts=struct('azimuthSegments',single(32),'checkpointFolder',sourceFolder,'evidenceSource',path);
view=exportRocketAxisymmetricStudy(item,result,folder,opts);
snapshot=load(view.snapshotFile);data=snapshot.data;visual=snapshot.visualMesh;
assert(isequaln(data.numericalSolver,rocketAxisymmetricSolverContract(item))&&isequaln(data.initialization,item.initialization));
assert(isequaln(snapshot.item,item)&&isequaln(snapshot.result,result)&& ...
    view.checkpointValidation.acceptedEndpointsReplayed==result.acceptedSteps&& ...
    view.checkpointValidation.immutableChainProvided&&view.spatialResidualReplayPassed&& ...
    ~view.physicalSolution&&~view.eligibleForSurfaceCoverage&&view.newProductionSolvedAreaM2==0);
identity=data.researchIdentity;
% 受控链必须逐字保留实际离散策略和原始审计，防止显示时混为旧链或改算法标签。
if isfield(item,'caseFactory')
    assert(strcmp(identity.caseFactory,item.caseFactory)&& ...
        strcmp(identity.wallVelocityGradientPolicy,item.inputContract.wallVelocityGradientPolicy)&& ...
        strcmp(identity.viscousWorkVelocityPolicy,item.inputContract.viscousWorkVelocityPolicy));
    for auditName={'wallVelocityGradient','viscousWorkVelocity'}
        key=auditName{1};
        assert(isfield(data,key)==isfield(result.finalInfo,key));
        if isfield(data,key),assert(isequaln(data.(key),result.finalInfo.(key)));end
    end
    if isfield(item.boundary,'pressureSourceQuadrature')
        assert(strcmp(identity.pressureSourceQuadraturePolicy,item.boundary.pressureSourceQuadrature)&& ...
            isequaln(data.pressureSourceQuadrature,result.finalInfo.pressureSourceQuadrature));
    end
    if isfield(item.boundary,'interiorConvectiveFlux')
        assert(strcmp(identity.interiorConvectiveFluxPolicy,item.boundary.interiorConvectiveFlux)&& ...
            isequaln(data.interiorConvectiveFlux,result.finalInfo.interiorConvectiveFlux));
        assert(~contains(result.finalInfo.method,'first-order HLLC')||strcmp(item.boundary.interiorConvectiveFlux,'hllc'));
    end
    if isfield(item.inputContract,'gridPolicy')
        assert(strcmp(identity.gridPolicy,item.inputContract.gridPolicy)&& ...
            isequaln(identity.effectiveCollarOptions,item.inputContract.effectiveCollarOptions));
        if isfield(item.domain,'collar'),assert(isequaln(data.normalCollarGeometry,item.domain.collar));end
    end
end
assert(identity.acceptedSteps==result.acceptedSteps&&identity.steadyConverged==result.steadyConverged&& ...
    ~identity.physicalSolution&&~identity.experimentalValidation&&~identity.pseudoTimeIsPhysical&& ...
    identity.physicalTimeAdvancedS==0&&identity.nativeOriginalTriangleProjectionCompleted&& ...
    ~identity.displayGeometryIsOriginalNativeTriangles);
w=result.finalInfo.wall;b=data.wallBinding;map=b.sourceWallIndex;area=data.faceArea;
q=field(data,'heatFlux');p=field(data,'pressure');tw=field(data,'wallInput');
gasProxy=field(data,'wallAdjacentGas');
ownerQ=result.Q(w.ownerCell,:);rho=ownerQ(:,1);velocity=ownerQ(:,2:3)./rho;
temperature=(ownerQ(:,4)./rho-sum(velocity.^2,2)/2)/item.gas.cv;
assert(isequal(gasProxy.values,temperature(map))&&all(gasProxy.status==1)&& ...
    isequal(data.wallBinding.sourceWallAdjacentGasProxyK,temperature(map))&& ...
    contains(gasProxy.methodNames{1},'不是壁温、外缘温度或恢复温度'));
assert(all(abs(q.values*1000-w.heatFluxIntoWall(map))<=256*eps*max(abs(w.heatFluxIntoWall(map)),realmin))&& ...
    all(abs(p.values*1000-w.pressure(map))<=256*eps*abs(w.pressure(map)))&& ...
    isequal(tw.values,w.temperatureInput(map))&&all(q.status==1)&&all(p.status==1)&&all(tw.status==2)&& ...
    isequal(sign(q.values),sign(w.heatFluxIntoWall(map))),'RocketViewTest:Fields','signed源场量或给定/计算状态不同。');
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    f=field(data,key{1});assert(all(isnan(f.values))&&all(f.status==0)&&~isempty(f.methodNames));
end
for k=1:numel(data.fields)
    f=data.fields(k);assert(~isempty(f.methodNames)&&~isempty(f.scope)&& ...
        numel(f.values)==size(visual.F,1)&&all(f.methodIndex==0));
end
assert(isequal(sort(unique(b.sourceFaceID)),sort(item.domain.wallBoundaryFaceID))&& ...
    view.audit.geometry.closedTwoManifold&&view.audit.geometry.axisPoleCount==2&& ...
    ~view.audit.fictitiousEndcapsAdded&&view.audit.fullBaseAreaM2==item.domain.fullBaseAreaM2);
edges=sort([visual.F(:,[1 2]);visual.F(:,[2 3]);visual.F(:,[3 1])],2);
[~,~,eid]=unique(edges,'rows');assert(all(accumarray(eid,1)==2));
trueHeat=2*pi*sum(w.areaPerRadian.*w.heatFluxIntoWall);
assert(view.audit.revolvedWallHeatW==trueHeat&& ...
    view.audit.displayFacetedHeatIntegralW==sum(area.*w.heatFluxIntoWall(map))&& ...
    ~view.audit.signedHeatFluxWasRenormalized);

%% 二、完整原生三角面保守曲面积投影；与细分旋转显示分开核对
projection=snapshot.projection;native=data.sourceSurfaceProjection;
[present,wallIndex]=ismember(item.domain.wallBoundaryFaceID,w.faceID);assert(all(present));
assert(isequal(native.vertices,item.surfaceMesh.V)&&isequal(native.faces,item.surfaceMesh.F-1)&& ...
    native.indexBase==0&&native.originalSourceFacesPreserved&&~native.candidateMappingOnly&& ...
    ~native.pressureForceOrMomentConservationEstablished);
assert(isequal(native.pressurePa,projection.weights*w.pressure(wallIndex))&& ...
    isequal(native.heatFluxIntoWallWm2,projection.weights*w.heatFluxIntoWall(wallIndex))&& ...
    isequal(native.wallAdjacentGasProxyK,projection.weights*temperature(wallIndex))&& ...
    all(native.wallAdjacentGasProxyStatus==1)&&~native.gasTemperatureProjectionProvesEnergyConservation&& ...
    isequal(native.wallInputK,projection.weights*w.temperatureInput(wallIndex))&& ...
    all(native.pressureStatus==1)&&all(native.heatFluxStatus==1)&&all(native.wallInputStatus==2)&& ...
    all(isnan(native.wallSolvedK))&&all(isnan(native.edgeTemperatureK))&&all(isnan(native.recoveryTemperatureK)));
assert(abs(native.completeCurvedHeatW-trueHeat)<=2e-11*sum(2*pi*w.areaPerRadian.*abs(w.heatFluxIntoWall)));
reconstructed=sparse(native.sourceFaceIndex+1,native.wallSegmentIndex+1,native.weight, ...
    size(item.surfaceMesh.F,1),numel(item.domain.wallBoundaryFaceID));
assert(isequal(reconstructed,projection.weights));

%% 三、每个非法输入必须在写页面前拒绝；不以新页覆盖任何失败身份
refusals=struct('label',{},'identifier',{});
bad=result;bad.acceptedSteps=0;reject('zero-step-initial-guess',item,bad,opts,'');
bad=result;bad.physicalSolution=true;reject('false-physical-qualification',item,bad,opts,'RocketAxisView:Identity');
bad=result;bad.initialQ(1,1)=bad.initialQ(1,1)*(1+1e-6);reject('changed-initial-state',item,bad,opts,'RocketAxisView:Identity');
% 选真实最大体积单元，使篡改能被全局积分解析；鼻尖极小格的改变量可能低于全域求和舍入。
[~,massProbeCell]=max(item.mesh.volumePerRadian);
bad=result;bad.Q(massProbeCell,1)=bad.Q(massProbeCell,1)*(1+1e-6);reject('changed-final-state',item,bad,opts,'RocketAxisView:Ledger');
bad=result;bad.finalInfo.wall.faceID=bad.finalInfo.wall.faceID(2:end);reject('missing-wall-info',item,bad,opts,'RocketAxisView:Checkpoint');
bad=result;bad.history.newtonTarget(end)=bad.history.newtonTarget(end)*2;reject('forged-newton-target',item,bad,opts,'RocketAxisView:Checkpoint');
bad=result;bad.time=1;reject('pseudo-time-as-physical-time',item,bad,opts,'RocketAxisView:Identity');
bad=item;bad.boundary.jMin.temperature=bad.boundary.jMin.temperature+1;reject('changed-wall-temperature',bad,result,opts,'RocketAxisView:Case');
if isfield(item,'caseFactory')
    bad=item;bad.caseFactory='unregisteredFactory';reject('unregistered-case-factory',bad,result,opts,'RocketAxisView:Case');
    bad=item;bad.inputContract.viscousWorkVelocityPolicy='legacy-wls';
    if strcmp(item.inputContract.viscousWorkVelocityPolicy,'legacy-wls'),bad.inputContract.viscousWorkVelocityPolicy='arithmetic';end
    reject('forged-work-velocity-policy',bad,result,opts,'RocketAxisView:Case');
    if isfield(item.inputContract,'gridPolicy')
        bad=item;bad.inputContract.gridPolicy='forged-grid';reject('forged-grid-policy',bad,result,opts,'RocketAxisView:Case');
    end
end
if result.acceptedSteps>1
    noChain=rmfield(opts,'checkpointFolder');reject('multistep-without-immutable-chain',item,result,noChain,'RocketAxisView:Checkpoint');
end
one=load(fullfile(sourceFolder,'steps','step_000001.mat'),'segment');singleResult=one.segment;
singleOptions=struct('azimuthSegments',8);
bad=singleResult;bad.history.newtonTarget(end)=bad.history.newtonTarget(end)*2;
reject('single-step-forged-newton-target',item,bad,singleOptions,'RocketAxisView:BE');
totalScale=sum(item.mesh.volumePerRadian)*item.residualScale;
bad=singleResult;bad.accumulatedBoundaryIntegral(1)=bad.accumulatedBoundaryIntegral(1)+1e-11*totalScale(1);
bad.globalBalanceRelative=max(abs(bad.finalConservedIntegral-bad.initialConservedIntegral+ ...
    bad.accumulatedBoundaryIntegral-bad.accumulatedSourceIntegral)./totalScale);
bad.history.balance(end)=bad.globalBalanceRelative;
reject('single-step-summary-ledger-with-false-flux-integral',item,bad,singleOptions,'RocketAxisView:BE');
bad=singleResult;bad.history.stepDefect(end,1)=bad.history.stepDefect(end,1)+1e-11*totalScale(1);
reject('single-step-forged-subtolerance-ledger',item,bad,singleOptions,'RocketAxisView:Checkpoint');
assert(strcmp(sourceSHA,fileSHA(path)),'RocketViewTest:Source','只读导出改变了原证据。');
report=struct('status','passed-complete-rocket-replayed-offline-view-and-native-projection', ...
    'view',view,'refusals',refusals,'sourceEvidenceSHA256',sourceSHA, ...
    'acceptedEndpointsReplayed',result.acceptedSteps,'displayFaceCount',size(visual.F,1), ...
    'nativeFaceCount',size(item.surfaceMesh.F,1),'signedHeatPreserved',true,'unknownTemperaturesPreserved',true, ...
    'browserVisualOrWebGLRuntimeVerified',false,'physicalSolution',false,'experimentalValidation',false, ...
    'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'view_test_evidence.mat'),'report');
fid=fopen(fullfile(folder,'view_test_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Rocket research view: %d real endpoints, %d display/%d native faces, %d refusals; visual runtime not verified.\n', ...
    result.acceptedSteps,size(visual.F,1),size(item.surfaceMesh.F,1),numel(refusals));
    function reject(label,c,r,o,identifier)
        target=fullfile(folder,['rejected_' label]);
        assert(~exist(target,'dir'),'RocketViewTest:Overwrite','拒绝测试目录已存在。');
        try,exportRocketAxisymmetricStudy(c,r,target,o);catch exception
            assert(isempty(identifier)||strcmp(exception.identifier,identifier), ...
                'RocketViewTest:Refusal','%s得到非预期拒绝%s。',label,exception.identifier);
            assert(~exist(target,'dir'),'RocketViewTest:Refusal','拒绝前不应写页面。');
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('RocketViewTest:Refusal','%s非法输入未拒绝。',label);
    end
end

function f=field(data,key)
f=data.fields(strcmp({data.fields.key},key));assert(isscalar(f));
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
