function report=testAxisymmetricRiemannBoundary(folder,sourceFile)
%TESTAXISYMMETRICRIEMANNBOUNDARY 无旋四量储库HLLE边界及真实出口声速转换校核。
% 模块一：正反向、亚/超声速均匀储库，核对度量、共享通量和四量账。
% 模块二：缺失/多量/非法状态必须拒绝；不暗补旋转速度或外部物性。
% 模块三：读取真实TMR第一次35步终态，仅改变出口一个owner的法向速度，
% 保持rho、T及切向速度，检查Mach=1两侧有效性、通量连续及旧边界拒绝ID。
% 模块四：保存原始证据和来源哈希；这是数值边界校核，无新增物理解面积。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','riemann_boundary');end
if nargin<2
    sourceFile=fullfile(root,'outputs','v63_validation','tmr_implicit_first','tmr_implicit_budget.mat');
    if ~exist(sourceFile,'file')
        sourceFile=fullfile(root,'validation','v63','tmr_implicit_first','tmr_implicit_budget.mat');
    end
end
assert(exist(sourceFile,'file')==2,'AxisFVTest:RiemannSource','缺少实际35步TMR终态；不能用均匀场替代真实转换校核。');
if ~exist(folder,'dir'),mkdir(folder);end

%% 一、同状态储库在正反向、亚声速和超声速均匀场上保持真实有限体积恒等式
[X,R]=ndgrid([.1 .05 0],[.03 .04 .06]);mesh=buildAxisymmetricFVGeometry(X,R);
gas=axisymmetricFVGas(struct('gamma',1.4,'R',287.05,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]));
velocities=[-1100 1100 -110 110];uniformRecords=struct([]);uniformInfos=cell(size(velocities));
metric=max(abs(mesh.incidence*mesh.faces.areaVectorPerRadian-[zeros(mesh.cellCount,1),mesh.planarArea])./mesh.planarArea,[],'all');
assert(metric<1e-11,'AxisFVTest:RiemannMetric','均匀场试验的轴对称度量恒等式不成立。');
for j=1:numel(velocities)
    primitive=[1 velocities(j) 0 300];Q=axisymmetricFVConservative(repmat(primitive,mesh.cellCount,1),gas);
    far=struct('type','riemannReservoir','primitive',primitive);
    boundary=struct('iMin',far,'iMax',far,'jMin',far,'jMax',far);
    [rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);P=axisymmetricFVState(Q,gas);
    assert(strcmp(info.boundaryRiemannSolver,'HLLE')&&contains(info.method,'open-face HLLE'), ...
        'AxisFVTest:RiemannMethod','开放储库必须明确记录统一HLLE数值方法。');
    scale=[1,max(P(:,6)),max(P(:,6)),max(Q(:,4))];
    change=max(abs(rate)./scale./(info.convectiveFrequency+info.diffusiveFrequency),[],'all');
    assert(change<2e-12,'AxisFVTest:RiemannUniform','同状态储库改变了均匀解。');
    assert(all(info.faceBoundaryKind(mesh.boundaryFaces)==7),'AxisFVTest:RiemannKind','新储库必须具有独立边界身份7。');
    physicalFlux=physical(Q(mesh.faces.owner,:),P(mesh.faces.owner,:),mesh.faces.normal);
    fluxScale=max(1,max(abs(physicalFlux.*mesh.faces.areaPerRadian),[],'all'));
    convectiveError=max(abs(info.faceConvectiveFluxPerRadian-physicalFlux.*mesh.faces.areaPerRadian),[],'all')/fluxScale;
    fluxError=max(abs(info.faceFluxIntegralPerRadian-physicalFlux.*mesh.faces.areaPerRadian),[],'all')/fluxScale;
    assert(fluxError<2e-12,'AxisFVTest:RiemannUniformFlux','同状态储库的真实总通量不等于均匀物理通量。');
    assert(convectiveError<2e-12&&max(abs(info.faceViscousFluxPerRadian),[],'all')/fluxScale<2e-12);
    assert(max(abs(info.faceFluxIntegralPerRadian-(info.faceConvectiveFluxPerRadian-info.faceViscousFluxPerRadian)),[],'all')/fluxScale<2e-12);
    speed=abs(sum(P(mesh.faces.owner,2:3).*mesh.faces.normal,2))+P(mesh.faces.owner,6);
    expectedFrequency=(mesh.cellFaceMultiplicity*(mesh.faces.areaPerRadian.*speed))./mesh.volumePerRadian;
    frequencyError=max(abs(info.convectiveFrequency-expectedFrequency)./expectedFrequency);
    assert(frequencyError<2e-12,'AxisFVTest:RiemannFrequency','均匀HLLE波速频率不匹配实际两侧状态。');
    assert(max(abs(info.faceWaveSpeed-speed)./speed)<2e-12);
    expectedIncoming=info.faceConvectiveFluxPerRadian(mesh.boundaryFaces,1)<0;
    assert(isequal(info.faceIncoming(mesh.boundaryFaces),expectedIncoming)&& ...
        isequal(info.faceNormalZeroGradient(mesh.boundaryFaces),~expectedIncoming));
    balance=max(abs(info.conservationDefect)./max(1,abs(info.sourceIntegral)));
    assert(balance<1e-12&&~info.physicalSolution&&~info.eligibleForSurfaceCoverage);
    item=struct('velocityX',velocities(j),'normalMachMagnitude',abs(velocities(j))/P(1,6), ...
        'maximumLocalStepChange',change,'metricRelative',metric,'integratedFaceFluxRelative',fluxError, ...
        'convectiveFrequencyRelative',frequencyError,'convectiveFluxRelative',convectiveError,'fourQuantityDefectRelative',balance, ...
        'boundaryFlux',info.netBoundaryFlux,'axisSource',info.sourceIntegral);
    if isempty(uniformRecords),uniformRecords=item;else,uniformRecords(end+1)=item;end %#ok<AGROW>
    uniformInfos{j}=info;
end

%% 一附、黏性策略按实际HLLE质量方向，验证内点速度反向和恰零质量的分支
directionRecords=struct([]);directionInfos=cell(1,3);
directionInputs={ [1 -30 0 300], [3 30 0 450], [1 0 0 300] };
directionExterior={ [3 100 0 450], [1 -100 0 300], [1 0 0 300] };
for j=1:3
    cellPrimitive=directionInputs{j};exterior=directionExterior{j};
    Qdirection=axisymmetricFVConservative(repmat(cellPrimitive,mesh.cellCount,1),gas);
    default=struct('type','riemannReservoir','primitive',cellPrimitive);
    directedBoundary=struct('iMin',default,'iMax',struct('type','riemannReservoir','primitive',exterior),'jMin',default,'jMax',default);
    [~,directionInfo]=axisymmetricNavierStokesResidual(Qdirection,mesh,gas,directedBoundary);
    faceIDs=mesh.boundaryFaces(mesh.faces.boundarySide(mesh.boundaryFaces)==2);testFace=faceIDs(1);
    normal=mesh.faces.normal(testFace,:);L=axisymmetricFVState(Qdirection(mesh.faces.owner(testFace),:),gas);
    Uexterior=axisymmetricFVConservative(exterior,gas);Rstate=axisymmetricFVState(Uexterior,gas);
    [Fc,waves,middle]=referenceHLLE(Qdirection(mesh.faces.owner(testFace),:),Uexterior,L,Rstate,normal);
    integratedReference=Fc*mesh.faces.areaPerRadian(testFace);
    assert(max(abs(directionInfo.faceConvectiveFluxPerRadian(testFace,:)-integratedReference)./max(1,abs(integratedReference)))<2e-12);
    mass=integratedReference(1);incoming=mass<0;
    assert(directionInfo.faceIncoming(testFace)==incoming&&directionInfo.faceNormalZeroGradient(testFace)==~incoming);
    assert(max(abs(directionInfo.faceRiemannExteriorState(testFace,:)-Rstate)./max(1,abs(Rstate)))<2e-12);
    if incoming,expectedTrace=Rstate;else,expectedTrace=L;end
    traceError=max(abs(directionInfo.faceViscousGradientTrace(testFace,:)-expectedTrace)./max(1,abs(expectedTrace)));
    assert(traceError<2e-12,'AxisFVTest:RiemannViscousTrace','黏性迹不符合实际HLLE质量方向。');
    cellNormalVelocity=dot(cellPrimitive(2:3),normal);
    if j<3
        assert(sign(mass)==-sign(cellNormalVelocity),'AxisFVTest:RiemannMassDirection','非均匀储库试验未覆盖质量通量与内点法向速度反向。');
    else
        assert(mass==0&&~incoming&&directionInfo.faceNormalZeroGradient(testFace),'AxisFVTest:RiemannZeroMass','恰零质量必须使用出流零法向梯度策略。');
    end
    assert(middle.used&&middle.density>0&&middle.internalEnergyPerVolume>0, ...
        'AxisFVTest:RiemannIntermediate','质量方向试验必须实际取用正密度/正内能的HLL中间态。');
    expectedSpeed=max(abs(waves));
    assert(abs(directionInfo.faceWaveSpeed(testFace)-expectedSpeed)/expectedSpeed<2e-12);
    item=struct('cellPrimitive',cellPrimitive,'reservoirPrimitive',exterior,'ownerNormalVelocity',cellNormalVelocity, ...
        'massFluxPerRadian',mass,'actualMassIsIncoming',incoming,'viscousTraceRelative',traceError, ...
        'hlleWaveSpeeds',waves,'hllIntermediateAudit',middle,'convectiveFluxIntegralPerRadian',integratedReference);
    if isempty(directionRecords),directionRecords=item;else,directionRecords(end+1)=item;end %#ok<AGROW>
    directionInfos{j}=directionInfo;
end

%% 二、四量储库是明确合同；缺失、多量旋转输入、非物理及温域外状态均不得暗补
refusals=struct('label',{},'identifier',{});goodBoundary=boundary;
invalidBoundary=goodBoundary;invalidBoundary.iMax=struct('type','riemannReservoir');
reject('missing-primitive',@()axisymmetricNavierStokesResidual(Q,mesh,gas,invalidBoundary));
invalidValues={ [1 110 0 0 300], [1 110 0 1 300], [0 110 0 300], [-1 110 0 300], ...
    [1 NaN 0 300], [1 110 0 79], [1 110 0 2001], [1 110 300] };
invalidLabels={'five-variable-zero-swirl-contract','five-variable-nonzero-swirl-contract','zero-density','negative-density', ...
    'nonfinite-velocity','below-temperature-range','above-temperature-range','incomplete-four-state'};
for j=1:numel(invalidValues)
    invalidBoundary=goodBoundary;invalidBoundary.iMax=struct('type','riemannReservoir','primitive',invalidValues{j});
    reject(invalidLabels{j},@()axisymmetricNavierStokesResidual(Q,mesh,gas,invalidBoundary));
end

%% 三、实际TMR受阻终态上的出口跨声速：保存旧特征guard，独立比较新HLLE数值通量
source=load(sourceFile,'implicit','mesh','gas','boundary','scale');
assert(all(isfield(source,{'implicit','mesh','gas','boundary','scale'}))&& ...
    source.mesh.cellCount==1000&&source.implicit.acceptedSteps==35&& ...
    ~source.implicit.physicalSolution&&~source.implicit.eligibleForSurfaceCoverage, ...
    'AxisFVTest:RiemannSource','必须使用实际35步受阻终态，而不是构造场或已验收物理解。');
tmrMesh=source.mesh;tmrGas=axisymmetricFVGas(source.gas);originalBoundary=source.boundary;
originalQ=source.implicit.Q;originalP=axisymmetricFVState(originalQ,tmrGas);
faces=tmrMesh.faces;exitFaces=tmrMesh.boundaryFaces(faces.boundarySide(tmrMesh.boundaryFaces)==2);
normalMach=sum(originalP(faces.owner(exitFaces),2:3).*faces.normal(exitFaces,:),2)./originalP(faces.owner(exitFaces),6);
[~,nearest]=min(abs(normalMach-1));face=exitFaces(nearest);cellID=faces.owner(face);normal=faces.normal(face,:);
assert(face==81&&cellID==40&&abs(normalMach(nearest)-1)<1e-8, ...
    'AxisFVTest:RiemannSource','真实受阻终态的近声速出口身份已改变，须重新审查来源。');
oldReservoir=originalBoundary.iMax.reservoirPrimitive;
assert(isnumeric(oldReservoir)&&numel(oldReservoir)==5&&oldReservoir(4)==0, ...
    'AxisFVTest:RiemannSource','真实试验需要显式冻结的无旋旧储库。');
reservoir=oldReservoir([1 2 3 5]);newBoundary=originalBoundary;
newBoundary.iMax=struct('type','riemannReservoir','primitive',reservoir);
offsets=[1e-4 1e-6 1e-8 -1e-8 -1e-6 -1e-4];probeRecords=struct([]);
probeQ=cell(size(offsets));probeInfo=cell(size(offsets));legacyInfo=cell(size(offsets));
convectiveExpected=zeros(numel(offsets),4);totalFlux=zeros(numel(offsets),4);
reservoirState=axisymmetricFVState(axisymmetricFVConservative(reservoir,tmrGas),tmrGas);
for j=1:numel(offsets)
    primitive=originalP(:,1:4);un=dot(primitive(cellID,2:3),normal);a=originalP(cellID,6);
    primitive(cellID,2:3)=primitive(cellID,2:3)+((1+offsets(j))*a-un)*normal;
    trial=axisymmetricFVConservative(primitive,tmrGas);P=axisymmetricFVState(trial,tmrGas);
    assert(max(abs(P(:,1)-originalP(:,1)))==0&& ...
        max(abs(P(:,4)-originalP(:,4)))<=128*eps(max(originalP(:,4))), ...
        'AxisFVTest:RiemannPerturbation','声速探测只能改变一个owner的法向速度，rho/T不可改变。');
    assert(max(abs(P(setdiff(1:tmrMesh.cellCount,cellID),:)-originalP(setdiff(1:tmrMesh.cellCount,cellID),:)),[],'all')< ...
        128*eps(max(abs(originalP),[],'all')));
    [rate,newInfo]=axisymmetricNavierStokesResidual(trial,tmrMesh,tmrGas,newBoundary);
    assert(all(isfinite(rate),'all')&&all(newInfo.faceBoundaryKind(exitFaces)==7)&& ...
        ~newInfo.physicalSolution&&~newInfo.eligibleForSurfaceCoverage);
    conservation=conservationAudit(rate,newInfo,tmrMesh);
    balance=conservation.maximumRelativeDefect;
    assert(conservation.normalizedPassed&&conservation.roundoffPassed&&conservation.reportedDefectAgrees, ...
        'AxisFVTest:RiemannBalance','真实转换探测的四量账不符合逐量物理尺度或显式浮点求和界。');
    % 负对照故意在rate中漏掉柱坐标源，保持边界和账本来源不变；两个验收均须拒绝。
    rateWithoutAxisSource=rate-newInfo.cellSourceIntegral./tmrMesh.volumePerRadian;
    missingSource=conservationAudit(rateWithoutAxisSource,newInfo,tmrMesh);
    assert(~missingSource.normalizedPassed&&~missingSource.roundoffPassed&& ...
        any(newInfo.sourceIntegral~=0), ...
        'AxisFVTest:RiemannBalanceSensitivity','漏掉真实柱坐标源的负对照未被守恒审计拒绝。');
    [expected,waves,middle]=referenceHLLE(trial(cellID,:),axisymmetricFVConservative(reservoir,tmrGas), ...
        P(cellID,:),reservoirState,normal);
    convectiveExpected(j,:)=expected*faces.areaPerRadian(face);totalFlux(j,:)=newInfo.faceFluxIntegralPerRadian(face,:);
    convectiveError=max(abs(newInfo.faceConvectiveFluxPerRadian(face,:)-convectiveExpected(j,:))./max(1,abs(convectiveExpected(j,:))));
    assert(convectiveError<2e-12&&strcmp(newInfo.boundaryRiemannSolver,'HLLE'), ...
        'AxisFVTest:RiemannNumericalFlux','真实转换处未使用统一开放面HLLE数值通量。');
    assert(max(abs(newInfo.faceFluxIntegralPerRadian(face,:)-(newInfo.faceConvectiveFluxPerRadian(face,:)-newInfo.faceViscousFluxPerRadian(face,:)))./max(1,abs(totalFlux(j,:))))<2e-12);
    massIsIncoming=convectiveExpected(j,1)<0;
    assert(newInfo.faceIncoming(face)==massIsIncoming&&newInfo.faceNormalZeroGradient(face)==~massIsIncoming);
    if offsets(j)<0
        assert(middle.used&&middle.density>0&&middle.internalEnergyPerVolume>0, ...
            'AxisFVTest:RiemannIntermediate','声速以下必须实际取用正密度/正内能HLL中间态。');
    else
        assert(~middle.used,'AxisFVTest:RiemannIntermediate','两波均向外时不得假装取用未使用的中间态。');
    end
    expectedSpeed=max(abs(waves));
    assert(abs(newInfo.faceWaveSpeed(face)-expectedSpeed)/expectedSpeed<2e-12, ...
        'AxisFVTest:RiemannWaves','真实储库频率没有包括两侧HLLE实际波速。');
    legacyIdentifier='';legacyMessage='';legacySucceeded=false;legacyAgreement=[];
    try
        [~,oldInfo]=axisymmetricNavierStokesResidual(trial,tmrMesh,tmrGas,originalBoundary);
        legacySucceeded=true;legacyInfo{j}=oldInfo;
    catch exception
        assert(startsWith(exception.identifier,'RocketAeroHeat:NSDraft:Characteristic'), ...
            'AxisFVTest:RiemannLegacyError','旧边界出现未分类代码错误：%s。',exception.identifier);
        legacyIdentifier=exception.identifier;legacyMessage=exception.message;
    end
    if offsets(j)<0
        assert(~legacySucceeded&&strcmp(legacyIdentifier,'RocketAeroHeat:NSDraft:CharacteristicRegime'), ...
            'AxisFVTest:RiemannLegacyGuard','冻结旧特征跨分支guard必须继续明确拒绝，不能删保护伪装修复。');
    end
    if offsets(j)>0
        assert(legacySucceeded,'AxisFVTest:RiemannLegacySupersonic','声速以上原特征超声速出流应可计算。');
        legacyAgreement=max(abs(oldInfo.faceFluxIntegralPerRadian(face,:)-totalFlux(j,:))./ ...
            max(1,abs(oldInfo.faceFluxIntegralPerRadian(face,:))));
        assert(legacyAgreement<2e-12,'AxisFVTest:RiemannSupersonic','新储库在此超声速出流面与原QL物理边界不一致。');
        direct=physical(trial(cellID,:),P(cellID,:),normal)*faces.areaPerRadian(face);
        assert(max(abs(convectiveExpected(j,:)-direct)./max(1,abs(direct)))<2e-12, ...
            'AxisFVTest:RiemannSupersonic','实际两侧波均向外时，HLLE应退化为QL物理通量。');
    end
    item=struct('normalMach',sum(P(cellID,2:3).*normal)/P(cellID,6), ...
        'requestedOffset',offsets(j),'newOperatorSucceeded',true,'boundaryKind',newInfo.faceBoundaryKind(face), ...
        'originalNSResidual',max(abs(rate)./source.scale,[],'all'),'fourQuantityDefectRelative',balance, ...
        'conservationAudit',conservation,'omittedAxisSourceNegativeControl',missingSource, ...
        'hlleReferenceFluxIntegralPerRadian',convectiveExpected(j,:),'hlleWaveSpeeds',waves, ...
        'hllIntermediateAudit',middle,'hlleFluxRelative',convectiveError, ...
        'totalFluxIntegralPerRadian',totalFlux(j,:),'legacySucceeded',legacySucceeded, ...
        'legacyIdentifier',legacyIdentifier,'legacyMessage',legacyMessage,'legacySupersonicFluxRelative',legacyAgreement);
    if isempty(probeRecords),probeRecords=item;else,probeRecords(end+1)=item;end %#ok<AGROW>
    probeQ{j}=trial;probeInfo{j}=newInfo;
end
% 以rho、a及真实r加权面积建立维度一致尺度；检查相同扰动两侧的通量差随扰动缩小。
rho=originalP(cellID,1);a=originalP(cellID,6);
fluxScale=faces.areaPerRadian(face)*rho*[a a^2 a^2 a^3];pairs=[1 6;2 5;3 4];
totalGaps=zeros(1,3);referenceGaps=totalGaps;
for j=1:3
    totalGaps(j)=max(abs(totalFlux(pairs(j,1),:)-totalFlux(pairs(j,2),:))./fluxScale);
    referenceGaps(j)=max(abs(convectiveExpected(pairs(j,1),:)-convectiveExpected(pairs(j,2),:))./fluxScale);
end
assert(all(totalGaps<200*[1e-4 1e-6 1e-8])&&all(referenceGaps<200*[1e-4 1e-6 1e-8])&& ...
    all(totalGaps(2:3)<.03*totalGaps(1:2))&&all(referenceGaps(2:3)<.03*referenceGaps(1:2)), ...
    'AxisFVTest:RiemannContinuity','实际出口Mach=1两侧通量没有随扰动连续收敛。');

%% 四、保存数值身份、来源和真实终态；不把开放边界校核算为区域物理解
bindings=sourceBindings(root);
report=struct('status','passed-four-variable-HLLE-reservoir-uniform-and-actual-TMR-sonic-probes', ...
    'uniformRecords',uniformRecords,'massDirectionRecords',directionRecords,'refusals',refusals,'rejections',numel(refusals), ...
    'sourceEvidenceSHA256',fileSHA(sourceFile),'sourceAcceptedSteps',source.implicit.acceptedSteps, ...
    'sourceFace',face,'sourceCell',cellID,'sourceNormalMach',normalMach(nearest), ...
    'prescribedReservoirPrimitive',reservoir,'probeRecords',probeRecords, ...
    'continuityOffsets',[1e-4 1e-6 1e-8],'totalFluxScaledGaps',totalGaps, ...
    'independentHLLEScaledGaps',referenceGaps,'sourceBindings',bindings, ...
    'method','Explicit 4-variable exterior reservoir; uniform open-face HLLE two-wave flux, positive actually used HLL intermediate state; interior HLLC unchanged; viscous inflow reservoir/outflow zero normal gradient', ...
    'scope','Uniform traces and six probes of actual unconverged 35-step TMR state; not steady inlet/domain/turbulence/experimental acceptance', ...
    'steadyConverged',false,'solverRun',false,'actualBoundaryLayerInletValidated',false,'experimentalValidation',false, ...
    'regionalFlowValidated',false,'methodEligibleForProduction',false,'fullSurfaceSolved',false, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'riemann_boundary_evidence.mat'),'report','mesh','gas','uniformInfos','directionInfos', ...
    'tmrMesh','tmrGas','originalQ','originalBoundary','newBoundary','probeQ','probeInfo','legacyInfo','convectiveExpected','totalFlux');
fid=fopen(fullfile(folder,'riemann_boundary_summary.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Axisymmetric Riemann reservoir passed: 4 uniform traces, %d refusals, actual TMR face%d Mach1 +/- 1e-4/1e-6/1e-8; no production area.\n',numel(refusals),face);

    function reject(label,action)
        try
            action();
        catch exception
            assert(startsWith(exception.identifier,{'AxisFV:','MATLAB:expected','MATLAB:invalidType','MATLAB:validators:'}), ...
                'AxisFVTest:RiemannUnexpectedRefusal','拒绝%s时出现未分类错误：%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier); %#ok<AGROW>
            return;
        end
        error('AxisFVTest:RiemannMissingRefusal','非法储库%s未被拒绝。',label);
    end
end

function audit=conservationAudit(rate,info,mesh)
% 每分量保留真实单位：kg/(s rad)、N/rad、N/rad、W/rad，不用常数1跨单位归一化。
% 对空间恒等式sum(V*rate)+sum(F_boundary)-sum(S_cell)=0直接重算。
weightedRate=rate.*mesh.volumePerRadian;boundaryFlux=info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:);
cellSource=info.cellSourceIntegral;
rateIntegral=sum(weightedRate,1);boundaryIntegral=sum(boundaryFlux,1);sourceIntegral=sum(cellSource,1);
defect=rateIntegral+boundaryIntegral-sourceIntegral;
rateAbsolute=sum(abs(weightedRate),1);boundaryAbsolute=sum(abs(boundaryFlux),1);sourceAbsolute=sum(abs(cellSource),1);
physicalScale=rateAbsolute+boundaryAbsolute+sourceAbsolute;
zeroScale=physicalScale==0;relative=zeros(1,4);positive=~zeroScale;
relative(positive)=abs(defect(positive))./physicalScale(positive);
relative(zeroScale&defect~=0)=Inf;

% 显式保守舍入界包含所有共享面按incidence次数贡献，内部面计两次。
% k覆盖单元散度/源加减、除V后乘回V、全域/边界/源求和和末次加减的运算路径。
% 此处eps采用2u的机器间距，比标准unit roundoff u更保守；不改变1e-11相对合同。
multiplicity=full(sum(abs(mesh.incidence),1))';
faceAbsolute=sum(abs(info.faceFluxIntegralPerRadian).*multiplicity,1);
maximumCellFaces=max(full(sum(spones(mesh.incidence),2)));
k=2*mesh.cellCount+numel(mesh.boundaryFaces)+maximumCellFaces+12;
machineSpacing=eps('double');assert(k*machineSpacing<1,'AxisFVTest:RiemannRoundoff','舍入界k过大。');
gammaK=k*machineSpacing/(1-k*machineSpacing);
allAbsolute=rateAbsolute+boundaryAbsolute+sourceAbsolute+faceAbsolute;
roundoffBound=gammaK*allAbsolute;
roundoffPassed=all(abs(defect)<=roundoffBound);
reportedDifference=info.conservationDefect-defect;
reportedDefectAgrees=all(abs(reportedDifference)<=roundoffBound);
audit=struct('quantityOrder',{{'mass','axialMomentum','radialMomentum','totalEnergy'}}, ...
    'units',{{'kg/(s rad)','N/rad','N/rad','W/rad'}},'rawDefect',defect, ...
    'reportedRawDefect',info.conservationDefect,'reportedMinusRecomputed',reportedDifference, ...
    'rateVolumeIntegral',rateIntegral,'boundaryFluxIntegral',boundaryIntegral,'sourceIntegral',sourceIntegral, ...
    'absoluteRateVolumeSum',rateAbsolute,'absoluteBoundaryFluxSum',boundaryAbsolute, ...
    'absoluteCellSourceSum',sourceAbsolute,'physicalSummationScale',physicalScale,'zeroScaleComponents',zeroScale, ...
    'relativeDefect',relative,'maximumRelativeDefect',max(relative),'relativeTolerance',1e-11, ...
    'normalizedPassed',all(relative<1e-11)&&all(defect(zeroScale)==0), ...
    'absoluteAllFaceIncidenceSum',faceAbsolute,'allAbsoluteRoundoffContributions',allAbsolute, ...
    'maximumCellFaceCount',maximumCellFaces,'roundoffOperationBoundK',k,'machineSpacing',machineSpacing, ...
    'gammaK',gammaK,'absoluteRoundoffBound',roundoffBound,'roundoffPassed',roundoffPassed, ...
    'reportedDefectAgrees',reportedDefectAgrees);
end

function F=physical(U,P,n)
% 独立物理通量只用于均匀/全特征出流核对；不代替数值储库通量。
un=sum(P(:,2:3).*n,2);F=[P(:,1).*un,U(:,2:3).*un+P(:,5).*n,(U(:,4)+P(:,5)).*un];
end

function [F,waves,middle]=referenceHLLE(UL,UR,L,R,n)
% 独立标量HLLE：两侧声学特征夹逼，不计算/裁剪接触星压，不按星压切换格式。
% 实际跨零波扇才取用HLL中间态；直接核对该态rho和rhoE-rhou^2/(2rho)。
uL=dot(L(2:3),n);uR=dot(R(2:3),n);SL=min(uL-L(6),uR-R(6));SR=max(uL+L(6),uR+R(6));
waves=[SL SR];FL=physical(UL,L,n);FR=physical(UR,R,n);
middle=struct('used',false,'state',nan(1,4),'density',NaN,'internalEnergyPerVolume',NaN);
assert(all(isfinite(waves))&&SR>SL,'AxisFVTest:RiemannReference','独立HLLE波序无效。');
if SL>=0
    F=FL;
elseif SR<=0
    F=FR;
else
    intermediate=(SR*UR-SL*UL-(FR-FL))/(SR-SL);
    density=intermediate(1);energy=intermediate(4)-sum(intermediate(2:3).^2)/(2*density);
    assert(all(isfinite(intermediate))&&density>0&&isfinite(energy)&&energy>0, ...
        'AxisFVTest:RiemannReference','独立HLLE实际中间态密度或内能非正，不得裁值通过。');
    F=(SR*FL-SL*FR+SL*SR*(UR-UL))/(SR-SL);
    % HLL星区左右两条Rankine-Hugoniot关系须给同一实际数值通量。
    reconstructedLeft=FL+SL*(intermediate-UL);reconstructedRight=FR+SR*(intermediate-UR);
    assert(max(abs(F-reconstructedLeft)./max(1,abs(F)))<2e-12&& ...
        max(abs(F-reconstructedRight)./max(1,abs(F)))<2e-12, ...
        'AxisFVTest:RiemannReference','独立HLL中间态与左右守恒跃迁通量不一致。');
    middle=struct('used',true,'state',intermediate,'density',density,'internalEnergyPerVolume',energy);
end
assert(all(isfinite(F)),'AxisFVTest:RiemannReference','独立HLLE参考产生非有限通量。');
end

function bindings=sourceBindings(root)
% 记录空间算子、状态合同及本测试的源码，防止把历史通过套到后来修改。
paths={'regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVState.m','regionalFlow/axisymmetric/axisymmetricFVConservative.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVGas.m','internalFlow/resolveNSCharacteristicBoundary3D.m', ...
    'tests/testAxisymmetricRiemannBoundary.m'};
bindings=struct([]);
for j=1:numel(paths)
    item=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
end
end

function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0,'AxisFVTest:RiemannSource','无法读取来源文件。');
cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
