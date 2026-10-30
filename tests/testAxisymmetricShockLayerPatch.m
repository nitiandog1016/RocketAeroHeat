function report=testAxisymmetricShockLayerPatch(folder,sourceFile,level,coreFile)
%TESTAXISYMMETRICSHOCKLAYERPATCH 正则耦合层的驻点、加密、物性及范围验证。
% 模块一：同一真实自由流，两种物性和7/9/11/13个径向节点；不读取旧壁面场。
% 模块二：独立正常激波驻止极限、统一能量/熵和连续式切向条件。
% 模块三：积分加密、外侧截断、刚体/尺度变换与变曲率曲面，保留候选身份。
% 模块四：smoke只形成两项诊断，不得把其状态当完整测试已通过。
if nargin<3,level='complete';end
assert(any(strcmp(level,{'complete','smoke'})));
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupEngineeringShockLayer;
assert(~exist(folder,'dir'),'ShockPatch:Overwrite','必须写入新验证目录。');mkdir(folder);
source=load(sourceFile,'fs','geometry','seed');fs=source.fs;geometry=source.geometry;
upstream=struct('pressure',double(fs.pressure),'temperature',double(fs.temperature),'density',double(fs.density), ...
    'velocityVector',double(fs.flowVelocityVector(:)),'method','冻结原19752面M5/30km/alpha20均匀上游；未读取旧壁面场');
ex=upstream.velocityVector/norm(upstream.velocityVector);ey=[0;1;0];ey=ey-ex*(ex'*ey);ey=ey/norm(ey);
k=14.535;shape=struct('origin',source.seed.noseCenter-(source.seed.referenceRadiusM+.008)*ex, ...
    'frame',[ex ey cross(ex,ey)],'K',k*eye(2),'profileB',-k);
modes={'caloricallyPerfect','thermallyPerfect'};nodes=[7 9 11 13 15 17 19];eta=linspace(1,0,401)';
if strcmp(level,'smoke'),nodes=7;eta=linspace(1,0,21)';end
cases=cell(2,numel(nodes));limits=cell(2,1);thermalChecks=cell(2,1);jumpChecks=cell(2,1);
coreReuse=struct('used',false,'file','','note','本次实际重新计算主矩阵');
if nargin<4||isempty(coreFile)
for mode=1:2
    gas=createSurfaceFlowGas(modes{mode},double(fs.R),1.4);
    limits{mode}=axisLimit(upstream,gas,shape);
    % 温度根包含原声明端点和1000K分段点，检验已收敛根不会再次被移动。
    temperatures=[80;300;999.999;1000;1000.001;1400;2000];pressures=linspace(1000,80000,numel(temperatures))';
    enthalpy=zeros(size(temperatures));entropy=enthalpy;
    for j=1:numel(temperatures)
        s=surfaceFlowThermodynamicState(temperatures(j),pressures(j),gas);enthalpy(j)=s.enthalpy;entropy(j)=s.entropy;
    end
    recovered=shockLayerStateFromHS(enthalpy,entropy,gas);
    assert(max(abs(recovered.temperature-temperatures))<2e-9);
    assert(max(abs(recovered.pressure./pressures-1))<2e-11);
    thermalChecks{mode}=struct('gas',gas,'temperatures',temperatures,'pressures',pressures,'state',recovered);
    a1=sqrt(limits{mode}.upstreamGamma*gas.R*upstream.temperature);speeds=a1*[1.02;1.5;3;5];
    actual=shockLayerNormalJumpBatch(upstream,gas,speeds);reference=cell(4,1);
    plane=struct('position',zeros(3,1),'firstDerivatives',[zeros(1,2);eye(2)],'secondDerivatives',zeros(3,3));
    for j=1:4
        local=upstream;local.velocityVector=[speeds(j);0;0];reference{j}=shockSurfaceJumpState(local,gas,plane,1);
        assert(abs(reference{j}.pressure/actual.pressure(j)-1)<2e-11);
        assert(abs(reference{j}.temperature/actual.temperature(j)-1)<2e-11);
        assert(abs(reference{j}.entropy-actual.entropy(j))<1e-8);
    end
    jumpChecks{mode}=struct('gas',gas,'speeds',speeds,'actual',actual,'reference',{reference});
    for j=1:numel(nodes)
        patch=solveAxisymmetricShockLayerPatch(upstream,gas,shape,.032,eta,struct('radialNodeCount',nodes(j)));
        assert(max(abs(patch.energyResidual),[],'all')<2e-12);
        assert(max(abs(patch.entropyResidual),[],'all')<1e-8);
        assert(max(abs(patch.massConstraintMS),[],'all')<2e-7);
        assert(abs(patch.pressurePa(end,1)/limits{mode}.pressurePa-1)<2e-9);
        assert(abs(patch.temperatureK(end,1)-limits{mode}.temperatureK)<2e-7);
        assert(abs(patch.normalSpeedMS(end,1))<1e-7&&all(patch.streamwiseSpeedMS(:,1)==0));
        assert(~patch.fullEulerSolved&&~patch.bodyMatched&&~patch.productionQualified&&patch.newProductionSolvedAreaM2==0);
        cases{mode,j}=patch;
        fprintf('COUPLED_PATCH mode=%s N=%d axisP=%.12g axisV=%.6g n0=%.12g mass=%.6g calls=%d\n', ...
            modes{mode},nodes(j),patch.pressurePa(end,1),patch.normalSpeedMS(end,1),patch.normalDistanceM(end,1), ...
            max(abs(patch.massConstraintMS),[],'all'),patch.integrationRHSCalls);
    end
end
else
    assert(strcmp(level,'complete'),'ShockPatch:CoreReuse','只允许完整范围复核复用已完成八组主矩阵。');
    prior=load(coreFile);assert(isequal(size(prior.cases),[2 numel(nodes)])&& ...
        isequaln(prior.upstream,upstream)&&isequaln(prior.shape,shape),'ShockPatch:CoreReuse','原主矩阵或输入身份不符。');
    cases=prior.cases;limits=prior.limits;thermalChecks=prior.thermalChecks;jumpChecks=prior.jumpChecks;
    for mode=1:2,for j=1:numel(nodes)
        assert(cases{mode,j}.radialNodeCount==nodes(j)&&isequal(cases{mode,j}.eta,eta)&& ...
            ~cases{mode,j}.productionQualified&&max(abs(cases{mode,j}.massConstraintMS),[],'all')<2e-7);
    end,end
    coreReuse=struct('used',true,'file',coreFile,'note','从首轮失败前检查点逐值复用八组已通过主矩阵；未重新计算');
end
if strcmp(level,'smoke')
    report=struct('status','draft-two-mode-smoke-only','completeValidation',false,'productionQualified',false,'newProductionSolvedAreaM2',0);
    save(fullfile(folder,'axisymmetric_patch_smoke.mat'),'cases','limits','thermalChecks','jumpChecks','shape','upstream','report','-v7');
    disp(report);return;
end
% 在额外范围试验前保存已完成的主矩阵，任何后续失败仍保留原计算证据。
if ~coreReuse.used
    save(fullfile(folder,'core_before_scope_checks.mat'),'cases','limits','thermalChecks','jumpChecks','shape','upstream','-v7');
end

%% 二、积分和径向加密、有限外侧截断；数值敏感性不等于真实体形已匹配
queries=[0 .005 .015 .02]'.^2;convergence=zeros(2,numel(nodes)-1);outerChecks=cell(2,2);integrationRefined=cell(2,1);
parabolic=cell(2,1);outerDifferences=nan(2,2);refinementDifferences=zeros(2,1);outerRefusals=cell(2,2);
for mode=1:2
    base=cases{mode,end};gas=base.gas;
    for j=1:numel(nodes)-1
        a=sample(cases{mode,j},queries,'normalDistanceM');b=sample(cases{mode,j+1},queries,'normalDistanceM');
        convergence(mode,j)=max(abs(a-b),[],'all');
    end
    assert(all(convergence(mode,2:3)<convergence(mode,1:2)),'ShockPatch:RadialRefinement','原7至13节点加密没有改善共同位置厚度。');
    assert(convergence(mode,end)<1e-7,'ShockPatch:RadialRefinement','本测试域厚度还未分辨至0.1微米。');
    integrationRefined{mode}=solveAxisymmetricShockLayerPatch(upstream,gas,shape,.032,eta, ...
        struct('radialNodeCount',nodes(end),'relativeTolerance',5e-11,'maximumEtaStep',.0015));
    refinementDifferences(mode)=max(abs(base.normalDistanceM-integrationRefined{mode}.normalDistanceM),[],'all');
    assert(refinementDifferences(mode)<2e-9,'ShockPatch:IntegrationRefinement','积分加密厚度变化超过2纳米。');
    for j=1:2
        radii=[.028 .036];checkpoint=fullfile(folder,sprintf('outer_mode%d_radius%d_rejection.mat',mode,j));
        try
            outerChecks{mode,j}=solveAxisymmetricShockLayerPatch(upstream,gas,shape,radii(j),eta, ...
                struct('radialNodeCount',nodes(end),'failureCheckpointFile',checkpoint));
            outerDifferences(mode,j)=max(abs(sample(outerChecks{mode,j},queries,'normalDistanceM')-sample(base,queries,'normalDistanceM')),[],'all');
            fprintf('OUTER_SCOPE mode=%d radius=%.3f solved=1\n',mode,radii(j));
        catch ex
            if ~strcmp(ex.identifier,'ShockPatch:Sonic'),rethrow(ex);end
            failed=load(checkpoint,'failure');
            outerRefusals{mode,j}=struct('identifier',ex.identifier,'checkpointFile',checkpoint, ...
                'trialEta',failed.failure.eta,'minimumTrialDeterminant',min(failed.failure.determinant), ...
                'acceptedFlowSonicLocationProven',false,'resultAvailable',false);
            fprintf('OUTER_SCOPE mode=%d radius=%.3f solved=0 eta=%.9g det=%.9g\n',mode,radii(j), ...
                failed.failure.eta,min(failed.failure.determinant));
        end
    end
    other=shape;other.profileB=0;
    parabolic{mode}=solveAxisymmetricShockLayerPatch(upstream,gas,other,.028,eta,struct('radialNodeCount',11));
    assert(max(parabolic{mode}.shockGeometry.meridionalCurvaturePerM)-min(parabolic{mode}.shockGeometry.meridionalCurvaturePerM)>1);
end
save(fullfile(folder,'scope_before_transform_checks.mat'),'integrationRefined','outerChecks','parabolic', ...
    'queries','convergence','outerDifferences','refinementDifferences','outerRefusals','-v7');

%% 三、刚体和长度缩放；所有标量及SI速度应保持不变
base=cases{2,end};angle=.37;rotation=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];shift=[.7;-.2;.4];
moved=shape;moved.origin=rotation*shape.origin+shift;moved.frame=rotation*shape.frame;
up=upstream;up.velocityVector=rotation*up.velocityVector;
transformed=solveAxisymmetricShockLayerPatch(up,base.gas,moved,.032,eta,struct('radialNodeCount',nodes(end)));
assert(max(abs(transformed.normalDistanceM-base.normalDistanceM),[],'all')<2e-9);
for j=1:3
    expected=shift(j)+sum(base.position.*reshape(rotation(j,:),1,1,3),3);
    assert(max(abs(transformed.position(:,:,j)-expected),[],'all')<2e-9);
end
scale=2.3;scaledShape=shape;scaledShape.origin=shape.origin*scale;scaledShape.K=shape.K/scale;scaledShape.profileB=shape.profileB/scale;
scaled=solveAxisymmetricShockLayerPatch(upstream,base.gas,scaledShape,.032*scale,eta,struct('radialNodeCount',nodes(end)));
save(fullfile(folder,'transforms_before_assertions.mat'),'transformed','rotation','shift','scale','scaled','-v7');
scaleDistanceDifference=max(abs(scaled.normalDistanceM/scale-base.normalDistanceM),[],'all');
scalePressureDifference=max(abs(scaled.pressurePa./base.pressurePa-1),[],'all');
scaleCheck=struct('scale',scale,'distanceDifferenceM',scaleDistanceDifference, ...
    'pressureRelativeDifference',scalePressureDifference,'distanceGateM',2e-9,'pressureRelativeGate',2e-9, ...
    'passed',scaleDistanceDifference<2e-9&&scalePressureDifference<2e-9, ...
    'note','保留原缩放门槛，失败作为未收口数值敏感性记录，不将后续诊断运行完成视为全范围通过');
fprintf('SCALE_SCOPE passed=%d distanceM=%.9g pressureRelative=%.9g\n',scaleCheck.passed,scaleDistanceDifference,scalePressureDifference);

%% 四、真实原面元几何只作距离与速度接口诊断，不赋予壁面压力/热流资格
angles=[0 2*pi/3 4*pi/3];radii=[.012 .022 .032];coordinates=zeros(9,2);q=0;
shockPositions=zeros(9,3);normals=shockPositions;predicted=zeros(9,1);normalRatio=nan(9,1);terminalVelocity=nan(9,3);
for radius=radii
    for angle=angles
        q=q+1;coordinates(q,:)=[radius*cos(angle),radius*sin(angle)];
        g=axisymmetricShockGraph(shape,radius^2);er=shape.frame(:,2)*cos(angle)+shape.frame(:,3)*sin(angle);
        shockPositions(q,:)=(shape.origin+shape.frame(:,1)*g.axialPositionM+radius*er)';
        normal=(shape.frame(:,1)*g.inverseMeridionalMetric-radius*g.azimuthalCurvaturePerM*er);normals(q,:)=normal';
        predicted(q)=sample(base,radius^2,'normalDistanceM');
        u=sample(base,radius^2,'streamwiseSpeedMS');v=sample(base,radius^2,'normalSpeedMS');
        tangent=radius*g.azimuthalCurvaturePerM*shape.frame(:,1)+g.inverseMeridionalMetric*er;
        terminalVelocity(q,:)=(u*tangent+v*normal)';
    end
end
mapping=intersectShockNormalsWithBody(double(geometry.V),double(geometry.F),shockPositions,normals);
residual=nan(9,1);residual(mapping.hit)=predicted(mapping.hit)-mapping.distanceM(mapping.hit);
for j=find(mapping.hit)'
    normalRatio(j)=abs(dot(terminalVelocity(j,:),mapping.bodyNormal(j,:)))/norm(terminalVelocity(j,:));
end
realBody=struct('coordinates',coordinates,'mapping',mapping,'predictedNormalDistanceM',predicted, ...
    'residualM',residual,'candidateVelocityMS',terminalVelocity,'velocityNormalFraction',normalRatio, ...
    'bodyMatchToleranceM',1e-5,'fitPerformed',false,'candidateAxisymmetryIsAdditionalAssumption',true, ...
    'flowValuesAssignedToProductionBody',false,'productionQualified',false);

%% 五、输入拒绝和资格身份，不能把拒绝改写为默认自由流或补零
refused={};bad=shape;bad.K(2,2)=1.1*k;
refused{end+1}=reject(@()solveAxisymmetricShockLayerPatch(upstream,base.gas,bad,.032,eta),'ShockPatch:Axisymmetry');
bad=shape;bad.frame(:,[1 3])=-bad.frame(:,[1 3]);
refused{end+1}=reject(@()solveAxisymmetricShockLayerPatch(upstream,base.gas,bad,.032,eta),'ShockPatch:Frame');
refused{end+1}=reject(@()solveAxisymmetricShockLayerPatch(upstream,base.gas,shape,1.1/k,eta),'ShockPatch:Graph');
refused{end+1}=reject(@()solveAxisymmetricShockLayerPatch(upstream,base.gas,shape,.99/k,eta),'ShockPatch:NormalMach');
refused{end+1}=reject(@()solveAxisymmetricShockLayerPatch(upstream,base.gas,shape,.032,[1;.4;.5;0]),'ShockPatch:MassCoordinate');
badGas=base.gas;badGas.temperatureRangeK=[80 1000];
refused{end+1}=reject(@()solveAxisymmetricShockLayerPatch(upstream,badGas,shape,.032,eta),'ShockPatch:Temperature');
report=struct('status','tested-coupled-axisymmetric-core-with-explicit-outer-domain-rejection', ...
    'completeValidation',~any(cellfun(@(x)~isempty(x),outerRefusals),'all')&&scaleCheck.passed,'primaryMatrixValidated',true, ...
    'allPlannedScopeProbesExecuted',true,'coreReuse',coreReuse,'outerDomainRefusals',{outerRefusals}, ...
    'thermalModes',2,'radialNodeCounts',nodes,'scaleCheck',scaleCheck,'radialRefinementThicknessDifferencesM',convergence, ...
    'integrationRefinementThicknessDifferencesM',refinementDifferences,'outerCutoffThicknessDifferencesM',outerDifferences, ...
    'outerCutoffScope','仅在r<=0.02m的同位置比较0.028/0.032/0.036m外侧截断；不是严格误差界或外域独立性证明', ...
    'originalFaceCount',size(geometry.F,1),'realBodyHitCount',nnz(mapping.hit),'realBodyMaximumResidualM',max(abs(residual)), ...
    'realBodyMaximumNormalVelocityFraction',max(normalRatio),'realBodyFitPerformed',false,'rejections',{refused}, ...
    'experimentalValidation',false,'fullEulerSolved',false,'productionQualified',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'axisymmetric_patch_evidence.mat'),'cases','limits','thermalChecks','jumpChecks','shape','upstream','nodes', ...
    'integrationRefined','outerChecks','parabolic','transformed','rotation','shift','scale','scaled','queries', ...
    'geometry','fs','realBody','coreReuse','outerRefusals','scaleCheck','report','-v7');
fid=fopen(fullfile(folder,'axisymmetric_patch_report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));disp(report);
end

function limit=axisLimit(upstream,gas,shape)
% 模块用途：独立调用正常激波入口后再以标量二分求同熵驻止态，不调用耦合层端点。
sample=ellipticConicShockSample(shape,[0;0]);jump=shockSurfaceJumpState(upstream,gas,sample,1);
first=surfaceFlowThermodynamicState(upstream.temperature,upstream.pressure,gas);
H=first.enthalpy+norm(upstream.velocityVector)^2/2;lo=gas.temperatureRangeK(1);hi=gas.temperatureRangeK(2);
for j=1:55
    T=(lo+hi)/2;s=surfaceFlowThermodynamicState(T,jump.pressure,gas);
    if s.enthalpy<H,lo=T;else,hi=T;end
end
T=(lo+hi)/2;s=surfaceFlowThermodynamicState(T,jump.pressure,gas);p=jump.pressure*exp((s.entropy-jump.entropy)/gas.R);
limit=struct('temperatureK',T,'pressurePa',p,'shock',jump,'upstreamGamma',first.gammaLocal, ...
    'method','原正常激波RH入口 + 独立标量二分同H0/熵驻止极限；不是实验基准');
end

function value=sample(patch,queries,name)
% 模块用途：同一z=r^2中的Chebyshev重心插值，仅比较已算域内的端点。
z=patch.radiusSquaredM2;values=patch.(name);N=numel(z);weights=(-1).^(0:N-1)';weights([1 end])=.5*weights([1 end]);
value=zeros(size(queries));
for j=1:numel(queries)
    [distance,index]=min(abs(z-queries(j)));
    if distance<16*eps*max(z),value(j)=values(end,index);else
        w=weights./(queries(j)-z);value(j)=(values(end,:)*w)/sum(w);
    end
end
end
function identifier=reject(action,expected)
identifier='';try,action();catch ex,identifier=ex.identifier;end
assert(strcmp(identifier,expected),'ShockPatch:TestRejection','实际拒绝%s，期望%s。',identifier,expected);
end
