function report=testAxisymmetricWallThermalGradient(folder,legacyEvidenceFile)
%TESTAXISYMMETRICWALLTHERMALGRADIENT 斜壁Dirichlet热梯度的解析、守恒和真实状态校核。
% 模块一：在强非正交的真斜墙上，用独立体积积分的线性守恒场检验冷热墙符号。
% 模块二：函数墙温显式提供沿t=[-nr,nx]的导数，核对独立解析梯度及法向Fourier通量。
% 模块三：非线性同一Q对照只允许墙热通量变化；检查四/五量账和完整Jacobian方向。
% 模块四：冻结v65稳态默认分支逐值重放，并对真实负热流两面只作新旧空间诊断。
% 本测试没有推进新稳态；修正后的单次空间结果不是实验精度或生产面积验收。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v66_validation','wall_thermal_gradient');end
if nargin<2
    legacyEvidenceFile=fullfile(root,'outputs','v65_validation','tmr_sa_continue1','tmr_sa_evidence.mat');
    if ~exist(legacyEvidenceFile,'file'),legacyEvidenceFile=fullfile(root,'validation','v65','tmr_sa_continue1','tmr_sa_evidence.mat');end
end
assert(exist(legacyEvidenceFile,'file')==2,'AxisFVTest:WallThermalLegacy','缺少冻结v65真实稳态，不能跳过旧分支核对。');
if ~exist(folder,'dir'),mkdir(folder);end
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[100 1500], ...
    'transport',struct('model','sutherland','referenceViscosity',.05,'referenceTemperature',300,'sutherlandTemperature',110.4)));
[s,eta]=ndgrid([0 .02 .05 .08],[0 .006 .017 .03]);a=1.8;m=.3;ri=.03;thickness=.03;
mesh=buildAxisymmetricFVGeometry(s+a*eta,ri+m*s+eta);f=mesh.faces;
wallFaces=find(f.boundarySide==3|f.boundarySide==4);inner=f.boundarySide(wallFaces)==3;
normal=f.normal(wallFaces,:);tangent=[-normal(:,2) normal(:,1)];
skew=abs(sum(f.delta(wallFaces,:).*tangent,2)./sum(f.delta(wallFaces,:).*normal,2));
assert(min(skew)>3,'AxisFVTest:WallThermalExcitation','斜网格没有真正激活强非正交位移。');
refusals=struct('label',{},'identifier',{});

%% 一、常温斜墙的线性法向场；rhoE线性时，体质心点值是真正保守FV均值
Tin=300;Tout=500;gradient=[-m 1]*(Tout-Tin)/(thickness*(1-a*m));
offset=Tin-(Tout-Tin)*ri/(thickness*(1-a*m));rho=1.2;
temperature=@(p)offset+p(:,1)*gradient(1)+p(:,2)*gradient(2);
[Q,quadrature]=linearFVMeans(mesh,gas,rho,offset,gradient);
trace=@(p)[rho 0 0 temperature(p)];edge=struct('type','prescribedTrace','primitive',trace,'nuTilde',.1);
boundary=struct('iMin',edge,'iMax',edge, ...
    'jMin',struct('type','wall','temperature',Tin,'velocity',[0 0],'thermalGradient','dirichlet-tangent'), ...
    'jMax',struct('type','wall','temperature',Tout,'velocity',[0 0],'thermalGradient','dirichlet-tangent'));
[rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);w=info.wall;
expectedHeat=-independentK(w.temperatureInput,gas).*sum(normal.*gradient,2);
normalHeatRelative=max(abs(w.heatFluxIntoWall-expectedHeat)./abs(expectedHeat));
assert(normalHeatRelative<2e-10&&all(w.heatFluxIntoWall(inner)>0)&&all(w.heatFluxIntoWall(~inner)<0), ...
    'AxisFVTest:WallThermalSign','冷墙加热或热墙冷却的符号/独立Fourier量值不正确。');
assert(isequal(w.faceID,wallFaces)&&isfield(info,'wallThermalGradient'));
constantAudit=info.wallThermalGradient;
gradientError=max(abs(constantAudit.gradientXR-gradient),[],'all')/norm(gradient);
assert(gradientError<2e-10&&all(constantAudit.tangentialGradientInputKm==0)&&~constantAudit.cellStateClipped);
constantConstraints=constraints(constantAudit,mesh);
constantLedger=ledger(rate,info,mesh);

%% 二、非常量墙温：边界导数是显式输入，沿墙方向与几何t同号
linearGradient=[-500 4000];linearOffset=250;
linearT=@(p)linearOffset+p(:,1)*linearGradient(1)+p(:,2)*linearGradient(2);
[linearQ,linearQuadrature]=linearFVMeans(mesh,gas,rho,linearOffset,linearGradient);
linearEdge=struct('type','prescribedTrace','primitive',@(p)[rho 0 0 linearT(p)],'nuTilde',.1);
innerTangent=[1 m]/sqrt(1+m^2);outerTangent=-innerTangent;
linearBoundary=struct('iMin',linearEdge,'iMax',linearEdge, ...
    'jMin',struct('type','wall','temperature',linearT,'velocity',[0 0], ...
        'thermalGradient','dirichlet-tangent','temperatureTangentialGradient',dot(linearGradient,innerTangent)), ...
    'jMax',struct('type','wall','temperature',linearT,'velocity',[0 0], ...
        'thermalGradient','dirichlet-tangent','temperatureTangentialGradient',@(p)dot(linearGradient,outerTangent)));
[linearRate,linearInfo]=axisymmetricNavierStokesResidual(linearQ,mesh,gas,linearBoundary);
linearAudit=linearInfo.wallThermalGradient;linearWall=linearInfo.wall;
linearExpected=-independentK(linearT(linearWall.centroid),gas).*sum(normal.*linearGradient,2);
functionHeatRelative=max(abs(linearWall.heatFluxIntoWall-linearExpected)./abs(linearExpected));
functionGradientRelative=max(abs(linearAudit.gradientXR-linearGradient),[],'all')/norm(linearGradient);
assert(functionHeatRelative<2e-10&&functionGradientRelative<2e-10&& ...
    all(linearWall.heatFluxIntoWall(inner)>0)&&all(linearWall.heatFluxIntoWall(~inner)<0), ...
    'AxisFVTest:WallThermalFunction','显式非常量Tw切向导数或解析法向通量不匹配。');
functionConstraints=constraints(linearAudit,mesh);functionLedger=ledger(linearRate,linearInfo,mesh);
% 直接helper合同：非墙保持输入原梯度逐值不变；墙解满足两个独立约束。
[oldX,oldR,Twall]=legacyTemperatureGradient(linearQ,mesh,gas,linearBoundary,wallFaces);
Plinear=axisymmetricFVState(linearQ,gas);
[directX,directR,directAudit]=axisymmetricWallThermalGradient(mesh,Plinear(:,4),linearBoundary,wallFaces,Twall,oldX,oldR);
nonwall=true(numel(f.owner),1);nonwall(wallFaces)=false;
assert(isequal(directX(nonwall),oldX(nonwall))&&isequal(directR(nonwall),oldR(nonwall))&& ...
    isequal(directAudit.gradientXR,linearAudit.gradientXR)&& ...
    isequal(directAudit.normalGradientKm,linearAudit.normalGradientKm)&& ...
    max(abs([directX(wallFaces) directR(wallFaces)]-linearGradient),[],'all')/norm(linearGradient)<2e-10, ...
    'AxisFVTest:WallThermalHelper','helper改动了非墙面，或与完整NS使用的墙梯度不同。');
% 独立旧WLS模板不复刻NS边界EOS往返的ulp；新约束解不依赖旧梯度，单独逐值核对。

%% 三、非线性状态只比较同一空间算子；不以新墙热正值冒充新稳态
cellPoint=mesh.volumeCentroid;logicalEta=(cellPoint(:,2)-ri-m*cellPoint(:,1))/(1-a*m);
logicalS=(cellPoint(:,1)-a*(cellPoint(:,2)-ri))/(1-a*m);
Tnonlinear=Tin+(Tout-Tin)*logicalEta/thickness+ ...
    120*sin(pi*logicalEta/thickness).*sin(pi*logicalS/.08);
Pnonlinear=[rho*ones(mesh.cellCount,1),zeros(mesh.cellCount,2),Tnonlinear];
Qnonlinear=axisymmetricFVConservative(Pnonlinear,gas);Qunchanged=Qnonlinear;
nonlinearT=@(p)nonlinearTemperature(p,ri,m,a,thickness,Tin,Tout);
nonlinearEdge=struct('type','prescribedTrace','primitive',@(p)[rho 0 0 nonlinearT(p)],'nuTilde',.1);
nonlinearBoundary=boundary;nonlinearBoundary.iMin=nonlinearEdge;nonlinearBoundary.iMax=nonlinearEdge;
legacyBoundary=removePolicy(nonlinearBoundary);
[legacyRate,legacyInfo]=axisymmetricNavierStokesResidual(Qnonlinear,mesh,gas,legacyBoundary);
[newRate,newInfo]=axisymmetricNavierStokesResidual(Qnonlinear,mesh,gas,nonlinearBoundary);
assert(isequal(Qnonlinear,Qunchanged)&&~isfield(legacyInfo,'wallThermalGradient'));
nonlinearAudit=newInfo.wallThermalGradient;nonlinearConstraints=constraints(nonlinearAudit,mesh);
oldTangential=max(abs(nonlinearAudit.oldTangentialGradientKm));
newTangential=max(abs(sum(nonlinearAudit.gradientXR.*nonlinearAudit.tangentXR,2)));
assert(oldTangential>100&&newTangential<1e-9*oldTangential, ...
    'AxisFVTest:WallThermalNegativeControl','非线性旧重构没有明显违反常Tw切向约束，或新梯度仍违反。');
fourQuantityChange=assertOnlyWallHeatChanged(legacyRate,legacyInfo,newRate,newInfo,mesh,wallFaces);
legacyLedger=ledger(legacyRate,legacyInfo,mesh);newLedger=ledger(newRate,newInfo,mesh);
N=.1+.01*sin((1:mesh.cellCount)');Q5=[Qnonlinear rho*N];Q5unchanged=Q5;
[legacyRate5,legacyInfo5]=axisymmetricRANSResidual(Q5,mesh,gas,legacyBoundary);
[newRate5,newInfo5]=axisymmetricRANSResidual(Q5,mesh,gas,nonlinearBoundary);
fiveQuantityChange=assertOnlyWallHeatChanged(legacyRate5,legacyInfo5,newRate5,newInfo5,mesh,wallFaces);
assert(isequal(Q5,Q5unchanged)&&isequaln(legacyInfo5.sa,newInfo5.sa)&& ...
    isequal(legacyRate5(:,5),newRate5(:,5))&&all(newInfo5.sa.faceMuT(wallFaces)==0), ...
    'AxisFVTest:WallThermalSA','改变热梯度同时篡改了第五量输运或使用墙湍流加热因子。');
legacyLedger5=ledger(legacyRate5,legacyInfo5,mesh);newLedger5=ledger(newRate5,newInfo5,mesh);
jacobian4=directionJacobian(Qnonlinear,mesh,gas,nonlinearBoundary,[rho 20 20 2000]);
jacobian5=directionJacobian(Q5,mesh,gas,nonlinearBoundary,[rho 20 20 2000 .2]);

%% 四、非法或不完整的墙温导数不能暗中回退，也不能把几何错误当有效导数
bad=boundary;bad.jMin.thermalGradient='normal-only-unreviewed';
reject('unknown-policy',@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad),'AxisFV:WallThermalBoundary');
bad=linearBoundary;bad.jMin=rmfield(bad.jMin,'temperatureTangentialGradient');
reject('function-Tw-missing-tangent',@()axisymmetricNavierStokesResidual(linearQ,mesh,gas,bad),'AxisFV:WallThermalBoundary');
bad=boundary;bad.jMin.temperatureTangentialGradient=1;
reject('constant-Tw-nonzero-tangent',@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad),'AxisFV:WallThermalBoundary');
bad=legacyBoundary;bad.jMin.temperatureTangentialGradient=0;
reject('orphan-tangent',@()axisymmetricNavierStokesResidual(Qnonlinear,mesh,gas,bad),'AxisFV:WallThermalBoundary');
bad=boundary;bad.iMin.thermalGradient='dirichlet-tangent';
reject('nonwall-policy',@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad),'AxisFV:WallThermalBoundary');
for badValue={NaN,Inf,[0 0],complex(1,1),@(p)NaN,@(p)[0 0]}
    bad=linearBoundary;bad.jMin.temperatureTangentialGradient=badValue{1};
    reject('invalid-tangent-value',@()axisymmetricNavierStokesResidual(linearQ,mesh,gas,bad),'AxisFV:WallThermalBoundary');
end
reject('wrong-cell-temperature-size',@()axisymmetricWallThermalGradient(mesh,Plinear(1:end-1,4),linearBoundary,wallFaces,Twall,oldX,oldR),'AxisFV:WallThermalGeometry');
badMesh=mesh;badMesh.faces.delta(wallFaces(1),:)=-badMesh.faces.delta(wallFaces(1),:);
reject('nonpositive-wall-normal-distance',@()axisymmetricWallThermalGradient(badMesh,Plinear(:,4),linearBoundary,wallFaces,Twall,oldX,oldR),'AxisFV:WallThermalGeometry');
badX=oldX;badX(wallFaces(1))=NaN;
reject('nonfinite-input-gradient',@()axisymmetricWallThermalGradient(mesh,Plinear(:,4),linearBoundary,wallFaces,Twall,badX,oldR),'AxisFV:WallThermalState');

%% 五、真实v65负热流位置的只读空间对照，保留冻结父状态和完整旧结果
legacySHA=fileSHA(legacyEvidenceFile);frozen=load(legacyEvidenceFile,'implicit','mesh','gas','boundary');
assert(frozen.implicit.steadyConverged&&frozen.implicit.residualNorm<=2e-6&&size(frozen.implicit.Q,2)==5);
[frozenRate,frozenInfo]=axisymmetricRANSResidual(frozen.implicit.Q,frozen.mesh,frozen.gas,frozen.boundary);
assert(isequaln(frozenInfo,frozen.implicit.finalInfo)&& ...
    max(abs(frozenRate)./frozen.implicit.residualScale,[],'all')==frozen.implicit.residualNorm, ...
    'AxisFVTest:WallThermalLegacyReplay','默认省略策略不能完整逐值重放v65真实终态。');
correctedBoundary=frozen.boundary;correctedBoundary.jMin.thermalGradient='dirichlet-tangent';
frozenQ=frozen.implicit.Q;[correctedRate,correctedInfo]=axisymmetricRANSResidual(frozenQ,frozen.mesh,frozen.gas,correctedBoundary);
assert(isequal(frozenQ,frozen.implicit.Q)&&strcmp(legacySHA,fileSHA(legacyEvidenceFile)));
frozenFaces=frozenInfo.wall.faceID;selected=[61;63];[present,indices]=ismember(selected,frozenFaces);assert(all(present));
newAudit=correctedInfo.wallThermalGradient;[present,ga]=ismember(selected,newAudit.faceID);assert(all(present));
beforeHeat=frozenInfo.wall.heatFluxIntoWall(indices);afterHeat=correctedInfo.wall.heatFluxIntoWall(indices);
assert(all(beforeHeat<0)&&all(afterHeat>0)&& ...
    all(frozenInfo.primitive(frozenInfo.wall.ownerCell(indices),4)>frozenInfo.wall.temperatureInput(indices)), ...
    'AxisFVTest:WallThermalRealExcitation','冻结真实负墙热或同Q导数空间对照不符；不能以改数据替代验证。');
realChange=assertOnlyWallHeatChanged(frozenRate,frozenInfo,correctedRate,correctedInfo,frozen.mesh,frozenFaces);
realConstraints=constraints(newAudit,frozen.mesh);realLedger=ledger(correctedRate,correctedInfo,frozen.mesh);
realDiagnostic=struct('sourceSHA256',legacySHA,'completeLegacyInfoReplay',true,'faceID',selected, ...
    'ownerCell',frozenInfo.wall.ownerCell(indices),'ownerTemperatureK',frozenInfo.primitive(frozenInfo.wall.ownerCell(indices),4), ...
    'prescribedWallTemperatureK',frozenInfo.wall.temperatureInput(indices),'oldSignedHeatWm2',beforeHeat, ...
    'sameStateNewSignedHeatWm2',afterHeat,'oldTangentialGradientKm',newAudit.oldTangentialGradientKm(ga), ...
    'newTangentialGradientKm',sum(newAudit.gradientXR(ga,:).*newAudit.tangentXR(ga,:),2), ...
    'oldSteadyOriginalResidual',frozen.implicit.residualNorm, ...
    'newSameStateOriginalResidual',max(abs(correctedRate)./frozen.implicit.residualScale,[],'all'), ...
    'newSteadySolvePerformed',false,'newStateIsSteady',false,'experimentalValidation',false);

report=struct('status','passed-explicit-sloping-wall-Dirichlet-gradient-spatial-contract', ...
    'cellCount',mesh.cellCount,'minimumWallSkewRatio',min(skew), ...
    'constantWallHeatRelative',normalHeatRelative,'constantGradientRelative',gradientError, ...
    'functionWallHeatRelative',functionHeatRelative,'functionGradientRelative',functionGradientRelative, ...
    'linearConservativeFVQuadrature',quadrature,'functionConservativeFVQuadrature',linearQuadrature, ...
    'constantConstraints',constantConstraints,'functionConstraints',functionConstraints, ...
    'constantLedger',constantLedger,'functionLedger',functionLedger, ...
    'oldNonlinearTangentialMaximumKm',oldTangential,'newNonlinearTangentialMaximumKm',newTangential, ...
    'nonlinearConstraints',nonlinearConstraints,'fourQuantityChange',fourQuantityChange, ...
    'fiveQuantityChange',fiveQuantityChange,'legacyLedger',legacyLedger,'newLedger',newLedger, ...
    'legacyLedger5',legacyLedger5,'newLedger5',newLedger5,'jacobian4',jacobian4,'jacobian5',jacobian5, ...
    'refusals',refusals,'realDiagnostic',realDiagnostic,'realChange',realChange, ...
    'realConstraints',realConstraints,'realLedger',realLedger,'sourceBindings',sourceBindings(root), ...
    'timeStepsTaken',0,'physicalSolution',false,'experimentalValidation',false, ...
    'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'wall_thermal_gradient_evidence.mat'),'report','gas','mesh','Q','boundary','rate','info', ...
    'linearQ','linearBoundary','linearRate','linearInfo','directAudit','oldX','oldR','directX','directR', ...
    'Qnonlinear','legacyBoundary','nonlinearBoundary','legacyRate','legacyInfo','newRate','newInfo', ...
    'Q5','legacyRate5','legacyInfo5','newRate5','newInfo5','frozenRate','frozenInfo','correctedRate','correctedInfo','correctedBoundary');
fid=fopen(fullfile(folder,'wall_thermal_gradient_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Wall thermal gradient passed: true sloping FV means, two signs, function Tw, 4/5 ledger/Jacobian, %d refusals, exact v65 replay and read-only negative-q comparison.\n',numel(refusals));
    function reject(label,action,expected)
        try,action();catch exception
            assert(strcmp(exception.identifier,expected),'AxisFVTest:WallThermalUnexpectedRefusal','%s期望%s，实际%s。',label,expected,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFVTest:WallThermalMissingRefusal','非法合同%s未拒绝。',label);
    end
end

function result=constraints(audit,mesh)
% 模块用途：独立检验物理切向导数和cell-to-wall温差两条约束，不复算被测法向公式。
assert(strcmp(audit.kind,'explicit-Dirichlet-tangent-wall-temperature-gradient-v1')&&~audit.cellStateClipped);
g=audit.gradientXR;f=audit.faceID;d=mesh.faces.delta(f,:);
tangentResidual=sum(g.*audit.tangentXR,2)-audit.tangentialGradientInputKm;
temperatureResidual=sum(g.*d,2)-(audit.wallTemperatureInputK-audit.cellTemperatureK);
tangentBudget=256*eps*(sum(abs(g.*audit.tangentXR),2)+abs(audit.tangentialGradientInputKm));
temperatureBudget=256*eps*(sum(abs(g.*d),2)+abs(audit.wallTemperatureInputK)+abs(audit.cellTemperatureK));
assert(all(abs(tangentResidual)<=tangentBudget)&&all(abs(temperatureResidual)<=temperatureBudget), ...
    'AxisFVTest:WallThermalConstraint','显式墙梯度未同时满足切向Dirichlet导数与真实位移温差。');
result=struct('maximumTangentialResidualKm',max(abs(tangentResidual)), ...
    'maximumDisplacementTemperatureResidualK',max(abs(temperatureResidual)), ...
    'maximumTangentBudgetKm',max(tangentBudget),'maximumTemperatureBudgetK',max(temperatureBudget),'passed',true);
end

function record=assertOnlyWallHeatChanged(oldRate,old,newRate,new,mesh,wallFaces)
% 模块用途：同Q只允许显式墙的能量通量改变，不改状态/剪切/本构系数/SA源或非墙面。
nq=size(oldRate,2);notEnergy=[1:3 5:nq];nonwall=true(numel(mesh.faces.owner),1);nonwall(wallFaces)=false;
assert(isequal(oldRate(:,notEnergy),newRate(:,notEnergy))&& ...
    isequal(old.faceFluxIntegralPerRadian(:,notEnergy),new.faceFluxIntegralPerRadian(:,notEnergy))&& ...
    isequal(old.faceFluxIntegralPerRadian(nonwall,4),new.faceFluxIntegralPerRadian(nonwall,4))&& ...
    isequal(old.cellSourceIntegral,new.cellSourceIntegral)&&isequal(old.primitive,new.primitive)&& ...
    isequal(old.wall.viscousTraction,new.wall.viscousTraction)&&isequaln(old.transport,new.transport)&& ...
    isequal(old.convectiveFrequency,new.convectiveFrequency)&&isequal(old.diffusiveFrequency,new.diffusiveFrequency), ...
    'AxisFVTest:WallThermalScope','显式墙热策略改变了非墙能量、其他守恒量、系数、状态或频率。');
change=new.faceFluxIntegralPerRadian-old.faceFluxIntegralPerRadian;
wallEnergy=(new.wall.heatFluxIntoWall-old.wall.heatFluxIntoWall).*new.wall.areaPerRadian;
heatDifference=max(abs(change(wallFaces,4)-wallEnergy));
assert(heatDifference<=128*eps*max(sum(abs([change(wallFaces,4) wallEnergy]),2))&&any(wallEnergy~=0), ...
    'AxisFVTest:WallThermalEnergy','墙面能量变化不是原signed Fourier热流变化，或测试未激活。');
record=struct('unchangedNonEnergyComponents',notEnergy,'nonwallEnergyFluxUnchanged',true, ...
    'stateCoefficientsAndSAUnchanged',true,'wallHeatEnergyDifferenceWPerRad',heatDifference, ...
    'netWallPowerChangeWPerRad',sum(wallEnergy),'rateOnlyEnergyChanged',true);
end

function result=directionJacobian(Q,mesh,gas,boundary,scale)
% 模块用途：实际有色完整Jacobian对全分量方向差分，另核对不同单位积分导数账。
nq=size(Q,2);n=mesh.cellCount;U=Q./scale;
[J,details]=axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,struct('differenceStep',2e-6,'differenceScheme','central'));
assert(details.success&&details.columnsDifferentiated==n*nq);
v=reshape(sin((1:n*nq)'*.73),n,nq);h=7e-7;
if nq==4
    [rp,ip]=axisymmetricNavierStokesResidual((U+h*v).*scale,mesh,gas,boundary);
    [rm,im]=axisymmetricNavierStokesResidual((U-h*v).*scale,mesh,gas,boundary);
else
    [rp,ip]=axisymmetricRANSResidual((U+h*v).*scale,mesh,gas,boundary);
    [rm,im]=axisymmetricRANSResidual((U-h*v).*scale,mesh,gas,boundary);
end
action=reshape(J*v(:),n,nq);direct=(rp-rm)./scale/(2*h);
relative=max(abs(action-direct),[],'all')/max(abs(direct),[],'all');
assert(relative<1e-5,'AxisFVTest:WallThermalJacobian','显式墙梯度没有进入完整耦合Jacobian。');
integral=sum(action.*scale.*mesh.volumePerRadian,1);
expected=(-(ip.netBoundaryFlux-im.netBoundaryFlux)+(ip.sourceIntegral-im.sourceIntegral))/(2*h);
arithmetic=sum(abs(action.*scale.*mesh.volumePerRadian),1)+abs(expected);
defect=integral-expected;budget=1e-5*arithmetic+512*eps*arithmetic;
assert(all(abs(defect)<=budget),'AxisFVTest:WallThermalJacobianLedger','完整方向导数的边界/源账不闭合。');
result=struct('quantityCount',nq,'directionRelative',relative,'integratedDerivativeDefect',defect, ...
    'integratedDerivativeBudget',budget,'columnCount',details.columnsDifferentiated,'residualEvaluations',details.residualEvaluations);
end

function result=ledger(rate,info,mesh)
% 模块用途：每个有单位守恒量独立预算，并包括共享面相消的舍入误差。
part=rate.*mesh.volumePerRadian;bf=info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:);
source=info.cellSourceIntegral;defect=sum(part,1)+sum(bf,1)-sum(source,1);
physical=sum(abs(part),1)+sum(abs(bf),1)+sum(abs(source),1);
multiplicity=full(sum(mesh.cellFaceMultiplicity,1))';
arithmetic=physical+sum(abs(info.faceFluxIntegralPerRadian).*multiplicity,1);
k=2*mesh.cellCount+numel(mesh.boundaryFaces)+max(full(sum(mesh.cellFaceMultiplicity,2)))+12;
gammaK=k*eps/(1-k*eps);budget=2e-12*physical+gammaK*arithmetic;
assert(all(abs(defect)<=budget),'AxisFVTest:WallThermalLedger','边界/源账未闭合，不能把径向动量误当无源量。');
result=struct('quantityCount',size(rate,2),'defect',defect,'physicalScale',physical,'roundoffBudget',gammaK*arithmetic, ...
    'combinedBudget',budget,'passed',true,'perRadian',true,'unitSeparated',true);
end

function [Q,record]=linearFVMeans(mesh,gas,rho,offset,gradient)
% 模块用途：独立Duffy/Gauss积分每个真实四边形的r权守恒均值，不读取生产梯度。
[q,w]=gauss4();Q=zeros(mesh.cellCount,4);volume=zeros(mesh.cellCount,1);
for cellID=1:mesh.cellCount
    vertices=mesh.vertices(mesh.cells(cellID,:),:);integral=zeros(1,4);
    for tri=[1 2 3;1 3 4]'
        A=vertices(tri(1),:);B=vertices(tri(2),:);C=vertices(tri(3),:);detJ=abs(det([B-A;C-A]));
        for i=1:4
            for j=1:4
                u=q(i);v=q(j);point=A+u*(B-A)+(1-u)*v*(C-A);
                weight=w(i)*w(j)*detJ*(1-u)*point(2);
                T=offset+dot(gradient,point);integral=integral+weight*[rho 0 0 rho*gas.cv*T];
                volume(cellID)=volume(cellID)+weight;
            end
        end
    end
    Q(cellID,:)=integral/volume(cellID);
end
exact=axisymmetricFVConservative([rho*ones(mesh.cellCount,1),zeros(mesh.cellCount,2),offset+mesh.volumeCentroid*gradient'],gas);
relative=max(abs(Q(:,[1 4])-exact(:,[1 4]))./abs(exact(:,[1 4])),[],'all');
volumeRelative=max(abs(volume-mesh.volumePerRadian)./mesh.volumePerRadian);
assert(relative<5e-13&&volumeRelative<5e-13,'AxisFVTest:WallThermalFVMean','独立r权保守积分与线性场真FV均值不符。');
record=struct('method','independent two-triangle Duffy 4x4 Gauss integral of r*Q / integral(r)', ...
    'relativeToLinearExactMean',relative,'volumeRelative',volumeRelative,'realConservativeFVMeans',true);
end

function [x,w]=gauss4()
% 模块用途：冻结四点Gauss-Legendre结点与权重，映射到[0,1]。
x=([-.8611363115940526;-.3399810435848563;.3399810435848563;.8611363115940526]+1)/2;
w=[.3478548451374539;.6521451548625461;.6521451548625461;.3478548451374539]/2;
end

function [Gx,Gr,Twall]=legacyTemperatureGradient(Q,mesh,gas,boundary,wallFaces)
% 模块用途：只为负对照取旧WLS/位移约束迹，独立解析验收不依赖本重构。
P=axisymmetricFVState(Q,gas);f=mesh.faces;bf=mesh.boundaryFaces;names={'iMin','iMax','jMin','jMax'};
trace=P(f.owner(bf),4);Twall=zeros(numel(wallFaces),1);
for k=1:numel(bf)
    face=bf(k);bc=boundary.(names{f.boundarySide(face)});
    if strcmp(bc.type,'wall')
        value=bc.temperature;if isa(value,'function_handle'),value=value(f.centroid(face,:));end
        trace(k)=value;Twall(wallFaces==face)=value;
    else
        value=bc.primitive;if isa(value,'function_handle'),value=value(f.centroid(face,:));end
        trace(k)=value(4);
    end
end
sample=[P(:,4);trace];gx=mesh.gradientX*sample;gr=mesh.gradientR*sample;
Gx=gx(f.owner);Gr=gr(f.owner);inside=f.neighbour>0;
Gx(inside)=.5*(gx(f.owner(inside))+gx(f.neighbour(inside)));Gr(inside)=.5*(gr(f.owner(inside))+gr(f.neighbour(inside)));
right=P(f.owner,4);right(inside)=P(f.neighbour(inside),4);right(bf)=trace;
correction=(right-P(f.owner,4)-Gx.*f.delta(:,1)-Gr.*f.delta(:,2))./sum(f.delta.^2,2);
Gx=Gx+correction.*f.delta(:,1);Gr=Gr+correction.*f.delta(:,2);
end

function boundary=removePolicy(boundary)
% 模块用途：生成完全省略新选项的旧合同，确保兼容性不是显式legacy伪装。
for name={'iMin','iMax','jMin','jMax'}
    bc=boundary.(name{1});
    for key={'thermalGradient','temperatureTangentialGradient'}
        if isfield(bc,key{1}),bc=rmfield(bc,key{1});end
    end
    boundary.(name{1})=bc;
end
end

function T=nonlinearTemperature(point,ri,m,a,h,Tin,Tout)
% 模块用途：常温墙内的非线性光滑温度，墙切向真导数严格为0。
eta=(point(2)-ri-m*point(1))/(1-a*m);s=(point(1)-a*(point(2)-ri))/(1-a*m);
T=Tin+(Tout-Tin)*eta/h+120*sin(pi*eta/h)*sin(pi*s/.08);
end

function k=independentK(T,gas)
% 模块用途：从已固定Sutherland常数独立求墙k，不调用生产本构helper。
tr=gas.transport;mu=tr.referenceViscosity*(T/tr.referenceTemperature).^1.5.* ...
    (tr.referenceTemperature+tr.sutherlandTemperature)./(T+tr.sutherlandTemperature);k=gas.cp*mu/gas.Pr;
end

function records=sourceBindings(root)
% 模块用途：将本测试和实际执行的关键算子绑定到完整文件SHA。
names={'tests/testAxisymmetricWallThermalGradient.m','regionalFlow/axisymmetric/axisymmetricWallThermalGradient.m', ...
    'regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m','regionalFlow/axisymmetric/axisymmetricRANSResidual.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVJacobian.m'};records=struct('path',{},'sha256',{});
for k=1:numel(names),records(k)=struct('path',names{k},'sha256',fileSHA(fullfile(root,strrep(names{k},'/',filesep))));end
end

function digest=fileSHA(path)
% 模块用途：只读输入检查点和完整源码，不改写旧版证据。
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
