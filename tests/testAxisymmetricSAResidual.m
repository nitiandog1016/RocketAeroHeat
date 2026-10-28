function report=testAxisymmetricSAResidual(folder)
%TESTAXISYMMETRICSARESIDUAL 五量SA共享通量、完整源、墙热与保守账的定向校核。
% 模块一：N=0严格退化原四量NS；常N只检查内部/开放面质量一致性。
% 模块二：线性rho/N内部梯度有独立解析值，逐项复算SA源及密度修正。
% 模块三：每弧度五量分别记账，故意漏密度项/重复轴源均必须被拒绝。
% 模块四：实际出流/零质量边界的零法向扩散、墙muT=0及非法状态拒绝。
% 仅为空间算子数值校核，不是湍流稳态、实验精度或生产面积验收。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v65_validation','sa_residual');end
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[200 800], ...
    'transport',struct('model','sutherland','referenceViscosity',.05,'referenceTemperature',300,'sutherlandTemperature',110.4)));
refusals=struct('label',{},'identifier',{});

%% 一、零SA变量的旧NS精确退化；真实温度梯度使墙热不是零值占位
[X,R]=ndgrid([0 .02 .05 .08],linspace(.03,.06,5));
mesh0=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));c0=mesh0.volumeCentroid;
Q0=axisymmetricFVConservative([ones(mesh0.cellCount,1),-2+10*c0(:,2),zeros(mesh0.cellCount,1),300+3000*c0(:,2)],gas);
per=struct('type','periodic');boundary0=struct('iMin',per,'iMax',per, ...
    'jMin',struct('type','wall','temperature',390,'velocity',[-1.7 0]), ...
    'jMax',struct('type','wall','temperature',480,'velocity',[-1.4 0]));
[laminarRate,laminarInfo]=axisymmetricNavierStokesResidual(Q0,mesh0,gas,boundary0);
[zeroRate,zeroInfo]=axisymmetricRANSResidual([Q0 zeros(mesh0.cellCount,1)],mesh0,gas,boundary0);
assert(isequal(zeroRate(:,1:4),laminarRate)&&all(zeroRate(:,5)==0)&& ...
    isequal(zeroInfo.faceFluxIntegralPerRadian(:,1:4),laminarInfo.faceFluxIntegralPerRadian)&& ...
    isequal(zeroInfo.cellSourceIntegral(:,1:4),laminarInfo.cellSourceIntegral)&& ...
    all(zeroInfo.sa.cellMuT==0)&&all(zeroInfo.sa.faceMuT==0), ...
    'AxisSATest:ZeroReduction','N=0没有精确退化旧NS，或出现伪SA源/涡黏性。');
zeroLedger=ledger(zeroRate,zeroInfo,mesh0,true);

%% 二、常N在实际内部/开放面逐面等于N乘同一NS质量通量；墙N=0另行保留
[X,R]=ndgrid(linspace(0,.08,5),linspace(.03,.06,5));mesh=buildAxisymmetricFVGeometry(X,R);n=mesh.cellCount;
t=(1:n)';Pconstant=[1+.025*sin(t),-2+.1*cos(.7*t),.015*sin(.4*t),450+8*cos(.3*t)];
constantQ4=axisymmetricFVConservative(Pconstant,gas);constantN=.27;
flatQ4=axisymmetricFVConservative(repmat([1 -2 0 450],n,1),gas);
open=struct('type','riemannReservoir','primitive',[1 -2 0 450],'nuTilde',constantN);
boundary=struct('iMin',open,'iMax',open, ...
    'jMin',struct('type','wall','temperature',450,'velocity',[-2 0]),'jMax',open);
[constantRate,constantInfo]=axisymmetricRANSResidual([constantQ4 constantQ4(:,1)*constantN],mesh,gas,boundary);
f=mesh.faces;inside=f.neighbour>0;bf=mesh.boundaryFaces;
wall=f.boundarySide==3;eligible=inside|(~inside&~wall);
massConsistent=constantInfo.sa.faceAdvectiveFluxPerRadian(eligible)-constantN*constantInfo.faceConvectiveFluxPerRadian(eligible,1);
massConsistencyBudget=16*eps*max(abs(constantInfo.sa.faceAdvectiveFluxPerRadian(eligible)));
assert(max(abs(massConsistent))<=massConsistencyBudget, ...
    'AxisSATest:MassConsistency','SA共享对流未逐面使用NS实际质量通量。');
% 独立owner物理质量通量故意不使用Riemann：非均匀四量必须使此错误替代明显失败。
ownerPhysicalMass=Pconstant(f.owner,1).*sum(Pconstant(f.owner,2:3).*f.normal,2).*f.areaPerRadian;
riemannMass=constantInfo.faceConvectiveFluxPerRadian(:,1);
massRiemannDifference=max(abs(riemannMass(eligible)-ownerPhysicalMass(eligible)));
massRiemannRelative=massRiemannDifference/max(abs(riemannMass(eligible)));
wrongAdvective=constantN*ownerPhysicalMass;
wrongMassConsistencyDefect=wrongAdvective(eligible)-constantN*riemannMass(eligible);
assert(massRiemannRelative>1e-3&&max(abs(wrongMassConsistencyDefect))>1e6*massConsistencyBudget&& ...
    max(abs(riemannMass(inside)-ownerPhysicalMass(inside)))>1e-3*max(abs(riemannMass(inside))), ...
    'AxisSATest:MassConsistencyNegativeControl','非均匀常N测试未发现另算owner物理质量通量的错误mN。');
assert(all(constantInfo.sa.boundaryTraceNuTilde(f.boundarySide(bf)==3)==0)&& ...
    any(constantInfo.faceIncoming(bf))&&any(~constantInfo.faceIncoming(bf)&~wall(bf)), ...
    'AxisSATest:BoundaryExcitation','常N测试缺少真实入/出流，或墙N未严格为0。');
constantLedger=ledger(constantRate,constantInfo,mesh,true);

%% 三、非均匀rho/N与真实体积求积；内部格线性场梯度由解析斜率独立给定
c=mesh.volumeCentroid;rho=1+.2*c(:,1)+.3*c(:,2);N=.15+.04*c(:,1)+.8*c(:,2);
P=[rho,-3+.4*c(:,1)+2*c(:,2),zeros(n,1),400+200*c(:,2)];Q4=axisymmetricFVConservative(P,gas);Q=[Q4 rho.*N];
primitive=@(point)[1+.2*point(1)+.3*point(2),-3+.4*point(1)+2*point(2),0,400+200*point(2)];
inN=@(point).15+.04*point(1)+.8*point(2);
edge=struct('type','riemannReservoir','primitive',primitive,'nuTilde',inN);
nonuniformBoundary=struct('iMin',edge,'iMax',edge, ...
    'jMin',struct('type','wall','temperature',406,'velocity',@(point)[-3+.4*point(1)+.06 0]), ...
    'jMax',struct('type','wall','temperature',412,'velocity',@(point)[-3+.4*point(1)+.12 0]));
[rate,info]=axisymmetricRANSResidual(Q,mesh,gas,nonuniformBoundary);sa=info.sa;
[ix,jr]=ind2sub(mesh.shape,(1:n)');interior=ix>1&ix<mesh.shape(1)&jr>1&jr<mesh.shape(2);
assert(nnz(interior)>=4&&max(abs(sa.cellDensityGradientX(interior)-.2))<2e-12&& ...
    max(abs(sa.cellDensityGradientR(interior)-.3))<2e-12&& ...
    max(abs(sa.cellNuTildeGradientX(interior)-.04))<2e-12&& ...
    max(abs(sa.cellNuTildeGradientR(interior)-.8))<2e-12, ...
    'AxisSATest:AnalyticGradient','内部真实线性rho/N梯度不满足独立解析斜率。');
Omega=abs(info.cellVelocityGradientX(:,2)-info.cellVelocityGradientR(:,1));
independent=localModel(rho,P(:,4),N,Omega,sa.wallDistance,gas);
sigma=2/3;densityExpected=-(independent.nu+N)/sigma.* ...
    (sa.cellDensityGradientX.*sa.cellNuTildeGradientX+sa.cellDensityGradientR.*sa.cellNuTildeGradientR);
analyticDensity=-(independent.nu(interior)+N(interior))/sigma*(.2*.04+.3*.8);
assert(relative(sa.sourceDensityCorrection(interior),analyticDensity)<2e-11&&any(abs(analyticDensity)>1e-3), ...
    'AxisSATest:DensityCorrection','完整密度修正未满足独立线性场公式，或激励太小。');
production=rho.*independent.production;destruction=rho.*independent.destruction;
gradSquare=rho*(.622/sigma).*(sa.cellNuTildeGradientX.^2+sa.cellNuTildeGradientR.^2);
fullSource=production-destruction+gradSquare+densityExpected;
sourceRelative=[relative(sa.sourceProduction,production),relative(sa.sourceDestruction,destruction), ...
    relative(sa.sourceGradientSquare,gradSquare),relative(sa.sourceDensityCorrection,densityExpected), ...
    relative(sa.cellSourceIntegral,fullSource.*mesh.volumePerRadian)];
assert(max(sourceRelative)<3e-12,'AxisSATest:SourceFormula','SA完整源、符号或真实体积求积不满足独立公式。');
assert(relative(sa.diffusionCoefficient,(independentMu(info.transport.faceTemperature,gas)+sa.faceDensity.*sa.faceNuTilde)/sigma)<3e-14&& ...
    any(abs(sa.diffusionCoefficient-sa.effectiveFaceMu/sigma)>1e-4), ...
    'AxisSATest:DiffusionCoefficient','SA扩散系数错误混成有效应力黏度。');
nonuniformLedger=ledger(rate,info,mesh,true);
missingDensity=rate;missingDensity(:,5)=missingDensity(:,5)-densityExpected;
missingDensityLedger=ledger(missingDensity,info,mesh,false);
duplicatedAxis=rate;duplicatedAxis(:,5)=duplicatedAxis(:,5)+(independent.nu+N)/sigma.*sa.cellNuTildeGradientR./c(:,2);
duplicateAxisLedger=ledger(duplicatedAxis,info,mesh,false);
assert(~missingDensityLedger.passed&&abs(missingDensityLedger.defect(5))>100*missingDensityLedger.budget(5)&& ...
    ~duplicateAxisLedger.passed&&abs(duplicateAxisLedger.defect(5))>100*duplicateAxisLedger.budget(5), ...
    'AxisSATest:LedgerNegativeControl','五量账未能发现漏密度修正或重复轴源。');
[sameStateLaminar,laminarWall]=axisymmetricNavierStokesResidual(Q4,mesh,gas,nonuniformBoundary); %#ok<ASGLU>
w=info.wall;wallMu=independentMu(w.temperatureInput,gas);
assert(sa.wallMuTIsZero&&~sa.wallHeatEnhancementFactorApplied&&all(sa.faceMuT(w.faceID)==0)&& ...
    isequal(w.heatFluxIntoWall,laminarWall.wall.heatFluxIntoWall)&& ...
    isequal(w.viscousTraction,laminarWall.wall.viscousTraction)&& ...
    relative(w.viscosity,wallMu)<3e-14&&relative(w.thermalConductivity,wallMu*gas.cp/gas.Pr)<3e-14&& ...
    any(abs(w.heatFluxIntoWall)>1e-3), ...
    'AxisSATest:MolecularWall','墙面出现涡黏性或人工热增强，或Tw处分子系数错误。');

%% 四、出流及恰零质量开放面的SA零法向梯度；负状态不得裁值后接受
NzeroMass=.2+.1*c(:,1)+2*c(:,2);
[zeroMassRate,zeroMassInfo]=axisymmetricRANSResidual([flatQ4 NzeroMass],mesh,gas,boundary); %#ok<ASGLU>
z=zeroMassInfo.sa;normalGradient=z.faceNuTildeGradientX.*f.normal(:,1)+z.faceNuTildeGradientR.*f.normal(:,2);
zeroMassFaces=bf(f.boundarySide(bf)==4);
outFaces=bf(zeroMassInfo.sa.faceMassFluxPerRadian(bf)>0);
assert(all(zeroMassInfo.sa.faceMassFluxPerRadian(zeroMassFaces)==0)&&~isempty(outFaces)&& ...
    all(z.faceNormalZeroGradient([zeroMassFaces;outFaces]))&& ...
    max(abs(normalGradient([zeroMassFaces;outFaces])))<1e-12&& ...
    max(abs(z.faceDiffusiveFluxPerRadian([zeroMassFaces;outFaces])))<1e-13, ...
    'AxisSATest:OutflowDiffusion','真实出流/零质量边界没有保留SA零法向扩散。');
bad=Q;bad(1,5)=-1e-12;reject('negative-W',@()axisymmetricRANSResidual(bad,mesh,gas,nonuniformBoundary),'AxisSA:State');
bad=Q;bad(1,5)=NaN;reject('nonfinite-W',@()axisymmetricRANSResidual(bad,mesh,gas,nonuniformBoundary),'AxisSA:State');
badBoundary=boundary;badBoundary.iMax=rmfield(badBoundary.iMax,'nuTilde');
reject('missing-actual-inflow-N',@()axisymmetricRANSResidual([constantQ4 constantQ4(:,1)*constantN],mesh,gas,badBoundary),'AxisSA:Boundary');
badBoundary=boundary;badBoundary.iMax.nuTilde=-.1;
reject('negative-inflow-N',@()axisymmetricRANSResidual([constantQ4 constantQ4(:,1)*constantN],mesh,gas,badBoundary),'AxisSA:Boundary');
assert(~info.physicalSolution&&~info.eligibleForSurfaceCoverage&&~info.turbulenceValidated&&info.turbulenceImplemented);
report=struct('status','passed-SA-five-variable-spatial-flux-source-wall-boundary-ledger', ...
    'cellCount',n,'zeroReductionExact',true,'constantNInternalOpenMassFluxDefect',max(abs(massConsistent)), ...
    'constantNMassConsistencyBudget',massConsistencyBudget,'RiemannVsOwnerMassDifference',massRiemannDifference, ...
    'RiemannVsOwnerMassRelative',massRiemannRelative,'wrongOwnerMassTimesNDefect',max(abs(wrongMassConsistencyDefect)), ...
    'wrongOwnerPhysicalMassNegativeControlRejected',true, ...
    'independentSourceRelative',sourceRelative,'analyticDensityCorrection',analyticDensity, ...
    'zeroLedger',zeroLedger,'constantLedger',constantLedger,'nonuniformLedger',nonuniformLedger, ...
    'omittedDensityCorrectionNegativeControl',missingDensityLedger,'duplicateAxisSourceNegativeControl',duplicateAxisLedger, ...
    'zeroMassOpenFaceCount',numel(zeroMassFaces),'outflowFaceCount',numel(outFaces),'refusals',refusals, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'sa_residual_evidence.mat'),'report','gas','mesh','Q','rate','info','boundary','nonuniformBoundary', ...
    'mesh0','Q0','boundary0','zeroRate','zeroInfo','laminarRate','laminarInfo','constantQ4','constantN','constantRate','constantInfo', ...
    'independent','fullSource','zeroMassRate','zeroMassInfo','NzeroMass','Pconstant','flatQ4', ...
    'ownerPhysicalMass','wrongAdvective','wrongMassConsistencyDefect');
fid=fopen(fullfile(folder,'sa_residual_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('SA residual passed: exact N0, shared mN, source/density, five-quantity ledger, wall and %d refusals.\n',numel(refusals));
    function reject(label,action,expected)
        try,action();catch e
            assert(strcmp(e.identifier,expected),'AxisSATest:UnexpectedRefusal','%s期望%s，实际%s。',label,expected,e.identifier);
            refusals(end+1)=struct('label',label,'identifier',e.identifier);return; %#ok<AGROW>
        end
        error('AxisSATest:MissingRefusal','非法输入%s未拒绝。',label);
    end
end

function result=ledger(rate,info,mesh,required)
% 模块用途：五个不同单位的量分别预算，内部共享面的大数相消也计入舍入预算。
change=sum(rate.*mesh.volumePerRadian,1);boundary=sum(info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:),1);
source=sum(info.cellSourceIntegral,1);defect=change+boundary-source;
arithmetic=sum(abs(rate.*mesh.volumePerRadian),1)+sum(abs(info.faceFluxIntegralPerRadian),1)+sum(abs(info.cellSourceIntegral),1);
physical=sum(abs(info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:)),1)+sum(abs(info.cellSourceIntegral),1);
budget=256*eps*arithmetic+2e-12*physical;passed=all(abs(defect)<=budget);
assert(~required||passed,'AxisSATest:Ledger','每弧度五量边界/源账不闭合，未用统一max(1)掩盖不同单位。');
result=struct('defect',defect,'budget',budget,'passed',passed,'unitSeparated',true,'perRadian',true);
end

function mu=independentMu(T,gas)
tr=gas.transport;mu=tr.referenceViscosity*(T/tr.referenceTemperature).^1.5.* ...
    (tr.referenceTemperature+tr.sutherlandTemperature)./(T+tr.sutherlandTemperature);
end

function value=relative(a,b)
den=max(abs(b),[],'all');if den==0,value=max(abs(a),[],'all');else,value=max(abs(a-b),[],'all')/den;end
end

function s=localModel(rho,T,N,Omega,d,gas)
% 模块用途：由冻结方程直接独立复算局部生产/破坏，不调用生产SAModel辅助函数。
mu=independentMu(T,gas);nu=mu./rho;chi=N./nu;fv1=chi.^3./(chi.^3+7.1^3);
fv2=1-chi./(1+chi.*fv1);bar=N.*fv2./(.41^2*d.^2);strain=Omega+bar;
k=bar<-.7*Omega;strain(k)=Omega(k)+Omega(k).*(.7^2*Omega(k)+.9*bar(k))./((.9-2*.7)*Omega(k)-bar(k));
r=10*ones(size(N));k=strain>0;r(k)=min(N(k)./(strain(k)*.41^2.*d(k).^2),10);
g=r+.3*(r.^6-r);fw=g.*((1+2^6)./(g.^6+2^6)).^(1/6);
cw1=.1355/.41^2+(1+.622)/(2/3);
s=struct('nu',nu,'production',.1355*strain.*N,'destruction',cw1*fw.*(N./d).^2);
end
