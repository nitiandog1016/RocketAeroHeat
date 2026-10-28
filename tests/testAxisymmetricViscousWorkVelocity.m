function report=testAxisymmetricViscousWorkVelocity(folder,v68File,v67File)
%TESTAXISYMMETRICVISCOUSWORKVELOCITY 内面一阶黏性功速度的独立数值合同。
% 模块一：用真实r权保守均值激活强偏斜网格，同Q只允许内部面能量通量改变。
% 模块二：独立重算旧WLS功迹、新算术平均、真实牵引功与配对动量/内能恒等式。
% 模块三：四/五量账、逐列Jacobian、方向导数和真实小BE步分别验收。
% 模块四：冻结v68/v67默认完整info逐值回放；保留真实初态+26km/s旧迹诊断。
% 本模块不证明一般偏斜网格正性、二阶线性精度或真实箭体稳态/实验资格。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v69_validation','viscous_work_velocity');end
if nargin<2,v68File=pick({fullfile(root,'outputs','v68_validation','actual_rocket','rocket_axisymmetric_evidence.mat'), ...
        fullfile(root,'validation','v68','actual_rocket','rocket_axisymmetric_evidence.mat')});end
if nargin<3,v67File=pick({fullfile(root,'outputs','v67_validation','both','tmr_sensitivity_evidence.mat'), ...
        fullfile(root,'validation','v67','both','tmr_sensitivity_evidence.mat')});end
v68SHA=fileSHA(v68File);v67SHA=fileSHA(v67File);
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.015,'Pr',.72,'temperatureRange',[100 1500], ...
    'transport',struct('model','sutherland','referenceViscosity',.015,'referenceTemperature',300,'sutherlandTemperature',110.4)));
[s,eta]=ndgrid([0 .02 .05 .08],[0 .006 .017 .03]);a=1.8;m=.3;ri=.03;
mesh=buildAxisymmetricFVGeometry(s+a*eta,ri+m*s+eta);f=mesh.faces;
reservoir=struct('type','riemannReservoir','primitive',[1.2 8 1 340],'nuTilde',.007);
boundary=struct('iMin',reservoir,'iMax',reservoir, ...
    'jMin',struct('type','wall','temperature',300,'velocity',[0 0],'thermalGradient','dirichlet-tangent'), ...
    'jMax',reservoir);
arithmetic=boundary;arithmetic.viscousWorkVelocity='arithmetic';
[Q,quadrature]=manufacturedFVMeans(mesh,gas,a,m,ri);
Q5=[Q,repmat(1.2*.007,mesh.cellCount,1)];
[oldRate,oldInfo]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);
[newRate,newInfo]=axisymmetricNavierStokesResidual(Q,mesh,gas,arithmetic);
[oldRate5,oldInfo5]=axisymmetricRANSResidual(Q5,mesh,gas,boundary);
[newRate5,newInfo5]=axisymmetricRANSResidual(Q5,mesh,gas,arithmetic);
same4=sameQContract(oldRate,oldInfo,newRate,newInfo,mesh);
same5=sameQContract(oldRate5,oldInfo5,newRate5,newInfo5,mesh);
pair4=pairedWork(newInfo,mesh);pair5=pairedWork(newInfo5,mesh);
ledger4=ledger(newRate,newInfo,mesh);ledger5=ledger(newRate5,newInfo5,mesh);

%% 一、常速度内部迹必须严格保持常值，即使给定静止墙使旧WLS梯度非零
uniformQ=axisymmetricFVConservative(repmat([1.2 35 0 350],mesh.cellCount,1),gas);
[~,uniformInfo]=axisymmetricNavierStokesResidual(uniformQ,mesh,gas,arithmetic);
u=uniformInfo.viscousWorkVelocity;
assert(all(u.velocityXR==[35 0],'all')&&all(u.ownerVelocityProxyXR==u.neighbourVelocityProxyXR,'all'), ...
    'WorkVelocityTest:Constant','内部均匀速度的功迹没有保持常量。');
constantPair=pairedWork(uniformInfo,mesh);
assert(constantPair.maximumOwnerInternalEnergyPower==0&&constantPair.maximumNeighbourInternalEnergyPower==0, ...
    'WorkVelocityTest:Constant','同速两单元之间产生了非零局部配对内能功。');
explicitLegacy=boundary;explicitLegacy.viscousWorkVelocity='legacy-wls';
[legacyRate,legacyInfo]=axisymmetricNavierStokesResidual(Q,mesh,gas,explicitLegacy);
assert(isequal(legacyRate,oldRate)&&isequaln(legacyInfo,oldInfo)&&~isfield(oldInfo,'viscousWorkVelocity'));

%% 二、四/五变量完整Jacobian与真实BE；原正性/温域和非线性接受门槛不变
scale=[1.2 50 50 max(Q(:,4))];scale5=[scale .02];
jacobian4=serialJacobian(Q,mesh,gas,arithmetic,scale);
jacobian5=serialJacobian(Q5,mesh,gas,arithmetic,scale5);
be4=oneBE(Q,mesh,gas,arithmetic,scale);be5=oneBE(Q5,mesh,gas,arithmetic,scale5);

%% 三、真实完整箭体与TMR冻结默认算子完整逐值重放，不因研究选择改旧info
frozen68=load(v68File,'caseDefinition','implicit','initialInfo');item=frozen68.caseDefinition;
[r68,i68]=axisymmetricRANSResidual(frozen68.implicit.Q,item.mesh,item.gas,item.boundary);
assert(isequaln(i68,frozen68.implicit.finalInfo)&& ...
    max(abs(r68)./frozen68.implicit.residualScale,[],'all')==frozen68.implicit.residualNorm);
frozen67=load(v67File,'implicit','mesh','gas','boundary');
[r67,i67]=axisymmetricRANSResidual(frozen67.implicit.Q,frozen67.mesh,frozen67.gas,frozen67.boundary);
assert(isequaln(i67,frozen67.implicit.finalInfo)&& ...
    max(abs(r67)./frozen67.implicit.residualScale,[],'all')==frozen67.implicit.residualNorm);
rocketBoundary=item.boundary;rocketBoundary.viscousWorkVelocity='arithmetic';
[initialRate,initialInfo]=axisymmetricRANSResidual(item.Q,item.mesh,item.gas,rocketBoundary);
[originalRate,originalInfo]=axisymmetricRANSResidual(item.Q,item.mesh,item.gas,item.boundary);
assert(isequaln(originalInfo,frozen68.initialInfo));
rocketContract=sameQContract(originalRate,originalInfo,initialRate,initialInfo,item.mesh);
audit=initialInfo.viscousWorkVelocity;k=find(audit.faceID==214);assert(isscalar(k));
assert(audit.oldVelocityXR(k,1)>2e4&&audit.velocityXR(k,1)<0&& ...
    isequal(audit.velocityXR(k,:),.5*(audit.ownerVelocityProxyXR(k,:)+audit.neighbourVelocityProxyXR(k,:)))&& ...
    isequal(audit.ownerVelocityProxyXR(k,:),audit.neighbourVelocityProxyXR(k,:))&&~audit.rawVelocityWasClipped, ...
    'WorkVelocityTest:Rocket','真实旧无界迹未保留，或新一阶迹没有保持两侧同一自由流。');
rocketProbe=struct('faceID',214,'ownerCell',audit.ownerCell(k),'neighbourCell',audit.neighbourCell(k), ...
    'oldUnlimitedWLSVelocityXR',audit.oldVelocityXR(k,:),'newArithmeticVelocityXR',audit.velocityXR(k,:), ...
    'ownerVelocityXR',audit.ownerVelocityProxyXR(k,:),'neighbourVelocityXR',audit.neighbourVelocityProxyXR(k,:), ...
    'sourceOldExtrapolationPreserved',true,'rawVelocityClipped',false,'notANewSteadyRocketSolution',true);

%% 四、错误策略、错误范围、非有限输入或假共享面显式拒绝
refusals=struct('label',{},'identifier',{});
bad=boundary;bad.viscousWorkVelocity='bounded';reject('unknown-policy',@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad),'AxisFV:ViscousWorkVelocityBoundary');
bad=boundary;bad.jMin.viscousWorkVelocity='arithmetic';reject('policy-on-wall-side',@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad),'AxisFV:ViscousWorkVelocityBoundary');
bad=boundary;bad.iMin.viscousWorkVelocity='arithmetic';reject('policy-on-open-side',@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad),'AxisFV:ViscousWorkVelocityBoundary');
bad=boundary;bad.viscousWorkVelocity={'arithmetic'};reject('policy-not-text',@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad),'AxisFV:ViscousWorkVelocityBoundary');
faceVelocity=oldFaceVelocity(oldInfo,mesh);traction=zeros(numel(f.owner),2);ids=find(f.areaPerRadian>0);
traction(ids,:)=oldInfo.faceViscousFluxPerRadian(ids,2:3)./f.areaPerRadian(ids);
badVelocity=faceVelocity;badVelocity(1,1)=NaN;
reject('nonfinite-work-state',@()axisymmetricViscousWorkVelocity(mesh,oldInfo.primitive(:,2:3),arithmetic,badVelocity,traction),'AxisFV:ViscousWorkVelocityState');
badMesh=mesh;first=find(f.neighbour>0,1);badMesh.faces.areaPerRadian(first)=0;
reject('zero-area-interior-face',@()axisymmetricViscousWorkVelocity(badMesh,oldInfo.primitive(:,2:3),arithmetic,faceVelocity,traction),'AxisFV:ViscousWorkVelocityGeometry');
assert(strcmp(v68SHA,fileSHA(v68File))&&strcmp(v67SHA,fileSHA(v67File)),'WorkVelocityTest:Frozen','只读测试改动了冻结证据。');
report=struct('status','passed-explicit-arithmetic-interior-viscous-work-contract', ...
    'cellCount',mesh.cellCount,'quadrature',quadrature,'sameQ4',same4,'sameQ5',same5,'pair4',pair4,'pair5',pair5, ...
    'constantVelocityPair',constantPair,'ledger4',ledger4,'ledger5',ledger5, ...
    'jacobian4',jacobian4,'jacobian5',jacobian5,'be4',be4,'be5',be5,'rocketSameQ',rocketContract,'rocketProbe',rocketProbe, ...
    'defaultV68CompleteInfoReplayed',true,'defaultV67CompleteInfoReplayed',true,'refusals',refusals, ...
    'frozenV68SHA256',v68SHA,'frozenV67SHA256',v67SHA, ...
    'arithmeticIsSecondOrderLinearReconstruction',false,'generalSkewPositivityEstablished',false, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'viscous_work_evidence.mat'),'report','mesh','gas','Q','Q5','boundary','arithmetic', ...
    'oldInfo','newInfo','oldInfo5','newInfo5','initialInfo','originalInfo');
fid=fopen(fullfile(folder,'viscous_work_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Viscous work velocity passed: same-Q 4/5, pair identities, uniform trace, Jacobian/BE, %d refusals and frozen v68/v67 complete replay.\n',numel(refusals));
    function reject(label,action,identifier)
        try,action();catch exception
            assert(strcmp(exception.identifier,identifier),'WorkVelocityTest:Refusal','%s得到错误拒绝%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('WorkVelocityTest:Refusal','%s错误输入未拒绝。',label);
    end
end

function record=sameQContract(oldRate,old,newRate,new,mesh)
% 模块用途：用实际牵引积分和独立旧迹预测唯一能量改变量，不能只检查“有变化”。
nq=size(oldRate,2);kept=setdiff(1:nq,4);f=mesh.faces;ids=find(f.neighbour>0);bf=mesh.boundaryFaces;
assert(isequal(oldRate(:,kept),newRate(:,kept))&&isequal(old.primitive,new.primitive)&& ...
    isequal(old.cellSourceIntegral,new.cellSourceIntegral)&& ...
    isequal(old.faceFluxIntegralPerRadian(:,kept),new.faceFluxIntegralPerRadian(:,kept))&& ...
    isequal(old.faceFluxIntegralPerRadian(bf,:),new.faceFluxIntegralPerRadian(bf,:))&& ...
    isequaln(old.wall,new.wall)&&isequaln(old.transport,new.transport)&& ...
    isequal(old.cellVelocityGradientX,new.cellVelocityGradientX)&&isequal(old.cellVelocityGradientR,new.cellVelocityGradientR)&& ...
    isequal(old.sourceHoopStress,new.sourceHoopStress)&&isequal(old.netBoundaryFlux,new.netBoundaryFlux));
if nq==5,assert(isequaln(old.sa,new.sa));end
u=oldFaceVelocity(old,mesh);arithmetic=.5*(old.primitive(f.owner(ids),2:3)+old.primitive(f.neighbour(ids),2:3));
a=new.viscousWorkVelocity;
assert(strcmp(a.kind,'explicit-arithmetic-interior-viscous-work-velocity-v1')&&isequal(a.faceID,ids)&& ...
    isequal(a.ownerCell,f.owner(ids))&&isequal(a.neighbourCell,f.neighbour(ids))&& ...
    isequal(a.facePointXR,f.centroid(ids,:))&&isequal(a.areaPerRadianM2,f.areaPerRadian(ids))&& ...
    isequal(a.oldVelocityXR,u(ids,:))&&isequal(a.velocityXR,arithmetic)&& ...
    isequal(a.boundaryVelocityXR,u(bf,:))&&isequal(a.boundaryWorkVelocityXR,u(bf,:))&& ...
    ~a.boundaryVelocityChanged&&~a.tractionChanged&&~a.cellStateClipped&&~a.rawVelocityWasClipped&& ...
    ~a.secondOrderLinearReconstructionClaimed&&~a.generalSkewGridPositivityProved);
momentum=old.faceViscousFluxPerRadian(ids,2:3);
deltaWork=sum((arithmetic-u(ids,:)).*momentum,2);
oldWork=sum(u(ids,:).*momentum,2);newWork=sum(arithmetic.*momentum,2);
deltaFlux=zeros(numel(f.owner),1);deltaFlux(ids)=-deltaWork;
actualFlux=new.faceFluxIntegralPerRadian(:,4)-old.faceFluxIntegralPerRadian(:,4);
budget=256*eps*(abs(old.faceFluxIntegralPerRadian(:,4))+abs(new.faceFluxIntegralPerRadian(:,4))+abs(deltaFlux));
budget(ids)=budget(ids)+256*eps*(abs(oldWork)+abs(newWork)+abs(deltaWork));
assert(all(abs(actualFlux-deltaFlux)<=budget+realmin),'WorkVelocityTest:Energy','能量通量变化不等于原牵引乘功速度差。');
expected=-mesh.incidence*deltaFlux./mesh.volumePerRadian;actual=newRate(:,4)-oldRate(:,4);
rateBudget=(abs(mesh.incidence)*budget+256*eps*(abs(mesh.incidence)*(abs(old.faceFluxIntegralPerRadian(:,4))+ ...
    abs(new.faceFluxIntegralPerRadian(:,4)))+abs(old.cellSourceIntegral(:,4))))./mesh.volumePerRadian;
assert(all(abs(actual-expected)<=rateBudget+realmin)&&any(actual~=0),'WorkVelocityTest:Rate','完整能量残差变化不符或未激活。');
record=struct('quantityCount',nq,'onlyInteriorEnergyChanged',true,'wallHeatAndPowerUnchanged',true,'SAAndCellGradientsUnchanged',true, ...
    'maximumEnergyFluxIdentityDefect',max(abs(actualFlux-deltaFlux)), ...
    'maximumEnergyRateIdentityDefect',max(abs(actual-expected)),'noStateOrGradientClipping',true);
end

function record=pairedWork(info,mesh)
% 模块用途：两单元同一共享牵引的内能功各为半个Δu·牵引，不推定该量普遍非负。
f=mesh.faces;ids=find(f.neighbour>0);u0=info.primitive(f.owner(ids),2:3);u1=info.primitive(f.neighbour(ids),2:3);
uw=info.viscousWorkVelocity.velocityXR;momentum=info.faceViscousFluxPerRadian(ids,2:3);
left=sum((uw-u0).*momentum,2);right=sum((u1-uw).*momentum,2);half=.5*sum((u1-u0).*momentum,2);
% 算术平均与随后减法也会形成舍入；同速/近同速不能只按最终速度差估算。
budget=256*eps*(sum((abs(uw)+abs(u0)).*abs(momentum),2)+ ...
    sum((abs(u1)+abs(uw)).*abs(momentum),2)+sum((abs(u1)+abs(u0)).*abs(momentum),2));
assert(all(abs(left-half)<=budget+realmin)&&all(abs(right-half)<=budget+realmin));
record=struct('interiorFaces',numel(ids),'equalHalfPairIdentityMaximumDefect',max([abs(left-half);abs(right-half)]), ...
    'maximumOwnerInternalEnergyPower',max(abs(left)),'maximumNeighbourInternalEnergyPower',max(abs(right)), ...
    'minimumPairDiscreteInternalPower',min(2*half),'negativePairCount',sum(half<0), ...
    'generalSkewPositivityEstablished',false,'units','W/rad');
end

function value=oldFaceVelocity(info,mesh)
% 模块用途：从完整原始量和单元梯度独立重建旧功速度，不读取新helper的old副本。
f=mesh.faces;o=f.owner;j=f.neighbour;inside=j>0;P=info.primitive;
value=P(o,2:3);gx=info.cellVelocityGradientX;gr=info.cellVelocityGradientR;
dl=f.centroid(inside,:)-mesh.volumeCentroid(o(inside),:);
dr=f.centroid(inside,:)-f.neighbourShift(inside,:)-mesh.volumeCentroid(j(inside),:);
value(inside,:)=.5*(P(o(inside),2:3)+gx(o(inside),:).*dl(:,1)+gr(o(inside),:).*dl(:,2)+ ...
    P(j(inside),2:3)+gx(j(inside),:).*dr(:,1)+gr(j(inside),:).*dr(:,2));
value(mesh.boundaryFaces,:)=info.faceViscousGradientTrace(mesh.boundaryFaces,2:3);
end

function record=serialJacobian(Q,mesh,gas,boundary,scale)
% 模块用途：完整有色Jacobian与独立逐列中央差分对照，再检查全分量方向及积分账。
n=mesh.cellCount;nq=size(Q,2);U=Q./scale;
[J,details]=axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,struct('differenceStep',2e-6,'differenceScheme','central'));
assert(details.success&&details.columnsDifferentiated==numel(Q));dense=zeros(n*nq);h=2e-6*max(1,abs(U(:)));
for k=1:numel(Q)
    plus=U;minus=U;plus(k)=plus(k)+h(k);minus(k)=minus(k)-h(k);
    rp=spatial(plus.*scale,mesh,gas,boundary)./scale;rm=spatial(minus.*scale,mesh,gas,boundary)./scale;
    dense(:,k)=(rp(:)-rm(:))/(2*h(k));
end
relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');assert(relative<2e-6);
v=reshape(sin((1:numel(Q))'*.63),size(Q));h=7e-7;
[rp,ip]=spatial((U+h*v).*scale,mesh,gas,boundary);[rm,im]=spatial((U-h*v).*scale,mesh,gas,boundary);
action=reshape(J*v(:),size(Q));direct=(rp-rm)./scale/(2*h);
direction=max(abs(action-direct),[],'all')/max(abs(direct),[],'all');assert(direction<1e-5);
integral=sum(action.*scale.*mesh.volumePerRadian,1);expected=(-(ip.netBoundaryFlux-im.netBoundaryFlux)+(ip.sourceIntegral-im.sourceIntegral))/(2*h);
arithmetic=sum(abs(action.*scale.*mesh.volumePerRadian),1)+abs(expected);
assert(all(abs(integral-expected)<=1e-5*arithmetic+512*eps*arithmetic));
record=struct('quantityCount',nq,'coloredVsSerialRelative',relative,'directionRelative',direction, ...
    'allColumnsChecked',numel(Q),'integratedDerivativeDefect',integral-expected);
end

function record=oneBE(Q,mesh,gas,boundary,scale)
% 模块用途：实际完整耦合BE一步；原绝对/相对Newton与五量账接受，不伪造稳态。
opts=struct('maxSteps',1,'initialPseudoCFL',10,'maximumPseudoCFL',10,'maximumPseudoStep',2e-8, ...
    'residualScale',scale,'maxNewtonIterations',18,'jacobianRefreshEvery',1);
result=advanceAxisymmetricSteady(Q,mesh,gas,boundary,opts);
assert(result.acceptedSteps==1&&~result.steadyConverged&&result.history.dt(end)==opts.maximumPseudoStep);
[rate,info]=spatial(result.Q,mesh,gas,boundary);dt=result.history.dt(end);
defect=max(abs((result.Q-Q-dt*rate)./scale),[],'all');
assert(defect<=result.history.newtonTarget(end)+512*eps*max(abs(Q./scale),[],'all')&&isequaln(info,result.finalInfo));
raw=sum((result.Q-Q).*mesh.volumePerRadian,1)+dt*(info.netBoundaryFlux-info.sourceIntegral);
relative=abs(raw)./(sum(mesh.volumePerRadian)*scale);assert(all(relative<=result.options.conservationTolerance));
assert(any(result.Q~=Q,'all')&&all(isfinite(result.Q),'all'));
record=struct('quantityCount',size(Q,2),'acceptedSteps',1,'steadyConverged',false,'actualDt',dt, ...
    'BEStateDefect',defect,'actualNewtonTarget',result.history.newtonTarget(end),'ledgerRelative',relative, ...
    'solverResidualEvaluations',result.residualEvaluations,'stateWasClipped',false,'physicalTimeAdvanced',0);
end

function record=ledger(rate,info,mesh)
% 模块用途：每个有单位量分开传播共享面、源和归并舍入，不用混合单位固定阈值。
part=rate.*mesh.volumePerRadian;source=info.cellSourceIntegral;flux=info.faceFluxIntegralPerRadian;
bf=flux(mesh.boundaryFaces,:);defect=sum(part,1)+sum(bf,1)-sum(source,1);
multiplicity=full(sum(mesh.cellFaceMultiplicity,1))';
arithmetic=sum(abs(part),1)+sum(abs(bf),1)+sum(abs(source),1)+sum(abs(flux).*multiplicity,1);
operationCount=8*mesh.cellCount+4*size(flux,1)+32;gamma=operationCount*eps/(1-operationCount*eps);
budget=gamma*arithmetic;assert(all(abs(defect)<=budget+realmin));
record=struct('quantityCount',size(rate,2),'rawDefect',defect,'formationBudget',budget,'unitSeparated',true,'perRadian',true);
end

function [Q,record]=manufacturedFVMeans(mesh,gas,a,m,ri)
% 模块用途：独立Duffy/Gauss积分r*Q；动能含二次/四次项，不能用体心原始量冒充FV均值。
x=([-.8611363115940526;-.3399810435848563;.3399810435848563;.8611363115940526]+1)/2;
w=[.3478548451374539;.6521451548625461;.6521451548625461;.3478548451374539]/2;
Q=zeros(mesh.cellCount,4);vol=zeros(mesh.cellCount,1);
for cellID=1:mesh.cellCount
    vertices=mesh.vertices(mesh.cells(cellID,:),:);integral=zeros(1,4);
    for tri=[1 2 3;1 3 4]'
        A=vertices(tri(1),:);B=vertices(tri(2),:);C=vertices(tri(3),:);detJ=abs(det([B-A;C-A]));
        for i=1:4
            for j=1:4
                u=x(i);v=x(j);p=A+u*(B-A)+(1-u)*v*(C-A);
                weight=w(i)*w(j)*detJ*(1-u)*p(2);eta=(p(2)-ri-m*p(1))/(1-a*m);s=p(1)-a*eta;
                ux=7+300*eta+1000*eta^2+80*s*eta;ur=.15*ux+40*eta^2;T=330+300*eta+100*s*eta;rho=1.2;
                integral=integral+weight*[rho,rho*ux,rho*ur,rho*gas.cv*T+.5*rho*(ux^2+ur^2)];
                vol(cellID)=vol(cellID)+weight;
            end
        end
    end
    Q(cellID,:)=integral/vol(cellID);
end
relative=max(abs(vol-mesh.volumePerRadian)./mesh.volumePerRadian);assert(relative<5e-13);
record=struct('method','independent two-triangle 4x4 Duffy/Gauss integral of r*Q', ...
    'volumeRelative',relative,'trueConservativeFVMeans',true,'EOSProxyNotCalledExactMeanTemperature',true);
end

function [rate,info]=spatial(Q,mesh,gas,boundary)
if size(Q,2)==5,[rate,info]=axisymmetricRANSResidual(Q,mesh,gas,boundary);else,[rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);end
end
function path=pick(paths)
path='';for k=1:numel(paths),if exist(paths{k},'file'),path=paths{k};break;end,end
assert(~isempty(path),'WorkVelocityTest:Evidence','缺少冻结真实证据，不能跳过旧路径回放。');
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
