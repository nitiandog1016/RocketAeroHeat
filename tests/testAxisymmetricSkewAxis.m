function report=testAxisymmetricSkewAxis(folder,legacyEvidenceFile)
%TESTAXISYMMETRICSKEWAXIS 真偏斜轴邻近的独立几何、奇偶性与完整方程校核。
% 模块一：真实r权保守均值用三角Duffy/Gauss积分，不把一般EOS点代理称精确均值。
% 模块二：独立2x2正规方程核对镜像采样；区分真实轴坐标、系数代理和原始外推。
% 模块三：检查完整NS/SA、密度梯度修正、稀疏Jacobian和真实小步BE账。
% 模块四：冻结v67默认结果逐值重放；坏策略/假轴/坏实际状态拒绝，未用负外推不裁值。
% 本测试只开通轴边界数值合同，不认证真实箭体、分离或新增生产面积。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v68_validation','skew_axis');end
if nargin<2
    legacyEvidenceFile=fullfile(root,'outputs','v67_validation','both','tmr_sensitivity_evidence.mat');
    if ~exist(legacyEvidenceFile,'file'),legacyEvidenceFile=fullfile(root,'validation','v67','both','tmr_sensitivity_evidence.mat');end
end
assert(exist(legacyEvidenceFile,'file')==2,'AxisTest:Legacy','必须提供冻结v67真实五量稳态，不跳过默认回放。');
if ~exist(folder,'dir'),mkdir(folder);end
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[40 2000]));
variable=gas;variable.transport=struct('model','sutherland','referenceViscosity',.05, ...
    'referenceTemperature',300,'sutherlandTemperature',110.4);variable=axisymmetricFVGas(variable);
[s,r]=ndgrid([0 .018 .045 .08],[0 .004 .012 .025]);shear=1.8;
mesh=buildAxisymmetricFVGeometry(s+shear*r,r);originalMesh=mesh;f=mesh.faces;axisIDs=find(f.boundarySide==3);
ratio=abs(f.centroid(axisIDs,1)-mesh.volumeCentroid(f.owner(axisIDs),1))./mesh.volumeCentroid(f.owner(axisIDs),2);
assert(min(ratio)>1,'AxisTest:Excitation','没有真正激活偏斜轴面与r权体质心位移。');
policy=struct('type','axis','symmetryGradient','mirror-parity');
refusals=struct('label',{},'identifier',{});

%% 一、偶温度线性场的真保守均值、独立WLS及零面积轴系数代理EOS
rho=1.2;T0=350;gT=30;
evenPrimitive=@(p)[rho+zeros(size(p,1),1),zeros(size(p,1),2),T0+gT*p(:,1)];
[evenQ,evenQuadrature]=trueMeans(mesh,gas,evenPrimitive);
boundary=allTraceBoundary(evenPrimitive,policy,.1);
legacy=boundary;legacy.jMin=struct('type','axis');
reject('skew-axis-without-explicit-policy',@()axisymmetricNavierStokesResidual(evenQ,mesh,gas,legacy),{'AxisFV:AxisSkew'});
[evenRate,evenInfo]=axisymmetricNavierStokesResidual(evenQ,mesh,gas,boundary);
evenP=axisymmetricFVState(evenQ,gas);evenA=evenInfo.axisSymmetry;
evenSamples=nsBoundarySamples(mesh,gas,boundary,evenP);
[oracleX,oracleR,oracleGeometry]=independentWLS(mesh,evenP(:,[2 3 4 5]),evenSamples,[1 -1 1 1]);
near(evenInfo.cellVelocityGradientX,oracleX(:,1:2),2e-10,2e-10,'NS全部速度X梯度');
near(evenInfo.cellVelocityGradientR,oracleR(:,1:2),2e-10,2e-10,'NS全部速度r梯度');
owners=f.owner(axisIDs);near(evenA.cellGradientXAtAxisOwner,oracleX(owners,:),2e-10,2e-10,'axis owner X梯度');
near(evenA.cellGradientRAtAxisOwner,oracleR(owners,:),2e-10,2e-10,'axis owner r梯度');
near(evenA.cellGradientXAtAxisOwner(:,3),gT+zeros(size(owners)),2e-10,2e-10,'真均值偶T轴向导数');
near(evenA.cellGradientRAtAxisOwner(:,3),zeros(size(owners)),2e-10,2e-10,'真均值偶T轴法向导数');
assertAxisGeometry(evenA,mesh,axisIDs);
expectedEven=evenPrimitive(f.centroid(axisIDs,:));
expectedProxy=evenP(owners,1:4);expectedProxy(:,3)=0;
near(evenA.coefficientPrimitiveRhoUxUrT,expectedProxy,0,0,'零面积轴owner投影系数代理');
near(evenA.coefficientPressurePa,expectedProxy(:,1).*gas.R.*expectedProxy(:,4),0,0,'系数代理EOS压力');
near(evenA.rawTangentialPrimitiveRhoUxUrT,expectedEven,5e-11,5e-11,'独立仿射场的未用原始切向外推');
assertAxisProxyIdentity(evenA);
assert(all(evenInfo.faceFluxIntegralPerRadian(axisIDs,:)==0,'all')&&isequaln(mesh,originalMesh), ...
    'AxisTest:PhysicalGeometry','虚采样改变了真实网格或非零轴面积通量。');
evenLedger=ledger(evenRate,evenInfo,mesh);

%% 二、奇ur线性场：选取总能量为常量，避免一般均值EOS造成伪解析温度
beta=4;oddT0=500;
oddPrimitive=@(p)[rho+zeros(size(p,1),1),zeros(size(p,1),1),beta*p(:,2), ...
    oddT0-beta^2*p(:,2).^2/(2*gas.cv)];
[oddQ,oddQuadrature]=trueMeans(mesh,gas,oddPrimitive);
oddBoundary=allTraceBoundary(oddPrimitive,policy,.1);
[oddRate,oddInfo]=axisymmetricNavierStokesResidual(oddQ,mesh,gas,oddBoundary);
oddP=axisymmetricFVState(oddQ,gas);oddA=oddInfo.axisSymmetry;
near(oddP(:,3),beta*mesh.volumeCentroid(:,2),2e-12,2e-12,'真均值奇ur');
near(oddP(:,4),oddT0-beta^2*mesh.volumeCentroid(:,2).^2/(2*gas.cv),3e-12,3e-12,'特定能量场的EOS温度');
near(oddInfo.cellVelocityGradientX(:,2),zeros(mesh.cellCount,1),1e-9,1e-9,'奇ur切向梯度');
near(oddInfo.cellVelocityGradientR(:,2),beta+zeros(mesh.cellCount,1),1e-9,1e-9,'奇ur径向梯度');
near(oddA.faceGradientX(:,2),zeros(size(owners)),1e-10,1e-10,'轴面ur切向正则条件');
near(oddA.faceGradientR(:,2),beta+zeros(size(owners)),1e-10,1e-10,'轴面ur径向正则极限');
near(oddA.hoopLimit,beta+zeros(size(owners)),1e-10,1e-10,'轴面ur/r正则极限');
near(oddInfo.sourceHoopStress,2*gas.mu*beta/3+zeros(mesh.cellCount,1),1e-9,1e-9,'独立圆周黏性应力');
assert(all(oddA.coefficientPrimitiveRhoUxUrT(:,3)==0)&&all(oddInfo.faceFluxIntegralPerRadian(axisIDs,:)==0,'all'));
oddLedger=ledger(oddRate,oddInfo,mesh);

%% 三、均匀轴向流的几何平衡；SA的零N不是给定非零湍流稳态
uniform=@(p)[.1+zeros(size(p,1),1),-12+zeros(size(p,1),1),zeros(size(p,1),1),350+zeros(size(p,1),1)];
uniformBoundary=allTraceBoundary(uniform,policy,0);
uniformBoundary.jMax=struct('type','wall','temperature',350,'velocity',[-12 0]);
uniformQ=axisymmetricFVConservative(uniform(mesh.volumeCentroid),variable);
[uniformRate,uniformInfo]=axisymmetricNavierStokesResidual(uniformQ,mesh,variable,uniformBoundary);
uniformBalance=uniformResidualBudget(uniformRate,uniformInfo,mesh,variable,uniformBoundary);
[uniformRate5,uniformInfo5]=axisymmetricRANSResidual([uniformQ zeros(mesh.cellCount,1)],mesh,variable,uniformBoundary);
uniformBalance5=uniformResidualBudget(uniformRate5,uniformInfo5,mesh,variable,uniformBoundary);
assert(all(uniformRate5(:,5)==0)&&all(uniformInfo5.sa.sourceDensityCorrection==0));

% 实际C域使用iMin/iMax两个轴端，不能只证明jMin上的剪切网格。
[angle,radius]=ndgrid([0 .25 .8 1.5 2.3 2.9 pi],[.05 .065 .09]);
sphereX=cos(angle).*radius;sphereR=sin(angle).*radius;sphereR([1 end],:)=0;
sphereMesh=buildAxisymmetricFVGeometry(sphereX,sphereR);
sphereEdge=struct('type','prescribedTrace','primitive',[1 0 0 311],'nuTilde',0);
sphereBoundary=struct('iMin',policy,'iMax',policy, ...
    'jMin',struct('type','wall','temperature',311,'velocity',[0 0]),'jMax',sphereEdge);
sphereQ=axisymmetricFVConservative(repmat([1 0 0 311],sphereMesh.cellCount,1),variable);
[sphereRate,sphereInfo]=axisymmetricNavierStokesResidual(sphereQ,sphereMesh,variable,sphereBoundary);
[sphereRate5,sphereInfo5]=axisymmetricRANSResidual([sphereQ zeros(sphereMesh.cellCount,1)],sphereMesh,variable,sphereBoundary);
sphereIDs=find(ismember(sphereMesh.faces.boundarySide,[1 2]));
assertAxisGeometry(sphereInfo.axisSymmetry,sphereMesh,sphereIDs);
assert(all(sphereInfo5.faceFluxIntegralPerRadian(sphereIDs,:)==0,'all'));
twoAxisSphere=struct('caseKind','static-sphere-exterior-C-domain-axis-contract-only', ...
    'cellCount',sphereMesh.cellCount,'axisFaceCount',numel(sphereIDs), ...
    'fourQuantityBalance',uniformResidualBudget(sphereRate,sphereInfo,sphereMesh,variable,sphereBoundary), ...
    'fiveQuantityBalance',uniformResidualBudget(sphereRate5,sphereInfo5,sphereMesh,variable,sphereBoundary), ...
    'notHypersonicSphereHeatingValidation',true);

%% 四、完整非均匀五量：独立rho/N梯度和密度修正，不能从保存梯度反算输入
c=mesh.volumeCentroid;t=(1:mesh.cellCount)';
P=[1.1+.025*sin(.7*t),-2+.15*cos(.4*t),.06*c(:,2).*sin(t),440+22*cos(.6*t)];
N=.15+.012*sin(.3*t);Q=axisymmetricFVConservative(P,variable);Q5=[Q P(:,1).*N];
edge=struct('type','prescribedTrace','primitive',[1.1 -2 0 440],'nuTilde',.15);
nonlinearBoundary=struct('iMin',edge,'iMax',edge,'jMin',policy, ...
    'jMax',struct('type','wall','temperature',450,'velocity',[0 0]));
[rate5,info5]=axisymmetricRANSResidual(Q5,mesh,variable,nonlinearBoundary);P5=axisymmetricFVState(Q5(:,1:4),variable);
[rhoSamples,nSamples]=saBoundarySamples(mesh,variable,nonlinearBoundary,P5,N);
[saOracleX,saOracleR]=independentWLS(mesh,[P5(:,1) N],[rhoSamples nSamples],[1 1]);
near(info5.sa.cellDensityGradientX,saOracleX(:,1),3e-10,3e-10,'SA rho X梯度');
near(info5.sa.cellDensityGradientR,saOracleR(:,1),3e-10,3e-10,'SA rho r梯度');
near(info5.sa.cellNuTildeGradientX,saOracleX(:,2),3e-10,3e-10,'SA N X梯度');
near(info5.sa.cellNuTildeGradientR,saOracleR(:,2),3e-10,3e-10,'SA N r梯度');
mu=independentMu(P5(:,4),variable);nu=mu./P5(:,1);sigma=2/3;
expectedDensity=-(nu+N)/sigma.*(saOracleX(:,1).*saOracleX(:,2)+saOracleR(:,1).*saOracleR(:,2));
near(info5.sa.sourceDensityCorrection,expectedDensity,2e-10,2e-10,'完整SA密度梯度修正');
assert(max(abs(expectedDensity))>1e-4,'AxisTest:SAExcitation','密度修正未被激活。');
sourceScale=sum(abs([info5.sa.sourceProduction info5.sa.sourceDestruction info5.sa.sourceGradientSquare expectedDensity]),2);
near(info5.sa.sourceRatePerVolume,info5.sa.sourceProduction-info5.sa.sourceDestruction+ ...
    info5.sa.sourceGradientSquare+expectedDensity,1e-10,1e-10,'SA完整源');
fiveLedger=ledger(rate5,info5,mesh);
scale4=[1 20 20 2000];scale5=[scale4 .5];
jacobian4=serialJacobian(Q,mesh,variable,nonlinearBoundary,scale4);
jacobian5=serialJacobian(Q5,mesh,variable,nonlinearBoundary,scale5);
be4=oneBE(Q,mesh,variable,nonlinearBoundary,scale4);
be5=oneBE(Q5,mesh,variable,nonlinearBoundary,scale5);

%% 五、拒绝坏实际合同；未参与方程的负原始切向外推只保留诊断，不能裁值冒充迹
bad=boundary;bad.jMin.symmetryGradient='unreviewed-normal';
reject('unknown-explicit-axis-policy',@()axisymmetricNavierStokesResidual(evenQ,mesh,gas,bad),{'AxisFV:AxisSymmetryBoundary'});
bad=boundary;bad.iMin.symmetryGradient='mirror-parity';
reject('axis-policy-on-nonaxis-boundary',@()axisymmetricNavierStokesResidual(evenQ,mesh,gas,bad),{'AxisFV:AxisSymmetryBoundary'});
badMesh=mesh;badMesh.faces.centroid(axisIDs(1),2)=1e-6;
reject('zero-area-face-at-finite-radius',@()axisymmetricAxisSymmetry(badMesh,boundary),{'AxisFV:AxisSymmetryGeometry'});
badMesh=mesh;badMesh.volumeCentroid(owners(1),2)=0;
reject('zero-radius-owner-is-not-finite-volume-centroid',@()axisymmetricAxisSymmetry(badMesh,boundary),{'AxisFV:AxisSymmetryGeometry'});
% 所有实际cell/非轴边界均正温，但零面积轴原始外推为负；必须保留负诊断。
[ss,rr]=ndgrid([0 .02 .05],[0 .01 .02]);steepMesh=buildAxisymmetricFVGeometry(ss+10*rr,rr);
steepPrimitive=@(p)[ones(size(p,1),1),zeros(size(p,1),2),-500+1e4*p(:,1)];
steepQ=axisymmetricFVConservative(steepPrimitive(steepMesh.volumeCentroid),gas);steepInput=steepQ;
steepBoundary=allTraceBoundary(steepPrimitive,policy,.1);
[steepRate,steepInfo]=axisymmetricNavierStokesResidual(steepQ,steepMesh,gas,steepBoundary);
steepA=steepInfo.axisSymmetry;assertAxisProxyIdentity(steepA);
assert(all(isfinite(steepRate),'all')&&min(steepA.rawTangentialPrimitiveRhoUxUrT(:,4))<0&& ...
    any(~steepA.rawTangentialStateAdmissible)&&all(steepA.coefficientPrimitiveRhoUxUrT(:,4)>0)&& ...
    isequal(steepQ,steepInput),'AxisTest:Clipping','未用负外推被裁值或污染了实际输入/系数。');
steepP=axisymmetricFVState(steepQ,gas);
near(steepA.coefficientPrimitiveRhoUxUrT(:,4),steepP(steepA.ownerCell,4),0,0,'负原始T外推仍使用owner系数代理');
negativeN=negativeScalarTraceFixture(mesh,variable,policy);
[negativeNRate,negativeNInfo]=axisymmetricRANSResidual(negativeN.Q,mesh,variable,negativeN.boundary);
negativeNA=negativeNInfo.sa.axisSymmetry;assertAxisProxyIdentity(negativeNA);
assert(all(isfinite(negativeNRate),'all')&&min(negativeNA.rawTangentialDensityNuTilde(:,2))<0&& ...
    any(~negativeNA.rawTangentialStateAdmissible)&&all(negativeNA.coefficientDensityNuTilde(:,2)>0));
near(min(negativeNA.rawTangentialDensityNuTilde(:,2)),negativeN.independentMinimumRawTangentialN, ...
    2e-10,2e-10,'独立负原始N外推不裁值');
badQ=evenQ;badQ(1,1)=-1;
reject('actual-negative-cell-density',@()axisymmetricNavierStokesResidual(badQ,mesh,gas,boundary),{'AxisFV:State'});
badN=negativeN.Q;badN(1,5)=-1;
reject('actual-negative-cell-N',@()axisymmetricRANSResidual(badN,mesh,variable,negativeN.boundary),{'AxisSA:State'});
unusedNegativeExtrapolation=struct('rawMinimumTemperature',min(steepA.rawTangentialPrimitiveRhoUxUrT(:,4)), ...
    'rawMinimumNuTilde',min(negativeNA.rawTangentialDensityNuTilde(:,2)), ...
    'rawTangentialStateUsed',false,'coefficientStateIsPhysicalFaceSolution',false, ...
    'positiveCellStatesUnchanged',true,'rawValuesNotClipped',true);

%% 六、冻结v67真实稳态默认完整info重放；研究身份与生产资格保持分开
legacySHA=fileSHA(legacyEvidenceFile);frozen=load(legacyEvidenceFile,'implicit','mesh','gas','boundary');
[frozenRate,frozenInfo]=axisymmetricRANSResidual(frozen.implicit.Q,frozen.mesh,frozen.gas,frozen.boundary);
assert(isequaln(frozenInfo,frozen.implicit.finalInfo)&& ...
    max(abs(frozenRate)./frozen.implicit.residualScale,[],'all')==frozen.implicit.residualNorm&& ...
    strcmp(legacySHA,fileSHA(legacyEvidenceFile)),'AxisTest:Legacy','新轴策略破坏v67默认完整重放或改了证据。');
report=struct('status','passed-explicit-skew-axis-complete-NS-SA-contract', ...
    'cellCount',mesh.cellCount,'minimumAxisSkewRatio',min(ratio),'evenQuadrature',evenQuadrature, ...
    'oddQuadrature',oddQuadrature,'independentMirrorGeometry',oracleGeometry,'evenLedger',evenLedger, ...
    'oddLedger',oddLedger,'uniformBalance',uniformBalance,'uniformBalance5',uniformBalance5, ...
    'twoAxisSphere',twoAxisSphere, ...
    'sourceDensityCorrectionExcited',true,'sourceDensityCorrectionScale',sourceScale, ...
    'fiveLedger',fiveLedger,'jacobian4',jacobian4,'jacobian5',jacobian5,'be4',be4,'be5',be5, ...
    'refusals',refusals,'unusedNegativeExtrapolation',unusedNegativeExtrapolation, ...
    'legacyEvidenceSHA256',legacySHA,'legacyCompleteInfoReplayed',true, ...
    'trueFVMeansUsedOnlyInExplicitlyControlledFields',true,'generalEOSProxyNotCalledExactMean',true, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'skew_axis_evidence.mat'),'report','gas','variable','mesh','boundary','evenQ','evenRate','evenInfo', ...
    'oddQ','oddBoundary','oddRate','oddInfo','uniformQ','uniformBoundary','uniformRate','uniformInfo', ...
    'uniformRate5','uniformInfo5','sphereMesh','sphereBoundary','sphereQ','sphereRate','sphereInfo','sphereRate5','sphereInfo5', ...
    'Q','Q5','N','nonlinearBoundary','rate5','info5','saOracleX','saOracleR','expectedDensity','negativeN', ...
    'steepMesh','steepQ','steepBoundary','steepRate','steepInfo','negativeNRate','negativeNInfo');
fid=fopen(fullfile(folder,'skew_axis_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Skew axis passed: real FV parity, independent WLS/EOS owner proxies/raw diagnostics/SA density source, 4/5 Jacobian/BE, %d refusals and complete v67 replay.\n',numel(refusals));
    function reject(label,action,expected)
        try,action();catch exception
            assert(ismember(exception.identifier,expected),'AxisTest:UnexpectedRefusal', ...
                '%s期望%s，实际%s。',label,strjoin(expected,','),exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('AxisTest:MissingRefusal','错误合同%s未被拒绝。',label);
    end
end

function boundary=allTraceBoundary(primitive,policy,N)
edge=struct('type','prescribedTrace','primitive',@(p)primitive(p),'nuTilde',N);
boundary=struct('iMin',edge,'iMax',edge,'jMin',policy,'jMax',edge);
end
function samples=nsBoundarySamples(mesh,gas,boundary,P)
bf=mesh.boundaryFaces;f=mesh.faces;names={'iMin','iMax','jMin','jMax'};samples=P(f.owner(bf),[2 3 4 5]);
for j=1:numel(bf)
    k=bf(j);bc=boundary.(names{f.boundarySide(k)});
    if strcmp(bc.type,'axis'),continue;end
    if strcmp(bc.type,'wall')
        T=bc.temperature;if isa(T,'function_handle'),T=T(f.centroid(k,:));end
        u=bc.velocity;if isa(u,'function_handle'),u=u(f.centroid(k,:));end
        samples(j,1:3)=[u(:)' T];
    else
        p=bc.primitive;if isa(p,'function_handle'),p=p(f.centroid(k,:));end
        samples(j,:)=[p(2:4) p(1)*gas.R*p(4)];
    end
end
end
function [rhoSamples,nSamples]=saBoundarySamples(mesh,gas,boundary,P,N)
bf=mesh.boundaryFaces;f=mesh.faces;names={'iMin','iMax','jMin','jMax'};
rhoSamples=P(f.owner(bf),1);nSamples=N(f.owner(bf));
for j=1:numel(bf)
    k=bf(j);bc=boundary.(names{f.boundarySide(k)});
    if strcmp(bc.type,'axis'),continue;end
    if strcmp(bc.type,'wall')
        rhoSamples(j)=P(f.owner(k),5)/(gas.R*bc.temperature);nSamples(j)=0;
    else
        p=bc.primitive;if isa(p,'function_handle'),p=p(f.centroid(k,:));end
        rhoSamples(j)=p(1);
        % 本受控测试全边界ux<0、ur=0，iMax为质量入流，iMin出流。
        if dot(p(2:3),f.normal(k,:))<0
            z=bc.nuTilde;if isa(z,'function_handle'),z=z(f.centroid(k,:));end;nSamples(j)=z;
        end
    end
end
end
function [gx,gr,record]=independentWLS(mesh,cellValues,boundaryValues,parity)
% 独立按几何邻接形成小正规方程，不读取任何生产gradient矩阵或保存梯度。
nc=mesh.cellCount;nv=size(cellValues,2);f=mesh.faces;bf=mesh.boundaryFaces;
axis=f.areaPerRadian(bf)==0&f.centroid(bf,2)==0;axisFaces=bf(axis);owner=f.owner(axisFaces);
points=f.centroid(bf,:);points(axis,:)=[mesh.volumeCentroid(owner,1),-mesh.volumeCentroid(owner,2)];
values=boundaryValues;values(axis,:)=cellValues(owner,:).*parity;
gx=zeros(nc,nv);gr=gx;determinants=zeros(nc,1);
for c=1:nc
    ids=find(f.owner==c|f.neighbour==c);M=zeros(2);rhs=zeros(2,nv);
    for id=ids'
        if f.neighbour(id)>0
            if f.owner(id)==c,other=f.neighbour(id);shift=f.neighbourShift(id,:);else,other=f.owner(id);shift=-f.neighbourShift(id,:);end
            delta=mesh.volumeCentroid(other,:)+shift-mesh.volumeCentroid(c,:);value=cellValues(other,:);
        else
            index=find(bf==id);delta=points(index,:)-mesh.volumeCentroid(c,:);value=values(index,:);
        end
        weight=1/sum(delta.^2);M=M+weight*(delta'*delta);rhs=rhs+weight*delta'*(value-cellValues(c,:));
    end
    determinants(c)=det(M);g=M\rhs;gx(c,:)=g(1,:);gr(c,:)=g(2,:);
end
record=struct('method','independent per-cell 2x2 normal equations with physical neighbours and mirrored parity samples', ...
    'minimumNormalMatrixDeterminant',min(determinants),'productionGradientMatrixUsedAsInput',false);
end
function assertAxisGeometry(a,mesh,ids)
assert(isequal(a.faceID,ids)&&isequal(a.ownerCell,mesh.faces.owner(ids)));
near(a.physicalPointXR,mesh.faces.centroid(ids,:),0,0,'真实轴面位置');
c=mesh.volumeCentroid(mesh.faces.owner(ids),:);near(a.ghostPointXR,[c(:,1),-c(:,2)],0,0,'独立镜像虚采样位置');
assert(all(a.physicalPointXR(:,2)==0)&&all(a.ghostPointXR(:,2)<0));
end
function assertAxisProxyIdentity(a)
assert(strcmp(a.coefficientStateKind,'owner-even-projection-at-zero-measure-axis-not-a-solved-face-trace')&& ...
    ~a.coefficientStateIsPhysicalFaceSolution&&~a.physicalTraceReconstructed&&~a.rawTangentialStateUsed&& ...
    ~a.stateClipped&&~a.physicalFaceGeometryChanged,'AxisTest:AxisProxy','零面积轴系数代理被误标为物理面解或裁过状态。');
end
function record=serialJacobian(Q,mesh,gas,boundary,scale)
n=mesh.cellCount;nq=size(Q,2);U=Q./scale;
[J,details,pattern]=axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,struct('differenceStep',2e-6,'differenceScheme','central'));
assert(details.success&&details.columnsDifferentiated==numel(Q));dense=zeros(n*nq);h=2e-6*max(1,abs(U(:)));
for k=1:numel(Q)
    plus=U;minus=U;plus(k)=plus(k)+h(k);minus(k)=minus(k)-h(k);
    rp=spatial(plus.*scale,mesh,gas,boundary)./scale;rm=spatial(minus.*scale,mesh,gas,boundary)./scale;
    dense(:,k)=(rp(:)-rm(:))/(2*h(k));
end
relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');assert(relative<2e-6,'AxisTest:Jacobian','完整有色与独立逐列Jacobian不符。');
outside=~full(kron(ones(nq),pattern.cellPattern));outsideMaximum=max([0;abs(dense(outside))]);
assert(outsideMaximum<1e-6,'AxisTest:JacobianStencil','镜像轴模板依赖超出既有保守稀疏支持。');
v=reshape(sin((1:numel(Q))'*.63),size(Q));epsilon=7e-7;
[rp,ip]=spatial((U+epsilon*v).*scale,mesh,gas,boundary);[rm,im]=spatial((U-epsilon*v).*scale,mesh,gas,boundary);
action=reshape(J*v(:),size(Q));direct=(rp-rm)./scale/(2*epsilon);
direction=max(abs(action-direct),[],'all')/max(abs(direct),[],'all');assert(direction<1e-5);
integral=sum(action.*scale.*mesh.volumePerRadian,1);expected=(-(ip.netBoundaryFlux-im.netBoundaryFlux)+(ip.sourceIntegral-im.sourceIntegral))/(2*epsilon);
arithmetic=sum(abs(action.*scale.*mesh.volumePerRadian),1)+abs(expected);
assert(all(abs(integral-expected)<=1e-5*arithmetic+512*eps*arithmetic));
record=struct('quantityCount',nq,'coloredVsSerialRelative',relative,'directionRelative',direction, ...
    'outsideOriginalStencilMaximum',outsideMaximum,'allColumnsChecked',numel(Q),'integratedDerivativeDefect',integral-expected);
end
function record=oneBE(Q,mesh,gas,boundary,scale)
opts=struct('maxSteps',1,'initialPseudoCFL',10,'maximumPseudoCFL',10,'maximumPseudoStep',2e-8, ...
    'residualScale',scale,'maxNewtonIterations',18,'jacobianRefreshEvery',1);
result=advanceAxisymmetricSteady(Q,mesh,gas,boundary,opts);
assert(result.acceptedSteps==1&&~result.steadyConverged&&result.history.dt(end)==opts.maximumPseudoStep);
[rate,info]=spatial(result.Q,mesh,gas,boundary);dt=result.history.dt(end);
defect=max(abs((result.Q-Q-dt*rate)./scale),[],'all');
assert(defect<=result.history.newtonTarget(end)+512*eps*max(abs(Q./scale),[],'all')&&isequaln(info,result.finalInfo));
ledgerRaw=sum((result.Q-Q).*mesh.volumePerRadian,1)+dt*(info.netBoundaryFlux-info.sourceIntegral);
relative=abs(ledgerRaw)./(sum(mesh.volumePerRadian)*scale);assert(all(relative<=result.options.conservationTolerance));
assert(any(result.Q~=Q,'all')&&all(isfinite(result.Q),'all'));
record=struct('quantityCount',size(Q,2),'acceptedSteps',result.acceptedSteps,'steadyConverged',false, ...
    'actualDt',dt,'BEStateDefect',defect,'actualNewtonTarget',result.history.newtonTarget(end), ...
    'fiveOrFourQuantityLedger',ledgerRaw,'ledgerRelative',relative,'solverResidualEvaluations',result.residualEvaluations, ...
    'stateWasClipped',false,'physicalTimeAdvanced',0);
end
function [rate,info]=spatial(Q,mesh,gas,boundary)
if size(Q,2)==5,[rate,info]=axisymmetricRANSResidual(Q,mesh,gas,boundary);else,[rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);end
end
function record=uniformResidualBudget(rate,info,mesh,gas,boundary)
% 常量温度的稀疏G*phi会发生抵消，不能再乘eps把微小导热当成无形成误差。
% 第一项只覆盖已形成的通量/源运算；第二项从有单位的WLS样本前向传播。
% 本专用校核仅用于这里的均匀轴向/静止场和普通墙温合同，不是任意状态误差估计。
P=info.primitive;f=mesh.faces;nc=mesh.cellCount;nf=numel(f.owner);bf=mesh.boundaryFaces;
assert(all(P==P(1,:),'all')&&P(1,3)==0,'AxisTest:UniformFixture','预算只能用于明确均匀且ur=0的受控状态。');
names={'iMin','iMax','jMin','jMax'};trace=P(f.owner(bf),:);
for b=1:numel(bf)
    k=bf(b);bc=boundary.(names{f.boundarySide(k)});
    switch bc.type
        case 'axis',trace(b,3)=0;
        case 'wall'
            T=bc.temperature;if isa(T,'function_handle'),T=T(f.centroid(k,:));end
            u=bc.velocity;if isa(u,'function_handle'),u=u(f.centroid(k,:));end
            trace(b,2:4)=[u(:)' T];trace(b,1)=trace(b,5)/(gas.R*T);
        case 'prescribedTrace'
            v=bc.primitive;if isa(v,'function_handle'),v=v(f.centroid(k,:));end
            trace(b,:)=axisymmetricFVState(axisymmetricFVConservative(v,gas),gas);
        otherwise,error('AxisTest:UniformFixture','均匀校核未建立此边界的形成误差合同。');
    end
end
trace=axisymmetricFVState(axisymmetricFVConservative(trace(:,1:4),gas),gas);
a=axisymmetricAxisSymmetry(mesh,boundary);assert(a.active);
sample=[P(:,[2 3 4 5]);trace(:,[2 3 4 5])];
sample(nc+a.boundaryIndex,:)=[P(a.ownerCell,2),-P(a.ownerCell,3),P(a.ownerCell,4:5)];
nominal=P(1,[2 3 4 5]);
% 直接差分正规方程对完全相同的常量样本严格零，独立确认模板没有物理常量缺陷。
[zeroX,zeroR]=independentWLS(mesh,repmat(nominal,nc,1),repmat(nominal,numel(bf),1),[1 -1 1 1]);
assert(all(zeroX==0,'all')&&all(zeroR==0,'all'),'AxisTest:UniformStencil','独立差分WLS未保持常量。');
% 每格d个邻接：稀疏对角汇总至多d-1次、乘加至多2(d+1)次；
% gamma_(4d+8)保守包围这两段运算。边界EOS往返的实际ulp差单列为输入差。
degree=full(sum(mesh.cellFaceMultiplicity,2));gamma=@(m)(m*eps)./(1-m*eps);g=gamma(4*degree+8);
deviation=abs(sample-nominal);
cellBX=abs(a.gradientX)*deviation+g.*(abs(a.gradientX)*abs(sample));
cellBR=abs(a.gradientR)*deviation+g.*(abs(a.gradientR)*abs(sample));
gx=a.gradientX*sample;gr=a.gradientR*sample;
assert(all(abs(gx)<=cellBX+realmin,'all')&&all(abs(gr)<=cellBR+realmin,'all'), ...
    'AxisTest:UniformStencil','常量WLS形成误差超过按样本数计算的前向界。');
near(info.cellVelocityGradientX,gx(:,1:2),0,0,'均匀预算真实速度X梯度');
near(info.cellVelocityGradientR,gr(:,1:2),0,0,'均匀预算真实速度r梯度');
inside=f.neighbour>0;owner=f.owner;other=f.neighbour;
BX=cellBX(owner,1:3);BR=cellBR(owner,1:3);
BX(inside,:)=.5*(cellBX(owner(inside),1:3)+cellBX(other(inside),1:3));
BR(inside,:)=.5*(cellBR(owner(inside),1:3)+cellBR(other(inside),1:3));
left=P(owner,2:4);right=left;right(inside,:)=P(other(inside),2:4);right(bf,:)=trace(:,2:4);
% 真实非正交位移约束的三项差、平方、除法和回加逐步有界；不加无量纲绝对阈值。
d=f.delta;correction=(abs(right-left)+abs(d(:,1)).*BX+abs(d(:,2)).*BR)./sum(d.^2,2);
BX=(BX+abs(d(:,1)).*correction)/(1-gamma(12));
BR=(BR+abs(d(:,2)).*correction)/(1-gamma(12));
af=a.faceID;ao=a.ownerCell;
BX(af,:)=[cellBX(ao,1),zeros(numel(ao),1),cellBX(ao,3)];BR(af,:)=0;
% ur样本严格零，因而圆周速度和ur/r也严格零；轴只保留偶量切向导数。
assert(all(gx(:,2)==0)&all(gr(:,2)==0));
mu=info.transport.faceMu;k=info.transport.faceK;
divB=BX(:,1)+BR(:,2);
txx=2*mu.*(BX(:,1)+divB/3);trr=2*mu.*(BR(:,2)+divB/3);txr=mu.*(BR(:,1)+BX(:,2));
n=abs(f.normal);traction=[txx.*n(:,1)+txr.*n(:,2),txr.*n(:,1)+trr.*n(:,2)];
velocity=abs(P(owner,2:3));
dl=f.centroid(inside,:)-mesh.volumeCentroid(owner(inside),:);
dr=f.centroid(inside,:)-f.neighbourShift(inside,:)-mesh.volumeCentroid(other(inside),:);
velocity(inside,:)=.5*(abs(P(owner(inside),2:3))+cellBX(owner(inside),1:2).*abs(dl(:,1))+ ...
    cellBR(owner(inside),1:2).*abs(dl(:,2))+abs(P(other(inside),2:3))+ ...
    cellBX(other(inside),1:2).*abs(dr(:,1))+cellBR(other(inside),1:2).*abs(dr(:,2)));
velocity(bf,:)=abs(trace(:,2:3));
heat=k.*(BX(:,3).*n(:,1)+BR(:,3).*n(:,2));
faceFormation=[zeros(nf,1),traction,sum(velocity.*traction,2)+heat].*f.areaPerRadian/(1-gamma(24));
faceFormation(af,:)=0;
delta=abs(mesh.planarCentroid-mesh.volumeCentroid);
sourceHoop=2*info.transport.sourceMu.*(cellBX(:,1)+cellBR(:,2))/3;
sourceFormation=zeros(nc,4);sourceFormation(:,3)= ...
    (cellBX(:,4).*delta(:,1)+cellBR(:,4).*delta(:,2)+sourceHoop).*mesh.planarArea/(1-gamma(16));
if size(rate,2)==5,faceFormation(:,5)=0;sourceFormation(:,5)=0;end
arithmetic=abs(mesh.incidence)*abs(info.faceFluxIntegralPerRadian)+abs(info.cellSourceIntegral);
arithmeticBudget=1024*eps*arithmetic./mesh.volumePerRadian;
formationBudget=(abs(mesh.incidence)*faceFormation+sourceFormation)./mesh.volumePerRadian;
budget=arithmeticBudget+formationBudget;
assert(all(abs(rate)<=budget+realmin,'all'),'AxisTest:Uniform','均匀场残差超过通量/源运算及WLS本构形成预算。');
record=struct('maximumAbsoluteResidual',max(abs(rate),[],1),'maximumRoundoffBudget',max(budget,[],1), ...
    'maximumFluxSourceArithmeticBudget',max(arithmeticBudget,[],1), ...
    'maximumWLSConstitutiveFormationBudget',max(formationBudget,[],1), ...
    'maximumCellGradientFormationX',max(cellBX,[],1),'maximumCellGradientFormationR',max(cellBR,[],1), ...
    'directDifferenceConstantWLSExactlyZero',true,'gradientGammaOperationCount',4*degree+8, ...
    'formationMethod','sample-count gamma; actual boundary EOS ulp differences; WLS to nonorthogonal stress/Fourier/source and true area/incidence/volume', ...
    'ledger',ledger(rate,info,mesh),'zeroNIsNotTurbulenceSteadyValidation',true);
end
function record=ledger(rate,info,mesh)
terms=rate.*mesh.volumePerRadian;source=info.cellSourceIntegral;
boundary=info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:);
defect=sum(terms,1)+sum(boundary,1)-sum(source,1);
arithmetic=sum(abs(terms),1)+sum(abs(source),1)+sum(abs(info.faceFluxIntegralPerRadian).*full(sum(mesh.cellFaceMultiplicity,1))',1);
budget=2048*eps*arithmetic+realmin;assert(all(abs(defect)<=budget),'AxisTest:Ledger','共享面/径向源账不闭合。');
record=struct('quantityCount',size(rate,2),'rawDefect',defect,'roundoffBudget',budget,'unitSeparated',true,'perRadian',true);
end
function [Q,record]=trueMeans(mesh,gas,primitive)
% 两个真三角形上的Duffy积分；直接积分rhoE，不对T或速度均值假做EOS。
x=([-.8611363115940526;-.3399810435848563;.3399810435848563;.8611363115940526]+1)/2;
w=[.3478548451374539;.6521451548625461;.6521451548625461;.3478548451374539]/2;
Q=zeros(mesh.cellCount,4);volume=zeros(mesh.cellCount,1);
for c=1:mesh.cellCount
    vertices=mesh.vertices(mesh.cells(c,:),:);integral=zeros(1,4);
    for tri=[1 2 3;1 3 4]'
        A=vertices(tri(1),:);B=vertices(tri(2),:);C=vertices(tri(3),:);jac=abs(det([B-A;C-A]));
        for i=1:4
            for j=1:4
                point=A+x(i)*(B-A)+(1-x(i))*x(j)*(C-A);weight=w(i)*w(j)*jac*(1-x(i))*point(2);
                p=primitive(point);q=[p(1),p(1)*p(2:3),p(1)*(gas.cv*p(4)+.5*sum(p(2:3).^2))];
                integral=integral+weight*q;volume(c)=volume(c)+weight;
            end
        end
    end
    Q(c,:)=integral/volume(c);
end
relative=max(abs(volume-mesh.volumePerRadian)./mesh.volumePerRadian);assert(relative<5e-13);
record=struct('kind','true-r-weighted-conservative-FV-means','method','independent two-triangle 4x4 Duffy/Gauss', ...
    'volumeRelative',relative,'primitivePointMeansNotAssumed',true);
end
function fixture=negativeScalarTraceFixture(mesh,gas,policy)
% 有限个正N状态找独立未用外推为负者，验证保留诊断且代理未假装求解值。
rho=ones(mesh.cellCount,1);P=[rho,-2*rho,zeros(size(rho)),350*rho];
edge=struct('type','prescribedTrace','primitive',[1 -2 0 350],'nuTilde',1e-5);
boundary=struct('iMin',edge,'iMax',edge,'jMin',policy,'jMax',struct('type','wall','temperature',350,'velocity',[-2 0]));
ids=find(mesh.faces.boundarySide==3);owner=mesh.faces.owner(ids);dx=mesh.faces.centroid(ids,1)-mesh.volumeCentroid(owner,1);
found=false;
for k=1:mesh.cellCount
    N=1e-5+zeros(mesh.cellCount,1);N(k)=.3;
    [rhoSamples,nSamples]=saBoundarySamples(mesh,gas,boundary,axisymmetricFVState(axisymmetricFVConservative(P,gas),gas),N);
    [gx,~]=independentWLS(mesh,[rho N],[rhoSamples nSamples],[1 1]);physical=N(owner)+gx(owner,2).*dx;
    if min(physical)<-1e-6,found=true;break;end
end
assert(found,'AxisTest:NegativeTraceExcitation','没有构造出正cell负raw外推，不可跳过代理语义测试。');
fixture=struct('Q',[axisymmetricFVConservative(P,gas) rho.*N],'boundary',boundary, ...
    'independentMinimumRawTangentialN',min(physical),'allCellNPositive',all(N>0),'chosenCell',k);
end
function mu=independentMu(T,gas)
t=gas.transport;mu=t.referenceViscosity*(T/t.referenceTemperature).^1.5.* ...
    ((t.referenceTemperature+t.sutherlandTemperature)./(T+t.sutherlandTemperature));
end
function near(a,b,relative,absolute,label)
budget=absolute+relative*max(abs(a),abs(b));
assert(all(abs(a-b)<=budget,'all'),'AxisTest:Value','%s不匹配，最大绝对差%.9g。',label,max(abs(a-b),[],'all'));
end
function digest=fileSHA(path)
f=fopen(path,'rb');assert(f>=0);guard=onCleanup(@()fclose(f));bytes=fread(f,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
