function report=testAxisymmetricVariableTransport(folder,legacyEvidenceFile)
%TESTAXISYMMETRICVARIABLETRANSPORT Sutherland变输运、墙热/应力及旧常物性逐值校核。
% 模块一：独立冻结温点、数组形状及非法温度/参数拒绝，不读取外部CFD解。
% 模块二：线性非等温场核对Tw处墙输运、内部线性迹、平面求积环向应力和稳定频率。
% 模块三：无源非等温轴向Couette的连续积分参考，真实四量FV均值作三档空间/墙热比较。
% 模块四：v63真实TMR常物性检查点完整info/rate逐值重放，防止变输运改坏旧分支。
% 本测试没有长时间推进；解析/数值校核不构成真实空气、湍流、实验或生产面积认证。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v64_validation','variable_transport');end
if nargin<2
    legacyEvidenceFile=fullfile(root,'outputs','v63_validation','tmr_hlle_continue1','tmr_implicit_budget.mat');
    if ~exist(legacyEvidenceFile,'file'),legacyEvidenceFile=fullfile(root,'validation','v63','tmr_hlle_continue1','tmr_implicit_budget.mat');end
end
assert(exist(legacyEvidenceFile,'file')==2,'AxisFVTest:TransportLegacy','缺少冻结v63常输运终态，不能跳过兼容性重放。');
if ~exist(folder,'dir'),mkdir(folder);end
rawGas=struct('gamma',1.4,'R',287.05,'mu',1.716e-5,'Pr',.72,'temperatureRange',[40 2000]);
rawGas.transport=struct('model','sutherland','referenceViscosity',1.716e-5,'referenceTemperature',273.15,'sutherlandTemperature',110.4);
gas=axisymmetricFVGas(rawGas);

%% 一、独立温点数表和形状；参考黏度与gas.mu身份不能各自漂移
temperatures=[80 150 273.15 300 600 1200 1800];
expectedMu=[5.47905297834979e-6 1.0285688602850156e-5 1.716e-5 1.8459162511975804e-5 ...
    3.0162086848898431e-5 4.6249330249785523e-5 5.8280317831058107e-5];
properties=axisymmetricFVTransport(temperatures,gas);
assert(isstruct(properties)&&all(isfield(properties,{'mu','k'}))&&isequal(size(properties.mu),size(temperatures))&& ...
    isequal(size(properties.k),size(temperatures))&&max(abs(properties.mu./expectedMu-1))<3e-14&& ...
    max(abs(properties.k./(expectedMu*gas.cp/gas.Pr)-1))<3e-14, ...
    'AxisFVTest:TransportTable','Sutherland独立温点、导热率或数组形状不匹配。');
column=axisymmetricFVTransport(temperatures',gas);assert(isequal(column.mu,properties.mu')&&isequal(column.k,properties.k'));
matrixTemperature=reshape([80 300 600 1200],2,2);matrix=axisymmetricFVTransport(matrixTemperature,gas);
assert(isequal(size(matrix.mu),[2 2])&&max(abs(matrix.mu-independentMu(matrixTemperature,gas)),[],'all')<1e-18);
constantGas=axisymmetricFVGas(rmfield(rawGas,'transport'));constant=axisymmetricFVTransport(temperatures,constantGas);
assert(isequal(constant.mu,repmat(constantGas.mu,size(temperatures)))&& ...
    isequal(constant.k,repmat(constantGas.k,size(temperatures))),'AxisFVTest:TransportConstant','旧常物性不应随T改变。');
refusals=struct('label',{},'identifier',{});
for value={39,2001,NaN,Inf,-1,complex(300,1),[]}
    reject('invalid-temperature',@()axisymmetricFVTransport(value{1},gas),'AxisFV:TransportTemperature');
end
bad=rawGas;bad.mu=1.01*bad.mu;reject('mismatched-reference-mu',@()axisymmetricFVGas(bad),'AxisFV:GasTransport');
bad=rawGas;bad.transport.model='unreviewed-model';reject('unknown-model',@()axisymmetricFVGas(bad),'AxisFV:GasTransport');
for key={'referenceViscosity','referenceTemperature','sutherlandTemperature'}
    bad=rawGas;bad.transport=rmfield(bad.transport,key{1});
    reject(['missing-' key{1}],@()axisymmetricFVGas(bad),'AxisFV:GasTransport');
end
bad=rawGas;bad.transport.referenceTemperature=0;reject('nonpositive-reference-temperature',@()axisymmetricFVGas(bad),'AxisFV:GasTransport');
bad=rawGas;bad.transport.referenceViscosity=NaN;reject('nonfinite-reference-viscosity',@()axisymmetricFVGas(bad),'AxisFV:GasTransport');
bad=rawGas;bad.transport.sutherlandTemperature=-110.4;reject('negative-sutherland-temperature',@()axisymmetricFVGas(bad),'AxisFV:GasTransport');

%% 二、线性原始量场使梯度独立可知，检验各位置的真实输运取样而非复刻梯度代码
ri=.03;ro=.06;T0=250;temperatureSlope=8000;u0=-20;shear=200;pressure=8e4;
[X,R]=ndgrid([0 .02 .04],linspace(ri,ro,9));mesh=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));
P=[pressure./(gas.R*(T0+temperatureSlope*mesh.volumeCentroid(:,2))), ...
    u0+shear*mesh.volumeCentroid(:,2),zeros(mesh.cellCount,1),T0+temperatureSlope*mesh.volumeCentroid(:,2)];
Q=axisymmetricFVConservative(P,gas);periodic=struct('type','periodic');
boundary=struct('iMin',periodic,'iMax',periodic, ...
    'jMin',struct('type','wall','temperature',T0+temperatureSlope*ri,'velocity',[u0+shear*ri 0]), ...
    'jMax',struct('type','wall','temperature',T0+temperatureSlope*ro,'velocity',[u0+shear*ro 0]));
[rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);
[constantRate,constantInfo]=axisymmetricNavierStokesResidual(Q,mesh,constantGas,boundary);
w=info.wall;wallMu=independentMu(w.temperatureInput,gas);wallK=wallMu*gas.cp/gas.Pr;
wallExpectedQ=-wallK*temperatureSlope.*w.normal(:,2);wallExpectedTraction=wallMu*shear.*w.normal(:,2);
wallHeatError=max(abs(w.heatFluxIntoWall-wallExpectedQ)./abs(wallExpectedQ));
wallTractionError=max(abs(w.viscousTraction(:,1)-wallExpectedTraction)./abs(wallExpectedTraction));
assert(wallHeatError<2e-10&&wallTractionError<2e-10&& ...
    max(abs(w.viscosity-wallMu)./wallMu)<3e-14&&max(abs(w.thermalConductivity-wallK)./wallK)<3e-14, ...
    'AxisFVTest:TransportWall','墙输运未在真实输入Tw取值，或线性场墙热/剪切不满足独立本构式。');
faceTemperature=T0+temperatureSlope*mesh.faces.centroid(:,2);
sourceTemperature=T0+temperatureSlope*mesh.planarCentroid(:,2);
assert(max(abs(info.transport.faceTemperature-faceTemperature))<2e-10&& ...
    max(abs(info.transport.faceMu-independentMu(faceTemperature,gas)))<2e-16&& ...
    max(abs(info.transport.faceK-independentMu(faceTemperature,gas)*gas.cp/gas.Pr))<2e-13&& ...
    max(abs(info.transport.sourceTemperature-sourceTemperature))<2e-10&& ...
    max(abs(info.transport.sourceMu-independentMu(sourceTemperature,gas)))<2e-16, ...
    'AxisFVTest:TransportSampling','面线性迹或平面求积温度没有使用各自对应输运系数。');
state=axisymmetricFVState(Q,gas);assert(max(abs(info.transport.cellMu-independentMu(state(:,4),gas)))<1e-18);
f=mesh.faces;left=state(f.owner,:);right=left;inside=f.neighbour>0;right(inside,:)=state(f.neighbour(inside),:);
bf=mesh.boundaryFaces;right(bf,4)=w.temperatureInput;
right(bf,1)=right(bf,5)./(gas.R*right(bf,4));
frequencyMu=max([independentMu(faceTemperature,gas),independentMu(left(:,4),gas),independentMu(right(:,4),gas)],[],2);
frequencyK=frequencyMu*gas.cp/gas.Pr;
assert(max(abs(info.transport.frequencyFaceMu-frequencyMu)./frequencyMu)<3e-14&& ...
    max(abs(info.transport.frequencyFaceK-frequencyK)./frequencyK)<3e-14);
nu=max(4*frequencyMu./(3*min(left(:,1),right(:,1))),frequencyK./(gas.cv*min(left(:,1),right(:,1))));
expectedFrequency=2*(mesh.cellFaceMultiplicity*(nu.*f.areaPerRadian./f.normalDistance))./mesh.volumePerRadian+ ...
    4*independentMu(state(:,4),gas)./(state(:,1).*mesh.volumeCentroid(:,2).^2);
frequencyRelative=max(abs(info.diffusiveFrequency-expectedFrequency)./expectedFrequency);
assert(frequencyRelative<3e-13&&all(info.diffusiveFrequency>constantInfo.diffusiveFrequency)&& ...
    isequal(info.convectiveFrequency,constantInfo.convectiveFrequency), ...
    'AxisFVTest:TransportFrequency','扩散频率未包含面/双侧输运上界或单元环向黏性，不能继续使用参考常mu。');
assert(~info.physicalSolution&&~info.eligibleForSurfaceCoverage&&contains(lower(info.temperatureModel),'sutherland'));
linearConservation=assertConservation(rate,info,mesh);

% ur=a*r的散度为2a；独立tau_theta=2mu(T_planar)*a/3，不能误用体质心温度。
radialGradient=3;[Xh,Rh]=ndgrid([0 .02 .04],linspace(ri,ro,7));hoopMesh=buildAxisymmetricFVGeometry(Xh+.13*Rh,Rh);
hp=[pressure./(gas.R*(T0+temperatureSlope*hoopMesh.volumeCentroid(:,2))), ...
    u0+shear*hoopMesh.volumeCentroid(:,2),radialGradient*hoopMesh.volumeCentroid(:,2),T0+temperatureSlope*hoopMesh.volumeCentroid(:,2)];
hq=axisymmetricFVConservative(hp,gas);
trace=@(point)[pressure/(gas.R*(T0+temperatureSlope*point(2))),u0+shear*point(2),radialGradient*point(2),T0+temperatureSlope*point(2)];
bc=struct('type','prescribedTrace','primitive',trace);hoopBoundary=struct('iMin',bc,'iMax',bc,'jMin',bc,'jMax',bc);
[hoopRate,hoopInfo]=axisymmetricNavierStokesResidual(hq,hoopMesh,gas,hoopBoundary);
expectedHoop=2*independentMu(T0+temperatureSlope*hoopMesh.planarCentroid(:,2),gas)*radialGradient/3;
hoopRelative=max(abs(hoopInfo.sourceHoopStress-expectedHoop)./expectedHoop);
assert(hoopRelative<2e-9,'AxisFVTest:TransportHoop','环向源应力未在平面求积温度取mu。');
hoopConservation=assertConservation(hoopRate,hoopInfo,hoopMesh);
% 故意在rate里漏掉全部柱坐标源，边界账保持不变；物理/舍入双预算均必须拒绝。
missingSourceRate=hoopRate-hoopInfo.cellSourceIntegral./hoopMesh.volumePerRadian;
missingSourceControl=assertConservation(missingSourceRate,hoopInfo,hoopMesh,false);
assert(~missingSourceControl.combinedBudgetPassed&&~missingSourceControl.roundoffPassed&& ...
    any(hoopInfo.sourceIntegral~=0),'AxisFVTest:TransportConservationSensitivity','漏柱坐标源负对照被账本尺度掩盖。');

%% 三、无源变输运Couette：连续积分参考与真正保守FV均值，三档仅计算空间及墙热
reference=couetteReference(gas,ri,ro,-100,0,300,800,1e5);levels=[12 24 48];couetteRecords=struct([]);couetteEvidence=cell(1,3);
for j=1:3
    [Xc,Rc]=ndgrid([0 .04],linspace(ri,ro,levels(j)+1));cm=buildAxisymmetricFVGeometry(Xc,Rc,struct('periodicI',true));
    cQ=couetteMeans(cm,reference,gas);cb=struct('iMin',periodic,'iMax',periodic, ...
        'jMin',struct('type','wall','temperature',reference.innerT,'velocity',[reference.innerU 0]), ...
        'jMax',struct('type','wall','temperature',reference.outerT,'velocity',[reference.outerU 0]));
    [cr,ci]=axisymmetricNavierStokesResidual(cQ,cm,gas,cb);cw=ci.wall;
    expectedShear=reference.Cu./cw.centroid(:,2).*cw.normal(:,2);
    expectedHeat=-reference.Cu*(gas.cp*reference.beta/gas.Pr-cw.velocityInput(:,1))./cw.centroid(:,2).*cw.normal(:,2);
    heatRMS=sqrt(mean(((cw.heatFluxIntoWall-expectedHeat)./expectedHeat).^2));
    shearRMS=sqrt(mean(((cw.viscousTraction(:,1)-expectedShear)./expectedShear).^2));
    currentWallTransport=axisymmetricFVTransport(cw.temperatureInput,gas);
    assert(any(expectedHeat>0)&&any(expectedHeat<0)&&all(isfinite(cr),'all')&& ...
        all(cw.viscosity==currentWallTransport.mu));
    couetteConservation=assertConservation(cr,ci,cm);
    item=struct('radialCells',levels(j),'wallHeatRelativeRMS',heatRMS,'wallShearRelativeRMS',shearRMS, ...
        'actualWallHeat',cw.heatFluxIntoWall,'continuousWallHeat',expectedHeat, ...
        'actualWallShear',cw.viscousTraction(:,1),'continuousWallShear',expectedShear, ...
        'sourceFreeContinuousReference',true,'realConservativeFVMeans',true,'timeStepsTaken',0, ...
        'conservationAudit',couetteConservation);
    if isempty(couetteRecords),couetteRecords=item;else,couetteRecords(end+1)=item;end %#ok<AGROW>
    couetteEvidence{j}=struct('mesh',cm,'Q',cQ,'boundary',cb,'rate',cr,'info',ci);
end
heatErrors=[couetteRecords.wallHeatRelativeRMS];shearErrors=[couetteRecords.wallShearRelativeRMS];
assert(all(heatErrors(2:end)<.8*heatErrors(1:end-1))&&all(shearErrors(2:end)<.8*shearErrors(1:end-1))&& ...
    heatErrors(end)<.04&&shearErrors(end)<.04,'AxisFVTest:TransportCouette','变输运非等温Couette墙热/剪切未向无源连续参考收敛。');

%% 四、冻结v63常物性真实终态逐值重放；不重新跑109步或改写既有检查点
legacySHA=fileSHA(legacyEvidenceFile);legacy=load(legacyEvidenceFile,'implicit','mesh','gas','boundary');
assert(~isfield(legacy.gas,'transport')&&legacy.implicit.steadyConverged&&legacy.implicit.residualNorm<=2e-6);
[legacyRate,legacyInfo]=axisymmetricNavierStokesResidual(legacy.implicit.Q,legacy.mesh,legacy.gas,legacy.boundary);
assert(isequaln(legacyInfo,legacy.implicit.finalInfo)&& ...
    isequal(legacyRate,(-legacy.mesh.incidence*legacy.implicit.finalInfo.faceFluxIntegralPerRadian+ ...
    legacy.implicit.finalInfo.cellSourceIntegral)./legacy.mesh.volumePerRadian)&& ...
    max(abs(legacyRate)./legacy.implicit.residualScale,[],'all')==legacy.implicit.residualNorm&& ...
    strcmp(legacySHA,fileSHA(legacyEvidenceFile)), ...
    'AxisFVTest:TransportLegacyReplay','旧常输运分支不能完整逐值重放，不可用新模型标注覆盖旧info。');
report=struct('status','passed-Sutherland-local-spatial-couette-and-constant-legacy-replay', ...
    'temperatureTable',temperatures,'independentViscosityTable',expectedMu,'refusals',refusals, ...
    'wallHeatRelative',wallHeatError,'wallTractionRelative',wallTractionError, ...
    'diffusiveFrequencyRelative',frequencyRelative,'hoopStressRelative',hoopRelative,'couetteRecords',couetteRecords, ...
    'linearConservationAudit',linearConservation,'hoopConservationAudit',hoopConservation, ...
    'omittedAxisSourceNegativeControl',missingSourceControl, ...
    'legacyEvidenceSHA256',legacySHA,'legacySteadyResidual',legacy.implicit.residualNorm, ...
    'legacyCompleteInfoReplay',true,'sourceBindings',sourceBindings(root),'timeStepsTaken',0, ...
    'couetteValidationKind','Source-free continuous integral reference, spatial/FV-mean wall heat and shear refinement; not an advanced NS steady solution', ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'variable_transport_evidence.mat'),'report','gas','mesh','Q','boundary','rate','info','constantRate','constantInfo', ...
    'hoopMesh','hq','hoopBoundary','hoopRate','hoopInfo','couetteEvidence','legacyRate','legacyInfo');
fid=fopen(fullfile(folder,'variable_transport_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Variable transport passed: Sutherland table, %d refusals, wall/source/frequency, source-free nonisothermal Couette spatial reference, complete old constant replay; no time/production claim.\n',numel(refusals));
    function reject(label,action,expectedID)
        try,action();catch exception
            assert(strcmp(exception.identifier,expectedID),'AxisFVTest:TransportUnexpectedRefusal','%s期望%s，实际%s。',label,expectedID,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('AxisFVTest:TransportMissingRefusal','非法合同%s未拒绝。',label);
    end
end

function mu=independentMu(T,gas)
% 独立Sutherland连续参考；不调用被测transport模块。
t=gas.transport;coefficient=t.referenceViscosity*(t.referenceTemperature+t.sutherlandTemperature)/t.referenceTemperature^1.5;
mu=coefficient*T.^1.5./(T+t.sutherlandTemperature);
end
function reference=couetteReference(gas,ri,ro,Ui,Uo,Ti,To,pressure)
% 连续无源方程r*mu*u'=Cu及r*(k*T'+mu*u*u')=CE导出T(u)二次式。
B=gas.Pr/(2*gas.cp);beta=(To-Ti+B*(Uo^2-Ui^2))/(Uo-Ui);alpha=Ti-beta*Ui+B*Ui^2;
T=@(u)alpha+beta*u-B*u.^2;
Cu=integral(@(u)independentMu(T(u),gas),Ui,Uo,'RelTol',1e-12,'AbsTol',1e-15)/log(ro/ri);
solution=ode45(@(~,u)Cu/independentMu(T(u),gas),[log(ri) log(ro)],Ui,odeset('RelTol',2e-12,'AbsTol',2e-12));
assert(abs(deval(solution,log(ro))-Uo)<1e-8,'AxisFVTest:TransportReference','独立无源Couette积分未到达外墙速度。');
reference=struct('ri',ri,'ro',ro,'innerU',Ui,'outerU',Uo,'innerT',Ti,'outerT',To,'pressure',pressure, ...
    'Cu',Cu,'beta',beta,'alpha',alpha,'B',B,'solution',solution);
end
function Q=couetteMeans(mesh,ref,gas)
% 对真实rho/rhou/rhoE作16点Gauss体积积分，不把EOS点值当控制体保守均值。
n=16;j=(1:n-1)';b=j./sqrt(4*j.^2-1);[vectors,D]=eig(diag(b,1)+diag(b,-1));[g,order]=sort(diag(D));weights=2*vectors(1,order)'.^2;
Q=zeros(mesh.cellCount,4);rNodes=mesh.R(1,:);
for c=1:mesh.cellCount
    lo=rNodes(c);hi=rNodes(c+1);r=(lo+hi)/2+(hi-lo)*g/2;u=deval(ref.solution,log(r))';
    T=ref.alpha+ref.beta*u-ref.B*u.^2;rho=ref.pressure./(gas.R*T);
    values=[rho,rho.*u,zeros(n,1),ref.pressure/(gas.gamma-1)+rho.*u.^2/2];
    Q(c,:)=sum(values.*(weights.*r),1)*(hi-lo)/(hi^2-lo^2);
end
end
function audit=assertConservation(rate,info,mesh,enforce)
% 逐分量有单位账：物理相对预算仍2e-12，另明确给出共享面消去的浮点求和界。
% 近零净质量/动量的相对缺陷可大，但只能在独立舍入界内接受；不引入跨单位常数1。
if nargin<4,enforce=true;end
contribution=rate.*mesh.volumePerRadian;bf=info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:);S=info.cellSourceIntegral;
rateIntegral=sum(contribution,1);boundaryIntegral=sum(bf,1);sourceIntegral=sum(S,1);
defect=rateIntegral+boundaryIntegral-sourceIntegral;
rateAbsolute=sum(abs(contribution),1);boundaryAbsolute=sum(abs(bf),1);sourceAbsolute=sum(abs(S),1);
scale=rateAbsolute+boundaryAbsolute+sourceAbsolute;
mask=scale>0;relative=zeros(1,4);relative(mask)=abs(defect(mask))./scale(mask);relative(~mask&defect~=0)=Inf;
% B记录面被控制体使用的次数：共享面两次，包含owner==neighbour的周期共享面。
% 相比先压缩的|C|，B不会把周期自配对的+1/-1在计算预算前消去。
multiplicity=full(sum(mesh.cellFaceMultiplicity,1))';
allFaceAbsolute=sum(abs(info.faceFluxIntegralPerRadian).*multiplicity,1);
maximumCellFaces=max(full(sum(mesh.cellFaceMultiplicity,2)));
k=2*mesh.cellCount+numel(mesh.boundaryFaces)+maximumCellFaces+12;
spacing=eps('double');assert(k*spacing<1,'AxisFVTest:TransportRoundoff','浮点求和界k过大。');
gammaK=k*spacing/(1-k*spacing);allAbsolute=scale+allFaceAbsolute;
roundoffBound=gammaK*allAbsolute;physicalRelativeTolerance=2e-12;
physicalBudget=physicalRelativeTolerance*scale;combinedBudget=physicalBudget+roundoffBound;
combinedPassed=all(abs(defect)<=combinedBudget);roundoffPassed=all(abs(defect)<=roundoffBound);
reportedDifference=info.conservationDefect-defect;reportedAgrees=all(abs(reportedDifference)<=roundoffBound);
audit=struct('quantityOrder',{{'mass','axialMomentum','radialMomentum','totalEnergy'}}, ...
    'quantityUnits',{{'kg/(s rad)','N/rad','N/rad','W/rad'}},'rawDefect',defect, ...
    'reportedRawDefect',info.conservationDefect,'reportedMinusRecomputed',reportedDifference, ...
    'rateVolumeIntegral',rateIntegral,'boundaryIntegral',boundaryIntegral,'sourceIntegral',sourceIntegral, ...
    'absoluteRateVolumeSum',rateAbsolute,'absoluteBoundaryFluxSum',boundaryAbsolute,'absoluteCellSourceSum',sourceAbsolute, ...
    'physicalSummationScale',scale,'zeroPhysicalScaleComponents',~mask,'rawRelativeDefect',relative, ...
    'physicalRelativeTolerance',physicalRelativeTolerance,'physicalRelativeBudget',physicalBudget, ...
    'physicalRelativeOnlyPassed',all(relative<=physicalRelativeTolerance), ...
    'allFaceMultiplicityAbsolute',allFaceAbsolute,'allAbsoluteRoundoffContributions',allAbsolute, ...
    'maximumCellFaceMultiplicity',maximumCellFaces,'roundoffOperationBoundK',k,'machineSpacing',spacing,'gammaK',gammaK, ...
    'absoluteRoundoffBound',roundoffBound,'combinedAbsoluteBudget',combinedBudget, ...
    'roundoffDominatedComponents',abs(defect)>physicalBudget&abs(defect)<=roundoffBound, ...
    'combinedBudgetPassed',combinedPassed,'roundoffPassed',roundoffPassed,'reportedDefectAgrees',reportedAgrees);
if enforce
    assert(combinedPassed&&roundoffPassed&&reportedAgrees,'AxisFVTest:TransportConservation', ...
        '变输运四量账不满足原2e-12物理相对预算加显式舍入界，或超过纯浮点求和界。');
end
end
function bindings=sourceBindings(root)
paths={'regionalFlow/axisymmetric/axisymmetricFVTransport.m','regionalFlow/axisymmetric/axisymmetricFVGas.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVTemperatureTrace.m', ...
    'regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m','tests/testAxisymmetricVariableTransport.m'};bindings=struct([]);
for j=1:numel(paths)
    item=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
