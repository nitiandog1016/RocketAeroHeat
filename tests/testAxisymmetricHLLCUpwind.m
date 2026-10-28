function report=testAxisymmetricHLLCUpwind(folder,evidenceFile)
%TESTAXISYMMETRICHLLCUPWIND 检查内部HLLC只对实际取用的混合波扇校验星态。
% 模块一：两单元正反全迎风强膨胀反例，未使用负星压不能阻断FL/FR。
% 模块二：独立固定真实TMR face52的两状态/法向小fixture，不依赖整场MAT。
% 模块三：正常混合波扇仍与原HLLC公式一致，实际取用的非正星态仍须拒绝。
% 模块四：本地或包内存在完整诊断证据时重放实际1000格失败Q1；缺失时明确未执行。
% 只验空间通量选择与保护，不是稳态、实验或实际箭体气动热验收。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','hllc_upwind');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v63_validation','tmr_explicit_hllc_diagnostic','tmr_explicit_hllc_evidence.mat');
    if ~exist(evidenceFile,'file')
        evidenceFile=fullfile(root,'validation','v63','tmr_explicit_hllc_diagnostic','tmr_explicit_hllc_evidence.mat');
    end
    if ~exist(evidenceFile,'file'),evidenceFile='';end
elseif ~isempty(evidenceFile)
    assert(exist(evidenceFile,'file')==2,'AxisFVTest:UpwindSource','明确给出的完整诊断来源不存在。');
end
if ~exist(folder,'dir'),mkdir(folder);end
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.01,'Pr',.72,'temperatureRange',[80 2000]));
[X,R]=ndgrid([0 .01 .02],[.03 .04]);mesh=buildAxisymmetricFVGeometry(X,R);
interior=find(mesh.faces.neighbour>0);assert(numel(interior)==1);face=interior(1);
assert(norm(mesh.faces.normal(face,:)-[1 0])<1e-14);
cases={ [1 40 0 300;1 200 0 300], [1 -200 0 300;1 -40 0 300], ...
    [1 1 .2 300;.125 -1 -.1 240], [.125 1 -.1 240;1 -1 .2 300] };
labels={'strong-all-rightward','strong-all-leftward','regular-left-star','regular-right-star'};
records=struct([]);caseQ=cell(1,4);caseInfo=cell(1,4);
for j=1:4
    primitive=cases{j};Q=axisymmetricFVConservative(primitive,gas);
    boundary=twoCellBoundary(primitive,.01,0,.03,false);
    [rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);
    P=axisymmetricFVState(Q,gas);
    [expected,math]=independentHLLC(Q(1,:),Q(2,:),P(1,:),P(2,:),[1 0]);
    flux=info.faceConvectiveFluxPerRadian(face,:);expected=expected*mesh.faces.areaPerRadian(face);
    relative=max(abs(flux-expected)./max(1,abs(expected)));
    assert(all(isfinite(rate),'all')&&relative<2e-12&&~info.physicalSolution&&~info.eligibleForSurfaceCoverage);
    if j==1
        assert(math.SL>0&&strcmp(math.actualFluxBranch,'FL')&&min(math.starPressure)<0, ...
            'AxisFVTest:UpwindFixture','强正向fixture未形成未用负星压与全向右波扇。');
    elseif j==2
        assert(math.SR<0&&strcmp(math.actualFluxBranch,'FR')&&min(math.starPressure)<0, ...
            'AxisFVTest:UpwindFixture','强反向fixture未形成未用负星压与全向左波扇。');
    elseif j==3
        assert(strcmp(math.actualFluxBranch,'left-star')&&all(math.starPressure>0)&&all(math.starInternalEnergy>0));
    else
        assert(strcmp(math.actualFluxBranch,'right-star')&&all(math.starPressure>0)&&all(math.starInternalEnergy>0));
    end
    item=struct('label',labels{j},'primitive',primitive,'mathematics',math,'convectiveFluxRelative',relative, ...
        'convectiveFluxIntegralPerRadian',flux,'expectedFluxIntegralPerRadian',expected);
    if isempty(records),records=item;else,records(end+1)=item;end %#ok<AGROW>
    caseQ{j}=Q;caseInfo{j}=info;
end

%% 二、实际TMR失败面的小fixture：保留两侧原始量、原斜法向及冻结气体
actualGas=axisymmetricFVGas(struct('gamma',1.4,'R',287.05,'mu',5.4794210031583635e-6,'Pr',.72,'temperatureRange',[40 2000]));
actualPrimitive=[.0046865086013313725 -527.09081364643737 10.216176274773794 531.97785012655959; ...
    .073936908680531405 -1069.4962400387171 385.02503986881430 294.24631939958880];
originalNormal=[-.9953794982811034 .09601903145553357];skew=originalNormal(2)/(-originalNormal(1));
[Xa,Ra]=ndgrid([.02 .01 0],[.03 .04]);actualMesh=buildAxisymmetricFVGeometry(Xa+skew*(Ra-.03),Ra);
actualFace=find(actualMesh.faces.neighbour>0);assert(numel(actualFace)==1);
assert(norm(actualMesh.faces.normal(actualFace,:)-originalNormal)<2e-14);
actualQ=axisymmetricFVConservative(actualPrimitive,actualGas);
actualBoundary=twoCellBoundary(actualPrimitive,.01,skew,.03,true);
[actualRate,actualInfo]=axisymmetricNavierStokesResidual(actualQ,actualMesh,actualGas,actualBoundary);
actualP=axisymmetricFVState(actualQ,actualGas);
[actualExpected,actualMath]=independentHLLC(actualQ(1,:),actualQ(2,:),actualP(1,:),actualP(2,:),originalNormal);
actualExpected=actualExpected*actualMesh.faces.areaPerRadian(actualFace);
actualRelative=max(abs(actualInfo.faceConvectiveFluxPerRadian(actualFace,:)-actualExpected)./max(1,abs(actualExpected)));
assert(actualRelative<2e-12&&all(isfinite(actualRate),'all')&&actualMath.SL>0&& ...
    min(actualMath.starPressure)<0&&strcmp(actualMath.actualFluxBranch,'FL'), ...
    'AxisFVTest:UpwindActualFixture','真实失败面的小fixture未保留原正向FL和未用负星压。');

%% 三、实际取用星态的负压混合波扇不能借迎风修复删除原正性保护
badPrimitive=[1 -200 0 300;1 200 0 300];badQ=axisymmetricFVConservative(badPrimitive,gas);
badP=axisymmetricFVState(badQ,gas);[~,badMath]=independentHLLC(badQ(1,:),badQ(2,:),badP(1,:),badP(2,:),[1 0]);
assert(badMath.SL<0&&badMath.SR>0&&min(badMath.starPressure)<0&& ...
    ismember(badMath.actualFluxBranch,{'left-star','right-star'}));
badBoundary=twoCellBoundary(badPrimitive,.01,0,.03,false);refused=false;refusalIdentifier='';
try
    axisymmetricNavierStokesResidual(badQ,mesh,gas,badBoundary);
catch exception
    refusalIdentifier=exception.identifier;
    assert(strcmp(refusalIdentifier,'AxisFV:HLLC'),'AxisFVTest:UpwindUnexpected','混合波扇出现非HLLC拒绝：%s。',refusalIdentifier);
    refused=true;
end
assert(refused,'AxisFVTest:UpwindGuard','实际混合波扇的负星压未被拒绝。');

%% 四、可选完整实际Q1空间重放：真实来源存在才声明该项已运行
actualReplay=struct('performed',false,'reason','No complete development MAT required by self-contained two-cell fixtures', ...
    'sourceEvidenceSHA256','','cellCount',[],'sourceStep',[],'sourceAttempt',[],'sourceFace',52, ...
    'sourceStage','','newOperatorSucceeded',false,'convectiveFluxRelative',[],'originalResidual',[]);
fullReplayQ=[];fullReplayInfo=[];fullReplayRate=[];
if ~isempty(evidenceFile)
    full=load(evidenceFile,'report','mesh','gas','boundary','scale','failureStageQ','failureDiagnostics');
    assert(full.mesh.cellCount==1000&&full.report.acceptedSteps==220&&full.report.failureStep==221&& ...
        strcmp(full.report.failureStage,'first-stage-spatial')&&strcmp(full.report.terminalIdentifier,'AxisFV:HLLC'), ...
        'AxisFVTest:UpwindSource','完整重放必须使用实际显式221步Q1失败记录。');
    fullReplayQ=full.failureStageQ;fullGas=axisymmetricFVGas(full.gas);
    [fullReplayRate,fullReplayInfo]=axisymmetricNavierStokesResidual(fullReplayQ,full.mesh,fullGas,full.boundary);
    fullP=axisymmetricFVState(fullReplayQ,fullGas);target=52;owner=full.mesh.faces.owner(target);neighbour=full.mesh.faces.neighbour(target);
    [expected,fullMath]=independentHLLC(fullReplayQ(owner,:),fullReplayQ(neighbour,:),fullP(owner,:),fullP(neighbour,:),full.mesh.faces.normal(target,:));
    expected=expected*full.mesh.faces.areaPerRadian(target);
    relative=max(abs(fullReplayInfo.faceConvectiveFluxPerRadian(target,:)-expected)./max(1,abs(expected)));
    assert(all(isfinite(fullReplayRate),'all')&&relative<2e-12&&fullMath.SL>0&&min(fullMath.starPressure)<0&& ...
        strcmp(fullMath.actualFluxBranch,'FL')&&~fullReplayInfo.physicalSolution&&~fullReplayInfo.eligibleForSurfaceCoverage);
    actualReplay=struct('performed',true,'reason','Replayed full recorded actual first RK stage through current spatial operator', ...
        'sourceEvidenceSHA256',fileSHA(evidenceFile),'cellCount',1000,'sourceStep',full.report.failureStep, ...
        'sourceAttempt',full.report.failureAttempt,'sourceFace',target,'sourceStage',full.report.failureStage, ...
        'newOperatorSucceeded',true,'convectiveFluxRelative',relative,'originalResidual',max(abs(fullReplayRate)./full.scale,[],'all'));
end
report=struct('status','passed-upwind-HLLC-flux-selection-and-retained-mixed-state-guard', ...
    'records',records,'actualFaceFixture',struct('originalNormal',originalNormal,'primitive',actualPrimitive, ...
        'mathematics',actualMath,'convectiveFluxRelative',actualRelative), ...
    'mixedNegativePressureFixture',badMath,'mixedStateRefused',refused,'refusalIdentifier',refusalIdentifier, ...
    'completeActualStageReplay',actualReplay,'sourceBindings',sourceBindings(root), ...
    'scope','Internal face flux selection and fixed actual-stage regression; no added HLLE fallback or positivity clipping', ...
    'steadyConverged',false,'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'hllc_upwind_evidence.mat'),'report','mesh','gas','caseQ','caseInfo','actualMesh','actualGas','actualQ','actualBoundary', ...
    'actualRate','actualInfo','badPrimitive','badMath','fullReplayQ','fullReplayInfo','fullReplayRate');
fid=fopen(fullfile(folder,'hllc_upwind_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Axisymmetric HLLC upwind passed: both raw directions, actual face52 fixture, both normal star branches, retained mixed refusal; full recorded stage replay=%d.\n',actualReplay.performed);
end

function boundary=twoCellBoundary(primitive,cutoff,skew,r0,descending)
% 所有开放面给出与相邻单元相同的明确储库，只让唯一内部面测试HLLC。
trace=@(point)piecewisePrimitive(point,primitive,cutoff,skew,r0,descending);
far=struct('type','riemannReservoir','primitive',trace);
boundary=struct('iMin',far,'iMax',far,'jMin',far,'jMax',far);
end

function p=piecewisePrimitive(point,primitive,cutoff,skew,r0,descending)
logicalX=point(1)-skew*(point(2)-r0);
if (descending&&logicalX>cutoff)||(~descending&&logicalX<cutoff),p=primitive(1,:);else,p=primitive(2,:);end
end

function [F,math]=independentHLLC(UL,UR,L,R,n)
% 完整保留计算出的未用星态供审计；参考通量只按真实波区选择FL/FR/对应星态。
unL=dot(L(2:3),n);unR=dot(R(2:3),n);SL=min(unL-L(6),unR-R(6));SR=max(unL+L(6),unR+R(6));
DL=L(1)*(SL-unL);DR=R(1)*(SR-unR);SM=(R(5)-L(5)+DL*unL-DR*unR)/(DL-DR);
rL=L(1)*(SL-unL)/(SL-SM);rR=R(1)*(SR-unR)/(SR-SM);
pL=L(5)+DL*(SM-unL);pR=R(5)+DR*(SM-unR);
USL=[rL,rL*(L(2:3)+(SM-unL)*n),rL*(UL(4)/L(1)+(SM-unL)*(SM+L(5)/DL))];
USR=[rR,rR*(R(2:3)+(SM-unR)*n),rR*(UR(4)/R(1)+(SM-unR)*(SM+R(5)/DR))];
eL=USL(4)-sum(USL(2:3).^2)/(2*rL);eR=USR(4)-sum(USR(2:3).^2)/(2*rR);
FL=[L(1)*unL,UL(2:3)*unL+L(5)*n,(UL(4)+L(5))*unL];
FR=[R(1)*unR,UR(2:3)*unR+R(5)*n,(UR(4)+R(5))*unR];
if SL>=0,F=FL;branch='FL';elseif SR<=0,F=FR;branch='FR';
elseif SM>=0,F=FL+SL*(USL-UL);branch='left-star';else,F=FR+SR*(USR-UR);branch='right-star';end
math=struct('SL',SL,'SM',SM,'SR',SR,'starDensity',[rL rR],'starPressure',[pL pR], ...
    'starInternalEnergy',[eL eR],'actualFluxBranch',branch,'starConservativeLeft',USL,'starConservativeRight',USR);
end

function bindings=sourceBindings(root)
paths={'regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m','regionalFlow/axisymmetric/axisymmetricFVConservative.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVState.m','regionalFlow/axisymmetric/buildAxisymmetricFVGeometry.m','tests/testAxisymmetricHLLCUpwind.m'};
bindings=struct([]);
for j=1:numel(paths)
    item=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
end
end

function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
