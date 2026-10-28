function report=testCutCellWlsDirection(folder)
%TESTCUTCELLWLSDIRECTION 特征边界、完整WLS黏性及全部五量方向的实际原NS校核。
% 模块一：覆盖超声速入/出流、亚声速出流、回流、压力出口与零面速度策略。
% 模块二：两个非仿射方向各作原NS中心/半步差分，核对梯度、面通量及完整散度。
% 模块三：显式核对功速度链式项、冻结边界差别和缓存拒绝，非物理精度验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'WlsDirectionTest:Overwrite','不得覆盖方向验证。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','newtonunit03','steady_newton_unit_evidence.mat'),'result','cut','gas','boundary');
Q=s.result.initialQ;baseInfo=s.result.initialInfo;cut=s.cut;gas=s.gas;boundary=s.boundary;scale=s.result.residualScale;clear s;
boundaryChecks=boundaryTest(gas);state=cutCellBoundaryStencilState(Q,cut,gas,boundary);cache=buildCutCellGradientOperator(cut,state);
n=size(Q,1);directions={reshape(sin((1:5*n)*.37),n,5),reshape(cos((1:5*n)*.51),n,5)};
measurements=struct([]);allEvidence=cell(2,1);calls=0;
for index=1:2
    directionU=directions{index};bd=cutCellBoundaryDirectionalState(Q,directionU,cut,gas,boundary,scale,state);
    [viscous,vd]=cutCellWlsViscousDirection(Q,directionU,baseInfo,cut,gas,boundary,scale,cache,state,bd);
    [euler,ed]=cutCellSelectedEulerDirection(Q,directionU,baseInfo,cut,gas,boundary,scale,state,bd);analytic=euler+viscous;
    frozen=cutCellWlsViscousDirection(Q,directionU,baseInfo,cut,gas,boundary,scale,cache,state,[]);
    halfResults=cell(2,1);errors=zeros(2,4);
    for refinement=1:2
        h=1e-6/max(abs(directionU),[],'all')/2^(refinement-1);
        [rp,ip]=cutCellNavierStokesResidual(Q+h*directionU.*scale,cut,gas,boundary);
        [rm,im]=cutCellNavierStokesResidual(Q-h*directionU.*scale,cut,gas,boundary);calls=calls+2;
        full=(rp-rm)./scale/(2*h);g=(ip.gradient-im.gradient)/(2*h);
        fv=(ip.viscousFaceFlux-im.viscousFaceFlux)/(2*h);fc=(ip.convectiveFaceFlux-im.convectiveFaceFlux)/(2*h);
        errors(refinement,:)=[relative(analytic,full),relative(vd.gradientDerivative,g),relative(vd.faceViscousFluxDerivative,fv),relative(ed.faceEulerFluxDerivative,fc)];
        assert(max(errors(refinement,:))<2e-6,'WlsDirectionTest:FD','实际NS方向/梯度/面通量没有达到明确差分标准。');
        halfResults{refinement}=struct('h',h,'ratePlus',rp,'rateMinus',rm,'gradientPlus',ip.gradient,'gradientMinus',im.gradient, ...
            'faceViscousPlus',ip.viscousFaceFlux,'faceViscousMinus',im.viscousFaceFlux,'faceEulerPlus',ip.convectiveFaceFlux,'faceEulerMinus',im.convectiveFaceFlux);
    end
    row=struct('fullNSRelativeDifference',errors(2,1),'gradientRelativeDifference',errors(2,2), ...
        'viscousFluxRelativeDifference',errors(2,3),'eulerFluxRelativeDifference',errors(2,4), ...
        'workVelocityChainNorm',norm(vd.workVelocityU(:)),'frozenBoundaryViscousDifference',norm(viscous(:)-frozen(:)));
    assert(row.workVelocityChainNorm>0&&row.frozenBoundaryViscousDifference>0);
    if isempty(measurements),measurements=row;else,measurements(end+1)=row;end %#ok<AGROW>
    allEvidence{index}=struct('directionU',directionU,'boundaryDirection',bd,'viscous',viscous,'viscousDetail',vd,'euler',euler,'eulerDetail',ed,'refinements',{halfResults});
end
bad=cache;bad.geometry.cells.volume(1)=bad.geometry.cells.volume(1)*1.01;rejected=false;
try,cutCellWlsViscousDirection(Q,directions{1},baseInfo,cut,gas,boundary,scale,bad,state,[]);catch ex,rejected=strcmp(ex.identifier,'CutDirection:Cache');end
assert(rejected,'WlsDirectionTest:Cache','不同几何缓存必须拒绝。');
report=struct('status','passed-characteristic-boundary-full-WLS-flux-and-complete-NS-direction-small-domain', ...
    'cellCount',n,'boundaryChecks',boundaryChecks,'directions',measurements,'actualOriginalNSEvaluations',calls, ...
    'cacheMismatchRejected',true,'oldNSChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'wls_direction_unit_evidence.mat'),'report','Q','baseInfo','cut','gas','boundary','scale','state','allEvidence','-v7');
fid=fopen(fullfile(folder,'wls_direction_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function error=relative(a,b)
error=max(abs(a-b),[],'all')/max(1,max(abs(b),[],'all'));
end
function report=boundaryTest(gas)
reservoir=[.8 -20 30 10 320];normal=[1 0 0];base=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',reservoir);
P=[1 800 20 -30 300;1 -800 20 -30 300;1 100 20 -30 300;1 -100 20 -30 300;1 100 20 -30 300;1 -100 20 -30 300;1 0 20 -30 300];
worst=0;branches=cell(7,1);pressure=0;zero=0;
for i=1:size(P,1)
    b=base;direction=[.08 30 -10 20 12];
    if i==2,b.reservoirPrimitive=[.8 -900 30 10 320];end
    if i==4||i==6,b.reservoirPrimitive=[.8 -300 30 10 320];end
    if i==5||i==6,b.mode='pressureOutlet';b.pressure=.9*gas.R*300;end
    if i==7,b.reservoirPrimitive=P(i,:);direction=[.1 0 2 -3 0];end
    [actual,tag]=cutCellCharacteristicBoundaryDirection(P(i,:),direction,b,normal,gas);h=1e-6;
    [plus,tp]=resolveNSCharacteristicBoundary3D(P(i,:)+h*direction,b,normal,gas);
    [minus,tm]=resolveNSCharacteristicBoundary3D(P(i,:)-h*direction,b,normal,gas);
    assert(strcmp(tp.branch,tag.branch)&&strcmp(tm.branch,tag.branch),'WlsDirectionTest:Branch','边界测试跨过了未声明的分支。');
    error=max(abs(actual-(plus-minus)/(2*h)))/max(1,max(abs(actual)));worst=max(worst,error);branches{i}=tag.branch;
    pressure=pressure+tag.pressureInputApplied;zero=zero+strcmp(tag.zeroTangentialPolicy,'interior-entropy-and-tangential-velocity');
end
assert(worst<2e-7&&numel(unique(branches))==4&&pressure>0&&zero==1);
report=struct('cases',7,'maximumRelativeDifference',worst,'pressureInputCases',pressure,'zeroSpeedPolicyCases',zero,'branches',{branches});
end
