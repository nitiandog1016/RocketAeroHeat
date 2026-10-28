function report=studyTmrAxisymmetricSteady(folder,implicitSteps,explicitSteps,boundaryMode)
%STUDYTMRAXISYMMETRICSTEADY 公开1000格同初态、同残差、固定工作量下的真实对照。
% 模块一：只读取v62自算初态/原网格合同，不读取外部CFD流场。
% 模块二：隐式仍用原2e-6稳态阈值；另列预先固定的1e-3初残差工作目标。
% 模块三：显式RK2与隐式全部残差调用分别计数，同目标都达到时才给调用量比。
% 模块四：常物性层流与公开湍流实验不等价，工作上限不算稳态或精度通过。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','tmr_steady');end
if nargin<2,implicitSteps=60;end
if nargin<3,explicitSteps=12000;end
if nargin<4,boundaryMode='legacyCharacteristic';end
assert(ismember(boundaryMode,{'legacyCharacteristic','riemannReservoir'}),'TmrSteady:Boundary','须明确选择既有特征边界或HLLE储库边界。');
validateattributes(implicitSteps,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(explicitSteps,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
source=fullfile(root,'outputs','v62_validation','tmr_startup','tmr_startup_evidence.mat');
if ~exist(source,'file'),source=fullfile(root,'validation','v62','tmr_startup','tmr_startup_evidence.mat');end
frozen=load(source,'result','mesh','gas','boundary','report');
assert(frozen.mesh.cellCount==1000&&~frozen.report.steadyConverged&&~frozen.report.physicalSolution);
mesh=frozen.mesh;gas=frozen.gas;boundary=frozen.boundary;
boundaryIdentity=frozen.report.sourceBoundaryIdentity;
if strcmp(boundaryMode,'riemannReservoir')
    names={'iMin','iMax','jMax'};
    for j=1:numel(names)
        old=boundary.(names{j});state=old.reservoirPrimitive;
        assert(isnumeric(state)&&numel(state)==5&&state(4)==0,'TmrSteady:Swirl','储库原始量合同不匹配。');
        boundary.(names{j})=struct('type','riemannReservoir','primitive',state([1 2 3 5]));
    end
    boundaryIdentity='Uniform public M/p0/T0 prescribed exterior states; HLLE finite-volume reservoir flux on all open faces; not exact WIND/NMF boundary reproduction';
end
Q=frozen.result.initialQ;scale=frozen.result.residualScale;
[rate,initialInfo]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);
initialResidual=max(abs(rate)./scale,[],'all');workTarget=1e-3*initialResidual;
if ~exist(folder,'dir'),mkdir(folder);end
opts=struct('maxSteps',implicitSteps,'initialPseudoCFL',20,'maximumPseudoCFL',1e7, ...
    'pseudoCFLGrowth',2,'residualScale',scale,'differenceScheme','forward', ...
    'steadyResidualTolerance',2e-6,'steadyRequiredSteps',3,'minimumSteadySteps',3);
clock=tic;implicit=advanceAxisymmetricSteady(Q,mesh,gas,boundary,opts);implicitSeconds=toc(clock);
save(fullfile(folder,'tmr_implicit_budget.mat'),'implicit','mesh','gas','boundary','Q','scale','initialResidual','workTarget','implicitSeconds');
explicit=[];explicitSeconds=[];explicitCalls=[];
if explicitSteps>0
    eopts=struct('maxSteps',explicitSteps,'convectiveCFL',.25,'diffusiveCFL',.25, ...
        'residualScale',scale,'steadyResidualTolerance',2e-6,'steadyRequiredSteps',3,'minimumSteadySteps',3);
    clock=tic;explicit=advanceAxisymmetricNavierStokes(Q,mesh,gas,boundary,eopts);explicitSeconds=toc(clock);
    explicitCalls=explicit.residualEvaluations;
    save(fullfile(folder,'tmr_explicit_budget.mat'),'explicit','explicitSeconds','explicitCalls');
end
indexI=find(implicit.history.residual<=workTarget,1);indexE=[];ratio=[];callsI=[];callsE=[];
if ~isempty(indexI),callsI=implicit.history.residualEvaluations(indexI);end
if ~isempty(explicit)
    indexE=find(explicit.history.residual<=workTarget,1);
    if ~isempty(indexE),callsE=explicit.history.residualEvaluations(indexE);end
end
if ~isempty(callsI)&&~isempty(callsE),ratio=callsE/callsI;end
[r,info]=axisymmetricNavierStokesResidual(implicit.Q,mesh,gas,boundary);
assert(isequaln(info,implicit.finalInfo)&&max(abs(r)./scale,[],'all')==implicit.residualNorm);
state=axisymmetricFVState(implicit.Q,gas);wall=info.wall;
report=struct('status','completed-budgeted-numerical-comparison-not-experimental-validation', ...
    'cellCount',1000,'initialOriginalResidual',initialResidual,'fixedResidualWorkTarget',workTarget, ...
    'requiredSteadyTolerance',2e-6,'implicitAcceptedSteps',implicit.acceptedSteps, ...
    'implicitRejectedSteps',implicit.rejectedAttempts,'implicitOriginalResidual',implicit.residualNorm, ...
    'implicitSteadyConverged',implicit.steadyConverged,'implicitStopReason',implicit.stopReason, ...
    'implicitTotalResidualEvaluations',implicit.residualEvaluations,'implicitSeconds',implicitSeconds, ...
    'implicitCallsAtWorkTarget',callsI,'explicitCallsAtWorkTarget',callsE, ...
    'sameTargetResidualEvaluationRatio',ratio,'explicitTotalResidualEvaluations',explicitCalls, ...
    'explicitSeconds',explicitSeconds,'timingsAreNotAtSameStoppingTarget',true, ...
    'implicitMinimumTemperature',min(state(:,4)),'implicitMaximumTemperature',max(state(:,4)), ...
    'signedWallHeatFluxRangeWm2',[min(wall.heatFluxIntoWall),max(wall.heatFluxIntoWall)], ...
    'BEGlobalBalanceRelative',implicit.globalBalanceRelative, ...
    'sourceEvidenceSHA256',fileSHA(source),'sourceBoundaryIdentity',boundaryIdentity,'boundaryMode',boundaryMode, ...
    'transportIdentity',frozen.report.transportIdentity,'pseudoTimeIsPhysical',false, ...
    'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
if ~isempty(explicit)
    report.explicitOriginalResidual=explicit.residualNorm;report.explicitSteadyConverged=explicit.steadyConverged;
end
save(fullfile(folder,'tmr_steady_comparison.mat'),'report','implicit','explicit','mesh','gas','boundary','initialInfo');
fid=fopen(fullfile(folder,'tmr_steady_comparison.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('TMR budgeted comparison: original %.8g -> implicit %.8g, steady=%d; work target %.8g; not experiment/production.\n', ...
    initialResidual,implicit.residualNorm,implicit.steadyConverged,workTarget);
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
