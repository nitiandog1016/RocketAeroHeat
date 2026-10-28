function report=testAxisymmetricSteadyAnnular(folder,evidenceFile)
%TESTAXISYMMETRICSTEADYANNULAR 原8/16格环隙的同初态/同原残差稳态与工作量对照。
% 模块一：复用v62实际RK2终态及独立连续参考，禁止重新定义稳态门槛。
% 模块二：隐式从原静止初态实际求解，核对Q、p、墙q/tau与闭域质量族。
% 模块三：统计所有Jacobian/拒绝/线搜索的完整残差调用，并分开报告历史耗时。
% 本例是人工气体无轴环隙，未认证空气、激波、真实箭体或任何新增生产面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','annular');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v62_validation','advance','axisymmetric_advance_evidence.mat');
    if ~exist(evidenceFile,'file'),evidenceFile=fullfile(root,'validation','v62','advance','axisymmetric_advance_evidence.mat');end
end
saved=load(evidenceFile,'report','solutions','references','gas','boundary');
assert(strcmp(saved.report.status,'passed-real-annular-NS-time-advance-and-analytic-steady-check')&& ...
    ~saved.report.physicalSolution&&~saved.report.experimentalValidation&&numel(saved.solutions)==2);
for j=1:numel(saved.report.sourceBindings)
    item=saved.report.sourceBindings(j);
    % v63残差新增明确的Riemann储库合同；旧环隙必须逐值重放旧算子输出，见下方。
    if ~ismember(item.path,{'regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m', ...
            'regionalFlow/axisymmetric/advanceAxisymmetricNavierStokes.m'})
        assert(strcmp(fileSHA(fullfile(root,item.path)),item.sha256),'v62其他对照源码改变。');
    end
end
if ~exist(folder,'dir'),mkdir(folder);end
cases=struct([]);solutions=cell(1,2);sourceSHA=fileSHA(evidenceFile);
for j=1:2
    baseline=saved.solutions{j}.result;mesh=saved.solutions{j}.mesh;exact=saved.references{j};
    [~,baselineReplay]=axisymmetricNavierStokesResidual(baseline.Q,mesh,saved.gas,saved.boundary);
    assert(isequaln(baselineReplay,baseline.finalInfo),'旧边界合同下的原始空间算子没有逐值重放。');
    assert(baseline.steadyConverged&&baseline.options.steadyResidualTolerance==2e-6&&baseline.options.steadyRequiredSteps==10);
    options=struct('maxSteps',80,'initialPseudoCFL',20,'maximumPseudoCFL',1e7,'pseudoCFLGrowth',2, ...
        'steadyResidualTolerance',baseline.options.steadyResidualTolerance, ...
        'steadyRequiredSteps',baseline.options.steadyRequiredSteps,'minimumSteadySteps',10, ...
        'residualScale',baseline.residualScale,'differenceScheme','central');
    clock=tic;result=advanceAxisymmetricSteady(baseline.initialQ,mesh,saved.gas,saved.boundary,options);seconds=toc(clock);
    save(fullfile(folder,sprintf('annular_n%02d_run.mat',mesh.cellCount)),'result','mesh','exact','seconds');
    assert(result.steadyConverged&&strcmp(result.stopReason,'steady-residual-converged'),'实际隐式环隙没有达到原稳态门槛。');
    assert(all(result.history.residual(end-9:end)<=2e-6)&&result.globalBalanceRelative<=optionsDefault('conservationTolerance'));
    scale=baseline.residualScale;V=mesh.volumePerRadian;
    qDifference=max(abs(result.Q-baseline.Q)./scale,[],'all');
    assert(qDifference<5e-7,'隐式改变了原空间离散的稳态。');
    [r,info]=axisymmetricNavierStokesResidual(result.Q,mesh,saved.gas,saved.boundary);
    assert(isequaln(info,result.finalInfo)&&max(abs(r)./scale,[],'all')<=2e-6);
    wall=info.wall;oldWall=baseline.finalInfo.wall;
    wallHeatDifference=max(abs(wall.heatFluxIntoWall./oldWall.heatFluxIntoWall-1));
    shearDifference=max(abs(wall.viscousTraction(:,1)./oldWall.viscousTraction(:,1)-1));
    assert(wallHeatDifference<3e-6&&shearDifference<3e-6,'墙热或牵引不对应同一稳态。');
    inner=wall.normal(:,2)<0;expected=zeros(size(inner));
    expected(inner)=exact.innerHeatFlux;expected(~inner)=exact.outerHeatFlux;
    qError=sqrt(mean((wall.heatFluxIntoWall./expected-1).^2));
    power=sum(wall.mechanicalPowerIntoFluidPerRadian);heat=sum(wall.heatFluxIntoWall.*wall.areaPerRadian);
    closure=abs(heat-power)/abs(power);
    assert(power>0&&all(wall.heatFluxIntoWall>0)&&closure<1e-5);
    mass=baseline.initialConservedIntegral(1);massDrift=abs(result.finalConservedIntegral(1)-mass)/mass;
    assert(massDrift<5e-10&&result.closedMassConstraintActive&&~result.pseudoTimeIsPhysical&&~isfield(result,'time'));
    % 冻结v62原RK2无拒绝时每步r1/rnew两次残差，初态一次；EOS检查不计。
    assert(baseline.rejectedAttempts==0);
    explicitCalls=1+2*baseline.acceptedSteps;
    assert(result.residualEvaluations<explicitCalls,'本环隙未取得残差调用量收益。');
    item=struct('radialCells',mesh.cellCount,'implicitAcceptedSteps',result.acceptedSteps, ...
        'implicitRejectedSteps',result.rejectedAttempts,'originalResidual',result.residualNorm, ...
        'implicitResidualEvaluations',result.residualEvaluations,'explicitResidualEvaluations',explicitCalls, ...
        'residualEvaluationRatio',explicitCalls/result.residualEvaluations,'implicitSeconds',seconds, ...
        'historicalExplicitSeconds',saved.report.cases(j).seconds,'timingsAreNotSameProcess',true, ...
        'QDifferenceScaledMaximum',qDifference,'wallHeatDifferenceRelative',wallHeatDifference, ...
        'wallShearDifferenceRelative',shearDifference,'analyticWallHeatRelativeRMS',qError, ...
        'massRelativeDrift',massDrift,'wallPowerRelativeClosure',closure, ...
        'BEGlobalBalanceRelative',result.globalBalanceRelative,'jacobianBuilds',result.jacobianBuilds, ...
        'linearSolves',result.linearSolves,'factorizations',result.factorizations);
    if isempty(cases),cases=item;else,cases(end+1)=item;end %#ok<AGROW>
    solutions{j}=struct('result',result,'mesh',mesh);
end
report=struct('status','passed-same-annular-steady-original-residual-and-work-count', ...
    'cases',cases,'v62EvidenceSHA256',sourceSHA,'testGasIsAir',false,'axisIncluded',false, ...
    'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'annular_steady_evidence.mat'),'report','solutions','saved');
fid=fopen(fullfile(folder,'annular_steady_summary.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Implicit annular same-steady passed: calls %d/%d versus RK2 %d/%d, no production coverage.\n', ...
    cases(1).implicitResidualEvaluations,cases(2).implicitResidualEvaluations,cases(1).explicitResidualEvaluations,cases(2).explicitResidualEvaluations);
end
function value=optionsDefault(name)
assert(strcmp(name,'conservationTolerance'));value=5e-10;
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
