function report=studyRocketWallVelocityTrialPaths(folder,sourceFolder)
%STUDYROCKETWALLVELOCITYTRIALPATHS 用同一真实Newton方向比较守恒直线与已有原语切向路径。
% 模块一：精确复现第8次端点的完整残差；第9次保守接受步作为逐值对照。
% 模块二：只改变候选坐标，不改空间NS、边界、壁温、物性或80—2000K门槛。
% 模块三：候选须全域合法并检查原RMS下降/五量账；比较不自动授予物理或面积资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
if nargin<2,sourceFolder=fullfile(root,'outputs','v83_validation','krylov03');end
assert(~exist(folder,'dir'),'WallTrial:Overwrite','不得覆盖真实路径对照。');mkdir(folder);
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
s=load(fullfile(sourceFolder,'wall_velocity_krylov_evidence.mat'),'result','gas','boundary');
assert(numel(s.result.trace)>=3,'WallTrial:Source','需要实际保存的第三次方向。');
t=s.result.trace{3};baseQ=t.initialQ;gas=s.gas;boundary=s.boundary;scale=s.result.residualScale;
assert(strcmp(cutCellWallVelocityPolicy(boundary.solidWall),'dirichlet-tangent'));
assert(isequal(gas.temperatureRange,[80 2000])&&s.result.options.linearRelativeTolerance==.2);
calls=0;timer=tic;[baseRate,baseInfo]=evaluate(baseQ);
assert(isequaln(baseRate,t.initialRate),'WallTrial:Source','没有精确复现保存方向的真实基态残差。');
baseSummary=measure(baseQ,baseRate,baseInfo);clear baseInfo;
specs={'conservative',1/64;'primitive-tangent',1/4;'primitive-tangent',1/8};cases=cell(3,1);summaries=struct([]);
for j=1:size(specs,1)
    coordinates=specs{j,1};alpha=specs{j,2};Q=rocketShellTrialState(baseQ,t.linearUpdateU,scale,alpha,gas,coordinates);
    identifier='';message='';info=[];rate=[];stats=struct([]);meritPassed=false;reproduced=false;
    try
        check(Q);[rate,info]=evaluate(Q);stats=measure(Q,rate,info);
        meritPassed=stats.rmsMerit<=(1-s.result.options.sufficientDecrease*alpha)*baseSummary.rmsMerit;
        if j==1
            assert(isequaln(Q,s.result.Q)&&isequaln(rate,s.result.rate)&&isequaln(info,s.result.info), ...
                'WallTrial:Control','守恒对照没有精确复现第9次完整Q/rate/info。');reproduced=true;
        end
    catch ex
        allowed={'CutNS:State','CutNS:InternalEnergy','CutNS:Temperature','CutExactEuler:State', ...
            'CutExactEuler:Temperature','CutExactEuler:InternalEnergy','CutExactEuler:Vacuum', ...
            'CutNS:BoundaryMach','RocketAeroHeat:NSDraft:CharacteristicState','RocketAeroHeat:NSDraft:CharacteristicSound', ...
            'RocketAeroHeat:NSDraft:CharacteristicRegime','RocketAeroHeat:NSDraft:CharacteristicSonic', ...
            'RocketAeroHeat:NSDraft:CharacteristicAmbiguous'};
        if ~ismember(ex.identifier,allowed),rethrow(ex);end
        identifier=ex.identifier;message=ex.message;
    end
    item=struct('coordinates',coordinates,'alpha',alpha,'completeOriginalNSPassed',isempty(identifier), ...
        'originalRmsDecreasePassed',meritPassed,'fullConservativeControlExactlyReproduced',reproduced, ...
        'errorIdentifier',identifier,'message',message,'metrics',stats,'candidatePromotedToSolution',false);
    if isempty(summaries),summaries=item;else,summaries(end+1)=item;end %#ok<AGROW>
    if ~isempty(info),info=pick(info);end
    cases{j}=struct('Q',Q,'rate',rate,'info',info,'summary',item);
    fprintf('WALL_TRIAL_PATH case=%d coordinates=%s alpha=%.9g complete=%d merit=%d error=%s\n', ...
        j,coordinates,alpha,isempty(identifier),meritPassed,identifier);
end
assert(summaries(1).fullConservativeControlExactlyReproduced);
report=struct('status','tested-same-saved-real-direction-complete-NS-trial-path-comparison', ...
    'sourceEvidenceFolder',sourceFolder,'sourceAcceptedEquationIteration',8,'sourceDirectionEquationIteration',9, ...
    'baseMetrics',baseSummary,'cases',summaries,'actualOriginalResidualEvaluations',calls,'elapsedSeconds',toc(timer), ...
    'spatialEquationChanged',false,'physicalBoundaryChanged',false,'originalTemperatureRangeKept',true, ...
    'candidatePromotedToSolution',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
directionU=t.linearUpdateU;options=s.result.options;
save(fullfile(folder,'wall_trial_path_evidence.mat'),'report','baseQ','baseRate','scale','directionU','cases','gas','boundary','options','-v7');
fid=fopen(fullfile(folder,'wall_trial_path_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function [rate,info]=evaluate(Q)
        calls=calls+1;[rate,info]=cutCellNavierStokesResidual(Q,d.domain,gas,boundary);
    end
    function check(Q)
        rho=Q(:,1);assert(all(isfinite(Q),'all')&&all(rho>0),'CutNS:State','候选密度或五量非法，不裁剪。');
        u=Q(:,2:4)./rho;internal=Q(:,5)-.5*rho.*sum(u.^2,2);
        assert(all(internal>0),'CutNS:InternalEnergy','候选内能非正，不裁剪。');
        T=internal./(rho*gas.R/(gas.gamma-1));assert(all(T>=80&T<=2000),'CutNS:Temperature','候选温度越过原物性域，不裁温。');
    end
    function stats=measure(Q,rate,info)
        f=d.domain.openFaces;w=d.domain.walls;mask=f.cellPlus==0;
        denominator=max(1,sum(abs(info.faceFlux(mask,:)).*f.area(mask),1)+sum(abs(info.wallFlux).*w.area,1));
        spatial=max(abs(info.balanceResidual)./denominator);assert(spatial<=s.result.options.spatialConservationTolerance);
        z=rate./scale;T=info.primitive(:,5);
        stats=struct('maximumResidual',max(abs(z),[],'all'),'rmsMerit',sqrt(mean(z.^2,'all')), ...
            'spatialConservationRelative',spatial,'steadyBoundaryFluxBalance',max(abs(info.netBoundaryFlux)./denominator), ...
            'cellTemperatureRangeK',[min(T),max(T)],'minimumDensityKgM3',min(Q(:,1)));
    end
end
function out=pick(info)
out=struct;
for key={'primitive','faceFlux','wallFlux','wallPressure','wallHeatFluxIntoWall','wallTemperatureInput','wallMechanicalPowerIntoFluid'}
    out.(key{1})=info.(key{1});
end
end
