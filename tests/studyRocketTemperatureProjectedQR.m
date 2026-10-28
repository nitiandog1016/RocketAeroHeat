function report=studyRocketTemperatureProjectedQR(folder,preparationFile,faceIDs)
%STUDYROCKETTEMPERATUREPROJECTEDQR 同当前状态、同128步预算比较自由与界面温度投影方向。
% 模块一：当前状态、完整C、W和penalty共同；只显式增加采样温度的一阶齐次方向限制。
% 模块二：两方向分别核对完整模型与新鲜原NS，均从同一真实前驱独立试算，不拼接并列历史。
% 模块三：保留全部拒绝/实际调用和研究身份，一阶投影不替代任何非线性温域护栏。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
folder=char(java.io.File(folder).getCanonicalPath());assert(~exist(folder,'dir'));mkdir(folder);
p=loadRocketPeakSystem(preparationFile);preparationSHA256=rocketStreamingFileSHA256(preparationFile);
a=load(p.artifacts.packed.file,'packed');b=load(p.artifacts.factor.file,'factor');m=load(p.artifacts.metric.file);
penalty=.1*norm(m.rhs,'fro')/p.metricNormalization;factor=buildRocketAugmentedPatchQR(a.packed,b.factor,m.metric,penalty);clear b;
c=load(p.artifacts.cache.file,'cache');rows=buildRocketInterfaceTemperatureRows(c.cache.Q,c.cache.geometry,c.cache.gas,c.cache.scale,faceIDs);clear c;
projector=buildRocketTemperatureProjector(factor,rows);
save(fullfile(folder,'temperature_projection_system.mat'),'factor','rows','projector','preparationFile','preparationSHA256','-v7');
[free,freeDetail,freeCoordinates]=solveRocketAugmentedQRDirection(a.packed,factor,m.rhs,m.metric);
[projected,projectedDetail,projectedCoordinates]=solveRocketTemperatureProjectedDirection(a.packed,factor,m.rhs,m.metric,projector,rows);clear a factor;
rhs=m.rhs;metric=m.metric;sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;geometryBinding=p.geometryBinding;clear m;
save(fullfile(folder,'temperature_projection_directions.mat'),'free','freeDetail','freeCoordinates','projected','projectedDetail','projectedCoordinates', ...
    'rhs','metric','rows','preparationFile','preparationSHA256','sourceBinding','sourceIdentity','geometryBinding','-v7');
c=load(p.artifacts.cache.file,'cache');[freeModel,freeModelDetail]=applyCutCellSelectedDirection(c.cache,free);
[projectedModel,projectedModelDetail]=applyCutCellSelectedDirection(c.cache,projected);clear c;
values={free,projected};models={freeModel,projectedModel};details={freeDetail,projectedDetail};branches=cell(1,2);names={'free','projected'};
clear free projected freeModel projectedModel;
for j=1:2
    branchFolder=fullfile(folder,names{j});mkdir(branchFolder);value=values{j};model=models{j};
    modelDefect=norm(rhs+model,'fro')/norm(rhs,'fro');detail=details{j};constraintDerivativeK=rows.H*value(:);
    branch=struct('method',names{j},'linearReport',detail,'completeSelectedModelRelativeDefect',modelDefect, ...
        'constraintDerivativeK',constraintDerivativeK,'probe',[],'main',[],'resultFolder','', ...
        'newProbeActualNS',0,'newMainActualNS',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    save(fullfile(branchFolder,'direction.mat'),'value','model','detail','sourceBinding','sourceIdentity','geometryBinding','preparationFile','preparationSHA256','-v7');
    if modelDefect<=.1
        probeFolder=fullfile(branchFolder,'probe');probe=createRocketNSDirectionProbe(probeFolder,p.sourceBinding.sourceFile,p.geometryBinding.file,p.geometryBinding.variable,value,'result');
        branch.probe=probe;branch.newProbeActualNS=probe.newActualNSEvaluations;
        if probe.originalNSDirectionTargetMet
            options=struct('linearEquation','steady-newton','normalEulerLinearization','exact-branch','jacobianAction','original-difference', ...
                'maximumIterations',1,'initialPseudoDt',1e-7,'maximumKrylovIterations',48,'linearRelativeTolerance',.1, ...
                'maximumStepRetries',0,'maximumLineSearch',16,'steadyInitialGuess', ...
                struct('file',fullfile(probeFolder,'direct_ns_direction_probe.mat'),'variant','ns-probe','scale',1));
            branch.resultFolder=fullfile(branchFolder,'krylov01');
            branch.main=continueRocketWallVelocityKrylov(branch.resultFolder,fileparts(p.sourceBinding.sourceFile),options);
            branch.newMainActualNS=branch.main.residualEvaluations;
        end
    end
    branches{j}=branch;
    fid=fopen(fullfile(branchFolder,'branch.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(branch));fclose(fid);
end
report=struct('status','tested-current-state-free-versus-temperature-row-projected-QR-bounded-NS-comparison', ...
    'sourceBinding',sourceBinding,'preparationFile',preparationFile,'preparationSHA256',preparationSHA256, ...
    'faceIDs',faceIDs,'branches',{branches},'sameInitialState',true,'nativeLSQRBudgetPerBranch',128, ...
    'sameOriginalCAndPrimitivePenalty',true,'additionalDirectionConstraintExplicit',true, ...
    'nonlinearTemperatureBoundGuaranteed',false,'parallelHistoriesMerged',false,'physicalSolution',false, ...
    'newSelectedModelApplications',2,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'temperature_projection_models.mat'),'report','freeModelDetail','projectedModelDetail','-v7');
fid=fopen(fullfile(folder,'temperature_projection_comparison.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report));fclose(fid);disp(report);
end
