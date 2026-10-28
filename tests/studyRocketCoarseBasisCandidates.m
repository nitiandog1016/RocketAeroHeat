function report=studyRocketCoarseBasisCandidates(folder,kind)
%STUDYROCKETCOARSEBASISCANDIDATES 在保存失败剩余上比较局部空间细分与原始变量耦合粗基。
% 模块一：同一原面/原CV/Q/缩放和真实原NS剩余比较，不按好看的模型标签选方法。
% 模块二：两种像都受256MiB预算约束，小域各两个方向直接对原NS中心差分。
% 模块三：保存细层作用重放、固定叠加和投影；不接受实际流场、不增加合格面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CoarseCandidates:Overwrite','不得覆盖粗基候选研究。');assert(ismember(kind,{'unit','real'}));mkdir(folder);
if strcmp(kind,'unit'),parentName='layercoarseunit01';else,parentName='layercoarserocket01';end
parentFile=fullfile(root,'outputs','v83_validation',parentName,'coarse_layer_evidence.mat');s=load(parentFile,'coarse');parent=s.coarse;clear s;
Q=parent.baseQ;cut=parent.space.geometry;gas=parent.gas;boundary=parent.boundary;scale=parent.scale;sourceFile='';sharedNS=0;
if strcmp(kind,'unit')
    [rate,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);sharedNS=1;rhs=rate./scale;
    remaining=reshape(sin((1:numel(Q))*.37),size(Q));
else
    sourceFile=fullfile(root,'outputs','v83_validation','krylov10','wall_velocity_krylov_evidence.mat');s=load(sourceFile,'result');
    assert(isequaln(Q,s.result.Q)&&isequaln(s.result.initialQ,s.result.Q));rhs=s.result.initialRate./scale;
    remaining=rhs+s.result.linearLog{end}.steadyDirectionEvidence.value;info=struct;
    for key={'primitive','gradient','viscousFaceFlux'},info.(key{1})=s.result.initialInfo.(key{1});end;clear s;
end
cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);clear info;
[parent,reuse]=loadRocketLeastSquaresCoarseOperator(parentFile,cache,parent.space,256*1024^2);
oldCorrection=solveRocketLeastSquaresCoarseSystem(parent,remaining);oldPredicted=remaining+applyCutCellSelectedDirection(cache,oldCorrection);
sharedModelApplications=reuse.reuseVerificationModelApplications+1;
coarse=buildRocketPrimitiveCoarseOperator(cache,parent.space,256*1024^2);
primitive=validateCandidate(coarse,'primitive',fullfile(folder,'primitive'));clear coarse;
options=struct;if strcmp(kind,'unit'),options.axialFractionThreshold=.5;end
space=buildRocketTargetedCoarseSpace(parent.space,options);
coarse=buildRocketRefinedCoarseOperator(cache,space,parent,256*1024^2);
spatial=validateCandidate(coarse,'spatial',fullfile(folder,'spatial'));clear coarse;
report=struct('status','tested-both-coarse-candidates-with-original-kernel-replay-and-no-flow-update','kind',kind,'primitive',primitive,'spatial',spatial, ...
    'parentCoarseUnknowns',parent.space.coarseUnknowns,'sharedModelApplications',sharedModelApplications,'sharedOriginalNSCalls',sharedNS, ...
    'actualOriginalNSEvaluations',sharedNS+primitive.validationOriginalNSCalls+spatial.validationOriginalNSCalls, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
fid=fopen(fullfile(folder,'coarse_candidates_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function result=validateCandidate(coarse,basis,destination)
        primitiveMode=strcmp(basis,'primitive');n=size(Q,1);ng=coarse.space.groupCount;calls=0;modelApplications=0;
        [correction,linear]=project(remaining);predicted=remaining+applyCutCellSelectedDirection(cache,correction);modelApplications=modelApplications+1;
        secondRemaining=remaining.*reshape(cos((1:numel(remaining))*.43),size(remaining));other=project(secondRemaining);combined=project(.7*remaining-1.3*secondRemaining);
        linearity=norm(combined-(.7*correction-1.3*other),'fro')/max(1,norm(combined,'fro'));assert(linearity<1e-8,'CoarseCandidates:Linearity','固定粗校正叠加失败。');
        coefficients=reshape(sin((1:5*ng)*.17),ng,5);
        if primitiveMode,direction=liftRocketPrimitiveCoarse(Q,gas,scale,coarse.space,coefficients);
        else,direction=coarse.space.prolongation*coefficients;end
        actual=-applyCutCellSelectedDirection(cache,direction);modelApplications=modelApplications+1;
        imageProduct=coarse.image*(coefficients(:)./coarse.columnLiftScale);imageError=norm(actual(:)-imageProduct)/max(1,norm(actual(:)));
        imageProbeCoefficients=coefficients;imageProbeDirection=direction;
        assert(imageError<1e-9,'CoarseCandidates:Image','粗像不符合当前完整空间作用。');
        validationDirections={};validationPlus={};validationMinus={};validationModel={};h=1e-6;fdError=[];
        if strcmp(kind,'unit')
            for j=1:2
                coefficients=reshape(sin((1:5*ng)*(.13+.09*j)),ng,5);
                if primitiveMode,direction=liftRocketPrimitiveCoarse(Q,gas,scale,coarse.space,coefficients);
                else,direction=coarse.space.prolongation*coefficients;end
                plus=cutCellNavierStokesResidual(Q+h*direction.*scale,cut,gas,boundary);
                minus=cutCellNavierStokesResidual(Q-h*direction.*scale,cut,gas,boundary);calls=calls+2;
                measured=(plus-minus)./scale/(2*h);model=applyCutCellSelectedDirection(cache,direction);modelApplications=modelApplications+1;
                error=norm(measured-model,'fro')/max(1,norm(measured,'fro'));assert(error<2e-6,'CoarseCandidates:NS','候选粗方向与原NS差分不一致。');
                fdError(end+1)=error;validationDirections{j}=direction;validationPlus{j}=plus;validationMinus{j}=minus;validationModel{j}=model; %#ok<AGROW>
            end
        end
        copied=0;if isfield(coarse,'copiedParentImageColumns'),copied=coarse.copiedParentImageColumns;end
        result=struct('basis',basis,'coarseUnknowns',5*ng,'coarseGroups',ng,'imageBytes',coarse.imageBytes,'imageNonzeros',coarse.imageNonzeros, ...
            'reciprocalConditionEstimate',coarse.reciprocalConditionEstimate,'originalRemainingRelativeDefect',norm(remaining,'fro')/norm(rhs,'fro'), ...
            'parentProjectedRelativeDefect',norm(oldPredicted,'fro')/norm(rhs,'fro'),'candidateProjectedRelativeDefect',norm(predicted,'fro')/norm(rhs,'fro'), ...
            'extraSquaredNormFractionRemoved',1-norm(predicted,'fro')^2/norm(oldPredicted,'fro')^2,'imageReplayRelativeDifference',imageError, ...
            'fixedLinearityRelativeDifference',linearity,'originalNSDirectionFDRelativeDifferences',fdError,'validationOriginalNSCalls',calls, ...
            'constructionModelApplications',coarse.constructionModelApplications,'copiedParentImageColumns',copied,'otherModelApplications',modelApplications, ...
            'buildSeconds',coarse.buildSeconds,'physicalSolution',false,'newProductionSolvedAreaM2',0);
        mkdir(destination);save(fullfile(destination,'coarse_candidate_evidence.mat'),'coarse','result','sourceFile','parentFile','rhs','remaining','oldCorrection','oldPredicted', ...
            'correction','predicted','linear','secondRemaining','other','combined','coefficients','direction','actual','imageProduct','imageProbeCoefficients','imageProbeDirection', ...
            'validationDirections','validationPlus','validationMinus','validationModel','h','-v7');
        fid2=fopen(fullfile(destination,'coarse_candidate_summary.json'),'w','n','UTF-8');assert(fid2>=0);guard2=onCleanup(@()fclose(fid2));fprintf(fid2,'%s',jsonencode(result));fprintf('%s\n',jsonencode(result));
        function [value,detail]=project(input)
            if primitiveMode,[value,detail]=solveRocketPrimitiveCoarseSystem(coarse,input);else,[value,detail]=solveRocketLeastSquaresCoarseSystem(coarse,input);end
        end
    end
end
