function report=testRocketRecycledGMRESKernel(folder)
%TESTROCKETRECYCLEDGMRESKERNEL 独立矩阵、退化/零右端、真实小域Jacobian及原NS方向验证。
% 模块一：同Arnoldi预算比较保留/不保留，不预先假设新方法一定更快或收敛。
% 模块二：保存小域完整V/Z/保留U/C，供Python独立重建当前像与最小残差。
% 模块三：小域完整模型矩阵由原算子逐列构造，组合方向另用原NS中心差分核对。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RecycleTest:Overwrite','不得覆盖增强Krylov研究。');mkdir(folder);
n=80;A=diag(logspace(-3,0,n))+diag(.08*ones(n-1,1),1);b=sin((1:n)'*.37);
o=struct('restartLength',12,'maximumCycles',4,'retainedDirections',4,'maximumArnoldiSteps',48, ...
    'absoluteTolerance',1e-9*norm(b),'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
[x,recycled]=solveRocketRecycledGMRES(@(v)A*v,@(v)v,b,o);
p=o;p.retainedDirections=0;[plainX,plain]=solveRocketRecycledGMRES(@(v)A*v,@(v)v,b,p);
assert(abs(norm(b-A*x)-recycled.absoluteResidual)<1e-10*max(1,norm(b))&& ...
    abs(norm(b-A*plainX)-plain.absoluteResidual)<1e-10*max(1,norm(b)));
[identityX,identity]=solveRocketRecycledGMRES(@(v)v,@(v)v,b,o);assert(identity.flag==0&&norm(identityX-b)<1e-12&&identity.arnoldiSteps==1);
[zeroX,zero]=solveRocketRecycledGMRES(@(v)v,@(v)v,zeros(n,1),o);assert(zero.flag==0&&zero.modelApplications==0&&all(zeroX==0));
bad=o;bad.maximumBasisBytes=1;caught=false;try,solveRocketRecycledGMRES(@(v)v,@(v)v,b,bad);catch ex,caught=strcmp(ex.identifier,'RecycledGMRES:Budget');end;assert(caught);
sourceFile=fullfile(root,'outputs','v83_validation','localunit02','local_coarse_evidence.mat');s=load(sourceFile,'coarse');c=s.coarse;
Q=c.baseQ;gas=c.gas;boundary=c.boundary;scale=c.scale;cut=c.space.geometry;clear c s;
[~,info]=cutCellNavierStokesResidual(Q,cut,gas,boundary);cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);
nu=numel(Q);unitA=zeros(nu);for col=1:nu,v=zeros(size(Q));v(col)=1;w=-applyCutCellSelectedDirection(cache,v);unitA(:,col)=w(:);end
graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-4,1);unitB=sin((1:nu)'*.29);unitOptions=o;
unitOptions.restartLength=3;unitOptions.maximumCycles=4;unitOptions.retainedDirections=2;unitOptions.maximumArnoldiSteps=12;unitOptions.absoluteTolerance=1e-10*norm(unitB);
preconditions=0;[unitX,unitReport]=solveRocketRecycledGMRES(@(v)unitA*v,@physicalPrecondition,unitB,unitOptions);
assert(norm(unitB-unitA*unitX)<=norm(unitB)*(1+1e-10)&&preconditions==unitReport.preconditionerApplications);
validationDirections={};validationPlus={};validationMinus={};validationModel={};errors=[];h=1e-6;
for j=1:2
    if j==1,direction=reshape(unitX,size(Q));else,direction=reshape(cos((1:nu)'*.31),size(Q));end
    direction=direction/max(abs(direction),[],'all');
    plus=cutCellNavierStokesResidual(Q+h*direction.*scale,cut,gas,boundary);minus=cutCellNavierStokesResidual(Q-h*direction.*scale,cut,gas,boundary);
    measured=(plus-minus)./scale/(2*h);model=reshape(-unitA*direction(:),size(Q));error=norm(measured-model,'fro')/max(1,norm(measured,'fro'));
    assert(error<2e-6,'RecycleTest:NS','当前小域完整Jacobian与原NS组合方向不符。');
    validationDirections{j}=direction;validationPlus{j}=plus;validationMinus{j}=minus;validationModel{j}=model;errors(end+1)=error; %#ok<AGROW>
end
report=struct('status','passed-bounded-SVD-recycled-kernel-matrix-orthogonality-degeneracy-budget-and-unit-NS-directions', ...
    'matrixPlainRelativeResidual',plain.relativeResidual,'matrixRecycledRelativeResidual',recycled.relativeResidual, ...
    'unitCurrentJacobianRelativeResidual',unitReport.relativeResidual,'unitOriginalNSDirectionFDErrors',errors, ...
    'actualOriginalNSEvaluations',5,'unitJacobianConstructionModelApplications',nu,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'recycled_kernel_evidence.mat'),'report','A','b','x','recycled','plainX','plain','identityX','identity','zero','unitA','unitB','unitX','unitReport', ...
    'Q','cut','gas','boundary','scale','sourceFile','validationDirections','validationPlus','validationMinus','validationModel','h','-v7');
fid=fopen(fullfile(folder,'recycled_kernel_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function out=physicalPrecondition(v)
        [value,~]=solveRocketShellGraphPreconditioner(graph,reshape(v,size(Q)));out=value(:);preconditions=preconditions+1;
    end
end
