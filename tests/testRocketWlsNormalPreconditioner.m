function report=testRocketWlsNormalPreconditioner(folder)
%TESTROCKETWLSNORMALPRECONDITIONER 核对完整WLS同线提取、五带分解和固定扫描。
% 模块一：小域逐列完整方向构造独立参考，验证同线五带而非两点近似。
% 模块二：原侧向耦合保留，独立稠密逐色求解对照前后代及线性叠加。
% 模块三：两组原NS方向差分和预算拒绝；本测试不推进真实箭体或增加面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'),fullfile(root,'tests'));setupCutCellFlow;
assert(~exist(folder,'dir'),'WlsNormalTest:Overwrite','不得覆盖五带研究证据。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','recycledkernel01','recycled_kernel_evidence.mat');
s=load(sourceFile,'Q','cut','gas','boundary','scale','unitA');Q=s.Q;cut=s.cut;gas=s.gas;boundary=s.boundary;scale=s.scale;unitA=s.unitA;clear s;
[rate,baseInfo]=cutCellNavierStokesResidual(Q,cut,gas,boundary);cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);
graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-4,1);limit=32*1024^2;
operator=buildRocketWlsNormalPreconditioner(graph,cache,limit);n=size(Q,1);nf=graph.sourceFaces;nl=graph.layers;N=numel(Q);
% 模块四：混合内部/入流/出流面及单元素尾批，防止嵌套函数改写外层批次数。
[batchReference,~]=assembleRocketWlsNormalCorrection(cache,limit);batchSizes=[1 3 17];batchErrors=zeros(size(batchSizes));
for j=1:numel(batchSizes)
    [candidate,~]=assembleRocketWlsNormalCorrection(cache,limit,batchSizes(j));
    batchErrors(j)=norm(candidate(:)-batchReference(:))/max(1,norm(batchReference(:)));assert(batchErrors(j)<2e-12,'WlsNormalTest:Batch','小批/尾批改变了同线装配。');
end
fullViscous=zeros(N);graphViscous=zeros(N);modelReferenceError=0;
for col=1:N
    basis=zeros(size(Q));basis(col)=1;[full,approx,detail]=splitRocketGraphDirection(cache,graph,basis);
    fullValue=full.momentumStress+full.heatConduction+full.workGradient+full.workVelocity;
    graphValue=approx.momentumStress+approx.heatConduction+approx.workGradient+approx.workVelocity;
    fullViscous(:,col)=fullValue(:);graphViscous(:,col)=graphValue(:);
    modelReferenceError=max(modelReferenceError,norm(detail.fullTotal(:)+unitA(:,col))/max(1,norm(unitA(:,col))));
end
cellIDs=repmat((1:n)',5,1);sameLine=mod(cellIDs-1,nf)==mod(cellIDs'-1,nf);difference=fullViscous-graphViscous;
sameLineReference=difference.*sameLine;oldMatrix=matrix(graph,false);newMatrix=matrix(operator,true);extracted=oldMatrix-newMatrix;
lineError=norm(extracted-sameLineReference,'fro')/max(1,norm(sameLineReference,'fro'));
assert(lineError<2e-10&&modelReferenceError<2e-13,'WlsNormalTest:Assembly','提取的完整WLS同线块或原完整作用不符。');
far=abs(floor((cellIDs-1)/nf)-floor((cellIDs'-1)/nf))>2;
assert(norm(sameLineReference.*far,'fro')<2e-12*max(1,norm(sameLineReference,'fro')),'WlsNormalTest:Bandwidth','同线真参考含五带之外的作用。');
rhs=reshape(sin((1:N)*.37),n,5);[update,scan]=solveRocketWlsNormalPreconditioner(operator,rhs);
expected=zeros(N,1);nc=numel(operator.colorGroups);order=[1:nc nc-1:-1:1];
for sweep=1:operator.fixedSymmetricSweeps,for color=order
    ids=operator.colorGroups{color}.cellIDs;rows=reshape(ids+(0:4)*n,[],1);others=setdiff((1:N)',rows);
    expected(rows)=newMatrix(rows,rows)\(rhs(rows)-newMatrix(rows,others)*expected(others));
end,end
scanError=norm(update(:)-expected)/max(1,norm(expected));applied=applyRocketWlsNormalOperator(operator,update);
matrixError=norm(applied(:)-newMatrix*update(:))/max(1,norm(applied(:)));
second=reshape(cos((1:N)*.23),n,5);other=solveRocketWlsNormalPreconditioner(operator,second);
combined=solveRocketWlsNormalPreconditioner(operator,.7*rhs-1.3*second);
linearity=norm(combined-.7*update+1.3*other,'fro')/max(1,norm(combined,'fro'));
assert(max([scanError matrixError linearity])<2e-11,'WlsNormalTest:Scan','五带扫描/矩阵作用/固定线性映射不符。');
badBudget=false;try,buildRocketWlsNormalPreconditioner(graph,cache,1);catch ex,badBudget=strcmp(ex.identifier,'WlsNormal:Budget');end;assert(badBudget);
h=1e-6;directions={update,rate./scale};fdPlus=cell(2,1);fdMinus=fdPlus;fdModels=fdPlus;fdErrors=zeros(2,1);
for j=1:2
    directions{j}=directions{j}/max(abs(directions{j}),[],'all');direction=directions{j};
    fdPlus{j}=cutCellNavierStokesResidual(Q+h*direction.*scale,cut,gas,boundary);fdMinus{j}=cutCellNavierStokesResidual(Q-h*direction.*scale,cut,gas,boundary);
    fdModels{j}=reshape(-unitA*direction(:),size(Q));measured=(fdPlus{j}-fdMinus{j})./scale/(2*h);
    fdErrors(j)=norm(measured-fdModels{j},'fro')/max(1,norm(measured,'fro'));assert(fdErrors(j)<2e-6,'WlsNormalTest:NS','原NS方向检查失败。');
end
report=struct('status','passed-current-WLS-same-line-extraction-five-band-dense-color-solve-linearity-and-original-NS-directions', ...
    'cellCount',n,'sourceFaces',nf,'layers',nl,'lineExtractionRelativeDifference',lineError,'currentFullReferenceRelativeDifference',modelReferenceError, ...
    'denseColoredScanRelativeDifference',scanError,'matrixActionRelativeDifference',matrixError,'linearSuperpositionRelativeDifference',linearity, ...
    'originalNSDirectionFDErrors',fdErrors,'batchSizes',batchSizes,'batchAssemblyRelativeDifferences',batchErrors, ...
    'budgetRejected',badBudget,'completeSelectedReferenceApplications',N, ...
    'actualOriginalNSEvaluations',5,'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'wls_normal_unit_evidence.mat'),'report','sourceFile','Q','cut','gas','boundary','scale','baseInfo','rate','unitA','graph','operator', ...
    'fullViscous','graphViscous','sameLineReference','oldMatrix','newMatrix','rhs','update','expected','scan','second','other','combined', ...
    'directions','fdPlus','fdMinus','fdModels','h','-v7');
fid=fopen(fullfile(folder,'wls_normal_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function M=matrix(op,second)
n=op.sourceFaces*op.layers;nf=op.sourceFaces;M=zeros(5*n);
for id=1:n
    rows=id+(0:4)*n;M(rows,rows)=op.diagonal(:,:,id);
    if id>nf,M(rows,rows-nf)=op.lower(:,:,id);end
    if id<=n-nf,M(rows,rows+nf)=op.upper(:,:,id);end
    if second
        if id>2*nf,M(rows,rows-2*nf)=op.secondLower(:,:,id);end
        if id<=n-2*nf,M(rows,rows+2*nf)=op.secondUpper(:,:,id);end
    end
end
for k=1:numel(op.sideMinus)
    rows=op.sideMinus(k)+(0:4)*n;cols=op.sidePlus(k)+(0:4)*n;
    M(rows,cols)=M(rows,cols)+op.minusCoupling(:,:,k);M(cols,rows)=M(cols,rows)+op.plusCoupling(:,:,k);
end
end
