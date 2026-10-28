function report=testRocketShellGraphPreconditioner(folder)
%TESTROCKETSHELLGRAPHPRECONDITIONER 核对完整侧向装配、固定线性扫描和原NS接入。
% 模块一：独立非线性两点替代通量逐未知量差分，保留所有内部Euler精确波系。
% 模块二：完整稠密矩阵逐色直接解作扫描参考；叠加检验固定线性映射。
% 模块三：调用原完整NS做一次真实Newton修正，并复现旧直接稳态完整输出。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'GraphTest:Overwrite','不得覆盖图预条件校核。');mkdir(folder);
s=load(fullfile(root,'outputs','v83_validation','newtonunit03','steady_newton_unit_evidence.mat'),'result','cut','gas','boundary');
cut=s.cut;gas=s.gas;boundary=s.boundary;prior=s.result;Q=prior.initialQ;scale=prior.residualScale;
tau=1e-4;operator=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,tau,2);
n=size(Q,1);M=matrix(operator);U=Q./scale;J=zeros(5*n);h=1e-6;
for id=1:n,for b=1:5
    plus=U;minus=U;step=h*max(1,abs(U(id,b)));plus(id,b)=plus(id,b)+step;minus(id,b)=minus(id,b)-step;
    rp=surrogate(plus.*scale,Q,cut,gas,boundary.solidWall.temperature)./scale;
    rm=surrogate(minus.*scale,Q,cut,gas,boundary.solidWall.temperature)./scale;
    J(:,(id-1)*5+b)=reshape(((rp-rm)/(2*step))',[],1);
end,end
fdMatrix=eye(5*n)/tau-J;assemblyError=max(abs(M-fdMatrix),[],'all')/max(1,max(abs(fdMatrix),[],'all'));
assert(assemblyError<5e-8,'GraphTest:Assembly','完整侧向图装配与独立替代通量差分不一致。');
rhs=reshape(sin((1:5*n)*.37),n,5);[update,scan]=solveRocketShellGraphPreconditioner(operator,rhs);
expected=zeros(5*n,1);rhsVector=reshape(rhs',[],1);nc=numel(operator.colorGroups);order=[1:nc nc-1:-1:1];
for sweep=1:operator.fixedSymmetricSweeps,for c=order
    ids=operator.colorGroups{c}.cellIDs;rows=reshape((5*(ids-1)+(1:5))',[],1);others=setdiff((1:5*n)',rows);
    expected(rows)=M(rows,rows)\(rhsVector(rows)-M(rows,others)*expected(others));
end,end
expected=reshape(expected,5,n)';scanError=max(abs(update-expected),[],'all')/max(1,max(abs(expected),[],'all'));
action=applyRocketShellGraphOperator(operator,update);matrixAction=reshape(M*reshape(update',[],1),5,n)';
actionError=max(abs(action-matrixAction),[],'all')/max(1,max(abs(matrixAction),[],'all'));
assert(scanError<1e-11&&actionError<1e-12,'GraphTest:Scan','逐色直接解或完整图作用核对失败。');
second=reshape(cos((1:5*n)*.51),n,5);u2=solveRocketShellGraphPreconditioner(operator,second);
combined=solveRocketShellGraphPreconditioner(operator,.7*rhs-1.3*second);linearExpected=.7*update-1.3*u2;
linearity=max(abs(combined-linearExpected),[],'all')/max(1,max(abs(linearExpected),[],'all'));
assert(linearity<1e-12,'GraphTest:Linearity','右预条件扫描不是固定线性映射。');
oldDirection=prior.linearLog{1}.steadyDirectionEvidence.directionU;
graphAction=oldDirection/tau-applyRocketShellGraphOperator(operator,oldDirection);
old=load(fullfile(root,'outputs','v83_validation','lateralunit01','lateral_coupling_evidence.mat'),'exactSideGraphAction');
priorActionError=max(abs(graphAction-old.exactSideGraphAction),[],'all')/max(1,max(abs(old.exactSideGraphAction),[],'all'));
assert(priorActionError<1e-11,'GraphTest:Direction','图装配没有复现此前已独立核对的固定方向。');

%% 模块四：实际NS和默认合同保持，失败预算仍由完整方程门槛拒绝。
options=prior.options;options.rightPreconditioner='normal-line';
normal=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(isequaln(untimed(normal),untimed(prior)),'GraphTest:Compatibility','显式旧法向线未逐字段复现先前直接稳态结果。');
options.rightPreconditioner='lateral-colored-line';options.graphSweeps=2;
result=solveRocketShellSteadyKrylov(Q,cut,gas,boundary,options);
assert(result.acceptedNonlinearIterations==1&&result.linearLog{1}.fullSystemTargetMet&&result.fixedLinearPreconditioner, ...
    'GraphTest:ActualNS','侧向图未通过小域原完整NS实际方向和非线性验收。');
rejected=0;
for kind=1:3
    bad=options;
    if kind==1,bad.rightPreconditioner='unknown';elseif kind==2,bad.normalEulerLinearization='rusanov-frozen';else,bad.rightPreconditioner='normal-line';end
    caught=false;try,solveRocketShellSteadyKrylov(Q,cut,gas,boundary,bad);catch ex,caught=strcmp(ex.identifier,'ShellKrylov:Preconditioner');end
    assert(caught);rejected=rejected+1;
end
report=struct('status','passed-full-small-domain-FD-graph-assembly-dense-color-solve-linearity-and-actual-NS', ...
    'cellCount',n,'colorCount',nc,'fullGraphFDRelativeDifference',assemblyError,'denseColoredScanRelativeDifference',scanError, ...
    'matrixActionRelativeDifference',actionError,'linearSuperpositionRelativeDifference',linearity, ...
    'priorIndependentDirectionRelativeDifference',priorActionError,'completePriorNormalResultExactlyReplayed',true, ...
    'fixedSymmetricSweeps',operator.fixedSymmetricSweeps,'invalidPolicyRejections',rejected, ...
    'actualAcceptedIterations',result.acceptedNonlinearIterations,'actualNSCalls',result.residualEvaluations, ...
    'actualSteadyLinearDefect',result.linearLog{1}.fullSystemRelativeDefect,'physicalSolution',false,'newProductionSolvedAreaM2',0, ...
    'engineProcessID',feature('getpid'));
save(fullfile(folder,'graph_preconditioner_evidence.mat'),'report','operator','M','fdMatrix','rhs','update','expected','scan','second','u2','combined', ...
    'graphAction','priorActionError','-v7');
save(fullfile(folder,'steady_newton_unit_evidence.mat'),'report','result','cut','gas','boundary','-v7');
fid=fopen(fullfile(folder,'graph_preconditioner_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function out=untimed(out)
out.history=rmfield(out.history,{'operatorBuildSeconds','linearSolveSeconds'});
end
function M=matrix(op)
nf=op.sourceFaces;nl=op.layers;n=nf*nl;M=zeros(5*n);
for id=1:n
    rows=(id-1)*5+(1:5);M(rows,rows)=op.diagonal(:,:,id);
    if id>nf,M(rows,rows-5*nf)=op.lower(:,:,id);end
    if id<=n-nf,M(rows,rows+5*nf)=op.upper(:,:,id);end
end
for k=1:numel(op.sideMinus)
    rows=(op.sideMinus(k)-1)*5+(1:5);cols=(op.sidePlus(k)-1)*5+(1:5);
    M(rows,cols)=M(rows,cols)+op.minusCoupling(:,:,k);M(cols,rows)=M(cols,rows)+op.plusCoupling(:,:,k);
end
end
function [P,p]=primitive(Q,gas)
rho=Q(:,1);u=Q(:,2:4)./rho;T=(Q(:,5)./rho-.5*sum(u.^2,2))/(gas.R/(gas.gamma-1));P=[rho,u,T];p=rho*gas.R.*T;
end
function value=flux(QL,QR,normal,distance,gas,wave,workVelocity,inside)
[PL,pL]=primitive(QL,gas);[PR,pR]=primitive(QR,gas);unL=sum(PL(:,2:4).*normal,2);unR=sum(PR(:,2:4).*normal,2);
if inside,value=cutCellExactEulerFlux(QL,QR,normal,gas);
else
    FL=[QL(:,1).*unL,QL(:,2:4).*unL+pL.*normal,(QL(:,5)+pL).*unL];
    FR=[QR(:,1).*unR,QR(:,2:4).*unR+pR.*normal,(QR(:,5)+pR).*unR];value=.5*(FL+FR)-.5*wave.*(QR-QL);
end
gradient=(PR(2:4)-PL(2:4))'*normal/distance;stress=gas.mu*(gradient+gradient'-(2/3)*trace(gradient)*eye(3));
traction=(stress*normal')';heat=gas.mu*gas.gamma*gas.R/(gas.gamma-1)/gas.Pr*(PR(5)-PL(5))/distance;
value=value-[0,traction,dot(workVelocity,traction)+heat];
end
function rate=surrogate(Q,base,cut,gas,Twall)
v=cut.cells.volume;n=numel(v);rate=zeros(n,5);f=cut.openFaces;w=cut.walls;[P,p]=primitive(Q,gas);[P0,~]=primitive(base,gas);
for j=1:numel(f.area)
    m=f.cellMinus(j);other=f.cellPlus(j);normal=f.normal(j,:);
    if other>0,qR=Q(other,:);pR=P0(other,:);distance=dot(cut.cells.centroid(other,:)-cut.cells.centroid(m,:),normal);
    else,qR=base(m,:);pR=P0(m,:);distance=dot(f.centroid(j,:)-cut.cells.centroid(m,:),normal);end
    pL=P0(m,:);wave=max(abs(dot(pL(2:4),normal))+sqrt(gas.gamma*gas.R*pL(5)),abs(dot(pR(2:4),normal))+sqrt(gas.gamma*gas.R*pR(5)));
    F=flux(Q(m,:),qR,normal,distance,gas,wave,.5*(pL(2:4)+pR(2:4)),other>0);
    rate(m,:)=rate(m,:)-f.area(j)/v(m)*F;if other>0,rate(other,:)=rate(other,:)+f.area(j)/v(other)*F;end
end
k=gas.mu*gas.gamma*gas.R/(gas.gamma-1)/gas.Pr;
for j=1:numel(w.area)
    id=w.ownerCell(j);normal=w.normal(j,:);distance=dot(w.centroid(j,:)-cut.cells.centroid(id,:),normal);
    gradient=-P(id,2:4)'*normal/distance;stress=gas.mu*(gradient+gradient'-(2/3)*trace(gradient)*eye(3));
    F=[0,p(id)*normal-(stress*normal')',k*(P(id,5)-Twall)/distance];rate(id,:)=rate(id,:)-w.area(j)/v(id)*F;
end
end
