function report=testRocketRegularizedParentDirection(folder)
%TESTROCKETREGULARIZEDPARENTDIRECTION 用两种真实小域排列和独立稠密解校核状态度量与增广方向。
% 模块一：度量对原始变量中心差分与完整稠密转置；各部件不以自比较代替独立参考。
% 模块二：完整增广问题对直接backslash，核对惩罚总方向及非零初猜右端。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);s=load(fullfile(root,'outputs','v83_validation','parentlsqrunit02','parent_least_squares_unit_evidence.mat'),'cases');
gas=struct('gamma',1.4,'R',287.05287);cases=cell(2,1);
for j=1:2
    a=s.cases{j};n=a.packed.cellCount;N=5*n;ids=(1:n)';rho=.01+.003*sin(ids).^2;u=[300+20*sin(ids) 15*cos(ids) -8*sin(ids/2)];T=300+10*cos(ids);cv=gas.R/(gas.gamma-1);
    Q=[rho rho.*u rho.*(cv*T+.5*sum(u.^2,2))];scale=[.01 3 3 3 1000];metric=buildRocketPrimitiveUpdateMetric(Q,scale,gas);
    W=zeros(N);for k=1:N,v=zeros(n,5);v(k)=1;z=applyRocketPrimitiveUpdateMetric(metric,v);W(:,k)=z(:);end
    v=reshape(sin((1:N)'),n,5);w=reshape(cos((1:N)'),n,5);direct=applyRocketPrimitiveUpdateMetric(metric,v);back=applyRocketPrimitiveUpdateMetric(metric,w,'transp');
    assert(norm(back(:)-W'*w(:))<1e-12&&norm(direct(:)-W*v(:))<1e-12);
    step=1e-5;P0=primitive(Q,gas);Pplus=primitive(Q+step*v.*scale,gas);Pminus=primitive(Q-step*v.*scale,gas);
    fd=(Pplus-Pminus)/(2*step);fd=fd./[P0(:,1) repmat(sqrt(gas.gamma*gas.R*P0(:,5)),1,3) P0(:,5)];
    fdDifference=norm(fd-direct,'fro')/norm(direct,'fro');assert(fdDifference<1e-8);
    initial=a.initial;rhs=a.rhs;options=struct('maximumIterations',900,'regularizationStrength',.1);
    [value,detail,coordinates]=solveRocketRegularizedParentDirection(a.packed,a.factor,rhs,initial,a.columnScale,metric,options);
    penalty=detail.penaltyInOriginalResidualUnits;A=[a.C;penalty*W];b=[rhs(:);zeros(N,1)];dense=A\b;
    difference=norm(value(:)-dense)/max(1,norm(dense));assert(difference<1e-8&&detail.adjointIdentityRelativeDifference<1e-12);
    assert(detail.totalDirectionPenalized&&~detail.originalNSAcceptanceClaimed&&detail.newActualNSEvaluations==0);
    cases{j}=struct('Q',Q,'scale',scale,'gas',gas,'metric',metric,'W',W,'fdDifference',fdDifference, ...
        'C',a.C,'M',a.M,'columnScale',a.columnScale,'rhs',rhs,'initial',initial,'value',value,'coordinates',coordinates, ...
        'detail',detail,'denseDifference',difference);
end
report=struct('status','passed-two-orderings-primitive-state-FD-full-adjoint-and-dense-augmented-total-direction', ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'regularized_parent_unit_evidence.mat'),'cases','report','-v7');fprintf('%s\n',jsonencode(report));
end
function P=primitive(Q,gas)
rho=Q(:,1);u=Q(:,2:4)./rho;P=[rho u (Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*gas.R/(gas.gamma-1))];
end
