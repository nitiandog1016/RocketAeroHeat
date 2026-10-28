function report=testAxisymmetricInteriorArithmetic(folder)
%TESTAXISYMMETRICINTERIORARITHMETIC 显式Euler双double与完整四/五量残差、Jacobian和BE耦合校核。
% 模块一：未选择算术时重放v74十个完整info；显式double只增加真实方法标签。
% 模块二：两种HLLE/四五量共享通量、原黏性/墙/源不变，SA质量必须完全一致。
% 模块三：有色Jacobian对独立逐列与真实小BE，保留原门槛与非法选项拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v75_validation','arithmetic_coupled');end
assert(~exist(folder,'dir'),'ArithmeticTest:Overwrite','不能覆盖测试。');
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[200 800], ...
    'transport',struct('model','sutherland','referenceViscosity',.05,'referenceTemperature',300,'sutherlandTemperature',110.4)));
[x,r]=ndgrid([0 .02 .05],[.03 .05 .07]);mesh=buildAxisymmetricFVGeometry(x+.07*(r-.03),r,struct('periodicI',true));
n=mesh.cellCount;t=(1:n)';P=[1+.02*sin(t),-2+.2*cos(t),.02*sin(t),500+20*cos(t)];
Q4=axisymmetricFVConservative(P,gas);Q5=[Q4,P(:,1).*(.3+.02*sin(t))];scale5=[1 20 20 2000 .5];
per=struct('type','periodic');bc=struct('iMin',per,'iMax',per, ...
    'jMin',struct('type','wall','temperature',350,'velocity',[-2 0]), ...
    'jMax',struct('type','wall','temperature',650,'velocity',[0 0]),'pressureSourceQuadrature','cell-constant');
cases=struct([]);evidence={};
for nq=[4 5]
    if nq==4,Q=Q4;else,Q=Q5;end;scale=scale5(1:nq);U=Q./scale;
    for policy={'hlle-davis','hlle-einfeldt'}
        ordinary=bc;ordinary.interiorConvectiveFlux=policy{1};[rate0,info0]=evaluate(Q,ordinary);
        explicit=ordinary;explicit.interiorFluxArithmetic='double';[same,tagged]=evaluate(Q,explicit);
        assert(isequal(rate0,same)&&isequal(info0.faceFluxIntegralPerRadian,tagged.faceFluxIntegralPerRadian));
        accurate=ordinary;accurate.interiorFluxArithmetic='two-double-euler';[rate,info]=evaluate(Q,accurate);
        assert(strcmp(info.interiorFluxArithmetic,'two-double-euler')&&contains(info.method,'Euler arithmetic=two-double-euler'));
        assert(isequaln(info.wall,info0.wall)&&isequal(info.cellSourceIntegral,info0.cellSourceIntegral)&& ...
            isequal(info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:),info0.faceFluxIntegralPerRadian(mesh.boundaryFaces,:)));
        expected=-mesh.incidence*(info.faceFluxIntegralPerRadian-info0.faceFluxIntegralPerRadian)./mesh.volumePerRadian;
        assert(max(abs(rate-rate0-expected),[],'all')<=4096*eps*max(1,max(abs(rate)+abs(rate0),[],'all')));
        if nq==5,assert(isequal(info.sa.faceMassFluxPerRadian,info.faceFluxIntegralPerRadian(:,1)));end
        opts=struct('differenceStep',2e-6,'differenceScheme','central');
        [J,jreport]=axisymmetricFVJacobianLocal(Q,mesh,gas,accurate,scale,opts);assert(jreport.success);
        [reference,~]=axisymmetricJacobianLocalReference(Q,gas,scale,opts.differenceStep);dense=zeros(n*nq);
        for col=1:n*nq
            plus=U;minus=U;h=opts.differenceStep*reference(col);plus(col)=plus(col)+h;minus(col)=minus(col)-h;
            rp=evaluate(plus.*scale,accurate)./scale;rm=evaluate(minus.*scale,accurate)./scale;
            dense(:,col)=(rp(:)-rm(:))/(2*h);
        end
        relative=max(abs(full(J)-dense),[],'all')/max(1,max(abs(dense),[],'all'));assert(relative<2e-7);
        solver=struct('maxSteps',1,'initialPseudoCFL',.01,'maximumPseudoCFL',.01,'maxStepRetries',0, ...
            'maxNewtonIterations',12,'jacobianRefreshEvery',1,'residualScale',scale);
        trial=advanceAxisymmetricSteadyLocal(Q,mesh,gas,accurate,solver);assert(trial.acceptedSteps==1&&trial.globalBalanceRelative<5e-10);
        row=struct('quantities',nq,'policy',policy{1},'coloredVsSerialRelative',relative, ...
            'acceptedSteps',trial.acceptedSteps,'BEConservationRelative',trial.globalBalanceRelative,'residualEvaluations',trial.residualEvaluations);
        if isempty(cases),cases=row;else,cases(end+1)=row;end %#ok<AGROW>
        evidence{end+1}=struct('Q',Q,'scale',scale,'boundary',accurate,'J',J,'dense',dense,'jreport',jreport,'trial',trial,'rate',rate,'info',info); %#ok<AGROW>
    end
end
source=fullfile(root,'outputs','v74_validation','local_Q23');d=load(fullfile(source,'rocket_axisymmetric_evidence.mat'),'caseDefinition','durable');
for k=1:10
    s=load(fullfile(source,d.durable.steps(k).relativePath),'segment');c=d.caseDefinition;
    [~,info]=axisymmetricRANSResidual(s.segment.Q,c.mesh,c.gas,c.boundary);assert(isequaln(info,s.segment.finalInfo));
end
refusals=0;bad=bc;bad.interiorConvectiveFlux='hlle-einfeldt';bad.interiorFluxArithmetic='unknown';reject(@()evaluate(Q5,bad));
bad.interiorFluxArithmetic='two-double-euler';bad.interiorConvectiveFlux='hllc';reject(@()evaluate(Q5,bad));
bad=bc;bad.jMin.interiorFluxArithmetic='two-double-euler';reject(@()evaluate(Q5,bad));
report=struct('status','passed-explicit-Euler-arithmetic-four-five-quantity-coupling-and-default-replay', ...
    'matlabEngineProcessID',feature('getpid'),'cases',cases,'legacyTenCompleteInfoReplayed',true, ...
    'refusals',refusals,'fullRANSHighPrecision',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'arithmetic_coupled_evidence.mat'),'report','evidence','mesh','gas','-v7');
fid=fopen(fullfile(folder,'arithmetic_coupled_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function [rate,info]=evaluate(state,boundary)
        if size(state,2)==5,[rate,info]=axisymmetricRANSResidual(state,mesh,gas,boundary);else,[rate,info]=axisymmetricNavierStokesResidual(state,mesh,gas,boundary);end
    end
    function reject(action)
        try,action();catch,refusals=refusals+1;return;end;error('ArithmeticTest:Refusal','非法选项未拒绝。');
    end
end
