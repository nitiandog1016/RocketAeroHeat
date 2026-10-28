function report=testAxisymmetricPressureSourceQuadrature(folder,oldFolder)
%TESTAXISYMMETRICPRESSURESOURCEQUADRATURE 独立积分、常压平衡、全耦合差分和真实旧端点回放。
% 模块一：凸偏斜轴网格与环带逐格解析线性压力积分，区分体/平面质心。
% 模块二：四/五量真实残差只改变径向压力源，源账、墙热及共享通量保持。
% 模块三：完整有色/串行Jacobian和真实小BE检验实际新方程，不以积分通过代替收敛。
% 模块四：只读逐值回放v70十个旧info，明确新源码不冒充旧链执行源或可原地续算。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v72_validation','pressure_quadrature');end
if nargin<2,oldFolder=fullfile(root,'outputs','v70_validation','normal_collar_continued');end
assert(~exist(folder,'dir'),'PressureQuadTest:Overwrite','不能覆盖已发生的测试证据。');
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[200 800], ...
    'transport',struct('model','sutherland','referenceViscosity',.05,'referenceTemperature',300,'sutherlandTemperature',110.4)));
quad=struct([]);grids={};orders=struct;
for lower=[0 .03]
    errors=[];
    for n=[4 8 16 32]
        [x,r]=ndgrid(linspace(0,.12,n+1),linspace(lower,.2,n+1));
        m=buildAxisymmetricFVGeometry(x+.2*r,r);p=300+.2*m.volumeCentroid(:,1)+.4*m.volumeCentroid(:,2);
        a=.2*ones(m.cellCount,1);b=.4*ones(m.cellCount,1);
        [linear,la]=axisymmetricPressureSourceQuadrature(m,p,a,b,'linear-planar');
        [constant,ca]=axisymmetricPressureSourceQuadrature(m,p,a,b,'cell-constant');
        exact=300+.2*m.planarCentroid(:,1)+.4*m.planarCentroid(:,2);
        [planarIntegral,weightedIntegral]=integrateAffine(m,[300 .2 .4]);
        % 从原顶点用独立Duffy/Gauss积分，不用保存质心反过来当积分真值。
        assert(max(abs(planarIntegral-exact.*m.planarArea))<3e-13&& ...
            max(abs(weightedIntegral./m.volumePerRadian-p))<3e-10);
        linearError=max(abs(linear-exact));assert(linearError<512*eps(400));
        assert(isequal(constant,p)&&isequal(ca.rawLinearPlanarPressurePa,la.rawLinearPlanarPressurePa)&& ...
            ~ca.adaptivePolicyFallback&&~ca.negativePressureClipped);
        error=sum(abs(constant.*m.planarArea-planarIntegral))/sum(m.planarArea);errors(end+1)=error; %#ok<AGROW>
        item=struct('minimumRadiusM',lower,'gridN',n,'cellCount',m.cellCount, ...
            'affineLinearMaximumErrorPa',linearError,'constantMeanAbsoluteIntegralErrorPa',error);
        if isempty(quad),quad=item;else,quad(end+1)=item;end %#ok<AGROW>
        if n==4,grids{end+1}=struct('mesh',m,'p',p,'gx',a,'gr',b,'linear',linear,'constant',constant, ...
                'exactPlanarIntegral',planarIntegral,'exactRWeightedIntegral',weightedIntegral);end %#ok<AGROW>
    end
    rates=log2(errors(1:end-1)./errors(2:end));assert(all(rates>1.5));
    if lower==0,orders.axis=rates;else,orders.annulus=rates;end
end

%% 二、常压在真实偏斜轴域保持；负外推不会被新策略裁成正值
m=grids{1}.mesh;restP=repmat([1 0 0 300],m.cellCount,1);
trace=struct('type','prescribedTrace','primitive',[1 0 0 300]);
bc=struct('iMin',trace,'iMax',trace,'jMin',struct('type','axis','symmetryGradient','mirror-parity'),'jMax',trace);
restQ=axisymmetricFVConservative(restP,gas);constantBC=bc;constantBC.pressureSourceQuadrature='cell-constant';
[rest,restInfo]=axisymmetricNavierStokesResidual(restQ,m,gas,constantBC);
roundBudget=2048*eps*max(abs(restInfo.faceFluxIntegralPerRadian),[],'all')/min(m.volumePerRadian);
assert(max(abs(rest),[],'all')<=roundBudget&&max(abs(rest),[],'all')<5e-8);
badGradient=1e6*ones(m.cellCount,1);one=ones(m.cellCount,1);
[positive,negativeAudit]=axisymmetricPressureSourceQuadrature(m,one,zeros(size(one)),badGradient,'cell-constant');
assert(isequal(positive,one)&&negativeAudit.rawNonpositiveCount>0&&~negativeAudit.adaptivePolicyFallback);
refusals=struct('label',{},'identifier',{});
reject('unbounded-linear-negative',@()axisymmetricPressureSourceQuadrature(m,one,zeros(size(one)),badGradient,'linear-planar'));
reject('unknown-policy',@()axisymmetricPressureSourceQuadrature(m,one,one,one,'clip'));
reject('nontext-policy',@()axisymmetricPressureSourceQuadrature(m,one,one,one,1));
reject('negative-cell-pressure',@()axisymmetricPressureSourceQuadrature(m,-one,one,one,'cell-constant'));
reject('complex-gradient',@()axisymmetricPressureSourceQuadrature(m,one,complex(one,one),one,'cell-constant'));
reject('nonfinite-gradient',@()axisymmetricPressureSourceQuadrature(m,one,Inf*one,one,'cell-constant'));
sideBad=bc;sideBad.jMin.pressureSourceQuadrature='cell-constant';
reject('silently-ignored-side-policy',@()axisymmetricNavierStokesResidual(restQ,m,gas,sideBad));

%% 三、同非均匀Q的四/五量真实方程，及两类全量Jacobian与BE
[x,r]=ndgrid([0 .02 .05 .08],linspace(.03,.06,5));
mesh=buildAxisymmetricFVGeometry(x+.07*(r-.03),r,struct('periodicI',true));n=mesh.cellCount;t=(1:n)';
P=[1+.02*sin(t),-2+.2*cos(t),.02*sin(.8*t),500+60*cos(.4*t)];
Q4=axisymmetricFVConservative(P,gas);Q5=[Q4,P(:,1).*(.3+.02*sin(.6*t))];
per=struct('type','periodic');legacyBC=struct('iMin',per,'iMax',per, ...
    'jMin',struct('type','wall','temperature',350,'velocity',[-2 0]), ...
    'jMax',struct('type','wall','temperature',650,'velocity',[0 0]));
newBC=legacyBC;newBC.pressureSourceQuadrature='cell-constant';
scale5=[1 20 20 2000 .5];records=struct([]);cases={};
for nq=[4 5]
    if nq==4,Q=Q4;else,Q=Q5;end
    scale=scale5(1:nq);U=Q./scale;
    [oldRate,oldInfo]=evaluate(Q,legacyBC);[newRate,newInfo]=evaluate(Q,newBC);
    assert(isequal(oldInfo.faceFluxIntegralPerRadian,newInfo.faceFluxIntegralPerRadian)&& ...
        isequaln(oldInfo.wall,newInfo.wall)&&isequal(oldInfo.sourceHoopStress,newInfo.sourceHoopStress)&& ...
        isequal(newInfo.cellPressureAtPlanarCentroid,newInfo.primitive(:,5)));
    others=setdiff(1:nq,3);assert(isequal(newRate(:,others),oldRate(:,others)));
    assert(max(abs(newRate(:,3)-oldRate(:,3)))>1e-3);
    ledger=sum(newRate.*mesh.volumePerRadian,1)+newInfo.netBoundaryFlux-newInfo.sourceIntegral;
    budget=2048*eps*max(1,max(abs(newInfo.netBoundaryFlux)+abs(newInfo.sourceIntegral)));
    assert(max(abs(ledger))<=budget);
    [J,details,pattern]=axisymmetricFVJacobian(Q,mesh,gas,newBC,scale,struct('differenceStep',2e-6));assert(details.success);
    dense=zeros(nq*n);
    for column=1:nq*n
        h=2e-6*max(1,abs(U(column)));up=U;um=U;up(column)=up(column)+h;um(column)=um(column)-h;
        rp=evaluate(up.*scale,newBC)./scale;rm=evaluate(um.*scale,newBC)./scale;
        dense(:,column)=(rp(:)-rm(:))/(2*h);
    end
    relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');assert(relative<3e-7);
    outside=~full(kron(ones(nq),pattern.cellPattern));assert(max([0;abs(dense(outside))])<1e-6);
    options=struct('maxSteps',1,'maximumPseudoStep',1e-8,'residualScale',scale,'steadyResidualTolerance',0,'jacobianRefreshEvery',1);
    be=advanceAxisymmetricSteady(Q,mesh,gas,newBC,options);assert(be.acceptedSteps==1&&~be.steadyConverged);
    [endRate,endInfo]=evaluate(be.Q,newBC);defect=be.Q./scale-U-be.history.dt(end)*(endRate./scale);
    assert(isequaln(endInfo,be.finalInfo)&&max(abs(defect),[],'all')<=be.history.newtonTarget(end)+1024*eps&& ...
        be.globalBalanceRelative<=5e-10);
    item=struct('quantities',nq,'jacobianColoredVsSerialRelative',relative,'ledger',ledger, ...
        'beResidualEvaluations',be.residualEvaluations,'beBalanceRelative',be.globalBalanceRelative);
    if isempty(records),records=item;else,records(end+1)=item;end %#ok<AGROW>
    cases{end+1}=struct('Q',Q,'mesh',mesh,'boundary',newBC,'oldInfo',oldInfo,'newInfo',newInfo, ...
        'newRate',newRate,'oldRate',oldRate,'J',J,'serial',dense,'be',be); %#ok<AGROW>
end

%% 四、旧默认完整info逐值回放；只验证相容数值，不重绑旧源码或原地续算
saved=load(fullfile(oldFolder,'rocket_axisymmetric_evidence.mat'),'caseDefinition','durable','implicit');c=saved.caseDefinition;
assert(saved.implicit.acceptedSteps==10&&~isfield(c.boundary,'pressureSourceQuadrature'));
for k=1:10
    payload=load(fullfile(oldFolder,saved.durable.steps(k).relativePath),'segment');s=payload.segment;
    [rate,info]=axisymmetricRANSResidual(s.Q,c.mesh,c.gas,c.boundary);
    assert(isequaln(info,s.finalInfo)&&max(abs(rate)./s.residualScale,[],'all')==s.residualNorm);
end
report=struct('status','passed-explicit-pressure-source-quadrature-numerical-contract', ...
    'matlabEngineProcessID',feature('getpid'),'affineQuadrature',quad,'integralConvergenceOrders',orders, ...
    'skewAxisRestMaximum',max(abs(rest),[],'all'),'restRoundoffBudget',roundBudget, ...
    'onlyRadialPressureSourceChanged',true,'coupledCases',records,'refusals',refusals, ...
    'legacyTenCompleteInfoReplayed',true,'oldChainSourceRebound',false, ...
    'sourcePressureClipped',false,'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'pressure_quadrature_evidence.mat'),'report','grids','gas','cases','negativeAudit','restInfo','-v7');
fid=fopen(fullfile(folder,'pressure_quadrature_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function [rate,info]=evaluate(state,bound)
        if size(state,2)==5,[rate,info]=axisymmetricRANSResidual(state,mesh,gas,bound);
        else,[rate,info]=axisymmetricNavierStokesResidual(state,mesh,gas,bound);end
    end
    function reject(label,action)
        try,action();catch exception
            assert(startsWith(exception.identifier,'AxisFV:'),'PressureQuadTest:WrongRefusal','非合同错误%s。',exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('PressureQuadTest:NoRefusal','非法输入%s没有被拒绝。',label);
    end
end

function [areaIntegral,volumeIntegral]=integrateAffine(m,coefficient)
% 从原顶点三角化，以三点Gauss×三点Gauss Duffy积分p及r*p；不调用保存质心。
z=[.1127016653792583 .5 .8872983346207417];w=[5 8 5]/18;
areaIntegral=zeros(m.cellCount,1);volumeIntegral=areaIntegral;
for index=1:m.cellCount
    p=m.vertices(m.cells(index,:),:);
    for tri={[1 2 3],[1 3 4]}
        v=p(tri{1},:);a=v(2,:)-v(1,:);b=v(3,:)-v(1,:);determinant=abs(det([a;b]));
        for i=1:3
            for j=1:3
                point=v(1,:)+z(i)*a+(1-z(i))*z(j)*b;
                weight=w(i)*w(j)*(1-z(i))*determinant;
                pressure=coefficient(1)+coefficient(2)*point(1)+coefficient(3)*point(2);
                areaIntegral(index)=areaIntegral(index)+weight*pressure;
                volumeIntegral(index)=volumeIntegral(index)+weight*point(2)*pressure;
            end
        end
    end
end
end
