function report=testAxisymmetricLocalJacobian(folder)
%TESTAXISYMMETRICLOCALJACOBIAN 验证局部物理扰动的四/五量完整导数，而非镜像实现。
% 模块一：跨密度/SA量级及严格零SA，独立逐列残差差分检查有色支持与解码。
% 模块二：物理扰动不随任意归一化尺度变化；导数变换与局部方向割线分别检查。
% 模块三：非法输入明确拒绝；不推导一般正性或箭体热流精度。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v74_validation','local_unit');end
assert(~exist(folder,'dir'),'LocalJacobianTest:Overwrite','不能覆盖证据。');
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[200 800], ...
    'transport',struct('model','sutherland','referenceViscosity',.05,'referenceTemperature',300,'sutherlandTemperature',110.4)));
[x,r]=ndgrid([0 .02 .05 .08],linspace(.03,.06,5));mesh=buildAxisymmetricFVGeometry(x+.07*(r-.03),r,struct('periodicI',true));
n=mesh.cellCount;t=(1:n)';rho=10.^(-4*(1+sin(.7*t))/2);
P=[rho,-2+.2*cos(t),.02*sin(t),500+20*cos(.4*t)];Q4=axisymmetricFVConservative(P,gas);
W=rho.*10.^(-9+7*(1+cos(t))/2);W(1)=0;Q5=[Q4 W];
per=struct('type','periodic');wall=struct('type','wall','temperature',500,'velocity',[0 0]);
bc=struct('iMin',per,'iMax',per,'jMin',wall,'jMax',wall,'pressureSourceQuadrature','cell-constant','interiorConvectiveFlux','hlle-einfeldt');
records=struct([]);evidence={};h=2e-6;scale5=[1 20 20 2000 .5];
for nq=[4 5]
    if nq==4,Q=Q4;else,Q=Q5;end;scale=scale5(1:nq);opts=struct('differenceStep',h,'differenceScheme','central');
    [J,details,pattern]=axisymmetricFVJacobianLocal(Q,mesh,gas,bc,scale,opts);assert(details.success);
    % 独立物理尺度：声速解析式与rho/速度构造，不调用参考尺度函数。
    speed=max(sqrt(gas.gamma*gas.R*P(:,4)),hypot(P(:,2),P(:,3)));
    physical=[rho,rho.*speed,rho.*speed,abs(Q(:,4))];
    if nq==5
        mu=.05*(P(:,4)/300).^1.5*(300+110.4)./(P(:,4)+110.4);
        physical(:,5)=max(W,64*eps*mu/h);
    end
    ref=physical./scale;near(details.localScaling.physicalReference,physical,128*eps);
    base=evaluate(Q)./scale;dense=zeros(n*nq);U=Q./scale;
    for col=1:n*nq
        plus=Q;minus=Q;step=h*ref(col);plus(col)=(U(col)+step)*scale(ceil(col/n));
        hp=(plus(col)-Q(col))/scale(ceil(col/n));rp=evaluate(plus)./scale;
        if nq==5&&col==4*n+1
            dense(:,col)=(rp(:)-base(:))/hp;
        else
            minus(col)=(U(col)-step)*scale(ceil(col/n));hm=(Q(col)-minus(col))/scale(ceil(col/n));
            rm=evaluate(minus)./scale;weight=hm/(hp+hm);
            dense(:,col)=weight*(rp(:)-base(:))/hp+(1-weight)*(base(:)-rm(:))/hm;
        end
    end
    relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');assert(relative<2e-7);
    assert(all(details.columnSchemes~=0,'all')&&details.columnsDifferentiated==n*nq);
    if nq==5,assert(details.columnSchemes(1,5)==1&&isnan(details.acceptedNegativeNormalizedStep(1,5)));end
    % 用与逐列扰动不同的方向/步长校核，零SA方向设零以保留中央方向的可行性。
    direction=ref.*reshape(sin((1:n*nq)*.63),n,nq);if nq==5,direction(1,5)=0;end
    hd=2e-5;rp=evaluate((U+hd*direction).*scale)./scale;rm=evaluate((U-hd*direction).*scale)./scale;
    direct=(rp-rm)/(2*hd);action=reshape(J*direction(:),n,nq);
    directionError=max(abs(action-direct),[],'all')/max(1,max(abs(direct),[],'all'));assert(directionError<2e-5);
    factor=[.7 3 .1 2 4];otherScale=scale.*factor(1:nq);
    [~,other]=axisymmetricJacobianLocalReference(Q,gas,otherScale,h);near(other.physicalReference,physical,128*eps);
    [K,changed]=axisymmetricFVJacobianLocal(Q,mesh,gas,bc,otherScale,opts,pattern);assert(changed.success);
    transform=reshape(repmat(factor(1:nq),n,1),[],1);expected=(full(J).*transform')./transform;
    scaleError=max(abs(full(K)-expected),[],'all')/max(1,max(abs(expected),[],'all'));assert(scaleError<2e-5);
    row=struct('quantities',nq,'cellCount',n,'coloredVsIndependentSerialRelative',relative, ...
        'independentDirectionRelative',directionError,'scaleTransformRelative',scaleError, ...
        'centralForwardBackwardColumns',[nnz(details.columnSchemes==2),nnz(details.columnSchemes==1),nnz(details.columnSchemes==-1)], ...
        'residualEvaluations',details.residualEvaluations);
    if isempty(records),records=row;else,records(end+1)=row;end %#ok<AGROW>
    evidence{end+1}=struct('J',J,'dense',dense,'K',K,'Q',Q,'scale',scale,'details',details,'pattern',pattern); %#ok<AGROW>
end
refusals=0;bad=Q5;bad(1,5)=-1;reject(@()axisymmetricJacobianLocalReference(bad,gas,scale5,h));
reject(@()axisymmetricJacobianLocalReference(Q5,gas,scale5,0));badScale=scale5;badScale(2)=0;
reject(@()axisymmetricJacobianLocalReference(Q5,gas,badScale,h));
report=struct('status','passed-local-conservative-four-five-quantity-colored-Jacobian-contract', ...
    'matlabEngineProcessID',feature('getpid'),'cases',records,'refusals',refusals, ...
    'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'local_jacobian_evidence.mat'),'report','evidence','mesh','gas','bc','-v7');
fid=fopen(fullfile(folder,'local_jacobian_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function rate=evaluate(state)
        if size(state,2)==5,rate=axisymmetricRANSResidual(state,mesh,gas,bc);else,rate=axisymmetricNavierStokesResidual(state,mesh,gas,bc);end
    end
    function reject(action)
        try,action();catch,refusals=refusals+1;return;end;error('LocalJacobianTest:Refusal','非法参数未拒绝。');
    end
end
function near(a,b,tol)
assert(isequal(size(a),size(b))&&all(abs(a-b)<=tol*max(abs(b),realmin),'all'));
end
