function report=testAxisymmetricSAJacobian(folder)
%TESTAXISYMMETRICSAJACOBIAN 真五量耦合有色Jacobian的逐列、方向和边界右导数校核。
% 模块一：偏斜周期环隙中非均匀rho/u/T/N激活完整NS与SA交叉导数。
% 模块二：逐列中心差分独立检查25块布局、稀疏支持及整体方向导数。
% 模块三：专门验证W对动量/能量及rho/E对SA的交叉块，不能退化分裂更新。
% 模块四：N=0负扰动必须拒绝，交付明确右导数并记入失败调用。
% 不作长稳态求解或实验验收；每个尺度开跑前冻结，CFL不进入导数。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v65_validation','sa_jacobian');end
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[200 800], ...
    'transport',struct('model','sutherland','referenceViscosity',.05,'referenceTemperature',300,'sutherlandTemperature',110.4)));
[X,R]=ndgrid([0 .02 .05 .08],linspace(.03,.06,5));
mesh=buildAxisymmetricFVGeometry(X+.07*(R-.03),R,struct('periodicI',true));n=mesh.cellCount;t=(1:n)';
P=[1+.02*sin(t),-2+.2*cos(t),.02*sin(.8*t),500+60*cos(.4*t)];N=.3+.02*sin(.6*t);
Q=[axisymmetricFVConservative(P,gas),P(:,1).*N];scale=[1 20 20 2000 .5];U=Q./scale;
per=struct('type','periodic');boundary=struct('iMin',per,'iMax',per, ...
    'jMin',struct('type','wall','temperature',350,'velocity',[-2 0]), ...
    'jMax',struct('type','wall','temperature',650,'velocity',[0 0]));
[J,details,pattern]=axisymmetricFVJacobian(Q,mesh,gas,boundary,scale, ...
    struct('differenceStep',2e-6,'differenceScheme','central'));
assert(details.success&&isequal(size(J),[5*n 5*n])&&details.columnsDifferentiated==5*n);
dense=zeros(5*n);h=2e-6*max(1,abs(U(:)));
for column=1:5*n
    plus=U;minus=U;plus(column)=plus(column)+h(column);minus(column)=minus(column)-h(column);
    rp=axisymmetricRANSResidual(plus.*scale,mesh,gas,boundary)./scale;
    rm=axisymmetricRANSResidual(minus.*scale,mesh,gas,boundary)./scale;
    dense(:,column)=(rp(:)-rm(:))/(2*h(column));
end
relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');assert(relative<3e-7);
outside=~full(kron(ones(5),pattern.cellPattern));outsideMax=max([0;abs(dense(outside))]);
assert(outsideMax<1e-6,'AxisSATest:JacobianSupport','五量真实依赖超出有色稀疏支持。');
blockRelative=zeros(5);blockAmplitude=zeros(5);
for output=1:5
    for input=1:5
        rr=(output-1)*n+(1:n);cc=(input-1)*n+(1:n);reference=dense(rr,cc);
        blockAmplitude(output,input)=max(abs(reference),[],'all');
        if blockAmplitude(output,input)>1e-5
            blockRelative(output,input)=max(abs(full(J(rr,cc))-reference),[],'all')/blockAmplitude(output,input);
            assert(blockRelative(output,input)<2e-5,'AxisSATest:JacobianBlock','第%d输出/第%d输入块不匹配。',output,input);
        end
    end
end
assert(all(blockAmplitude(2:4,5)>1e-3)&&all(blockAmplitude(5,[1 4])>1e-3)&&blockAmplitude(1,5)==0, ...
    'AxisSATest:JacobianCouplingExcitation','W↔NS交叉耦合未激活，不能以分裂/代数muT代替五量求导。');
v=reshape(sin((1:5*n)'*.73),n,5);epsilon=1e-6;
[rp,ip]=axisymmetricRANSResidual((U+epsilon*v).*scale,mesh,gas,boundary);
[rm,im]=axisymmetricRANSResidual((U-epsilon*v).*scale,mesh,gas,boundary);
action=reshape(J*v(:),n,5);direct=(rp-rm)./scale/(2*epsilon);
directionRelative=max(abs(action-direct),[],'all')/max(abs(direct),[],'all');assert(directionRelative<7e-6);
integrated=sum(action.*scale.*mesh.volumePerRadian,1);
expected=(-(ip.netBoundaryFlux-im.netBoundaryFlux)+(ip.sourceIntegral-im.sourceIntegral))/(2*epsilon);
derivativeScale=sum(abs(action.*scale.*mesh.volumePerRadian),1)+abs(expected);
derivativeDefect=integrated-expected;
assert(all(abs(derivativeDefect)<=5e-6*derivativeScale+256*eps*derivativeScale), ...
    'AxisSATest:JacobianLedger','五量方向导数积分与边界/源导数账不闭合。');

%% 四、N=0在非负域边界使用右导数，不能裁负状态构造虚假的中央差分
zeroQ=Q;zeroQ(:,5)=0;zeroU=zeroQ./scale;
[Jzero,zeroDetails]=axisymmetricFVJacobian(zeroQ,mesh,gas,boundary,scale, ...
    struct('differenceStep',2e-6,'differenceScheme','central','maxPerturbationRetries',1),pattern);
assert(zeroDetails.success&&all(zeroDetails.columnSchemes(:,5)==1)&& ...
    zeroDetails.oneSidedBatches>0&&zeroDetails.columnsDifferentiated==5*n, ...
    'AxisSATest:ZeroJacobianScheme','N=0没有交付明确第五量右导数。');
base=axisymmetricRANSResidual(zeroQ,mesh,gas,boundary)./scale;zeroDense=zeros(5*n,n);
for cellID=1:n
    plus=zeroU;plus(cellID,5)=2e-6;
    positive=axisymmetricRANSResidual(plus.*scale,mesh,gas,boundary)./scale;
    delta=(plus(cellID,5)*scale(5)-zeroQ(cellID,5))/scale(5);
    zeroDense(:,cellID)=(positive(:)-base(:))/delta;
end
zeroRelative=max(abs(full(Jzero(:,4*n+(1:n)))-zeroDense),[],'all')/max(abs(zeroDense),[],'all');
assert(zeroRelative<3e-7,'AxisSATest:ZeroJacobianValue','N=0第五量右导数不匹配显式正向独立扰动。');
if zeroDetails.failedResidualEvaluations>0
    idx=strcmp(zeroDetails.failureIdentifiers,'AxisSA:State');
    assert(any(idx)&&sum(zeroDetails.failureCounts(idx))>0,'AxisSATest:ZeroJacobianRefusal','负N扰动失败未按真实状态拒绝记录。');
end
report=struct('status','passed-SA-complete-five-variable-Jacobian-and-zero-right-derivative', ...
    'cellCount',n,'coloredVsSerialRelative',relative,'blockRelative',blockRelative,'blockAmplitude',blockAmplitude, ...
    'outsideSupportMaximum',outsideMax,'directionRelative',directionRelative, ...
    'integratedDerivativeDefect',derivativeDefect,'integratedDerivativeScale',derivativeScale, ...
    'zeroNuTildeRightDerivativeRelative',zeroRelative,'colorCount',details.colorCount, ...
    'zeroFailedResidualEvaluations',zeroDetails.failedResidualEvaluations, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'sa_jacobian_evidence.mat'),'report','Q','gas','mesh','boundary','scale','J','dense','v','details','pattern', ...
    'zeroQ','Jzero','zeroDense','zeroDetails');
fid=fopen(fullfile(folder,'sa_jacobian_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('SA Jacobian passed: complete 25 blocks, serial %.3g, direction %.3g, N0-right %.3g.\n',relative,directionRelative,zeroRelative);
end
