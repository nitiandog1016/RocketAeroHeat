function report=testAxisymmetricVariableJacobian(folder)
%TESTAXISYMMETRICVARIABLEJACOBIAN 变mu/k的全耦合Jacobian与独立逐列/方向差分。
% 模块一：非均匀T、径向速度、偏斜周期网格，实际激活系数温度依赖。
% 模块二：冻结尺度，不差分CFL；逐列中心差分检查有色稀疏支持及布局。
% 模块三：保留有单位的积分导数账，不能把径向动量当成无源量。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v64_validation','variable_jacobian');end
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 800], ...
    'transport',struct('model','sutherland','referenceViscosity',.5,'referenceTemperature',300,'sutherlandTemperature',110.4));
[X,R]=ndgrid([0 .02 .05 .08],linspace(.03,.06,5));
mesh=buildAxisymmetricFVGeometry(X+.07*(R-.03),R,struct('periodicI',true));n=mesh.cellCount;t=(1:n)';
Q=axisymmetricFVConservative([1+.02*sin(t),-.5+.05*cos(t),.02*sin(.8*t),500+60*cos(.4*t)],gas);
scale=[1 20 20 2000];U=Q./scale;per=struct('type','periodic');
boundary=struct('iMin',per,'iMax',per, ...
    'jMin',struct('type','wall','temperature',350,'velocity',[-2 0]), ...
    'jMax',struct('type','wall','temperature',650,'velocity',[0 0]));
[J,details,pattern]=axisymmetricFVJacobian(Q,mesh,gas,boundary,scale, ...
    struct('differenceStep',2e-6,'differenceScheme','central'));
assert(details.success);
dense=zeros(4*n);h=2e-6*max(1,abs(U(:)));
for column=1:4*n
    plus=U;minus=U;plus(column)=plus(column)+h(column);minus(column)=minus(column)-h(column);
    rp=axisymmetricNavierStokesResidual(plus.*scale,mesh,gas,boundary)./scale;
    rm=axisymmetricNavierStokesResidual(minus.*scale,mesh,gas,boundary)./scale;
    dense(:,column)=(rp(:)-rm(:))/(2*h(column));
end
relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');assert(relative<2e-7);
outside=~full(kron(ones(4),pattern.cellPattern));assert(max([0;abs(dense(outside))])<1e-6);
v=reshape(sin((1:4*n)'*.73),n,4);epsilon=1e-6;
[rp,ip]=axisymmetricNavierStokesResidual((U+epsilon*v).*scale,mesh,gas,boundary);
[rm,im]=axisymmetricNavierStokesResidual((U-epsilon*v).*scale,mesh,gas,boundary);
action=reshape(J*v(:),n,4);direct=(rp-rm)./scale/(2*epsilon);
directionRelative=max(abs(action-direct),[],'all')/max(abs(direct),[],'all');assert(directionRelative<4e-6);
integrated=sum(action.*scale.*mesh.volumePerRadian,1);
expected=(-(ip.netBoundaryFlux-im.netBoundaryFlux)+(ip.sourceIntegral-im.sourceIntegral))/(2*epsilon);
ledgerRelative=max(abs(integrated-expected)./max(1,abs(expected)));assert(ledgerRelative<3e-6);
assert(max(ip.transport.faceMu)-min(ip.transport.faceMu)>.1&& ...
    max(abs(ip.transport.faceMu-im.transport.faceMu))>1e-8, ...
    'AxisFVTest:VariableJacobianExcitation','测试没有激活真实局部输运温度导数。');
report=struct('status','passed-variable-transport-full-Jacobian-support-and-direction', ...
    'cellCount',n,'coloredVsSerialRelative',relative,'directionRelative',directionRelative, ...
    'integratedDerivativeRelative',ledgerRelative,'colorCount',details.colorCount, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'variable_jacobian_evidence.mat'),'report','Q','gas','mesh','boundary','scale','J','dense','v','details','pattern');
fid=fopen(fullfile(folder,'variable_jacobian_summary.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Variable-transport Jacobian passed: serial %.3g, direction %.3g.\n',relative,directionRelative);
end
