function report=testAxisymmetricFVJacobian(folder)
%TESTAXISYMMETRICFVJACOBIAN 核对原残差稀疏有色Jacobian，非制造解求解。
% 模块一：有色列与逐列中心差分独立对照，检查两层梯度模板和分量布局。
% 模块二：随机方向及边界/环向源的积分导数独立对照，不把径向动量当无源。
% 模块三：周期偏斜几何、复用身份和非法扰动合同均显式检查。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','jacobian');end
if ~exist(folder,'dir'),mkdir(folder);end
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
[X,R]=ndgrid(linspace(.06,0,6),linspace(.03,.06,6));
mesh=buildAxisymmetricFVGeometry(X+.07*(R-.03),R,struct('periodicI',true));
n=mesh.cellCount;t=(1:n)';
primitive=[1+.03*sin(t),-1-.1*cos(.7*t),.01*sin(.9*t),300+.3*cos(t)];
Q=axisymmetricFVConservative(primitive,gas);scale=[1 25 25 800];U=Q./scale;
per=struct('type','periodic');wall=struct('type','wall','temperature',300,'velocity',[0 0]);
boundary=struct('iMin',per,'iMax',per,'jMin',wall,'jMax',wall);
opts=struct('differenceStep',2e-6,'differenceScheme','central');
[J,details,pattern]=axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,opts);
assert(details.success&&issparse(J)&&isequal(size(J),[4*n 4*n])&&all(isfinite(nonzeros(J))));
assert(details.residualEvaluations>0&&details.colorCount<=n);
dense=zeros(4*n);h=opts.differenceStep*max(1,abs(U(:)));
for column=1:4*n
    plus=U;minus=U;plus(column)=plus(column)+h(column);minus(column)=minus(column)-h(column);
    rp=axisymmetricNavierStokesResidual(plus.*scale,mesh,gas,boundary)./scale;
    rm=axisymmetricNavierStokesResidual(minus.*scale,mesh,gas,boundary)./scale;
    dense(:,column)=(rp(:)-rm(:))/(2*h(column));
end
relative=max(abs(full(J)-dense),[],'all')/max(abs(dense),[],'all');
assert(relative<2e-7,'有色Jacobian漏模板、布局或解码错误。');
outside=~full(kron(ones(4),pattern.cellPattern));
assert(max(abs(dense(outside)))<1e-6,'真实残差在保守pattern之外仍有导数。');
v=reshape(sin((1:4*n)'*.63),n,4);epsDirection=2e-6;
[rp,ip]=axisymmetricNavierStokesResidual((U+epsDirection*v).*scale,mesh,gas,boundary);
[rm,im]=axisymmetricNavierStokesResidual((U-epsDirection*v).*scale,mesh,gas,boundary);
action=reshape(J*v(:),n,4);direct=(rp-rm)./scale/(2*epsDirection);
actionRelative=max(abs(action-direct),[],'all')/max(abs(direct),[],'all');
assert(actionRelative<3e-6,'Jacobian作用与独立方向差分不符。');
integrated=sum(action.*scale.*mesh.volumePerRadian,1);
boundarySource=(-(ip.netBoundaryFlux-im.netBoundaryFlux)+(ip.sourceIntegral-im.sourceIntegral))/(2*epsDirection);
ledgerRelative=max(abs(integrated-boundarySource)./max(1,abs(boundarySource)));
assert(ledgerRelative<2e-6,'积分导数未保留边界与柱坐标源。');
[J2,repeat]=axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,opts,pattern);
assert(repeat.success&&max(abs(full(J-J2)),[],'all')<1e-8);
rejections=0;
bad=scale;bad(1)=0;reject(@()axisymmetricFVJacobian(Q,mesh,gas,boundary,bad,opts));
bad=opts;bad.differenceStep=-1;reject(@()axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,bad));
bad=opts;bad.differenceScheme='complex';reject(@()axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,bad));
bad=pattern;bad.cellPattern=speye(n);reject(@()axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,opts,bad));
bad=pattern;bad.colors=ones(n,1);reject(@()axisymmetricFVJacobian(Q,mesh,gas,boundary,scale,opts,bad));
report=struct('status','passed-colored-full-residual-Jacobian-only','cellCount',n, ...
    'coloredVsSerialRelative',relative,'directionRelative',actionRelative, ...
    'boundarySourceDerivativeRelative',ledgerRelative,'colorCount',details.colorCount, ...
    'residualEvaluations',details.residualEvaluations,'rejections',rejections, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'jacobian_evidence.mat'),'Q','scale','mesh','gas','boundary','J','dense','v','details','pattern','report');
fid=fopen(fullfile(folder,'jacobian_summary.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Axisymmetric full-residual colored Jacobian passed: %d cells/%d colors, serial %.3g, direction %.3g, %d refusals.\n',n,details.colorCount,relative,actionRelative,rejections);
    function reject(action)
        try,action();catch,rejections=rejections+1;return;end
        error('AxisSteadyTest:MissingRefusal','非法Jacobian合同没有拒绝。');
    end
end
