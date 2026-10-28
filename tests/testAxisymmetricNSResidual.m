function report=testAxisymmetricNSResidual(folder)
%TESTAXISYMMETRICNSRESIDUAL 常压静止、负X平移、非零hoop应力及完整源账。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v62_validation','residual');end
if ~exist(folder,'dir'),mkdir(folder);end
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
[X,R]=ndgrid([.2,.1,0],[0,.01,.03,.06]);
mesh=buildAxisymmetricFVGeometry(X,R);nc=mesh.cellCount;P=repmat([1 0 0 300],nc,1);
wall=@(u)struct('type','wall','temperature',300,'velocity',[u 0]);
trace=@(value)struct('type','prescribedTrace','primitive',value);
bc=struct('iMin',trace([1 0 0 300]),'iMax',trace([1 0 0 300]), ...
    'jMin',struct('type','axis'),'jMax',wall(0));
[rest,restInfo]=axisymmetricNavierStokesResidual(axisymmetricFVConservative(P,gas),mesh,gas,bc);
assert(max(abs(rest),[],'all')<2e-8,'常压静止场未保持，可能缺少径向压力几何源。');
assert(max(abs(restInfo.wall.heatFluxIntoWall))<1e-8); % 原始量/守恒量往返的浮点舍入

%% 真实负X整体平移；周期通量与同速墙不能产生虚假剪切或热
[X0,R0]=ndgrid([.1,0],[.03,.04,.06]);periodic=buildAxisymmetricFVGeometry(X0,R0,struct('periodicI',true));
pb=repmat([1 -1200 0 300],periodic.cellCount,1);per=struct('type','periodic');
periodBC=struct('iMin',per,'iMax',per,'jMin',wall(-1200),'jMax',wall(-1200));
[uniform,uniformInfo]=axisymmetricNavierStokesResidual(axisymmetricFVConservative(pb,gas),periodic,gas,periodBC);
assert(max(abs(uniform),[],'all')<2e-8&&max(abs(uniformInfo.wall.heatFluxIntoWall))<1e-7);
assert(max(abs(uniformInfo.wall.viscousTraction),[],'all')<1e-8);

%% 非零径向线性速度：两种mu相减，独立消去相同一阶对流截断误差
% 3D Cartesian真值为diag(0,a,a)，故tau_rr=tau_theta=2mu*a/3，
% 黏性力严格0，黏性功散度/r体积=4mu*a^2/3，不添加制造源。
a=.2;linearP=[ones(nc,1),-7*ones(nc,1),a*mesh.volumeCentroid(:,2),300*ones(nc,1)];
exact=@(point)[1 -7 a*point(2) 300];fullBC=struct('iMin',trace(exact),'iMax',trace(exact),'jMin',trace(exact),'jMax',trace(exact));
Q=axisymmetricFVConservative(linearP,gas);
[r1,i1]=axisymmetricNavierStokesResidual(Q,mesh,gas,fullBC);gas2=gas;gas2.mu=1;
[r2,i2]=axisymmetricNavierStokesResidual(Q,mesh,gas2,fullBC);difference=r2-r1;
assert(max(abs(i1.cellVelocityGradientR(:,2)-a))<1e-12);
assert(max(abs(i1.sourceHoopStress-2*gas.mu*a/3))<1e-12);
assert(max(abs(difference(:,2:3)),[],'all')<2e-8);
assert(max(abs(difference(:,4)-4*(gas2.mu-gas.mu)*a^2/3))<2e-7);
assert(max(abs(difference(:,1)))==0);
for current={i1,i2}
    info=current{1};assert(max(abs(info.conservationDefect)./max(1,abs(info.sourceIntegral)))<1e-11);
end

%% 完整非均匀四量共享通量与径向源账；不得将源误算为外边界通量
variable=linearP;variable(:,1)=1+.02*sin((1:nc)');variable(:,4)=300+2*cos((1:nc)');
[nonuniform,nonuniformInfo]=axisymmetricNavierStokesResidual(axisymmetricFVConservative(variable,gas),mesh,gas,fullBC);
ledger=sum(nonuniform.*mesh.volumePerRadian,1)+nonuniformInfo.netBoundaryFlux-nonuniformInfo.sourceIntegral;
assert(max(abs(ledger))<1e-10);

%% 无效物性、状态、轴/周期及周向速度身份须被拒绝
rejections=0;bad=Q;bad(1,1)=-1;reject(@()axisymmetricNavierStokesResidual(bad,mesh,gas,fullBC));
reject(@()axisymmetricNavierStokesResidual(complex(Q,1),mesh,gas,fullBC));
bad=gas;bad.mu=-1;reject(@()axisymmetricNavierStokesResidual(Q,mesh,bad,fullBC));
bad=gas;bad.k=1;reject(@()axisymmetricNavierStokesResidual(Q,mesh,bad,fullBC));
bad=mesh;bad.volumePerRadian(1)=NaN;reject(@()axisymmetricNavierStokesResidual(Q,bad,gas,fullBC));
bad=fullBC;bad.jMax=struct('type','axis');reject(@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad));
bad=fullBC;bad.iMin=per;reject(@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad));
bad=fullBC;bad.iMin=struct('type','supersonicInflow','primitive',[1 0 0 300]);reject(@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad));
bad=fullBC;bad.jMax=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',[1 0 0 1 300]);
reject(@()axisymmetricNavierStokesResidual(Q,mesh,gas,bad));
skew=buildAxisymmetricFVGeometry(X+.11*R,R);reject(@()axisymmetricNavierStokesResidual(Q,skew,gas,bc));
report=struct('status','passed-numerical-operator-only','restMaximum',max(abs(rest),[],'all'), ...
    'negativeXUniformMaximum',max(abs(uniform),[],'all'),'nonzeroHoopValue',2*gas.mu*a/3, ...
    'viscousHeatingDifferenceMaxError',max(abs(difference(:,4)-4*(gas2.mu-gas.mu)*a^2/3)), ...
    'nonuniformLedger',ledger,'rejections',rejections,'physicalSolution',false,'newSolvedAreaM2',0);
save(fullfile(folder,'residual_evidence.mat'),'mesh','Q','r1','r2','i1','i2','nonuniform','nonuniformInfo','periodic','uniform','uniformInfo','report');
fid=fopen(fullfile(folder,'residual_summary.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Axisymmetric NS residual passed: rest/negativeX, nonzero hoop stress and energy, shared source ledger, %d refusals.\n',rejections);
    function reject(action)
        try,action();catch,rejections=rejections+1;return;end
        error('AxisFVTest:MissingRefusal','无效合同未被拒绝。');
    end
end
