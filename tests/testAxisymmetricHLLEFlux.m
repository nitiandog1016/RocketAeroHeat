function report=testAxisymmetricHLLEFlux(folder)
%TESTAXISYMMETRICHLLEFLUX 实际四量通量的极限、旋转/迎风、正性及激波管解析收敛。
% 模块一：均匀、Sod/Lax、强激波、对向/双膨胀、近真空、接触/剪切和单向波区。
% 模块二：面反向守恒、二维坐标旋转、切向Galilean变换和独立标量HLL通量。
% 模块三：真实一阶Godunov-Sod演化与解析Riemann密度比较，逐步检查rho/内能。
% 本测试是Euler数值校核，不能证明NS/SA/几何源/隐式全局正性或箭体热流精度。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v73_validation','hlle_flux');end
assert(~exist(folder,'dir'),'HLLETest:Overwrite','不能覆盖实际测试证据。');
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.01,'Pr',.72,'temperatureRange',[1e-14 1e12]));
left=[1 3 .2 2;1 0 0 1;.445 .698 0 3.528;1 0 0 1000;1 2 0 .4;1 -2 0 .4; ...
    1 -10 0 1e-8;1 0 0 1;1 0 2 1;1 40 0 300;1 -200 0 300];
right=[left(1,:);.125 0 0 .1;.5 0 0 .571;1 0 0 .01;1 -2 0 .4;1 2 0 .4; ...
    1 10 0 1e-8;2 0 0 1;1 0 -2 1;1 200 0 300;1 -40 0 300];
labels={'uniform','Sod','Lax','strong-shock','head-on','double-rarefaction','near-vacuum', ...
    'stationary-contact','stationary-shear','all-right','all-left'};
QL=state(left);QR=state(right);normal=repmat([1 0],size(QL,1),1);records=struct([]);fixtures={};
for policy={'davis','einfeldt'}
    name=policy{1};[F,a]=axisymmetricHLLEFlux(QL,QR,normal,gas,name);
    assert(all(isfinite(F),'all')&&all(a.intermediateConservativeState(a.middleMask,1)>0)&& ...
        all(a.intermediateInternalEnergyPerVolume(a.middleMask)>0)&& ...
        a.sampleBranch(10)==1&&a.sampleBranch(11)==-1&& ...
        all(isnan(a.intermediateConservativeState(~a.middleMask,:)),'all'));
    P=axisymmetricFVState(QL,gas);physical=[P(:,1).*P(:,2),QL(:,2:3).*P(:,2)+P(:,5).*[ones(size(P,1),1) zeros(size(P,1),1)],(QL(:,4)+P(:,5)).*P(:,2)];
    assert(max(abs(F(1,:)-physical(1,:))./max(1,abs(physical(1,:))))<1e-12&&isequal(F(10,:),physical(10,:)));
    PR=axisymmetricFVState(QR,gas);rightFlux=[PR(:,1).*PR(:,2),QR(:,2:3).*PR(:,2)+PR(:,5).*[ones(size(PR,1),1) zeros(size(PR,1),1)],(QR(:,4)+PR(:,5)).*PR(:,2)];
    assert(isequal(F(11,:),rightFlux(11,:)));
    [reverse,~]=axisymmetricHLLEFlux(QR,QL,-normal,gas,name);
    antisym=max(abs(F+reverse)./max(1,abs(F)),[],'all');assert(antisym<3e-12);
    angle=.63;rotation=[cos(angle) -sin(angle);sin(angle) cos(angle)];
    rotatedL=QL;rotatedR=QR;rotatedL(:,2:3)=QL(:,2:3)*rotation';rotatedR(:,2:3)=QR(:,2:3)*rotation';
    transformed=axisymmetricHLLEFlux(rotatedL,rotatedR,normal*rotation',gas,name);
    expected=F;expected(:,2:3)=F(:,2:3)*rotation';
    rotateError=max(abs(transformed-expected)./max(1,abs(expected)),[],'all');assert(rotateError<3e-10);
    boost=[0 4];boostL=QL;boostR=QR;
    boostL(:,2:3)=QL(:,2:3)+QL(:,1).*boost;boostR(:,2:3)=QR(:,2:3)+QR(:,1).*boost;
    boostL(:,4)=QL(:,4)+sum(QL(:,2:3).*boost,2)+.5*QL(:,1)*sum(boost.^2);
    boostR(:,4)=QR(:,4)+sum(QR(:,2:3).*boost,2)+.5*QR(:,1)*sum(boost.^2);
    boosted=axisymmetricHLLEFlux(boostL,boostR,normal,gas,name);expectedBoost=F;
    expectedBoost(:,2:3)=F(:,2:3)+F(:,1).*boost;
    expectedBoost(:,4)=F(:,4)+sum(F(:,2:3).*boost,2)+.5*F(:,1)*sum(boost.^2);
    galilean=max(abs(boosted-expectedBoost)./max(1,abs(expectedBoost)),[],'all');assert(galilean<3e-10);
    scalar=zeros(size(F));
    for k=1:size(F,1)
        lo=a.waveSpeedsMPerS(k,1);hi=a.waveSpeedsMPerS(k,2);
        if lo>=0,scalar(k,:)=physical(k,:);elseif hi<=0,scalar(k,:)=rightFlux(k,:);
        else,scalar(k,:)=(hi*physical(k,:)-lo*rightFlux(k,:)+lo*hi*(QR(k,:)-QL(k,:)))/(hi-lo);end
    end
    scalarError=max(abs(F-scalar)./max(1,abs(F)),[],'all');assert(scalarError<3e-12);
    % 接触/剪切的耗散必须真实记录：新通量不能伪称保留了HLLC接触分辨能力。
    assert(abs(F(8,1))>0.1&&abs(F(9,3))>0.1);
    row=struct('wavePolicy',name,'antisymmetryRelative',antisym,'rotationRelative',rotateError, ...
        'tangentialGalileanRelative',galilean,'scalarReferenceRelative',scalarError, ...
        'stationaryContactArtificialMassFlux',F(8,1),'stationaryShearArtificialMomentumFlux',F(9,3), ...
        'minimumUsedIntermediateDensity',min(a.intermediateConservativeState(a.middleMask,1)), ...
        'minimumUsedIntermediateInternalEnergy',min(a.intermediateInternalEnergyPerVolume(a.middleMask)));
    if isempty(records),records=row;else,records(end+1)=row;end %#ok<AGROW>
    fixtures{end+1}=struct('policy',name,'flux',F,'audit',a); %#ok<AGROW>
end
refusals=struct('label',{},'identifier',{});
reject('unknown-waves',@()axisymmetricHLLEFlux(QL,QR,normal,gas,'guess'));
bad=normal;bad(1,:)=[2 0];reject('nonunit-normal',@()axisymmetricHLLEFlux(QL,QR,bad,gas,'einfeldt'));
bad=QL;bad(1,1)=-1;reject('negative-density',@()axisymmetricHLLEFlux(bad,QR,normal,gas,'einfeldt'));
bad=QL;bad(1,4)=-1;reject('negative-energy',@()axisymmetricHLLEFlux(bad,QR,normal,gas,'einfeldt'));
reject('mismatched-states',@()axisymmetricHLLEFlux(QL,QR(1,:),normal,gas,'einfeldt'));

%% 三、Sod真实时间推进；解析密度来自独立星压标量根和精确扇区关系
solutions={};tube=struct([]);
for n=[80 160]
    dx=1/n;x=((1:n)'-.5)*dx;initial=state([1 0 0 1;.125 0 0 .1]);
    Q=repmat(initial(1,:),n,1);Q(x>=.5,:)=repmat(initial(2,:),nnz(x>=.5),1);
    Q0=Q;t=0;net=zeros(1,4);count=0;minimumRho=Inf;minimumInternal=Inf;
    while t<.2
        [f,a]=axisymmetricHLLEFlux([initial(1,:);Q],[Q;initial(2,:)],repmat([1 0],n+1,1),gas,'einfeldt');
        dt=min(.4*dx/max(abs(a.waveSpeedsMPerS),[],'all'),.2-t);Q=Q-dt/dx*diff(f,1,1);
        internal=Q(:,4)-sum(Q(:,2:3).^2,2)./(2*Q(:,1));
        assert(all(isfinite(Q),'all')&&all(Q(:,1)>0&internal>0));
        minimumRho=min(minimumRho,min(Q(:,1)));minimumInternal=min(minimumInternal,min(internal));
        net=net+dt*(f(end,:)-f(1,:));t=t+dt;count=count+1;assert(count<1000);
    end
    [density,star]=sodDensity(x-.5,.2,gas.gamma);error=dx*sum(abs(Q(:,1)-density));
    ledger=dx*sum(Q-Q0,1)+net;assert(max(abs(ledger))<2e-13);
    row=struct('cells',n,'steps',count,'densityL1Error',error,'minimumDensity',minimumRho, ...
        'minimumInternalEnergy',minimumInternal,'conservationLedger',ledger,'starPressure',star(1),'starVelocity',star(2));
    if isempty(tube),tube=row;else,tube(end+1)=row;end %#ok<AGROW>
    solutions{end+1}=struct('Q',Q,'x',x,'exactDensity',density,'initialQ',Q0,'time',t); %#ok<AGROW>
end
order=log2(tube(1).densityL1Error/tube(2).densityL1Error);assert(order>.4);
assert(abs(tube(1).starPressure-.30313017805064685)<2e-12&&abs(tube(1).starVelocity-.9274526200489499)<2e-12);
report=struct('status','passed-explicit-HLLE-four-quantity-Riemann-and-Sod-numerical-contract', ...
    'matlabEngineProcessID',feature('getpid'),'caseLabels',{labels},'wavePolicies',records, ...
    'SodRuns',tube,'SodDensityL1Order',order,'refusals',refusals,'EulerBenchmarkOnly',true, ...
    'globalNSorSAorImplicitPositivityProved',false,'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'hlle_flux_evidence.mat'),'report','gas','QL','QR','normal','fixtures','solutions','-v7');
fid=fopen(fullfile(folder,'hlle_flux_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function Q=state(P)
        Q=axisymmetricFVConservative([P(:,1:3),P(:,4)./(P(:,1)*gas.R)],gas);
    end
    function reject(label,action)
        try,action();catch exception
            assert(startsWith(exception.identifier,'AxisFV:'),'HLLETest:Unexpected','意外错误%s。',exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('HLLETest:NoRefusal','非法输入%s未拒绝。',label);
    end
end

function [rho,star]=sodDensity(x,t,g)
% 仅Sod已知左稀疏/右激波支；不是隐藏的生产外流求解器。
al=sqrt(g);ar=sqrt(g*.1/.125);beta=(g-1)/(g+1);
fl=@(p)2*al/(g-1)*((p)^((g-1)/(2*g))-1);
fr=@(p)(p-.1)*sqrt(2/((g+1)*.125)/(p+beta*.1));
ps=fzero(@(p)fl(p)+fr(p),[.01 10]);us=(fr(ps)-fl(ps))/2;
rl=ps^(1/g);ratio=ps/.1;rr=.125*(ratio+beta)/(beta*ratio+1);
head=-al;tail=us-al*ps^((g-1)/(2*g));shock=ar*sqrt((g+1)/(2*g)*ratio+(g-1)/(2*g));
xi=x/t;rho=.125*ones(size(x));rho(xi<=head)=1;
fan=xi>head&xi<tail;rho(fan)=(2/(g+1)-(g-1)/(g+1)*xi(fan)/al).^(2/(g-1));
rho(xi>=tail&xi<=us)=rl;rho(xi>us&xi<=shock)=rr;star=[ps us];
end
