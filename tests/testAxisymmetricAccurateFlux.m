function report=testAxisymmetricAccurateFlux(folder)
%TESTAXISYMMETRICACCURATEFLUX 核对双double算术与同一HLLE模型，保存独立高精度参考所需原输入。
% 模块一：加法/乘法消去、除法和平方根余项，用明确解析恒等式核对。
% 模块二：十一类Riemann/近真空/单向波，两种波速、反法向与坐标旋转。
% 模块三：全部真实末期选中内部面逐一重算，原数据不改写；未认证完整RANS。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v75_validation','accurate_flux');end
assert(~exist(folder,'dir'),'AccurateFluxTest:Overwrite','不能覆盖测试证据。');
D=@(a)axisymmetricDoubleDouble('from',a);add=@(a,b)axisymmetricDoubleDouble('add',a,b);sub=@(a,b)axisymmetricDoubleDouble('sub',a,b);
mul=@(a,b)axisymmetricDoubleDouble('mul',a,b);div=@(a,b)axisymmetricDoubleDouble('div',a,b);
lost=sub(add(D(1e16),D(1)),D(1e16));assert(lost.hi==1&&lost.lo==0);
product=sub(mul(add(D(1),D(2^-52)),sub(D(1),D(2^-52))),D(1));assert(product.hi==-2^-104&&product.lo==0);
quotient=div(D(1),D(3));remainder=sub(mul(quotient,D(3)),D(1));assert(abs(remainder.hi+remainder.lo)<1e-31);
sqrt2=axisymmetricDoubleDouble('sqrt',D(2));remainder=sub(mul(sqrt2,sqrt2),D(2));assert(abs(remainder.hi+remainder.lo)<1e-30);
arithmetic=struct('lostUnit',lost,'cancelledProduct',product,'oneThird',quotient,'sqrt2',sqrt2);
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.01,'Pr',.72,'temperatureRange',[1e-14 1e12]));
left=[1 3 .2 2;1 0 0 1;.445 .698 0 3.528;1 0 0 1000;1 2 0 .4;1 -2 0 .4; ...
    1 -10 0 1e-8;1 0 0 1;1 0 2 1;1 40 0 300;1 -200 0 300];
right=[left(1,:);.125 0 0 .1;.5 0 0 .571;1 0 0 .01;1 -2 0 .4;1 2 0 .4; ...
    1 10 0 1e-8;2 0 0 1;1 0 -2 1;1 200 0 300;1 -40 0 300];
labels={'uniform','Sod','Lax','strong-shock','head-on','double-rarefaction','near-vacuum','stationary-contact','stationary-shear','all-right','all-left'};
QL=axisymmetricFVConservative([left(:,1:3),left(:,4)./(left(:,1)*gas.R)],gas);
QR=axisymmetricFVConservative([right(:,1:3),right(:,4)./(right(:,1)*gas.R)],gas);normal=repmat([1 0],size(QL,1),1);
policies=struct([]);benchmarks={};theta=.527;rotation=[cos(theta) -sin(theta);sin(theta) cos(theta)];
for policy={'davis','einfeldt'}
    name=policy{1};[flux,audit]=axisymmetricHLLEFluxAccurate(QL,QR,normal,gas,name);
    [reverse,~]=axisymmetricHLLEFluxAccurate(QR,QL,-normal,gas,name);
    error=max(abs(reverse+flux),[],'all')/max(1,max(abs(flux),[],'all'));assert(error<1e-13);
    ql=QL;qr=QR;ql(:,2:3)=QL(:,2:3)*rotation';qr(:,2:3)=QR(:,2:3)*rotation';
    rotated=axisymmetricHLLEFluxAccurate(ql,qr,normal*rotation',gas,name);expected=flux;expected(:,2:3)=flux(:,2:3)*rotation';
    rotationError=max(abs(rotated-expected),[],'all')/max(1,max(abs(expected),[],'all'));assert(rotationError<2e-12);
    assert(all(isnan(audit.intermediateConservativeState(~audit.middleMask,:)),'all')&& ...
        all(audit.intermediateConservativeState(audit.middleMask,1)>0)&&all(audit.intermediateInternalEnergyPerVolume(audit.middleMask)>0));
    row=struct('wavePolicy',name,'antisymmetryRelative',error,'rotationRelative',rotationError);
    if isempty(policies),policies=row;else,policies(end+1)=row;end %#ok<AGROW>
    benchmarks{end+1}=struct('wavePolicy',name,'QL',QL,'QR',QR,'normal',normal,'gas',gas,'flux',flux,'audit',audit); %#ok<AGROW>
end
path=fullfile(root,'outputs','v75_validation','late_euler','late_euler_witnesses.mat');d=load(path);real={};
for k=1:numel(d.witnesses)
    w=d.witnesses{k};f=w.faces;[found,order]=ismember(f.internalFaceID,f.faceID);assert(all(found));
    [flux,audit]=axisymmetricHLLEFluxAccurate(f.QL,f.QR,f.normal(order,:),d.c.gas,'einfeldt');
    real{end+1}=struct('alpha',w.case.alpha,'faceID',f.internalFaceID,'QL',f.QL,'QR',f.QR, ...
        'normal',f.normal(order,:),'gas',d.c.gas,'flux',flux,'audit',audit); %#ok<AGROW>
end
report=struct('status','passed-two-double-arithmetic-and-explicit-HLLE-basic-contract', ...
    'matlabEngineProcessID',feature('getpid'),'caseLabels',{labels},'policies',policies, ...
    'realSelectedCandidateCount',numel(real),'independentHighPrecisionComparisonPending',true, ...
    'fullRANSActuallyAdvanced',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'accurate_flux_evidence.mat'),'report','arithmetic','benchmarks','real','-v7');
fid=fopen(fullfile(folder,'accurate_flux_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
end
