function report=testAxisymmetricTemperatureTrace(folder,evidenceFile)
%TESTAXISYMMETRICTEMPERATURETRACE 本构温度局部限制器的线性、斜网格与真实TMR校核。
% 模块一：严格位于邻域界内的线性标量保持alpha=1，不裁单元状态或原梯度。
% 模块二：强斜率及周期自配对同时检查左右shift面，不能被平均温度掩盖越界迹。
% 模块三：实际v63终态在新输运下空间计算，6个域外raw面在真实局部界内恢复。
% 模块四：限制只供mu/k取样；墙面输入Tw保持精确，不声明完整NS正性或新物理解。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v64_validation','temperature_trace');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v63_validation','tmr_hlle_continue1','tmr_implicit_budget.mat');
    if ~exist(evidenceFile,'file'),evidenceFile=fullfile(root,'validation','v63','tmr_hlle_continue1','tmr_implicit_budget.mat');end
end
assert(exist(evidenceFile,'file')==2,'AxisFVTest:TemperatureTraceSource','缺少实际TMR冻结终态。');
if ~exist(folder,'dir'),mkdir(folder);end

%% 一、T=r的正线性标量，计算差与局部界差严格同源，避免舍入把线性迹误称过冲
[X,R]=ndgrid([0 .02 .05],[.03 .04 .06 .09]);linearMesh=buildAxisymmetricFVGeometry(X,R);
T=linearMesh.volumeCentroid(:,2);gx=zeros(size(T));gr=ones(size(T));
boundaryT=linearMesh.faces.centroid(linearMesh.boundaryFaces,2);
linear=axisymmetricFVTemperatureTrace(linearMesh,T,gx,gr,boundaryT);
assert(all(linear.cellLimiter==1)&&linear.limitedCellCount==0&& ...
    max(abs(linear.faceTemperature-linearMesh.faces.centroid(:,2)))<1e-15&& ...
    max(abs(linear.sourceTemperature-linearMesh.planarCentroid(:,2)))<1e-15, ...
    'AxisFVTest:TemperatureTraceLinear','处于局部界内的线性标量被限幅或不再线性精确。');
linearAudit=independentAudit(linearMesh,T,gx,gr,boundaryT,linear);

%% 二、强斜率斜网格；独立逐cell收集全部目标点，核对统一alpha与状态保持
[Xs,Rs]=ndgrid([0 .02 .05],[.03 .04 .06]);strongMesh=buildAxisymmetricFVGeometry(Xs+.7*Rs,Rs);
strongT=[300;360;420;330];strongGX=[1e5;-8e4;3e4;-2e4];strongGR=[5e4;5e4;-6e4;7e4];
strongBoundary=320+1000*strongMesh.faces.centroid(strongMesh.boundaryFaces,2);
copies={strongT,strongGX,strongGR,strongBoundary};
strong=axisymmetricFVTemperatureTrace(strongMesh,strongT,strongGX,strongGR,strongBoundary);
assert(isequal(copies,{strongT,strongGX,strongGR,strongBoundary})&&any(strong.cellLimiter<1)&& ...
    all(strong.cellLimiter>=0&strong.cellLimiter<=1)&&~strong.cellStateClipped&&~strong.originalViscousThermalGradientLimited);
strongAudit=independentAudit(strongMesh,strongT,strongGX,strongGR,strongBoundary,strong);

% 单轴向单元周期自配对：左增量-500给alpha<=.4，右shift增量+500给更严alpha<=.02。
[Xp,Rp]=ndgrid([0 .1],[.03 .06]);periodicMesh=buildAxisymmetricFVGeometry(Xp,Rp,struct('periodicI',true));
assert(periodicMesh.cellCount==1);pT=300;pGX=10000;pGR=0;
pBoundary=300*ones(numel(periodicMesh.boundaryFaces),1);side=periodicMesh.faces.boundarySide(periodicMesh.boundaryFaces);
pBoundary(side==3)=100;pBoundary(side==4)=310;
periodic=axisymmetricFVTemperatureTrace(periodicMesh,pT,pGX,pGR,pBoundary);
self=find(periodicMesh.faces.neighbour>0);assert(isscalar(self)&&periodicMesh.faces.owner(self)==periodicMesh.faces.neighbour(self));
expectedAlpha=.02*(1-16*eps);
assert(abs(periodic.cellLimiter-expectedAlpha)<1e-15&&periodic.cellLimiter<.1&& ...
    abs(periodic.rawLeftTemperature(self)+200)<1e-10&&abs(periodic.rawRightTemperature(self)-800)<1e-10&& ...
    periodic.leftTemperature(self)>=100&&periodic.rightTemperature(self)<=310&& ...
    abs(periodic.faceTemperature(self)-300)<1e-12, ...
    'AxisFVTest:TemperatureTracePeriodic','周期自配对的右shift目标没有参与统一alpha；平均迹不能掩盖两侧过冲。');
periodicAudit=independentAudit(periodicMesh,pT,pGX,pGR,pBoundary,periodic);
refusals=0;reject(@()axisymmetricFVTemperatureTrace(periodicMesh,-1,pGX,pGR,pBoundary));
reject(@()axisymmetricFVTemperatureTrace(periodicMesh,pT,NaN,pGR,pBoundary));
reject(@()axisymmetricFVTemperatureTrace(periodicMesh,pT,pGX,pGR,zeros(size(pBoundary))));
reject(@()axisymmetricFVTemperatureTrace(periodicMesh,[pT;pT],pGX,pGR,pBoundary));

%% 三、实际TMR同一Q只换显式输运；完整保存raw/limited/alpha并独立复算局部界
sourceSHA=fileSHA(evidenceFile);frozen=load(evidenceFile,'implicit','mesh','gas','boundary');
mesh=frozen.mesh;Q=frozen.implicit.Q;originalQ=Q;boundary=frozen.boundary;gas=frozen.gas;
for key={'cp','cv','k'},if isfield(gas,key{1}),gas=rmfield(gas,key{1});end,end
gas.mu=1.716e-5;gas.transport=struct('model','sutherland','referenceViscosity',1.716e-5, ...
    'referenceTemperature',273.15,'sutherlandTemperature',110.4);gas=axisymmetricFVGas(gas);
[rate,info]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);reconstruction=info.transport.reconstruction;
assert(isequal(Q,originalQ)&&isequal(info.primitive,axisymmetricFVState(originalQ,gas))&&all(isfinite(rate),'all')&& ...
    ~reconstruction.cellStateClipped&&~reconstruction.originalViscousThermalGradientLimited&& ...
    reconstruction.limitedCellCount>0&&all(reconstruction.cellLimiter>=0&reconstruction.cellLimiter<=1));
rawBad=find(reconstruction.rawFaceTemperature<gas.temperatureRange(1)|reconstruction.rawFaceTemperature>gas.temperatureRange(2));
assert(isequal(rawBad,[58;60;62;64;66;80]),'AxisFVTest:TemperatureTraceActual','实际旧Q的6个无限制域外面身份改变，需要重新审查。');
assert(all(reconstruction.faceTemperature>=gas.temperatureRange(1)&reconstruction.faceTemperature<=gas.temperatureRange(2))&& ...
    all(reconstruction.sourceTemperature>=gas.temperatureRange(1)&reconstruction.sourceTemperature<=gas.temperatureRange(2)));
wall=info.wall;assert(isequal(reconstruction.faceTemperature(wall.faceID),wall.temperatureInput)&& ...
    isequal(info.transport.faceTemperature,reconstruction.faceTemperature)&& ...
    isequal(info.transport.sourceTemperature,reconstruction.sourceTemperature));
sampleBoundary=info.faceViscousGradientTrace(mesh.boundaryFaces,4);sample=[info.primitive(:,4);sampleBoundary];
actualGX=mesh.gradientX*sample;actualGR=mesh.gradientR*sample;transportBoundary=sampleBoundary;
boundaryMap=zeros(numel(mesh.faces.owner),1);boundaryMap(mesh.boundaryFaces)=1:numel(mesh.boundaryFaces);
transportBoundary(boundaryMap(wall.faceID))=wall.temperatureInput;
actualAudit=independentAudit(mesh,info.primitive(:,4),actualGX,actualGR,transportBoundary,reconstruction);
assert(~info.physicalSolution&&~info.eligibleForSurfaceCoverage&&strcmp(sourceSHA,fileSHA(evidenceFile)));
report=struct('status','passed-local-constitutive-temperature-traces-and-actual-TMR-state', ...
    'linearAudit',linearAudit,'strongAudit',strongAudit,'periodicAudit',periodicAudit,'refusals',refusals, ...
    'sourceEvidenceSHA256',sourceSHA,'actualCellTemperatureRange',[min(info.primitive(:,4)) max(info.primitive(:,4))], ...
    'actualRawFaceRange',[min(reconstruction.rawFaceTemperature) max(reconstruction.rawFaceTemperature)], ...
    'actualLimitedFaceRange',[min(reconstruction.faceTemperature) max(reconstruction.faceTemperature)], ...
    'actualRawInvalidFaces',rawBad,'actualLimitedCellCount',reconstruction.limitedCellCount,'actualAudit',actualAudit, ...
    'cellQUnchanged',true,'wallInputTemperatureExact',true,'physicalSolution',false, ...
    'completeNSPositivityValidated',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0, ...
    'sourceBindings',sourceBindings(root));
save(fullfile(folder,'temperature_trace_evidence.mat'),'report','linear','linearMesh','strong','strongMesh', ...
    'strongT','strongGX','strongGR','strongBoundary','periodic','periodicMesh','mesh','Q','gas','boundary','rate','info');
fid=fopen(fullfile(folder,'temperature_trace_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Temperature trace passed: unchanged linear, strong local bounds, both periodic-self shifted sides, %d refusals; actual TMR 6 raw bad faces/ %d limited cells, Q unchanged; no production claim.\n',refusals,reconstruction.limitedCellCount);
    function reject(action)
        try,action();catch exception
            assert(strcmp(exception.identifier,'AxisFV:TemperatureTrace'),'AxisFVTest:TemperatureTraceRefusal','非法输入出现非预期拒绝。');
            refusals=refusals+1;return;
        end
        error('AxisFVTest:TemperatureTraceNoRefusal','非法重构输入没有被拒绝。');
    end
end

function audit=independentAudit(mesh,T,gx,gr,boundaryT,trace)
% 独立按cell遍历目标点及邻接边界；不用被测模块的向量targetCell/accumarray路径。
n=mesh.cellCount;f=mesh.faces;lower=T;upper=T;alpha=ones(n,1);bfIndex=zeros(numel(f.owner),1);bfIndex(mesh.boundaryFaces)=1:numel(mesh.boundaryFaces);
for c=1:n
    related=find(f.owner==c|f.neighbour==c);values=T(c);points=mesh.planarCentroid(c,:);
    for face=related'
        if f.owner(face)==c
            points(end+1,:)=f.centroid(face,:); %#ok<AGROW>
            if f.neighbour(face)>0,values(end+1)=T(f.neighbour(face));else,values(end+1)=boundaryT(bfIndex(face));end %#ok<AGROW>
        end
        if f.neighbour(face)==c
            points(end+1,:)=f.centroid(face,:)-f.neighbourShift(face,:);values(end+1)=T(f.owner(face)); %#ok<AGROW>
        end
    end
    lower(c)=min(values);upper(c)=max(values);d=points-mesh.volumeCentroid(c,:);increment=d*[gx(c);gr(c)];
    for q=1:numel(increment)
        if increment(q)>0,alpha(c)=min(alpha(c),(upper(c)-T(c))/increment(q));
        elseif increment(q)<0,alpha(c)=min(alpha(c),(lower(c)-T(c))/increment(q));end
    end
    if alpha(c)<1,alpha(c)=alpha(c)*(1-16*eps);end
end
assert(isequal(lower,trace.cellLower)&&isequal(upper,trace.cellUpper)&&max(abs(alpha-trace.cellLimiter))<3e-14);
owner=f.owner;inside=f.neighbour>0;other=f.neighbour(inside);
margin=256*eps(max([T;boundaryT]));
assert(all(trace.leftTemperature>=lower(owner)-margin&trace.leftTemperature<=upper(owner)+margin)&& ...
    all(trace.rightTemperature(inside)>=lower(other)-margin&trace.rightTemperature(inside)<=upper(other)+margin)&& ...
    all(trace.sourceTemperature>=lower-margin&trace.sourceTemperature<=upper+margin));
dL=f.centroid-mesh.volumeCentroid(owner,:);left=T(owner)+alpha(owner).*(gx(owner).*dL(:,1)+gr(owner).*dL(:,2));
right=left;dR=f.centroid(inside,:)-f.neighbourShift(inside,:)-mesh.volumeCentroid(other,:);
right(inside)=T(other)+alpha(other).*(gx(other).*dR(:,1)+gr(other).*dR(:,2));right(mesh.boundaryFaces)=boundaryT;
expectedFace=.5*(left+right);expectedFace(mesh.boundaryFaces)=boundaryT;
delta=mesh.planarCentroid-mesh.volumeCentroid;expectedSource=T+alpha.*(gx.*delta(:,1)+gr.*delta(:,2));
assert(max(abs(expectedFace-trace.faceTemperature))<=margin&&max(abs(expectedSource-trace.sourceTemperature))<=margin);
audit=struct('localBoundsRecomputed',true,'allOwnerNeighbourShiftAndSourceTargetsChecked',true, ...
    'maximumAlphaDifference',max(abs(alpha-trace.cellLimiter)),'limitedCellCount',nnz(alpha<1), ...
    'minimumAlpha',min(alpha),'maximumAlpha',max(alpha),'originalCellTemperaturesUnchanged',true);
end
function bindings=sourceBindings(root)
paths={'regionalFlow/axisymmetric/axisymmetricFVTemperatureTrace.m','regionalFlow/axisymmetric/axisymmetricNavierStokesResidual.m', ...
    'regionalFlow/axisymmetric/axisymmetricFVTransport.m','tests/testAxisymmetricTemperatureTrace.m'};bindings=struct([]);
for j=1:numel(paths)
    item=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
