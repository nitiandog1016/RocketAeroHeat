function report=testAxisymmetricSAModel(folder,evidenceFile)
%TESTAXISYMMETRICSAMODEL 独立代数公式、全部真实墙段及负输入拒绝校核。
% 模块一：逐标量直接幂公式复算固定SA候选，覆盖零场、两应变分支及r上限。
% 模块二：直墙解析距离和凸角端点距离，比较最近节点/同网格线的错误替代。
% 模块三：冻结真实TMR的1000体积质心与全部40壁段独立逐点逐段复算。
% 模块四：仅是局部函数与几何测试，没有推进SA或通过湍流/实验验收。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v65_validation','sa_model');end
if nargin<2
    evidenceFile=fullfile(root,'outputs','v64_validation','tmr_sutherland','tmr_sutherland_evidence.mat');
    if ~exist(evidenceFile,'file'),evidenceFile=fullfile(root,'validation','v64','tmr_sutherland','tmr_sutherland_evidence.mat');end
end
assert(exist(evidenceFile,'file')==2,'AxisFVTest:SAModelEvidence','缺少真实TMR冻结几何。');
if ~exist(folder,'dir'),mkdir(folder);end
gas=createTmrAxisymmetricAirGas;

%% 一、独立逐标量直接公式，不调用被测函数内部的稳定比值表达式
rho=1.2;T=311;mu=gas.transport.referenceViscosity*(T/gas.transport.referenceTemperature)^1.5* ...
    (gas.transport.referenceTemperature+gas.transport.sutherlandTemperature)/(T+gas.transport.sutherlandTemperature);
chi=[0 0 .5 5 5 20 1e3];Omega=[0 100 0 0 1e4 10 1e8];d=.003;
N=chi*(mu/rho);model=axisymmetricSAModel(rho,T,N,Omega,d,gas);
fields={'muT','production','destruction','nu','chi','fv1','fv2','sBar','sTilde','rSA','g','fw'};
reference=struct;for j=1:numel(fields),reference.(fields{j})=zeros(size(chi));end
maximumScaledDifference=0;
for i=1:numel(N)
    oracle=scalarReference(rho,mu,N(i),Omega(i),d);
    for j=1:numel(fields)
        key=fields{j};reference.(key)(i)=oracle.(key);difference=abs(model.(key)(i)-oracle.(key));
        denominator=max(abs([model.(key)(i) oracle.(key)]));
        if denominator>0,maximumScaledDifference=max(maximumScaledDifference,difference/denominator);
        else,assert(difference==0);end
    end
end
assert(maximumScaledDifference<2e-13,'AxisFVTest:SAModelFormula','直接公式与局部SA模型不一致。');
assert(all(model.muT(1:2)==0)&&all(model.production(1:2)==0)&&all(model.destruction(1:2)==0)&& ...
    model.rSA(1)==10&&model.rSA(2)==0&&model.sTilde(4)==0&&model.rSA(4)==10&& ...
    any(model.modifiedStrainBranch)&&any(~model.modifiedStrainBranch)&& ...
    all(isfinite(model.fw))&&~model.transportEquationSolved&&~model.physicalSolution&&~model.eligibleForSurfaceCoverage);
entry=axisymmetricSAModel(rho,T,5*mu/rho,0,d,gas);
assert(abs(entry.muT/mu-1.2942343413175514)<3e-14&& ...
    entry.constants.sigma==2/3&&entry.constants.cb2==.622&&entry.constants.turbulentPrandtl==.9&& ...
    entry.constants.cw1==.1355/.41^2+(1+.622)/(2/3)&& ...
    abs(entry.turbulentK/(gas.cp*entry.muT/.9)-1)<5e-15, ...
    'AxisFVTest:SAModelConstants','SA常数、nuTilde/nu=5身份或Prt不一致。');

% 以实际sBar建立AJS切换点，左右两侧均与直接公式相符且保持连续。
omega0=-entry.sBar/.7;nearOmega=omega0*[1-1e-9 1 1+1e-9];
near=axisymmetricSAModel(rho,T,5*mu/rho,nearOmega,d,gas);
for i=1:3
    oracle=scalarReference(rho,mu,5*mu/rho,nearOmega(i),d);
    assert(abs(near.sTilde(i)-oracle.sTilde)<1e-12*omega0);
end
assert(near.modifiedStrainBranch(1)&&~near.modifiedStrainBranch(3)&& ...
    max(near.sTilde)-min(near.sTilde)<1e-8*omega0);
large=axisymmetricSAModel(rho,T,1e100,1e100,1e100,gas);
assert(isfinite(large.fv1)&&large.fv1==1,'AxisFVTest:SAModelLargeChi','有限大chi不应仅因直接立方溢出而失效。');

%% 二、直墙解析式与凸肩端点，投影为几何操作而非裁物理场
boundary=plainBoundary;
[X,R]=ndgrid([0 .05 .1],[0 .002 .005]);R=.04+.3*X+R;
straightMesh=buildAxisymmetricFVGeometry(X,R);straight=axisymmetricSAWallDistance(straightMesh,boundary);
analytic=(straightMesh.volumeCentroid(:,2)-.3*straightMesh.volumeCentroid(:,1)-.04)/sqrt(1.09);
assert(max(abs(straight.distance-analytic))<1e-16&&all(straight.projectionFraction>0&straight.projectionFraction<1));
nodes=[X(:,1) R(:,1)];nearestNode=min(sqrt(sum((reshape(straightMesh.volumeCentroid,[],1,2)-reshape(nodes,1,[],2)).^2,3)),[],2);
assert(any(nearestNode>1.5*straight.distance),'AxisFVTest:SAWallNodeNegativeControl','最近节点负对照未区分真实壁段距离。');

[X,R]=ndgrid([-.03 0 .03],[0 .005 .02 .06]);R=.05-.7*abs(X)+R;
cornerMesh=buildAxisymmetricFVGeometry(X,R);corner=axisymmetricSAWallDistance(cornerMesh,boundary);
endpoint=corner.projectionFraction==0|corner.projectionFraction==1;
assert(any(endpoint)&&all(vecnorm(corner.point(endpoint,:)-[0 .05],2,2)<1e-15), ...
    'AxisFVTest:SAWallCorner','凸角端点没有参与最短壁距比较。');
cornerAudit=independentDistance(cornerMesh,corner);
multiBoundary=boundary;multiBoundary.iMin.type='wall';multi=axisymmetricSAWallDistance(cornerMesh,multiBoundary);
assert(numel(multi.wallFaceID)==numel(corner.wallFaceID)+size(X,2)-1&& ...
    all(multi.distance<=corner.distance+1e-16)&&any(multi.distance<corner.distance-1e-5));

% 周期平直环隙与无限直墙等价；起伏墙若不搜索镜像就可能漏掉跨缝最近点。
[Xp,Rp]=ndgrid([0 .02 .04],[.03 .04 .06]);periodicMesh=buildAxisymmetricFVGeometry(Xp,Rp,struct('periodicI',true));
periodicBoundary=boundary;periodicBoundary.iMin.type='periodic';periodicBoundary.iMax.type='periodic';periodicBoundary.jMax.type='wall';
periodic=axisymmetricSAWallDistance(periodicMesh,periodicBoundary);
assert(periodic.periodicGeometry&&max(abs(periodic.distance-min(periodicMesh.volumeCentroid(:,2)-.03, ...
    .06-periodicMesh.volumeCentroid(:,2))))<1e-16);

%% 三、真实TMR40段；另一实现从原始逻辑网格恢复，逐质心遍历全部段
beforeSHA=fileSHA(evidenceFile);frozen=load(evidenceFile,'mesh','boundary');mesh=frozen.mesh;
actual=axisymmetricSAWallDistance(mesh,frozen.boundary);
assert(mesh.cellCount==1000&&numel(actual.wallFaceID)==40&&all(actual.distance>0));
actualAudit=independentDistance(mesh,actual);
wallNodes=[mesh.X(:,1) mesh.R(:,1)];nearestNode=zeros(mesh.cellCount,1);
for i=1:mesh.cellCount,nearestNode(i)=min(vecnorm(wallNodes-mesh.volumeCentroid(i,:),2,2));end
assert(nnz(nearestNode>actual.distance+1e-10)>10&&strcmp(beforeSHA,fileSHA(evidenceFile)), ...
    'AxisFVTest:SAWallActualNegativeControl','真实TMR必须区别最近节点，冻结源不得改写。');

%% 四、非法局部状态和几何拒绝；不接受SA-neg、无墙或墙上零距离
refusals=0;
reject(@()axisymmetricSAModel(rho,T,-1,0,d,gas),'AxisFV:SAModelInput');
reject(@()axisymmetricSAModel(rho,T,1,-1,d,gas),'AxisFV:SAModelInput');
reject(@()axisymmetricSAModel(rho,T,1,0,0,gas),'AxisFV:SAModelInput');
reject(@()axisymmetricSAModel(0,T,1,0,d,gas),'AxisFV:SAModelInput');
reject(@()axisymmetricSAModel(rho,NaN,1,0,d,gas),'AxisFV:SAModelInput');
reject(@()axisymmetricSAModel(rho,T,[1 2],[1;2],d,gas),'AxisFV:SAModelInput');
reject(@()axisymmetricSAModel(rho,T,1,0,realmin,gas),'AxisFV:SAModelState');
reject(@()axisymmetricSAModel(rho,T,1e100,1e100,1e110,gas),'AxisFV:SAModelState');
reject(@()axisymmetricSAModel(rho,T,Inf,0,d,gas),'AxisFV:SAModelInput');
reject(@()axisymmetricSAModel(rho,10,1,0,d,gas),'AxisFV:TransportTemperature');
noWall=boundary;noWall.jMin.type='zeroGradient';
reject(@()axisymmetricSAWallDistance(straightMesh,noWall),'AxisFV:SAWallDistance');
bad=straightMesh;bad.volumeCentroid(1,:)=straight.wallStart(1,:);
reject(@()axisymmetricSAWallDistance(bad,boundary),'AxisFV:SAWallDistance');
bad=straightMesh;bad.faces.centroid(straight.wallFaceID(1),1)=bad.faces.centroid(straight.wallFaceID(1),1)+.001;
reject(@()axisymmetricSAWallDistance(bad,boundary),'AxisFV:SAWallDistance');
% 五节点周期起伏壁明确需要跨缝镜像搜索，不能套用有限段算法。
[Xw,Rw]=ndgrid([0 .01 .02 .03 .04],[.03 .04 .06]);Rw=Rw+.002*sin(2*pi*Xw/.04);
waveMesh=buildAxisymmetricFVGeometry(Xw,Rw,struct('periodicI',true));
reject(@()axisymmetricSAWallDistance(waveMesh,periodicBoundary),'AxisFV:SAWallDistance');
report=struct('status','passed-SA-local-formula-and-actual-shortest-wall-distance-only', ...
    'modelID',model.modelID,'maximumIndependentFormulaScaledDifference',maximumScaledDifference, ...
    'chi5MuTRatio',entry.muT/mu,'strainSwitchChecked',true,'largeFiniteChiChecked',true, ...
    'straightMaximumDistanceErrorM',max(abs(straight.distance-analytic)), ...
    'cornerEndpointCellCount',nnz(endpoint),'cornerAudit',cornerAudit,'actualAudit',actualAudit, ...
    'periodicStraightWallChecked',true,'periodicVaryingWallRejected',true, ...
    'actualCellCount',mesh.cellCount,'actualWallSegmentCount',numel(actual.wallFaceID), ...
    'actualNearestNodeWrongCellCount',nnz(nearestNode>actual.distance+1e-10), ...
    'actualDistanceRangeM',[min(actual.distance) max(actual.distance)],'refusals',refusals, ...
    'sourceEvidenceSHA256',beforeSHA,'sourceBindings',sourceBindings(root), ...
    'SAEquationSolved',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'sa_model_evidence.mat'),'report','model','reference','near','large', ...
    'straightMesh','straight','cornerMesh','corner','multi','periodicMesh','periodic','mesh','actual');
fid=fopen(fullfile(folder,'sa_model_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('SA local model/wall distance passed: independent %.6g; %d TMR cells/%d actual wall segments; %d refusals. No SA transport or production coverage claim.\n', ...
    maximumScaledDifference,mesh.cellCount,numel(actual.wallFaceID),refusals);
    function reject(action,identifier)
        try,action();catch exception
            assert(strcmp(exception.identifier,identifier),'AxisFVTest:SAModelRefusal', ...
                '拒绝身份不一致，实际%s，预期%s。',exception.identifier,identifier);
            refusals=refusals+1;return;
        end
        error('AxisFVTest:SAModelNoRefusal','非法输入没有被拒绝。');
    end
end

function r=scalarReference(rho,mu,N,omega,d)
% 独立直接幂公式逐标量复算，不调用被测的分段chi稳定表达式。
nu=mu/rho;chi=N/nu;fv1=chi^3/(chi^3+7.1^3);fv2=1-chi/(1+chi*fv1);
sBar=N*fv2/(.41^2*d^2);
if sBar>=-.7*omega,sTilde=omega+sBar;
else,sTilde=omega+omega*(.7^2*omega+.9*sBar)/((.9-2*.7)*omega-sBar);end
if sTilde==0,rSA=10;else,rSA=min(N/(sTilde*.41^2*d^2),10);end
g=rSA+.3*(rSA^6-rSA);fw=g*((1+2^6)/(g^6+2^6))^(1/6);
r=struct('muT',rho*N*fv1,'production',.1355*sTilde*N, ...
    'destruction',(.1355/.41^2+(1+.622)/(2/3))*fw*(N/d)^2, ...
    'nu',nu,'chi',chi,'fv1',fv1,'fv2',fv2,'sBar',sBar,'sTilde',sTilde,'rSA',rSA,'g',g,'fw',fw);
end

function boundary=plainBoundary
boundary=struct('iMin',struct('type','zeroGradient'),'iMax',struct('type','zeroGradient'), ...
    'jMin',struct('type','wall'),'jMax',struct('type','zeroGradient'));
end

function audit=independentDistance(mesh,result)
% 另一逐点算法：原始X/R的四条逻辑边建立段集，用夹角分类投影/两个端点。
segments=struct('a',{},'b',{},'faceID',{});names={'iMin','iMax','jMin','jMax'};
for side=1:4
    faceIDs=result.wallFaceID(mesh.faces.boundarySide(result.wallFaceID)==side);
    if isempty(faceIDs),continue;end
    if side==1,p=[mesh.X(1,:)' mesh.R(1,:)'];elseif side==2,p=[mesh.X(end,:)' mesh.R(end,:)'];
    elseif side==3,p=[mesh.X(:,1) mesh.R(:,1)];else,p=[mesh.X(:,end) mesh.R(:,end)];end
    for j=1:size(p,1)-1
        a=p(j,:);b=p(j+1,:);mid=(a+b)/2;radius=mid(2);
        if radius>0,centroid=((2*a(2)+b(2))*a+(a(2)+2*b(2))*b)/(6*radius);else,centroid=mid;end
        [errorM,index]=min(vecnorm(mesh.faces.centroid(faceIDs,:)-centroid,2,2));
        assert(errorM<1e-14&&strcmp(names{side},names{mesh.faces.boundarySide(faceIDs(index))}));
        segments(end+1)=struct('a',a,'b',b,'faceID',faceIDs(index)); %#ok<AGROW>
    end
end
maximumDistanceError=0;maximumPointError=0;nearestFaceTieDifferences=0;
for c=1:mesh.cellCount
    point=mesh.volumeCentroid(c,:);best=Inf;bestPoint=[NaN NaN];bestFace=Inf;
    for j=1:numel(segments)
        a=segments(j).a;b=segments(j).b;ab=b-a;
        if dot(point-a,ab)<=0,foot=a;
        elseif dot(point-b,ab)>=0,foot=b;
        else,foot=a+dot(point-a,ab)/dot(ab,ab)*ab;end
        candidate=norm(point-foot);
        if candidate<best||(candidate==best&&segments(j).faceID<bestFace)
            best=candidate;bestPoint=foot;bestFace=segments(j).faceID;
        end
    end
    maximumDistanceError=max(maximumDistanceError,abs(best-result.distance(c)));
    maximumPointError=max(maximumPointError,norm(bestPoint-result.point(c,:)));
    % 角点可以由两个真实段共同取得同一最短点；另按返回face重算以拒绝伪来源。
    source=find([segments.faceID]==result.faceID(c));assert(isscalar(source));
    a=segments(source).a;b=segments(source).b;ab=b-a;
    if dot(point-a,ab)<=0,sourceFoot=a;
    elseif dot(point-b,ab)>=0,sourceFoot=b;
    else,sourceFoot=a+dot(point-a,ab)/dot(ab,ab)*ab;end
    assert(norm(sourceFoot-result.point(c,:))<5e-15&&abs(norm(point-sourceFoot)-best)<5e-15, ...
        'AxisFVTest:SAWallFace','返回的源face不对应真实最短线段投影。');
    nearestFaceTieDifferences=nearestFaceTieDifferences+(result.faceID(c)~=bestFace);
end
assert(maximumDistanceError<5e-15&&maximumPointError<5e-15,'AxisFVTest:SAWallIndependent','独立最短段距/坐标不一致。');
audit=struct('cellCount',mesh.cellCount,'allLogicalWallSegmentCount',numel(segments), ...
    'maximumDistanceErrorM',maximumDistanceError,'maximumNearestPointErrorM',maximumPointError, ...
    'sourceFaceProjectionIndependentlyChecked',true,'nearestFaceTieDifferences',nearestFaceTieDifferences, ...
    'nearestNodeOrSameGridLineSubstitution',false,'allEndpointsChecked',true);
end
function bindings=sourceBindings(root)
paths={'regionalFlow/axisymmetric/axisymmetricSAModel.m','regionalFlow/axisymmetric/axisymmetricSAWallDistance.m', ...
    'tests/testAxisymmetricSAModel.m'};bindings=struct([]);
for j=1:numel(paths)
    item=struct('path',paths{j},'sha256',fileSHA(fullfile(root,paths{j})));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
end
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
