function report=testRocketAxisymmetricWallMapping(folder)
%TESTROCKETAXISYMMETRICWALLMAPPING 用独立一维积分检验原始三角面保守映射。
% 模块一：完整默认箭体及改变周向/底面/墙角细分的第二几何均保留原F/V。
% 模块二：独立由端环节点计数推导三角形周向宽度，再积分二次多项式。
% 模块三：正负制造热流的曲面总功率、常值保持及真实底面覆盖分别验收。
% 模块四：改动度量/源面绑定须拒绝；制造数据仅校核映射，不是流场/热流验证。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v68_validation','wall_mapping');end
options={struct,struct('geometry',struct('resolution',16,'protrusionEnabled',false), ...
    'domain',struct('radialLayers',4,'maxWallAngleStepDeg',3,'baseRadialSegments',3))};
cases=struct([]);
for c=1:numel(options)
    item=createRocketAxisymmetricCase(options{c});map=mapRocketAxisymmetricWall(item);
    d=item.domain;surface=item.surfaceMesh;F=surface.F;V=surface.V;
    nf=size(F,1);nw=numel(d.wallBoundaryFaceID);reference=zeros(nf,nw);closedRowArea=zeros(nf,1);
    for f=1:nf
        p=find(d.originalWallFaceBandID==surface.param.faceBandID(f));
        A=d.originalWallPolylineXR(p,:);B=d.originalWallPolylineXR(p+1,:);E=B-A;L=norm(E);
        ids=F(f,:);xr=[V(ids,1),item.geometryTruth.nodeBaselineRadius(ids)];
        atA=vecnorm(xr-A,2,2)<1e-12;atB=vecnorm(xr-B,2,2)<1e-12;
        assert(all(xor(atA,atB)));
        nonpole=item.geometryTruth.nodeBaselineRadius(ids)>0;
        angles=sort(mod(atan2(V(ids(nonpole),3),V(ids(nonpole),2)),2*pi));
        gaps=diff([angles;angles(1)+2*pi]);deltaTheta=2*pi-max(gaps);
        if nnz(nonpole)==2,c0=1;c1=0;
        elseif nnz(atA)==2,c0=1;c1=-1;
        else,assert(nnz(atB)==2);c0=0;c1=1;end
        r0=A(2);dr=B(2)-r0;
        closedRowArea(f)=L*deltaTheta*(r0*(c0+c1/2)+dr*(c0/2+c1/3));
        for j=find(d.wallSourceOriginalSegment==p)'
            lo=dot(d.wallPolylineXR(j,:)-A,E)/(L*L);
            hi=dot(d.wallPolylineXR(j+1,:)-A,E)/(L*L);
            reference(f,j)=L*deltaTheta*(r0*c0*(hi-lo)+ ...
                (r0*c1+dr*c0)*(hi^2-lo^2)/2+dr*c1*(hi^3-lo^3)/3);
        end
    end
    relative=max(abs(full(map.overlapAreaM2)-reference)./max(map.sourceCurvedAreaM2,realmin),[],'all');
    assert(relative<2e-11,'RocketMapTest:Independent','独立多项式积分与裁切面积矩阵不符。');
    assert(max(abs(map.sourceCurvedAreaM2-closedRowArea)./closedRowArea)<2e-11, ...
        'RocketMapTest:ClosedRow','整面闭式曲面积与墙段裁切分解不符。');
    assert(isequal(map.sourceVertices,V)&&isequal(map.sourceFaces,F)&& ...
        all(abs(full(sum(map.weights,2))-1)<2e-13)&& ...
        abs(sum(map.sourceCurvedAreaM2)-d.fullRevolvedWallAreaM2)<2e-12*d.fullRevolvedWallAreaM2);
    q=2e4*sin((1:nw)'*.7)-1e3;p=4e3+(1:nw)';
    mappedQ=map.weights*q;mappedP=map.weights*p;
    sourceIntegral=sum(map.sourceCurvedAreaM2.*mappedQ);wallIntegral=sum(map.wallAreaM2.*q);
    budget=2e-11*sum(map.wallAreaM2.*abs(q));
    assert(any(q<0)&&any(q>0)&&abs(sourceIntegral-wallIntegral)<=budget, ...
        'RocketMapTest:Conservation','正负制造场曲面功率不守恒：误差%.9g，预算%.9g。',sourceIntegral-wallIntegral,budget);
    % sparse内积和除以逐面曲面积各自有舍入；先传播真实权重和行和误差，不能裁映射值。
    operationCount=2*full(sum(spones(map.weights),2))+4;
    gamma=operationCount*eps./(1-operationCount*eps);
    pBudget=abs(full(sum(map.weights,2))-1)*max(abs(p))+gamma.*(abs(map.weights)*abs(p));
    assert(all(mappedP>=min(p)-pBudget&mappedP<=max(p)+pBudget), ...
        'RocketMapTest:Convexity','平均值超出独立逐行点积形成误差界。');
    base=surface.param.faceBandID==-1;
    assert(abs(sum(map.sourceCurvedAreaM2(base))-pi*item.geometryTruth.baseRadius^2)<2e-12&& ...
        all(map.sourceFaceIsPoleSector(base))&&all(map.sourceFaceIsPoleSector(surface.param.faceBandID==max(surface.param.faceBandID))));
    summary=struct('sourceFaces',nf,'wallSegments',nw,'overlapEntries',nnz(map.overlapAreaM2), ...
        'independentIntegralMaxRelativeError',relative,'signedManufacturedIntegralErrorW',sourceIntegral-wallIntegral, ...
        'curvedAreaM2',sum(map.sourceCurvedAreaM2),'facetAreaM2',sum(map.sourceDisplayFacetAreaM2), ...
        'originalFacesPreserved',true,'realPoleSectorsTested',true);
    if isempty(cases),cases=summary;else,cases(c)=summary;end %#ok<AGROW>
end
bad=item;bad.mesh.R(2,2)=bad.mesh.R(2,2)+1e-6;reject(bad,'RocketMap:Binding');
bad=item;bad.domain.wallSourceOriginalSegment(1)=2;reject(bad,'RocketMap:Binding');
bad=item;bad.surfaceMesh.F(1,:)=bad.surfaceMesh.F(1,[2 3 1]);reject(bad,'RocketMap:NativeSource');
report=struct('status','passed-native-source-face-curved-area-conservative-projection', ...
    'cases',cases,'invalidBindingRefusals',3,'manufacturedValuesArePhysicalData',false, ...
    'flowStateAdvanced',false,'physicalSolution',false,'experimentalValidation',false, ...
    'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'rocket_wall_mapping_evidence.mat'),'report','map');
fid=fopen(fullfile(folder,'rocket_wall_mapping_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Native source face projection: %d cases, independent parameter integrals and signed conservation passed.\n',numel(cases));
end

function reject(item,identifier)
% 模块用途：旧度量/源面绑定不能沿用到新几何。
try,mapRocketAxisymmetricWall(item);catch exception
    assert(strcmp(exception.identifier,identifier),'RocketMapTest:Refusal','非预期拒绝%s。',exception.identifier);return;
end
error('RocketMapTest:Refusal','错误几何绑定未拒绝。');
end
