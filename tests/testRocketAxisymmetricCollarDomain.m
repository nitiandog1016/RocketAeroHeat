function report=testRocketAxisymmetricCollarDomain(folder)
%TESTROCKETAXISYMMETRICCOLLARDOMAIN 完整箭体法向近壁候选的独立几何校核。
% 模块一：两套真实箭体核对完整墙/远场/轴与原面保留，不以几何通过冒称流场完成。
% 模块二：在真实三角分片上用Duffy/Gauss直接积分A、r、xr、r²，核对FV体积/质心。
% 模块三：核对实际法向偏移及墙距/切向偏移，比较v69热格对应墙的偏斜程度。
% 模块四：输入超域必须拒绝，不缩层、不移动给定墙、不裁原始控制体。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v70_validation','collar_domain');end
before=bindings(root);cases=cell(2,1);records=cell(2,1);
configs={struct('protrusionEnabled',false),struct('protrusionEnabled',false,'resolution',32, ...
    'bodyRadius',.1,'baseRadius',.06,'afterbodyLength',.18,'shoulderLength',.1,'noseRadius',.02,'coneAngleDeg',14)};
for k=1:numel(configs)
    [surface,truth]=buildRocketTestBody(configs{k});
    [mesh,domain]=buildRocketAxisymmetricCollarDomain(surface,truth);
    [legacy,original]=buildRocketAxisymmetricDomain(surface,truth,domain.options);
    assert(mesh.shape(2)==24&&isequal(mesh.X(:,1),domain.wallPolylineXR(:,1))&& ...
        isequal(mesh.R(:,1),domain.wallPolylineXR(:,2))&& ...
        isequal([mesh.X(:,end),mesh.R(:,end)],original.outerBoundaryPolylineXR)&& ...
        isequal(domain.originalWallPolylineXR,original.originalWallPolylineXR)&& ...
        isequal(domain.sourceVertices,surface.V)&&isequal(domain.sourceFaces,surface.F)&& ...
        isequal(domain.wallSourceTriangleIDs,original.wallSourceTriangleIDs)&& ...
        isequal(domain.wallSourceFaceBandID,original.wallSourceFaceBandID)&& ...
        isequal(domain.wallAreaPerRadianM2,original.wallAreaPerRadianM2));
    assert(~domain.physicalSolution&&~domain.eligibleForSurfaceCoverage&&~domain.nearWallResolutionValidated&& ...
        ~domain.farfieldAdequacyValidated&&~domain.collar.sourceFieldProjectionCompleted&& ...
        ~domain.collar.flowAdvancePerformed&&domain.newProductionSolvedAreaM2==0);
    axis=[domain.axisBoundaryFaceID{1};domain.axisBoundaryFaceID{2}];
    assert(all(mesh.faces.centroid(axis,2)==0)&&all(mesh.faces.areaPerRadian(axis)==0)&& ...
        all(mesh.volumePerRadian>0)&&all(mesh.faces.normalDistance>0)&& ...
        size(mesh.vertices,1)-numel(mesh.faces.owner)+mesh.cellCount==1&& ...
        numel(mesh.boundaryFaces)==2*sum(mesh.shape));
    m=moments(mesh.vertices,mesh.cells);center=m(:,3:4)./m(:,2);
    errA=max(abs(m(:,1)-mesh.planarArea)./m(:,1));errV=max(abs(m(:,2)-mesh.volumePerRadian)./m(:,2));
    errC=max(abs(center-mesh.volumeCentroid),[],'all');
    assert(errA<3e-11&&errV<3e-11&&errC<3e-11,'RocketCollarTest:Moments','独立r权积分不一致。');
    assert(abs(sum(m(:,2))-domain.fluidVolumePerRadianM3)<=3e-12*domain.fluidVolumePerRadianM3&& ...
        abs(domain.fluidVolumePerRadianM3-original.fluidVolumePerRadianM3)<=3e-12*original.fluidVolumePerRadianM3);
    wall=domain.wallPolylineXR;edge=diff(wall);out=[edge(:,2),-edge(:,1)]./vecnorm(edge,2,2);
    offsetError=0;scale=max(1,max(abs(mesh.vertices),[],'all'));
    for j=1:domain.collar.options.normalLayers
        points=[mesh.X(:,j+1),mesh.R(:,j+1)];delta=points-wall;
        offset=domain.collar.normalOffsetM(j+1);
        error=[sum(delta(1:end-1,:).*out,2)-offset;sum(delta(2:end,:).*out,2)-offset];
        offsetError=max(offsetError,max(abs(error)));
        assert(all(abs(error)<=1024*eps(scale)),'RocketCollarTest:Normal','实际节点未满足真实墙线的法向偏移。');
    end
    boundary=mesh.boundaryFaces;sample=[mesh.volumeCentroid;mesh.faces.centroid(boundary,:)];
    value=1+.2*sample(:,1)+.3*sample(:,2);
    grad=max([abs(mesh.gradientX*value-.2);abs(mesh.gradientR*value-.3)]);
    assert(grad<3e-9,'RocketCollarTest:WLS','真实几何WLS未再现仿射场。');
    oldSkew=domain.collar.legacyFirstCellSkewRatio;newSkew=domain.collar.firstCellSkewRatio;
    contraction=strcmp(domain.wallRegion,'contraction');afterbody=strcmp(domain.wallRegion,'afterbody');
    [~,tail]=min(domain.wallPolylineXR(find(afterbody),1));aid=find(afterbody);tail=aid(tail);
    assert(newSkew(tail)<oldSkew(tail),'RocketCollarTest:Tail','后体近底面墙的偏斜未得到实际改善。');
    records{k}=struct('cellCount',mesh.cellCount,'wallSegments',numel(domain.wallBoundaryFaceID), ...
        'relativeAreaError',errA,'relativeVolumeError',errV,'centroidErrorM',errC, ...
        'normalOffsetErrorM',offsetError,'affineGradientError',grad, ...
        'minimumWallDistanceM',min(domain.wallNormalDistanceM),'maximumWallDistanceM',max(domain.wallNormalDistanceM), ...
        'contractionMedianSkewBeforeAfter',[median(oldSkew(contraction)),median(newSkew(contraction))], ...
        'afterbodyMedianSkewBeforeAfter',[median(oldSkew(afterbody)),median(newSkew(afterbody))], ...
        'tailWallID',domain.wallBoundaryFaceID(tail),'tailSkewBeforeAfter',[oldSkew(tail),newSkew(tail)]);
    cases{k}=struct('surface',surface,'truth',truth,'mesh',mesh,'domain',domain,'legacy',legacy);
end

%% 四、非法近壁层/几何不以静默回退伪装成功
surface=cases{1}.surface;truth=cases{1}.truth;refusals=struct('label',{},'identifier',{});
reject('unknown-collar-option',@()buildRocketAxisymmetricCollarDomain(surface,truth,struct,struct('smoothAwayCorners',true)));
reject('no-normal-layer',@()buildRocketAxisymmetricCollarDomain(surface,truth,struct,struct('normalLayers',0)));
reject('no-outer-connection',@()buildRocketAxisymmetricCollarDomain(surface,truth,struct('radialLayers',4),struct('normalLayers',3)));
reject('negative-offset',@()buildRocketAxisymmetricCollarDomain(surface,truth,struct,struct('firstNormalSpacing',-1)));
reject('shrinking-offset',@()buildRocketAxisymmetricCollarDomain(surface,truth,struct,struct('normalGrowth',.5)));
reject('collar-outside-farfield',@()buildRocketAxisymmetricCollarDomain(surface,truth,struct,struct('firstNormalSpacing',100)));
[featured,featuredTruth]=buildRocketTestBody(struct('resolution',24));
reject('nonaxisymmetric-protrusion',@()buildRocketAxisymmetricCollarDomain(featured,featuredTruth));
assert(isequaln(before,bindings(root)),'RocketCollarTest:Source','测试期间源码改变。');
report=struct('status','passed-complete-rocket-normal-collar-geometry-only-unreleased', ...
    'matlabEngineProcessID',feature('getpid'),'cases',{records},'refusals',refusals,'sourceBindings',before, ...
    'sourceFieldProjectionPerformed',false,'flowResidualEvaluated',false,'flowAdvancePerformed',false, ...
    'wallNormalResolutionValidated',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'collar_domain_evidence.mat'),'report','cases');
fid=fopen(fullfile(folder,'collar_domain_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Normal collar geometry only: two complete rockets, positive convex FV/r-weighted moments/offsets, %d refusals; no flow claim.\n',numel(refusals));
    function reject(label,action)
        try,action();catch exception
            assert(~isempty(exception.identifier));refusals(end+1)=struct('label',label,'identifier',exception.identifier);return;
        end
        error('RocketCollarTest:Refusal','非法合同%s没有拒绝。',label);
    end
end

function m=moments(vertices,cells)
% 模块用途：两三角形上的直接Gauss积分，与构造器边界矩公式独立。
g=[(1-sqrt(3/5))/2,.5,(1+sqrt(3/5))/2];w=[5/18,4/9,5/18];m=zeros(size(cells,1),4);
for triplet={[1 2 3],[1 3 4]}
    ids=triplet{1};A=vertices(cells(:,ids(1)),:);E=vertices(cells(:,ids(2)),:)-A;D=vertices(cells(:,ids(3)),:)-A;
    jac=abs(E(:,1).*D(:,2)-E(:,2).*D(:,1));
    for i=1:3
        for j=1:3
            p=A+g(i)*E+(1-g(i))*g(j)*D;v=jac*(1-g(i))*w(i)*w(j);
            m=m+v.*[ones(size(v)),p(:,2),p(:,1).*p(:,2),p(:,2).^2];
        end
    end
end
end
function out=bindings(root)
paths={'buildRocketTestBody.m','regionalFlow/axisymmetric/buildAxisymmetricFVGeometry.m', ...
    'regionalFlow/axisymmetric/buildRocketAxisymmetricDomain.m','regionalFlow/axisymmetric/buildRocketAxisymmetricCollarDomain.m', ...
    'tests/testRocketAxisymmetricCollarDomain.m'};
out=struct('path',{},'sha256',{});
for k=1:numel(paths)
    fid=fopen(fullfile(root,paths{k}),'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));v=fread(fid,Inf,'*uint8');clear guard;
    h=java.security.MessageDigest.getInstance('SHA-256');h.update(v);
    out(k)=struct('path',paths{k},'sha256',lower(reshape(dec2hex(typecast(h.digest(),'uint8'),2)',1,[])));
end
end
