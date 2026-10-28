function report=testRocketWallResolutionCase(folder)
%TESTROCKETWALLRESOLUTIONCASE 核对实际完整箭体近壁网格与嵌套二分合同。
% 模块一：母线/原面/物理边界不变；实际节点偏移、体积和共享边重新计算。
% 模块二：所有父节点精确保留，体积/一阶矩分区及真实墙/轴/外域引用逐项检查。
% 模块三：初猜保留五量供体与非保守身份；几何通过不冒充物理解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v77_validation','factory_test');end
assert(~exist(folder,'dir'),'RocketWallResolutionTest:Overwrite','不能覆盖既有测试。');cases=cell(2,1);rows=struct([]);
for k=1:2
    kinds={'collar72','nested144'};c=createRocketAxisymmetricWallResolutionCase(struct('kind',kinds{k}));
    assert(isequaln(c,createRocketAxisymmetricWallResolutionCase(c.options))&& ...
        isequaln(c.mesh,buildAxisymmetricFVGeometry(c.mesh.X,c.mesh.R)));
    ni=c.mesh.shape(1);nr=c.mesh.shape(2);assert(ni==115&&nr==72*k&&c.mesh.cellCount==115*72*k);
    d=c.domain;m=c.mesh;ids=d.wallBoundaryFaceID;
    assert(isequal(m.faces.owner(ids),(1:ni)')&&isequal(m.faces.centroid(ids,:), ...
        casesWallCenters(d.wallPolylineXR))&&all(m.volumePerRadian>0)&& ...
        isequal(m.faces.areaPerRadian(ids),d.wallAreaPerRadianM2));
    axis=[d.axisBoundaryFaceID{1};d.axisBoundaryFaceID{2}];
    assert(all(m.faces.centroid(axis,2)==0)&&all(m.faces.areaPerRadian(axis)==0)&& ...
        numel(m.boundaryFaces)==2*sum(m.shape));
    normals=d.collar.wallOutwardNormalXR;offset=d.collar.normalOffsetM;
    offsetError=0;
    for j=1:numel(offset)
        delta=[m.X(:,j),m.R(:,j)]-d.wallPolylineXR;
        error=[sum(delta(1:end-1,:).*normals,2)-offset(j);sum(delta(2:end,:).*normals,2)-offset(j)];
        offsetError=max(offsetError,max(abs(error)));
    end
    assert(offsetError<2e-13&&min(d.wallNormalDistanceM)>0&&max(d.wallNormalDistanceM)<1.01e-6/k);
    assert(~c.initialization.mapping.conservativeRemappingClaimed&&~c.initialization.seed.externalCFDImported);
    if k==2
        p=cases{1};assert(isequal(m.X(:,1:2:end),p.mesh.X)&&isequal(m.R(:,1:2:end),p.mesh.R)&& ...
            isequaln(c.surfaceMesh,p.surfaceMesh)&&isequal(c.boundary,p.boundary));
        ab=d.radialRefinement.parentToChildren;vol=sum(m.volumePerRadian(ab),2);
        assert(max(abs(vol-p.mesh.volumePerRadian)./p.mesh.volumePerRadian)<3e-10);
    end
    a=c.initialization.mapping;inside=a.insideParentDomain;
    assert(all(inside)&&all(sum(a.donorWeights,2)==1)&&all(a.donorWeights>=0,'all'));
    row=struct('kind',kinds{k},'cellCount',m.cellCount,'radialLayers',nr, ...
        'nearWallNormalLayers',d.collar.options.normalLayers,'maximumNormalOffsetErrorM',offsetError, ...
        'actualWallDistanceM',[min(d.wallNormalDistanceM),max(d.wallNormalDistanceM)]);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    cases{k}=c;
end
reject(@()createRocketAxisymmetricWallResolutionCase(struct('kind','unknown')),'RocketWallResolution:Kind');
reject(@()createRocketAxisymmetricWallResolutionCase(struct('alpha',20)),'RocketWallResolution:Options');
reject(@()createRocketAxisymmetricWallResolutionCase(struct('parent72SHA256',repmat('0',1,64))),'RocketWallResolution:Parent');
report=struct('status','passed-real-complete-wall-grid-and-parent-node-preserving-radial-bisection', ...
    'cases',rows,'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'factory_test_evidence.mat'),'report','cases','-v7');
fid=fopen(fullfile(folder,'factory_test_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('%s\n',jsonencode(report));
end
function center=casesWallCenters(w)
% 模块用途：直线墙段r权面心的独立端点积分表达式，与几何装配无共享索引。
a=w(1:end-1,:);b=w(2:end,:);ra=a(:,2);rb=b(:,2);
center=((2*ra+rb).*a+(ra+2*rb).*b)./(3*(ra+rb));
end
function reject(action,id)
try,action();error('RocketWallResolutionTest:Refusal','应当拒绝。');catch ex,assert(strcmp(ex.identifier,id));end
end
