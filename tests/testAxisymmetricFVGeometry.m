function report=testAxisymmetricFVGeometry(folder)
%TESTAXISYMMETRICFVGEOMETRY 核对轴对称真矩、共享面、梯度及拒绝合同。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v62_validation','geometry');end
if ~exist(folder,'dir'),mkdir(folder);end
[X,R]=ndgrid([.2,.13,0],[0,.01,.03,.06]);mesh=buildAxisymmetricFVGeometry(X,R);
[dx,ra]=ndgrid(-diff(X(:,1)),R(1,1:end-1));[~,rb]=ndgrid(-diff(X(:,1)),R(1,2:end));
truthV=dx.*(rb.^2-ra.^2)/2;truthA=dx.*(rb-ra);
truthR=(2/3)*(rb.^3-ra.^3)./(rb.^2-ra.^2);
assert(max(abs(mesh.volumePerRadian-truthV(:)))<1e-18);
assert(max(abs(mesh.planarArea-truthA(:)))<1e-17);
assert(max(abs(mesh.volumeCentroid(:,2)-truthR(:)))<1e-16);
assert(abs(sum(mesh.volumePerRadian)-.2*.06^2/2)<1e-18);
assert(norm(mesh.incidence*mesh.faces.areaVectorPerRadian-[zeros(mesh.cellCount,1),mesh.planarArea],Inf)<1e-16);
pts=[mesh.volumeCentroid;mesh.faces.centroid(mesh.boundaryFaces,:)];field=[ones(size(pts,1),1),pts];
assert(max(abs(mesh.gradientX*field-repmat([0 1 0],mesh.cellCount,1)),[],'all')<1e-12);
assert(max(abs(mesh.gradientR*field-repmat([0 0 1],mesh.cellCount,1)),[],'all')<1e-12);

%% 偏斜贴体控制体也保留真质心和线性精确性，平移不改体积
Xs=X+.13*R+.3*R.^2;skew=buildAxisymmetricFVGeometry(Xs,R);
pts=[skew.volumeCentroid;skew.faces.centroid(skew.boundaryFaces,:)];
assert(max(abs(skew.gradientX*pts-repmat([1 0],skew.cellCount,1)),[],'all')<1e-12);
assert(max(abs(skew.gradientR*pts-repmat([0 1],skew.cellCount,1)),[],'all')<1e-12);
translated=buildAxisymmetricFVGeometry(Xs+.731,R);
assert(max(abs(translated.volumePerRadian-skew.volumePerRadian))<1e-17);
periodic=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));
assert(nnz(periodic.faces.periodic)==size(X,2)-1&&all(periodic.faces.boundarySide~=1&periodic.faces.boundarySide~=2));
assert(abs(sum(periodic.volumePerRadian)-sum(mesh.volumePerRadian))<1e-18);

%% 退化/翻折/负半径/非有限/未知选项不能进入空间算子
rejections=0;bad=R;bad(1,1)=-1;reject(@()buildAxisymmetricFVGeometry(X,bad));
bad=X;bad(1,1)=NaN;reject(@()buildAxisymmetricFVGeometry(bad,R));
reject(@()buildAxisymmetricFVGeometry(complex(X,1),R));reject(@()buildAxisymmetricFVGeometry(X,R(1:end-1,:)));
reject(@()buildAxisymmetricFVGeometry([0 0;0 0],[0 1;0 1]));
reject(@()buildAxisymmetricFVGeometry([0 1;1 0],[0 1;0 1]));
reject(@()buildAxisymmetricFVGeometry(X,R,struct('periodicI',1)));
reject(@()buildAxisymmetricFVGeometry(X,R,struct('unknown',true)));
bad=R;bad(end,:)=bad(end,:)+.002;reject(@()buildAxisymmetricFVGeometry(X,bad,struct('periodicI',true)));
report=struct('status','passed-geometric-identities-only','rejections',rejections, ...
    'perRadianVolume',sum(mesh.volumePerRadian),'maximumMetricClosure',mesh.metricClosureRelative, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false);
save(fullfile(folder,'geometry_evidence.mat'),'mesh','skew','periodic','report');
fid=fopen(fullfile(folder,'geometry_summary.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Axisymmetric geometry passed: true radial moments, metric closure, skew linear gradients, periodic topology and %d refusals.\n',rejections);
    function reject(action)
        try,action();catch,rejections=rejections+1;return;end
        error('AxisFVTest:MissingRefusal','无效几何未被拒绝。');
    end
end
