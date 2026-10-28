function report=prepareQuadraticCurvedFixtures(folder)
%PREPAREQUADRATICCURVEDFIXTURES 固定有向球面三角形、三档流体网格的精度输入。
% 模块用途：只改变Cartesian控制体分辨率，不把几何变化混入热流误差。
% 球面由320个固定三角面表示，验证对象是该多面曲面上的制造解，不是球体实验。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,folder=fullfile(root,'outputs','v59_validation','geometry');end
if ~isfolder(folder),mkdir(folder);end
maxNumCompThreads(1);addpath(fullfile(root,'internalFlow'),'-begin');setupCutCellFlow;
[V,F]=icosphere(2,.32,[.017 -.023 .011]);counts=[8 12 16];checks=cell(3,1);
for j=1:numel(counts)
    n=counts(j);h=1.2/n;grid=struct('origin',[-.6 -.6 -.6],'spacing',[h h h],'cellCounts',[n n n]);
    fprintf('CURVED_QUADRATIC_GEOMETRY_START grid=%d triangles=%d\n',n,size(F,1));started=tic;
    source=buildOriginalSurfaceCutCells(V,F,grid);
    [group,aggregation]=agglomerateFluidControlVolumes(source,struct('failOnUnresolved',true));
    identity=struct('name','fixed_320_triangle_spherical_polyhedron','gridCount',n,'spacing',h, ...
        'sourceFaces',size(F,1),'fluidGroups',numel(group.cells.volume),'wallPatches',numel(group.walls.area), ...
        'elapsedSeconds',toc(started),'physicalSolution',false,'eligibleForSurfaceCoverage',false);
    save(fullfile(folder,sprintf('sphere_grid%d.mat',n)),'source','group','aggregation','identity','-v7');
    checks{j}=identity;fprintf('CURVED_QUADRATIC_GEOMETRY_SAVED grid=%d groups=%d wall=%d elapsed=%.4g\n', ...
        n,identity.fluidGroups,identity.wallPatches,identity.elapsedSeconds);
end
report=struct('status','prepared-manufactured-curved-grid-inputs','checks',{checks}, ...
    'geometryFixedAcrossGrids',true,'heatComputed',false,'physicalSolution',false,'eligibleForSurfaceCoverage',false);
fid=fopen(fullfile(folder,'geometry_summary.json'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
function [V,F]=icosphere(level,radius,center)
% 模块用途：确定性固定曲面，边中点共享编号并明确修正为实体外法向。
t=(1+sqrt(5))/2;V=[-1 t 0;1 t 0;-1 -t 0;1 -t 0;0 -1 t;0 1 t;0 -1 -t;0 1 -t;t 0 -1;t 0 1;-t 0 -1;-t 0 1];
V=V./vecnorm(V,2,2);F=convhull(V);
for k=1:level
    m=size(F,1);edges=sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2);[uniqueEdges,~,where]=unique(edges,'rows');
    points=.5*(V(uniqueEdges(:,1),:)+V(uniqueEdges(:,2),:));points=points./vecnorm(points,2,2);
    index=size(V,1)+where;V=[V;points];ab=index(1:m);bc=index(m+(1:m));ca=index(2*m+(1:m));
    F=[F(:,1),ab,ca;ab,F(:,2),bc;ca,bc,F(:,3);ab,bc,ca];
end
normal=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
inside=sum(normal.*(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:)),2)<0;
F(inside,[2 3])=F(inside,[3 2]);V=radius*V+center;
end
