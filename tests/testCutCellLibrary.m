function report=testCutCellLibrary(folder)
%TESTCUTCELLLIBRARY 核对完整实验库目录解析、源迁移和冻结数值一致性。
% 模块用途：仅验证目录迁移，不重复几何构造或六步基础长回归。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1,folder=fullfile(root,'outputs','v59_validation','library');end
if ~isfolder(folder),mkdir(folder);end
maxNumCompThreads(1);addpath(fullfile(root,'internalFlow'),'-begin');library=setupCutCellFlow;
names={'partitionOrientedTriangles','buildOriginalSurfaceCutCells','agglomerateFluidControlVolumes', ...
    'cutCellNavierStokesResidual','cutCellBoundaryStencilState','buildCutCellDiffusionStencil', ...
    'advanceCutCellNavierStokes','mapCutCellWallFields'};resolved=struct;
for j=1:numel(names)
    actual=which(names{j});assert(strcmpi(actual,fullfile(library,[names{j},'.m'])), ...
        'CutLibrary:Shadow','%s没有从新完整库解析。',names{j});resolved.(names{j})=actual;
end
assert(strcmpi(which('resolveNSCharacteristicBoundary3D'),fullfile(root,'internalFlow','resolveNSCharacteristicBoundary3D.m')));
base=fullfile(root,'outputs','v58_validation');
if ~isfile(fullfile(base,'cut_cell_draft','cut_cell_geometry_results.mat')),base=fullfile(root,'validation','v58');end
d=load(fullfile(base,'cut_cell_draft','cut_cell_geometry_results.mat'),'actual');c=d.actual;
d=load(fullfile(base,'cut_cell_flow_draft','cut_cell_flow_evidence.mat'),'actualArbitrary','report');gas=d.report.testGas;
wall=struct('type','wall','temperature',303,'velocity',[0 0 0]);bc=struct;
for key={'xmin','xmax','ymin','ymax','zmin','zmax','solidWall'},bc.(key{1})=wall;end
[rate,info]=cutCellNavierStokesResidual(d.actualArbitrary.Q,c,gas,bc);
fluid=c.cells.volume>0;rateDifference=max(abs(rate(fluid,:)-d.actualArbitrary.rate(fluid,:)),[],'all');
assert(isequaln(rate,d.actualArbitrary.rate),'CutLibrary:Numeric','路径整理后实际任意场残差不再逐值相同。');
assert(isequaln(info.wallFlux,d.actualArbitrary.info.wallFlux));
mapped=mapCutCellWallFields(c,info,false);assert(~mapped.eligibleForSurfaceCoverage);
% 迁移后的边界helper及批量模板也必须实际解析到新库，而非只测文件存在。
a=load(fullfile(base,'cut_cell_advance_draft','cut_cell_advance_evidence.mat'),'actualRun','actualCut','air','actualBC');
b=cutCellBoundaryStencilState(a.actualRun.initialQ,a.actualCut,a.air,a.actualBC);
stencil=buildCutCellDiffusionStencil(a.actualCut,b);
reference=load(fullfile(base,'cut_cell_advance_draft','actual_stencil_evidence.mat'),'actualStencil');
assert(isequaln(stencil,reference.actualStencil),'CutLibrary:Stencil','迁移后完整初态模板不同。');
report=struct('status','tested-complete-experimental-library-migration','resolved',resolved, ...
    'actualRateMaximumDifference',rateDifference,'actualResidualBitwiseIdentical',true, ...
    'actualWallFluxBitwiseIdentical',true,'actualStencilBitwiseIdentical',true,'physicalSolution',false,'eligibleForSurfaceCoverage',false);
save(fullfile(folder,'library_integration.mat'),'report','mapped','-v7');
fid=fopen(fullfile(folder,'library_integration.json'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('CUT_LIBRARY_MIGRATION_PASS bitwise actual residual and original-wall flux; %d complete functions.\n',numel(names));
end
