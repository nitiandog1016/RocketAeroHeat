function report=testRocketBlockFlowOrdering(folder)
%TESTROCKETBLOCKFLOWORDERING 小域独立稠密C/确定性流向排序和块三角解核对。
% 模块一：从两种原排序重排到相同给定流向，保留原C全部系数与物理作用。
% 模块二：分解后的前后代与完整L/U参考一致；同投影CV按原编号打破并列。
% 模块三：错误质心、方向、预算和二次重排拒绝，所有测试无新J/NS。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'BlockFlowOrderTest:Overwrite','不得覆盖流向排序测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','blockiluunit01','block_ilu_unit_evidence.mat');s=load(sourceFile);C=s.C;
g=load(fullfile(root,'outputs','v83_validation','wlsnormalunit02','wls_normal_unit_evidence.mat'),'cut');centroid=g.cut.cells.centroid;clear g;
direction=[-1 .35 .07];expectedProjection=centroid*(direction'/norm(direction));[~,expectedOrder]=sortrows([expectedProjection (1:size(centroid,1))'],[1 2]);
cases=cell(2,1);errors=zeros(2,2);
for j=1:2
    original=s.cases{j}.packed;before=original;ordered=reorderRocketBlockGraphByFlow(original,centroid,direction);
    assert(isequal(ordered.order,expectedOrder)&&isequaln(original,before));
    assert(isequaln(ordered.values,original.values(:,:,ordered.flowReorderingReport.newToOldBlockPosition)));
    probe=s.cases{j}.probe;applied=applyRocketPackedGraph(ordered,probe);errors(j,1)=norm(applied(:)-C*probe(:))/max(1,norm(C*probe(:)));
    factor=factorRocketBlockILU0(ordered);[L,U]=matrices(factor);rhs=s.cases{j}.rhs;solution=solveRocketBlockILU0(factor,rhs);
    physical=solution(factor.order,:)';right=rhs(factor.order,:)';reference=U\(L\right(:));errors(j,2)=norm(physical(:)-reference)/max(1,norm(reference));assert(max(errors(j,:))<1e-9);
    cases{j}=struct('original',original,'ordered',ordered,'factor',factor,'L',L,'U',U,'rhs',rhs,'solution',solution,'probe',probe,'applied',applied);
end
assert(isequal(cases{1}.ordered.order,cases{2}.ordered.order)&&isequal(cases{1}.ordered.column,cases{2}.ordered.column));
tied=reorderRocketBlockGraphByFlow(original,zeros(size(centroid)),direction);assert(isequal(tied.order,(1:size(centroid,1))'));
rejected=false(6,1);
try,reorderRocketBlockGraphByFlow(original,centroid(1:end-1,:),direction);catch ex,rejected(1)=strcmp(ex.identifier,'BlockFlowOrder:Geometry');end
badCentroid=centroid;badCentroid(1)=NaN;try,reorderRocketBlockGraphByFlow(original,badCentroid,direction);catch ex,rejected(2)=strcmp(ex.identifier,'BlockFlowOrder:Geometry');end
try,reorderRocketBlockGraphByFlow(original,centroid,[0 0 0]);catch ex,rejected(3)=strcmp(ex.identifier,'BlockFlowOrder:Direction');end
bad=original;bad.options.maximumConstructionBytes=1;try,reorderRocketBlockGraphByFlow(bad,centroid,direction);catch ex,rejected(4)=strcmp(ex.identifier,'BlockFlowOrder:Budget');end
try,reorderRocketBlockGraphByFlow(ordered,centroid,direction);catch ex,rejected(5)=strcmp(ex.identifier,'BlockFlowOrder:Kind');end
bad=original;bad.options.maximumPersistentBytes=original.packingReport.persistentFactorStorageUpperBoundBytes+4*numel(original.column)-1;
try,reorderRocketBlockGraphByFlow(bad,centroid,direction);catch ex,rejected(6)=strcmp(ex.identifier,'BlockFlowOrder:Budget');end;assert(all(rejected));
report=struct('status','passed-small-flow-projection-ordering-all-C-values-triangular-solves-ties-and-six-rejections', ...
    'relativeErrors',errors,'rejected',rejected,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'block_flow_order_unit_evidence.mat'),'report','sourceFile','C','centroid','direction','cases','expectedOrder','tied','-v7');
fid=fopen(fullfile(folder,'block_flow_order_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function [L,U]=matrices(factor)
n=factor.cellCount;L=eye(5*n);U=zeros(5*n);
for row=1:n
    for entry=double(factor.rowPointer(row)):double(factor.rowPointer(row+1))-1
        col=double(factor.column(entry));rr=(row-1)*5+(1:5);cc=(col-1)*5+(1:5);
        if col<row,L(rr,cc)=factor.values(:,:,entry);else,U(rr,cc)=factor.values(:,:,entry);end
    end
end
end
