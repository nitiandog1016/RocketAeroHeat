function report=testRocketFullFGraph(folder)
%TESTROCKETFULLFGRAPH 两排序稠密参考核对全部F子块、已有局部逆及固定C值回代。
% 模块一：所有跨组F-F原系数精确保留，不以局部块近似替换完整A_FF。
% 模块二：子系统求解/物理Arnoldi与原C回代分开验收，零右端和预算单列。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FullFTest:Overwrite','不得覆盖完整F测试。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','spatialreductionunit01','spatial_reduction_unit_evidence.mat');s=load(parentFile,'cases');cases=cell(2,1);
for j=1:2
    source=s.cases{j};packed=source.packed;split=source.reduction.split;before=packed;
    options=struct('maximumGraphBytes',8*1024^2,'maximumConstructionBytes',64*1024^2);[graph,extraction]=extractRocketFullFGraph(packed,split,options);
    assert(isequaln(before,packed));A=densePacked(graph);assert(isequal(A,source.A));nf=split.fineCells;n=split.cellCount;
    fixedC=source.value(split.coarseOriginalCVs,:);bF=source.rhs(split.fineOriginalCVs,:)-applyRocketRectangularCoupling(source.reduction.B,fixedC);
    opt=struct('maximumIterations',5*nf,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
    [value,detail]=solveRocketFullFKrylov(graph,source.reduction.Ffactor,bF,opt);direct=A\bF(:);
    assert(detail.krylov.flag==0&&detail.relativeFullFResidual<1.1e-10&&norm(value(:)-direct)/max(1,norm(direct))<1e-8);
    complete=zeros(n,5);complete(split.coarseOriginalCVs,:)=fixedC;complete(split.fineOriginalCVs,:)=value;parentImage=applyRocketPackedGraph(packed,complete);
    parentResidual=norm(source.rhs-parentImage,'fro')/norm(source.rhs,'fro');assert(parentResidual<1e-8);
    compact=compactRocketFullFReport(detail);raw=jsonencode(compact);assert(~contains(raw,'"V":')&&~contains(raw,'"Z":')&&numel(raw)<128*1024);
    [zero,z]=solveRocketFullFKrylov(graph,source.reduction.Ffactor,zeros(nf,5),opt);assert(all(zero==0,'all')&&z.FApplications==0&&z.localFApplications==0);
    bad=options;bad.maximumGraphBytes=1;rejected=false;try,extractRocketFullFGraph(packed,split,bad);catch ex,rejected=strcmp(ex.identifier,'FullFGraph:Budget');end;assert(rejected);
    bad=opt;bad.maximumBasisBytes=1;rejected=false;try,solveRocketFullFKrylov(graph,source.reduction.Ffactor,bF,bad);catch ex,rejected=strcmp(ex.identifier,'RecycledGMRES:Budget');end;assert(rejected);
    cases{j}=struct('packed',packed,'split',split,'Ffactor',source.reduction.Ffactor,'graph',graph,'extraction',extraction,'A',A,'P',source.MF, ...
        'rhs',bF,'fixedC',fixedC,'value',value,'detail',detail,'parentRhs',source.rhs,'complete',complete,'parentImage',parentImage,'parentResidual',parentResidual);
end
report=struct('status','passed-two-orderings-exact-all-FF-extraction-full-F-solve-and-fixed-C-parent-reference', ...
    'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'full_F_unit_evidence.mat'),'report','cases','parentFile','-v7');fid=fopen(fullfile(folder,'full_F_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function C=densePacked(packed)
n=packed.cellCount;C=zeros(5*n);
for row=1:n,rr=packed.order(row)+(0:4)*n;for edge=double(packed.rowPointer(row)):double(packed.rowPointer(row+1))-1,cc=packed.order(packed.column(edge))+(0:4)*n;C(rr,cc)=packed.values(:,:,edge);end,end
end
