function [update,report]=solveRocketWlsNormalPreconditioner(operator,rhs)
%SOLVEROCKETWLSNORMALPRECONDITIONER 固定着色扫描，线内保留真实WLS五带耦合。
% 模块一：从零开始按固定次数对称着色，不按右端或残差提前终止。
% 模块二：前代含两个下带，回代含两个上带；同色各原面线互不耦合。
% 模块三：返回B的诊断残差，不授予原NS接受、物理稳态或覆盖资格。
nf=operator.sourceFaces;nl=operator.layers;n=nf*nl;nc=numel(operator.colorGroups);factor=operator.wlsNormalFactor;
assert(isequal(size(rhs),[n 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'WlsNormal:RHS','五量右端无效。');
assert(operator.fixedLinearPreconditioner,'WlsNormal:Policy','要求固定线性映射。');
update=zeros(n,5);order=[1:nc nc-1:-1:1];
for sweep=1:operator.fixedSymmetricSweeps
    for color=order
        g=operator.colorGroups{color};ids=g.cellIDs;count=numel(g.sourceFaceIDs);b=rhs(ids,:);
        b=subtract(b,g.minusSideIDs,g.minusLocalRows,operator.minusCoupling,operator.sidePlus);
        b=subtract(b,g.plusSideIDs,g.plusLocalRows,operator.plusCoupling,operator.sideMinus);
        solved=zeros(5,1,numel(ids));r=permute(b,[2 3 1]);
        for layer=1:nl
            local=(layer-1)*count+(1:count);globalIDs=ids(local);current=r(:,:,local);
            if layer>2,current=current-pagemtimes(operator.secondLower(:,:,globalIDs),solved(:,:,local-2*count));end
            if layer>1,current=current-pagemtimes(factor.effectiveLower1(:,:,globalIDs),solved(:,:,local-count));end
            solved(:,:,local)=pagemtimes(factor.inverse(:,:,globalIDs),current);
        end
        for layer=nl-1:-1:1
            local=(layer-1)*count+(1:count);globalIDs=ids(local);
            solved(:,:,local)=solved(:,:,local)-pagemtimes(factor.upper1(:,:,globalIDs),solved(:,:,local+count));
            if layer<nl-1,solved(:,:,local)=solved(:,:,local)-pagemtimes(factor.upper2(:,:,globalIDs),solved(:,:,local+2*count));end
        end
        update(ids,:)=reshape(solved,5,[])';
    end
end
assert(all(isfinite(update),'all'),'WlsNormal:Nonfinite','五带固定扫描非有限。');
defect=applyRocketWlsNormalOperator(operator,update)-rhs;
report=struct('status','applied-fixed-colored-WLS-normal-five-band-line-sweeps','graphRelativeResidual',norm(defect,'fro')/max(1,norm(rhs,'fro')), ...
    'relativeLinearDefect',max(abs(defect),[],'all')/max(1,max(abs(rhs),[],'all')),'fixedLinearMap',true, ...
    'colorCount',nc,'fixedSymmetricSweeps',operator.fixedSymmetricSweeps,'exactGraphSolveClaimed',false, ...
    'fullOriginalNSJacobian',false,'physicalTimeAdvanced',0,'newProductionSolvedAreaM2',0);
    function out=subtract(out,sideIDs,rows,blocks,neighbours)
        for first=1:4096:numel(sideIDs)
            pick=first:min(first+4095,numel(sideIDs));selected=sideIDs(pick);
            value=reshape(pagemtimes(blocks(:,:,selected),permute(update(neighbours(selected),:),[2 3 1])),5,[])';
            for a=1:5,out(:,a)=out(:,a)-accumarray(rows(pick),value(:,a),[size(out,1) 1]);end
        end
    end
end
