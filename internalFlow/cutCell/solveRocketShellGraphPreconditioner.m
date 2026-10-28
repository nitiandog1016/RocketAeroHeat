function [update,report]=solveRocketShellGraphPreconditioner(operator,rhs)
%SOLVEROCKETSHELLGRAPHPRECONDITIONER 固定次数对称着色法向线扫描，形成线性右预条件。
% 模块一：从零开始，逐色读取已更新侧向邻居，每条线用缓存块分解完整求解。
% 模块二：扫描次数固定、不按右端或容差早停；可用于现有固定右预条件GMRES。
% 模块三：返回真实图残差但不要求它是精确图解；最终方向仍须完整NS验收。
nf=operator.sourceFaces;nl=operator.layers;n=nf*nl;nc=numel(operator.colorGroups);
assert(isequal(size(rhs),[n 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'ShellGraph:RHS','侧向图右端无效。');
assert(operator.fixedLinearPreconditioner,'ShellGraph:Policy','固定GMRES要求明确固定线性预条件映射。');
update=zeros(n,5);order=[1:nc nc-1:-1:1];
for sweep=1:operator.fixedSymmetricSweeps
    for c=order
        g=operator.colorGroups{c};ids=g.cellIDs;count=numel(g.sourceFaceIDs);b=rhs(ids,:);
        b=subtract(b,g.minusSideIDs,g.minusLocalRows,operator.minusCoupling,operator.sidePlus);
        b=subtract(b,g.plusSideIDs,g.plusLocalRows,operator.plusCoupling,operator.sideMinus);
        solved=zeros(5,1,numel(ids));r=permute(b,[2 3 1]);
        for layer=1:nl
            local=(layer-1)*count+(1:count);globalIDs=ids(local);current=r(:,:,local);
            if layer>1,current=current-pagemtimes(operator.lower(:,:,globalIDs),solved(:,:,local-count));end
            solved(:,:,local)=pagemtimes(operator.lineInverse(:,:,globalIDs),current);
        end
        for layer=nl-1:-1:1
            local=(layer-1)*count+(1:count);globalIDs=ids(local);
            solved(:,:,local)=solved(:,:,local)-pagemtimes(operator.lineUpperFactor(:,:,globalIDs),solved(:,:,local+count));
        end
        update(ids,:)=reshape(solved,5,[])';
    end
end
assert(all(isfinite(update),'all'),'ShellGraph:Nonfinite','侧向图扫描非有限，不提交裁剪结果。');
defect=applyRocketShellGraphOperator(operator,update)-rhs;
report=struct('status','applied-fixed-linear-symmetric-colored-normal-line-graph-sweeps', ...
    'relativeLinearDefect',max(abs(defect),[],'all')/max(1,max(abs(rhs),[],'all')), ...
    'graphRelativeResidual',norm(defect(:))/max(1,norm(rhs(:))), ...
    'colorCount',nc,'fixedSymmetricSweeps',operator.fixedSymmetricSweeps,'fixedLinearMap',true, ...
    'exactGraphSolveClaimed',false,'physicalTimeAdvanced',0,'fullOriginalNSJacobian',false);
    function out=subtract(out,sideIDs,rows,blocks,neighbours)
        for first=1:4096:numel(sideIDs)
            pick=first:min(first+4095,numel(sideIDs));selected=sideIDs(pick);
            value=reshape(pagemtimes(blocks(:,:,selected),permute(update(neighbours(selected),:),[2 3 1])),5,[])';
            for a=1:5,out(:,a)=out(:,a)-accumarray(rows(pick),value(:,a),[size(out,1) 1]);end
        end
    end
end
