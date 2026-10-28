function value=applyRocketShellGraphOperator(operator,update)
%APPLYROCKETSHELLGRAPHOPERATOR 应用已装配的近似B=I/tau-Jgraph，核对预条件残差。
% 模块一：五量法向块与同层左右邻居均保留；不调用或替代正式NS。
% 模块二：按真实所有者和缩放后5x5块累加，结果仅用于预条件器诊断。
nf=operator.sourceFaces;nl=operator.layers;n=nf*nl;
assert(isequal(size(update),[n 5])&&all(isfinite(update),'all'),'ShellGraph:RHS','图算子输入无效。');
pages=permute(update,[2 3 1]);product=pagemtimes(operator.diagonal,pages);
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);
    if layer>1,product(:,:,ids)=product(:,:,ids)+pagemtimes(operator.lower(:,:,ids),pages(:,:,ids-nf));end
    if layer<nl,product(:,:,ids)=product(:,:,ids)+pagemtimes(operator.upper(:,:,ids),pages(:,:,ids+nf));end
end
value=reshape(product,5,n)';m=operator.sideMinus;p=operator.sidePlus;
for first=1:4096:numel(m)
    pick=first:min(first+4095,numel(m));left=reshape(pagemtimes(operator.minusCoupling(:,:,pick),pages(:,:,p(pick))),5,[])';
    right=reshape(pagemtimes(operator.plusCoupling(:,:,pick),pages(:,:,m(pick))),5,[])';
    [owners,~,group]=unique([m(pick);p(pick)]);
    for a=1:5,value(owners,a)=value(owners,a)+accumarray(group,[left(:,a);right(:,a)],[numel(owners) 1]);end
end
assert(all(isfinite(value),'all'),'ShellGraph:Nonfinite','图预条件算子作用非有限。');
end
