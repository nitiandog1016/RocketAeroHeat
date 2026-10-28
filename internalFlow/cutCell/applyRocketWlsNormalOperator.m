function value=applyRocketWlsNormalOperator(operator,update)
%APPLYROCKETWLSNORMALOPERATOR 应用含第二法向邻居的研究预条件矩阵B。
% 模块一：原图全部侧向块保留；额外两个法向带只属于预条件近似。
% 模块二：不能用该作用替代当前完整J或原NS总方向验收。
value=applyRocketShellGraphOperator(operator,update);nf=operator.sourceFaces;nl=operator.layers;
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);
    if layer>2,value(ids,:)=value(ids,:)+reshape(pagemtimes(operator.secondLower(:,:,ids),permute(update(ids-2*nf,:),[2 3 1])),5,[])';end
    if layer<nl-1,value(ids,:)=value(ids,:)+reshape(pagemtimes(operator.secondUpper(:,:,ids),permute(update(ids+2*nf,:),[2 3 1])),5,[])';end
end
assert(all(isfinite(value),'all'),'WlsNormal:Nonfinite','五带预条件作用非有限。');
end
