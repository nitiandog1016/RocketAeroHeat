function [update,report]=solveRocketShellLineSystem(operator,rhs)
%SOLVEROCKETSHELLLINESYSTEM 批量块Thomas求解每原面的完整五量法向线。
% 模块一：仅求近似预条件系统，法向线中保留5x5块耦合，不称物理时间步。
% 模块二：每页5x5由带主元稠密求解，随后以原线算子重新核对线性缺陷。
nf=operator.sourceFaces;nl=operator.layers;n=nf*nl;
assert(isequal(size(rhs),[n 5])&&all(isfinite(rhs),'all'),'ShellLine:RHS','线系统右端无效。');
D=operator.diagonal;L=operator.lower;U=operator.upper;
r=permute(rhs,[2 3 1]);factor=zeros(size(U));solved=zeros(size(r));
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);d=D(:,:,ids);b=r(:,:,ids);
    if layer>1
        previous=ids-nf;d=d-pagemtimes(L(:,:,ids),factor(:,:,previous));b=b-pagemtimes(L(:,:,ids),solved(:,:,previous));
    end
    combined=pagemldivide(d,cat(2,U(:,:,ids),b));
    assert(all(isfinite(combined),'all'),'ShellLine:Singular','法向线分解失败，不补对角或裁结果。');
    factor(:,:,ids)=combined(:,1:5,:);solved(:,:,ids)=combined(:,6,:);
end
for layer=nl-1:-1:1
    ids=(layer-1)*nf+(1:nf);solved(:,:,ids)=solved(:,:,ids)-pagemtimes(factor(:,:,ids),solved(:,:,ids+nf));
end
update=reshape(solved,5,n)';check=pagemtimes(D,solved);
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);
    if layer>1,check(:,:,ids)=check(:,:,ids)+pagemtimes(L(:,:,ids),solved(:,:,ids-nf));end
    if layer<nl,check(:,:,ids)=check(:,:,ids)+pagemtimes(U(:,:,ids),solved(:,:,ids+nf));end
end
check=reshape(check,5,n)';defect=max(abs(check-rhs),[],'all')/max(1,max(abs(rhs),[],'all'));
assert(defect<=1e-8,'ShellLine:LinearDefect','近似线系统未满足明确线性缺陷标准。');
report=struct('status','solved-independent-five-quantity-normal-lines','relativeLinearDefect',defect, ...
    'lineCount',nf,'unknownsPerLine',5*nl,'physicalTimeAdvanced',0,'fullOriginalNSJacobian',false);
end
