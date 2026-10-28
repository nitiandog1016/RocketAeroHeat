function value=solveRocketMomentumLines(factor,rhs)
%SOLVEROCKETMOMENTUMLINES 将全部三动量右端映射到原面法向列，再完整还原代数顺序。
% 模块一：不作全局长三角传播，各列缓存逆同时作用；原A的跨列耦合仍由外层作用检查。
% 模块二：固定线性映射，无残差早停、裁值或选择部分CV。
assert(strcmp(factor.kind,'fixed-three-momentum-original-normal-line-inverse'),'MomentumLines:Kind','需要已核对的法向列逆。');
n=factor.cellCount;assert(isequal(size(rhs),[n 3])&&isreal(rhs)&&all(isfinite(rhs),'all'),'MomentumLines:RHS','三动量右端无效。');
normal=rhs(factor.inverseOrder,:);pages=reshape(permute(reshape(normal,factor.sourceFaces,factor.layers,3),[3 2 1]),3*factor.layers,1,factor.sourceFaces);
solved=pagemtimes(factor.inverseLine,pages);original=reshape(permute(reshape(solved,3,factor.layers,factor.sourceFaces),[3 2 1]),n,3);
value=zeros(n,3);value(factor.inverseOrder,:)=original;assert(all(isfinite(value),'all'),'MomentumLines:Nonfinite','法向列作用非有限，拒绝使用。');
end
