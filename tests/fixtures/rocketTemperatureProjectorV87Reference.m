function projector=rocketTemperatureProjectorV87Reference(factor,rows)
%BUILDROCKETTEMPERATUREPROJECTOR 将采样温度齐次约束转为QR右坐标中的正交投影。
% 模块一：对每行H形成P转置乘H转置，再作薄QR，保留全部显式约束。
% 模块二：投影Pi=I-U*U转置，方向v=P*Pi*y满足H*v约为0；不构造百万维稠密Pi。
% 模块三：只限制迭代方向的一阶温度变化，不固定最终物理温度，不保证非线性合法。
assert(strcmp(factor.kind,'fixed-complete-augmented-image-patch-QR-right-map')&& ...
    strcmp(rows.kind,'current-interior-selected-Riemann-sample-temperature-tangent-rows')&& ...
    factor.cellCount==rows.cellCount,'TemperatureProjection:Input','当前QR映射与真实温度行不匹配。');
N=5*factor.cellCount;m=size(rows.H,1);assert(m>=1&&m<=16&&3*8*N*m<=448*1024^2, ...
    'TemperatureProjection:Budget','约束向量构造超声明预算，不自动减少约束。');
columns=zeros(N,m);shape=[factor.cellCount 5];
for j=1:m
    v=reshape(full(rows.H(j,:))',shape);v=applyRocketAugmentedPatchQR(factor,v,'transp');columns(:,j)=v(:);
end
norms=vecnorm(columns,2,1);assert(all(isfinite(norms)&norms>0),'TemperatureProjection:Rank','约束在右坐标中为零或非有限。');
[basis,R]=qr(columns./norms,0);pivots=abs(diag(R));
assert(all(pivots>1e-10),'TemperatureProjection:Rank','显式约束数值相关，拒绝而不截断约束。');
orthogonality=norm(basis'*basis-eye(m),Inf);assert(orthogonality<=1e-10,'TemperatureProjection:Orthogonality','投影基未通过正交检查。');
projector=struct('kind','fixed-current-temperature-row-nullspace-in-QR-right-coordinates','cellCount',factor.cellCount, ...
    'basis',basis,'constraintCount',m,'constraintColumnNorms',norms,'normalizedQRDiagonal',pivots, ...
    'basisOrthogonalityError',orthogonality,'faceIDs',rows.faceIDs,'homogeneousDerivativeTargetK',zeros(m,1), ...
    'physicalBoundaryConditionAdded',false,'nonlinearTemperatureBoundGuaranteed',false,'newProductionSolvedAreaM2',0);
end
