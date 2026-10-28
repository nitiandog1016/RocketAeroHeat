function projector=buildRocketTemperatureProjector(factor,rows)
%BUILDROCKETTEMPERATUREPROJECTOR 将采样温度齐次约束转为QR右坐标中的正交投影。
% 模块一：逐项完全相同的H行共用一个表示，保留全部原行与面号；不截断近似相关行。
% 模块二：投影Pi=I-U*U转置，方向v=P*Pi*y满足H*v约为0；不构造百万维稠密Pi。
% 模块三：只限制迭代方向的一阶温度变化，不固定最终物理温度，不保证非线性合法。
assert(strcmp(factor.kind,'fixed-complete-augmented-image-patch-QR-right-map')&& ...
    strcmp(rows.kind,'current-interior-selected-Riemann-sample-temperature-tangent-rows')&& ...
    factor.cellCount==rows.cellCount,'TemperatureProjection:Input','当前QR映射与真实温度行不匹配。');
N=5*factor.cellCount;m=size(rows.H,1);
assert(m>=1&&issparse(rows.H)&&size(rows.H,2)==N&&all(isfinite(nonzeros(rows.H))), ...
    'TemperatureProjection:Input','需要保留全部原面的有限稀疏温度行。');
% 模块四：精确重复是同一线性条件的多个面标签，并非新的独立方向。
representatives=zeros(0,1);originalToUnique=zeros(m,1);Ht=rows.H';
for j=1:m
    found=0;
    for k=1:numel(representatives)
        if isequal(Ht(:,j),Ht(:,representatives(k))),found=k;break;end
    end
    if found==0,representatives(end+1,1)=j;found=numel(representatives);end %#ok<AGROW>
    originalToUnique(j)=found;
end
q=numel(representatives);
assert(isequal(rows.H,rows.H(representatives(originalToUnique),:)), ...
    'TemperatureProjection:ExactDuplicates','原约束无法逐项精确重构，拒绝合并表示。');
clear Ht;
% 模块五：旧小规模路径逐值保留；超过原16行或稠密预算时按完整局部支撑存储。
% 扩展路径仍使用同一P转置、薄QR和1e-10秩门槛，不能靠删行或增加预算通过。
if m>16||3*8*N*q>448*1024^2
    projector=buildRocketCompactTemperatureProjector(factor,rows,representatives,originalToUnique);return;
end
columns=zeros(N,q);shape=[factor.cellCount 5];
for j=1:q
    v=reshape(full(rows.H(representatives(j),:))',shape);v=applyRocketAugmentedPatchQR(factor,v,'transp');columns(:,j)=v(:);
end
norms=vecnorm(columns,2,1);assert(all(isfinite(norms)&norms>0),'TemperatureProjection:Rank','约束在右坐标中为零或非有限。');
[basis,R]=qr(columns./norms,0);pivots=abs(diag(R));
assert(all(pivots>1e-10),'TemperatureProjection:Rank','显式约束数值相关，拒绝而不截断约束。');
orthogonality=norm(basis'*basis-eye(q),Inf);assert(orthogonality<=1e-10,'TemperatureProjection:Orthogonality','投影基未通过正交检查。');
projector=struct('kind','fixed-current-temperature-row-nullspace-in-QR-right-coordinates','cellCount',factor.cellCount, ...
    'basis',basis,'constraintCount',m,'constraintColumnNorms',norms(originalToUnique),'normalizedQRDiagonal',pivots, ...
    'basisOrthogonalityError',orthogonality,'faceIDs',rows.faceIDs,'homogeneousDerivativeTargetK',zeros(m,1), ...
    'physicalBoundaryConditionAdded',false,'nonlinearTemperatureBoundGuaranteed',false,'newProductionSolvedAreaM2',0, ...
    'independentConstraintCount',q,'exactRepresentativeRowIndices',representatives, ...
    'originalToUniqueRow',originalToUnique,'allOriginalRowsExactlyRepresented',true,'rankTruncationUsed',false);
end
