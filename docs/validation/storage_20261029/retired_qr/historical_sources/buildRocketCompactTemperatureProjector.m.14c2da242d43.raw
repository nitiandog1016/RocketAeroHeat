function projector=buildRocketCompactTemperatureProjector(factor,rows,representatives,originalToUnique)
%BUILDROCKETCOMPACTTEMPERATUREPROJECTOR 在完整局部支撑上表示原QR坐标投影。
% 模块一：H涉及的CV扩展到其所属完整QR组，不忽略组内任何五量耦合。
% 模块二：逐组形成P转置H转置，仅存支撑坐标；外部坐标上的基严格为零。
% 模块三：保留全部原行及精确重复映射，不截断相关行，不增加448MiB构造预算。
n=factor.cellCount;N=5*n;m=size(rows.H,1);q=numel(representatives);
[~,physicalColumns]=find(rows.H);touched=false(n,1);touched(mod(physicalColumns-1,n)+1)=true;
covered=false(n,1);pages=cell(numel(factor.buckets),1);
for k=1:numel(factor.buckets)
    b=factor.buckets{k};if isempty(b.patchIDs),pages{k}=zeros(0,1);continue;end
    pages{k}=find(any(reshape(touched(b.members(:)),size(b.members)),1));
    ids=double(b.members(:,pages{k}));
    assert(~any(covered(ids(:))),'TemperatureProjection:Support','QR组对约束支撑重复覆盖。');
    covered(ids(:))=true;
end
assert(all(covered(touched)),'TemperatureProjection:Support','温度行支撑未被当前QR分组完整覆盖。');
cells=find(covered);s=5*numel(cells);assert(q<=s,'TemperatureProjection:Rank','约束数超过完整支撑维数。');
% 上界包括四套支撑稠密面板、六套小矩阵、稀疏行副本、索引和逐组临时数组。
% 输入factor/rows单列为调用者持有，本报告不声称进程峰值内存。
bound=8*(4*s*q+6*q*q+8*N+32*m+64*s)+64*nnz(rows.H)+65536;
assert(bound<=448*1024^2,'TemperatureProjection:Budget', ...
    '保留全部原行的紧凑投影超过448MiB额外数组上界，不删约束。');
lookup=zeros(n,1,'uint32');lookup(cells)=uint32(1:numel(cells));
coordinateIDs=reshape(cells+(0:4)*n,[],1);columns=zeros(s,q);maximumPanelBytes=0;
for k=1:numel(factor.buckets)
    b=factor.buckets{k};
    for page=reshape(pages{k},1,[])
        members=double(b.members(:,page));
        original=reshape((members+(0:4)*n)',[],1);
        compact=reshape((double(lookup(members))+(0:4)*numel(cells))',[],1);
        raw=full(rows.H(representatives,original))';
        image=b.map(:,:,page)'*raw;columns(compact,:)=image;
        info=whos('raw','image','original','compact');maximumPanelBytes=max(maximumPanelBytes,sum([info.bytes]));
    end
end
norms=vecnorm(columns,2,1);assert(all(isfinite(norms)&norms>0),'TemperatureProjection:Rank','约束在右坐标中为零或非有限。');
[basis,R]=qr(columns./norms,0);pivots=abs(diag(R));
assert(all(pivots>1e-10),'TemperatureProjection:Rank','显式约束数值相关，拒绝而不截断约束。');
orthogonality=norm(basis'*basis-eye(q),Inf);
assert(orthogonality<=1e-10,'TemperatureProjection:Orthogonality','紧凑投影基未通过原正交检查。');
projector=struct('kind','fixed-current-temperature-row-nullspace-in-QR-right-coordinates','cellCount',n, ...
    'basis',basis,'constraintCount',m,'constraintColumnNorms',norms(originalToUnique),'normalizedQRDiagonal',pivots, ...
    'basisOrthogonalityError',orthogonality,'faceIDs',rows.faceIDs,'homogeneousDerivativeTargetK',zeros(m,1), ...
    'physicalBoundaryConditionAdded',false,'nonlinearTemperatureBoundGuaranteed',false,'newProductionSolvedAreaM2',0, ...
    'independentConstraintCount',q,'exactRepresentativeRowIndices',representatives, ...
    'originalToUniqueRow',originalToUnique,'allOriginalRowsExactlyRepresented',true,'rankTruncationUsed',false, ...
    'basisCoordinateIndices',coordinateIDs,'basisStorage','complete-patch-support-QR-coordinate-basis', ...
    'supportCellIDs',cells,'additionalConstructionArrayBoundBytes',bound,'originalConstructionBudgetBytes',448*1024^2, ...
    'basisArrayBytes',8*numel(basis),'equivalentFullBasisArrayBytes',8*N*q, ...
    'maximumTemporaryPatchPanelBytes',maximumPanelBytes,'processPeakMemoryClaimed',false);
end
