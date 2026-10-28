function operator=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,pseudoDt,sweeps)
%BUILDROCKETSHELLGRAPHPRECONDITIONER 装配保留侧向邻居的原面层域预条件器。
% 模块一：原法向线保留，全部内部Euler采用当前精确波系导数；正式NS不改变。
% 模块二：侧面左右5x5块批量保存；对角同步替换，不能只改非对角漏改本格。
% 模块三：原面邻接确定性着色，缓存法向线分解；固定对称扫描保持线性映射。
% 黏性仍为冻结两点近似，外边界仍为原储库近似，不称完整NS Jacobian。
validateattributes(sweeps,{'numeric'},{'real','finite','scalar','integer','positive'});
operator=buildRocketShellLineOperator(Q,cut,gas,boundary,scale,pseudoDt,'exact-branch');
scale=reshape(scale,1,5);nf=operator.sourceFaces;nl=operator.layers;n=nf*nl;
f=cut.openFaces;v=cut.cells.volume(:);sideIDs=find(f.cellPlus>0&f.kind==3);
m=double(f.cellMinus(sideIDs));p=double(f.cellPlus(sideIDs));ns=numel(sideIDs);
minusCoupling=zeros(5,5,ns);plusCoupling=minusCoupling;batch=2048;sampleCounts=zeros(6,1);
for first=1:batch:ns
    pick=first:min(first+batch-1,ns);ids=sideIDs(pick);mi=m(pick);pi=p(pick);
    normal=f.normal(ids,:);distance=sum((cut.cells.centroid(pi,:)-cut.cells.centroid(mi,:)).*normal,2);
    [DL,DR]=cutCellLineFaceLinearization(Q(mi,:),Q(pi,:),normal,distance,gas);
    [EL,ER,audit]=cutCellExactLineFaceLinearization(Q(mi,:),Q(pi,:),normal,distance,gas);
    sampleCounts=sampleCounts+accumarray(audit.sampleRegion,1,[6 1]);
    am=f.area(ids)./v(mi);ap=f.area(ids)./v(pi);[owners,~,group]=unique([mi;pi]);
    for a=1:5,for b=1:5
        ratio=scale(b)/scale(a);left=reshape(EL(a,b,:),[],1);right=reshape(ER(a,b,:),[],1);
        change=[am.*(left-reshape(DL(a,b,:),[],1));-ap.*(right-reshape(DR(a,b,:),[],1))]*ratio;
        operator.diagonal(a,b,owners)=operator.diagonal(a,b,owners)+reshape(accumarray(group,change,[numel(owners) 1]),1,1,[]);
        minusCoupling(a,b,pick)=reshape(am.*right*ratio,1,1,[]);
        plusCoupling(a,b,pick)=reshape(-ap.*left*ratio,1,1,[]);
    end,end
end
assert(all(isfinite(minusCoupling),'all')&&all(isfinite(plusCoupling),'all')&&all(isfinite(operator.diagonal),'all'), ...
    'ShellGraph:Nonfinite','侧向预条件块非有限，不裁剪或补对角。');

%% 模块四：只使用第一层原面邻接着色；每种颜色内不存在侧向耦合。
originalM=mod(m-1,nf)+1;originalP=mod(p-1,nf)+1;
assert(all(originalM~=originalP),'ShellGraph:Topology','侧面不能连接同一原面法向线。');
firstLayer=m<=nf&p<=nf;edges=unique(sort([originalM(firstLayer) originalP(firstLayer)],2),'rows');
adjacency=sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],1,nf,nf);color=zeros(nf,1);
for id=1:nf
    occupied=color(find(adjacency(id,:)));choice=1;while any(occupied==choice),choice=choice+1;end
    color(id)=choice;
end
assert(all(color(originalM)~=color(originalP)),'ShellGraph:Coloring','同色原面有侧向邻居，不能并行法向线扫描。');
groups=cell(max(color),1);
for c=1:numel(groups)
    faces=find(color==c);cells=reshape(faces+nf*(0:nl-1),[],1);map=zeros(n,1);map(cells)=1:numel(cells);
    left=find(color(originalM)==c);right=find(color(originalP)==c);
    groups{c}=struct('sourceFaceIDs',faces,'cellIDs',cells,'minusSideIDs',left,'plusSideIDs',right, ...
        'minusLocalRows',map(m(left)),'plusLocalRows',map(p(right)));
end

%% 模块五：缓存每条线的Schur块逆及上块，正式右预条件不按残差阈值早停。
inverse=zeros(5,5,n);upperFactor=inverse;identity=repmat(eye(5),1,1,nf);maximumFactorDefect=0;
for layer=1:nl
    ids=(layer-1)*nf+(1:nf);D=operator.diagonal(:,:,ids);
    if layer>1,D=D-pagemtimes(operator.lower(:,:,ids),upperFactor(:,:,ids-nf));end
    combined=pagemldivide(D,cat(2,identity,operator.upper(:,:,ids)));
    assert(all(isfinite(combined),'all'),'ShellGraph:Singular','侧向修正后的法向线分解失败，不补对角。');
    inverse(:,:,ids)=combined(:,1:5,:);upperFactor(:,:,ids)=combined(:,6:10,:);
    defect=max(abs(pagemtimes(D,inverse(:,:,ids))-identity),[],'all');maximumFactorDefect=max(maximumFactorDefect,defect);
end
assert(maximumFactorDefect<1e-8,'ShellGraph:FactorDefect','缓存法向线分解未达到明确数值标准。');
operator.kind='original-facet-exact-interior-Euler-frozen-viscous-colored-line-graph-preconditioner';
operator.sideFaceIDs=sideIDs;operator.sideMinus=m;operator.sidePlus=p;
operator.minusCoupling=minusCoupling;operator.plusCoupling=plusCoupling;
operator.color=color;operator.colorGroups=groups;operator.lineInverse=inverse;operator.lineUpperFactor=upperFactor;
operator.fixedSymmetricSweeps=sweeps;operator.fixedLinearPreconditioner=true;operator.maximumFactorDefect=maximumFactorDefect;
operator.sideEulerLinearization='exact-branch';operator.exactSideSampleRegionCounts=sampleCounts;
operator.sideAndExternalEulerLinearization='side-exact-branch-external-rusanov-frozen';
operator.internalViscousLinearization='frozen-two-point';operator.externalStateSensitivityIncluded=false;
end
