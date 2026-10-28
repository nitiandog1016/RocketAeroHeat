function cache=buildCutCellGradientOperator(cut,boundaryState)
%BUILDCELLGRADIENTOPERATOR 仅装配原面积加权WLS梯度及开放边界样本提升，不沿用旧墙应力矩阵。
% 模块一：约束次序、面积/位移平方权重及3x3求逆与原模板相同。
% 模块二：C作用于内部phi扰动；提升项允许特征面样本随当前状态变化。
% 模块三：缓存绑定完整几何和入流分支，不能在修改网格后继续套用。
v=cut.cells.volume(:);c=cut.cells.centroid;n=numel(v);f=cut.openFaces;w=cut.walls;
assert(all(v>0),'CutDirection:Geometry','当前方向模块只支持全部正体积控制体。');
inside=f.area>0&f.cellMinus>0&f.cellPlus>0;
outside=f.area>0&xor(f.cellMinus>0,f.cellPlus>0)&boundaryState.incoming;
capacity=4*nnz(inside)+nnz(outside)+nnz(w.ownerCell>0);
A=zeros(3,3,n);row=zeros(capacity,1);col=row;value=zeros(capacity,3);cursor=0;
for k=find(f.area>0)'
    m=f.cellMinus(k);p=f.cellPlus(k);
    if m>0&&p>0,add(m,p,c(p,:)-c(m,:),f.area(k));add(p,m,c(m,:)-c(p,:),f.area(k));
    elseif boundaryState.incoming(k),id=max(m,p);add(id,0,f.centroid(k,:)-c(id,:),f.area(k));end
end
for k=find(w.ownerCell>0)',id=w.ownerCell(k);add(id,0,w.centroid(k,:)-c(id,:),w.area(k));end
assert(cursor==capacity,'CutDirection:Constraint','梯度约束计数不符。');
H=cell(1,3);C=cell(1,3);inverse=zeros(3,3,n);
for d=1:3,H{d}=sparse(row,col,value(:,d),n,n);C{d}=sparse(n,n);end
for id=1:n
    assert(rcond(A(:,:,id))>1e-13,'CutDirection:GradientRank','原WLS约束退化。');inverse(:,:,id)=A(:,:,id)\eye(3);
end
for d=1:3,for e=1:3,C{d}=C{d}+spdiags(reshape(inverse(d,e,:),n,1),0,n,n)*H{e};end,end
ids=find(outside);owner=max(f.cellMinus(ids),f.cellPlus(ids));delta=f.centroid(ids,:)-c(owner,:);
weight=f.area(ids)./sum(delta.^2,2);coefficient=reshape(pagemtimes(inverse(:,:,owner),permute(weight.*delta,[2 3 1])),3,[])';
cache=struct('geometry',cut,'incoming',boundaryState.incoming,'wall',boundaryState.wall,'gradientOperator',{C}, ...
    'incomingFaceIDs',ids,'incomingOwner',owner,'boundaryCoefficientXYZ',coefficient, ...
    'oldStressOrHeatMatrixReused',false,'givenBoundarySampleSensitivitySupported',true,'fullNSJacobianClaimed',false);
    function add(id,other,delta,area)
        d2=dot(delta,delta);assert(isfinite(d2)&&d2>0,'CutDirection:Distance','WLS约束位移为零。');
        weight=area/d2;A(:,:,id)=A(:,:,id)+weight*(delta'*delta);
        cursor=cursor+1;row(cursor)=id;col(cursor)=id;value(cursor,:)=-weight*delta;
        if other>0,cursor=cursor+1;row(cursor)=id;col(cursor)=other;value(cursor,:)=weight*delta;end
    end
end
