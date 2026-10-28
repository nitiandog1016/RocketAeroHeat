function out=quadraticTemperatureFourierFlux(meanField,cut,conductivity,boundary,stencil)
%QUADRATICTEMPERATUREFOURIERFLUX 二次FV标量温度的共享Fourier通量实验草稿。
% 模块用途：在面质心求两侧重构梯度的同一平均，一面一通量保守组装。
% 墙梯度直接为g+H*delta；不能再套旧线性secant修正，否则破坏Cov一致性。
% 标量扩散不是完整NS温度闭合；没有时间推进、稳定/单调或物性/湍流证明。
validateattributes(conductivity,{'numeric'},{'real','finite','scalar','positive'});
r=evaluateQuadraticTemperature(meanField,cut,boundary,stencil);f=cut.openFaces;w=cut.walls;
n=numel(cut.cells.volume);rateIntegral=zeros(n,1);faceFlux=nan(numel(f.area),1);wallHeat=nan(numel(w.area),1);
wallKind=repmat({'inactive-no-domain-fluid'},numel(w.area),1);boundaryFlux=0;
for k=find(f.area>0)'
    m=f.cellMinus(k);p=f.cellPlus(k);point=f.centroid(k,:);normal=f.normal(k,:);
    if m>0&&p>0,G=.5*(atPoint(m,point)+atPoint(p,point));
    else
        id=max(m,p);if m>0,side=2;else,side=1;end;names={'xmin','xmax';'ymin','ymax';'zmin','zmax'};
        bc=boundary.(names{f.axis(k),side});G=atPoint(id,point);
        if strcmp(bc.type,'zeroNormalGradient'),G=G-dot(G,normal)*normal;end
    end
    q=-conductivity*dot(G,normal);
    assert(isfinite(q),'CutQuadratic:NonfiniteFlux','共享面Fourier通量溢出或非有限，未裁剪。');faceFlux(k)=q;
    if m>0,rateIntegral(m)=rateIntegral(m)-f.area(k)*q;else,boundaryFlux=boundaryFlux-f.area(k)*q;end
    if p>0,rateIntegral(p)=rateIntegral(p)+f.area(k)*q;else,boundaryFlux=boundaryFlux+f.area(k)*q;end
end
for k=find(w.ownerCell>0)'
    id=w.ownerCell(k);
    if strcmp(boundary.solidWall.type,'zeroNormalGradient'),q=0;wallKind{k}='prescribed-zero-normal-gradient';
    else,q=-conductivity*dot(atPoint(id,w.centroid(k,:)),w.normal(k,:));wallKind{k}='derived-from-declared-FV-temperature';end
    assert(isfinite(q),'CutQuadratic:NonfiniteFlux','墙片Fourier通量溢出或非有限，未裁剪。');
    wallHeat(k)=q;rateIntegral(id)=rateIntegral(id)-w.area(k)*q;boundaryFlux=boundaryFlux+w.area(k)*q;
end
fluid=cut.cells.volume>0;rate=nan(n,1);rate(fluid)=rateIntegral(fluid)./cut.cells.volume(fluid);
balance=sum(rateIntegral(fluid))+boundaryFlux;
assert(all(isfinite(rateIntegral(fluid)))&&all(isfinite(rate(fluid)))&&isfinite(boundaryFlux)&&isfinite(balance), ...
    'CutQuadratic:NonfiniteFlux','热功率、体积散度或全域热账溢出，未返回伪有限结果。');
out=struct('reconstruction',r,'faceHeatFluxPlusNormal',faceFlux,'wallHeatFluxIntoWall',wallHeat, ...
    'wallHeatValueKind',{wallKind},'thermalDiffusionRate',rate,'netBoundaryHeatPower',boundaryFlux, ...
    'thermalBalance',balance,'wallTemperatureInput',r.wallTemperatureInput, ...
    'wallTemperatureValueKind',{r.wallTemperatureValueKind}, ...
    'method','真实温度FV均值/Cov二次SVD重构 + 共享Fourier面通量', ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'timeIntegratorImplemented',false, ...
    'fullNavierStokesStabilityProven',false,'newSolvedSurfaceAreaM2',0);
    function value=atPoint(id,point)
        value=r.gradient(id,:)+(reshape(r.hessian(id,:,:),3,3)*(point-cut.cells.centroid(id,:))')';
    end
end
