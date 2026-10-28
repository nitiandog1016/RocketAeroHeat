function sa=axisymmetricSADiscretization(mesh,P,boundary,right,massFlux,gx,gr,gas,N,distance,temperature,axisSymmetry)
%AXISYMMETRICSADISCRETIZATION SA第五量的质量一致共享通量、扩散及完整可压缩源。
% 模块一：真实质量方向决定开放N迹；墙N=0，出流/零质量为owner且零法向扩散。
% 模块二：NS有效mu/k与SA扩散系数分开；正cell/边界一阶系数均值不外推N/rho。
% 模块三：cell体积质心求SA源，全部乘r加权V；密度梯度修正不可省略或重复加轴源。
% 显式第十二参量仅由NS的mirror-parity轴模板传入，缺省保持旧全部表达式和info。
% 本模块只构造实际残差，不宣称实验、湍流稳态或真实箭体区域已验收。
f=mesh.faces;o=f.owner;j=f.neighbour;inside=j>0;bf=mesh.boundaryFaces;nc=mesh.cellCount;nf=numel(o);
assert(isequal(size(N),[nc 1])&&all(isfinite(N)&N>=0),'AxisSA:State','nuTilde须有限非负，未裁值。');
assert(isequal(size(distance),[nc 1])&&all(isfinite(distance)&distance>0),'AxisSA:Distance','SA源须真实正壁距。');
rho=P(:,1);cellT=P(:,4);model=axisymmetricSAModel(rho,cellT,N,abs(gx(:,2)-gr(:,1)),distance,gas);
sigma=model.constants.sigma;Prt=model.constants.turbulentPrandtl;
Nright=N(o);Nright(inside)=N(j(inside));normalZero=false(nf,1);wall=false(nf,1);axisFace=false(nf,1);
names={'iMin','iMax','jMin','jMax'};
for b=1:numel(bf)
    k=bf(b);bc=boundary.(names{f.boundarySide(k)});
    if strcmp(bc.type,'wall'),Nright(k)=0;wall(k)=true;
    elseif strcmp(bc.type,'axis'),normalZero(k)=true;axisFace(k)=true;
    elseif massFlux(k)<0
        assert(isfield(bc,'nuTilde'),'AxisSA:Boundary','实际入流/回流必须明确给出nuTilde边界输入。');
        value=bc.nuTilde;if isa(value,'function_handle'),value=value(f.centroid(k,:));end
        assert(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>=0, ...
            'AxisSA:Boundary','给定入流nuTilde须有限非负，不能暗补0。');Nright(k)=value;
    else,normalZero(k)=true;
    end
end
mirror=nargin>=12&&axisSymmetry.active;
sample=[rho N;right(bf,1) Nright(bf)];
if mirror
    ai=axisSymmetry.boundaryIndex;ao=axisSymmetry.ownerCell;af=axisSymmetry.faceID;
    % 密度和N都是偶量；虚点用owner同值，不能把真实偏斜面迹冒充镜像样本。
    sample(nc+ai,:)=[rho(ao),N(ao)];
    sx=axisSymmetry.gradientX*sample;sr=axisSymmetry.gradientR*sample;
    % 切向外推仅诊断；零面积轴面系数固定使用owner偶量代理，不作限幅/回退。
    dx=axisSymmetry.physicalDisplacementXR(:,1);
    rawAxisRhoN=[rho(ao)+sx(ao,1).*dx,N(ao)+sx(ao,2).*dx];axisN=N(ao);
    assert(all(isfinite(axisN)&axisN>=0),'AxisSA:AxisTraceState','零面积轴面owner投影N非有限或为负，未裁值。');
    Nright(af)=axisN;
else
    sx=mesh.gradientX*sample;sr=mesh.gradientR*sample;
end
Gx=sx(o,2);Gr=sr(o,2);
Gx(inside)=.5*(sx(o(inside),2)+sx(j(inside),2));Gr(inside)=.5*(sr(o(inside),2)+sr(j(inside),2));
delta=f.delta;correction=(Nright-N(o)-Gx.*delta(:,1)-Gr.*delta(:,2))./sum(delta.^2,2);
Gx=Gx+correction.*delta(:,1);Gr=Gr+correction.*delta(:,2);
normal=Gx.*f.normal(:,1)+Gr.*f.normal(:,2);
Gx(normalZero)=Gx(normalZero)-normal(normalZero).*f.normal(normalZero,1);
Gr(normalZero)=Gr(normalZero)-normal(normalZero).*f.normal(normalZero,2);
if mirror
    Gx(af)=sx(ao,2);Gr(af)=0;
end
faceRho=rho(o);faceN=N(o);
faceRho(inside)=.5*(rho(o(inside))+rho(j(inside)));faceN(inside)=.5*(N(o(inside))+N(j(inside)));
faceRho(bf)=right(bf,1);faceN(bf)=Nright(bf);
lam=axisymmetricFVTransport(temperature.face,gas);sourceLam=axisymmetricFVTransport(temperature.source,gas);
cellLam=axisymmetricFVTransport(cellT,gas);
faceMuT=eddy(faceRho,faceN,lam.mu);sourceMuT=eddy(rho,N,sourceLam.mu);
effectiveFaceMu=lam.mu+faceMuT;effectiveFaceK=lam.k+gas.cp*faceMuT/Prt;
diffusion=(lam.mu+faceRho.*faceN)/sigma;
upwind=N(o);negative=massFlux<0;upwind(negative)=Nright(negative);
advective=massFlux.*upwind;
diffusive=diffusion.*(Gx.*f.normal(:,1)+Gr.*f.normal(:,2));
advective(axisFace)=0;diffusive(axisFace)=0;
production=rho.*model.production;destruction=rho.*model.destruction;
gradientSquare=rho*(model.constants.cb2/sigma).*(sx(:,2).^2+sr(:,2).^2);
densityCorrection=-(model.nu+N)/sigma.*(sx(:,1).*sx(:,2)+sr(:,1).*sr(:,2));
source=production-destruction+gradientSquare+densityCorrection;
assert(all(isfinite([source;advective;diffusive;effectiveFaceMu;effectiveFaceK])), ...
    'AxisSA:Arithmetic','SA通量或源非有限，未裁值。');
sa=struct('modelID',model.modelID,'nuTilde',N,'wallDistance',distance,'cellModel',model, ...
    'boundaryTraceNuTilde',Nright(bf),'faceNuTilde',faceN,'faceDensity',faceRho, ...
    'cellDensityGradientX',sx(:,1),'cellDensityGradientR',sr(:,1), ...
    'cellNuTildeGradientX',sx(:,2),'cellNuTildeGradientR',sr(:,2), ...
    'faceNuTildeGradientX',Gx,'faceNuTildeGradientR',Gr,'faceNormalZeroGradient',normalZero, ...
    'faceMassFluxPerRadian',massFlux.*f.areaPerRadian, ...
    'faceAdvectiveFluxPerRadian',advective.*f.areaPerRadian, ...
    'faceDiffusiveFluxPerRadian',diffusive.*f.areaPerRadian, ...
    'faceFluxIntegralPerRadian',(advective-diffusive).*f.areaPerRadian, ...
    'diffusionCoefficient',diffusion,'faceMuT',faceMuT,'cellMuT',model.muT,'sourceMuT',sourceMuT, ...
    'effectiveFaceMu',effectiveFaceMu,'effectiveFaceK',effectiveFaceK, ...
    'effectiveSourceMu',sourceLam.mu+sourceMuT,'effectiveCellMu',cellLam.mu+model.muT, ...
    'sourceProduction',production,'sourceDestruction',destruction,'sourceGradientSquare',gradientSquare, ...
    'sourceDensityCorrection',densityCorrection,'sourceRatePerVolume',source, ...
    'cellSourceIntegral',source.*mesh.volumePerRadian,'coefficientSampling','positive cell/actual-boundary arithmetic mean; no reconstructed N or rho', ...
    'sourceQuadrature','one-point volume-centroid, r-weighted volume; no duplicate scalar axis source', ...
    'wallMuTIsZero',all(faceMuT(wall)==0),'wallHeatEnhancementFactorApplied',false, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false);
if mirror
    sa.axisSymmetry=struct('kind',axisSymmetry.kind,'policy',axisSymmetry.policy, ...
        'faceID',af,'boundaryIndex',ai,'ownerCell',ao,'physicalPointXR',axisSymmetry.physicalPointXR, ...
        'ghostPointXR',axisSymmetry.ghostPointXR,'ghostDensityNuTilde',[rho(ao),N(ao)], ...
        'coefficientDensityNuTilde',[right(af,1),Nright(af)], ...
        'coefficientStateKind','owner-even-projection-at-zero-measure-axis-not-a-solved-face-trace', ...
        'coefficientStateIsPhysicalFaceSolution',false,'physicalTraceReconstructed',false, ...
        'rawTangentialDensityNuTilde',rawAxisRhoN,'rawTangentialStateUsed',false, ...
        'rawTangentialStateAdmissible',all(isfinite(rawAxisRhoN),2)&rawAxisRhoN(:,1)>0&rawAxisRhoN(:,2)>=0, ...
        'cellDensityGradientXAtAxisOwner',sx(ao,1),'cellDensityGradientRAtAxisOwner',sr(ao,1), ...
        'cellNuTildeGradientXAtAxisOwner',sx(ao,2),'cellNuTildeGradientRAtAxisOwner',sr(ao,2), ...
        'faceNuTildeGradientX',Gx(af),'faceNuTildeGradientR',Gr(af),'cellStencilUnchanged',true, ...
        'stateClipped',false,'physicalFaceGeometryChanged',false);
    sa.coefficientSampling=[sa.coefficientSampling '; opt-in zero-measure axis uses EOS-consistent owner rho and owner N proxies; raw tangential extrapolation diagnostic only; gradient samples use even mirrors'];
end
end
function muT=eddy(rho,N,mu)
% 模块用途：明确SA代数黏性；扩散中使用rho*N而非rho*N*fv1。
chi=rho.*N./mu;fv1=zeros(size(chi));small=chi<=7.1;
ratio=(chi(small)/7.1).^3;fv1(small)=ratio./(1+ratio);
ratio=(7.1./chi(~small)).^3;fv1(~small)=1./(1+ratio);muT=rho.*N.*fv1;
assert(all(isfinite(muT)&muT>=0),'AxisSA:Arithmetic','SA涡黏性非有限或负值。');
end
