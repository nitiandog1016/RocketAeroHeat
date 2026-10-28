function out=buildCutCellRiemannReservoirBoundaryState(Q,cut,gas,boundary)
%BUILDCUTCELLRIEMANNRESERVOIRBOUNDARYSTATE 按六个显式边界组批量求开放储库迹。
% 模块一：只处理type='riemannReservoir'，旧特征/墙/超声速边界保持原路径。
% 模块二：owner、外法向和axis/side沿用真实开放面，通量保存为向外方向。
% 模块三：状态和通量用于同一NS算子；不把储库输入或未收敛场计入覆盖。
f=cut.openFaces;names={'xmin','xmax';'ymin','ymax';'zmin','zmax'};
out=struct('faceID',zeros(0,1),'slot',zeros(0,1),'primitive',zeros(0,5),'outwardConvectiveFlux',zeros(0,5), ...
    'incoming',false(0,1),'sampleRegion',zeros(0,1),'reservoirPrimitive',zeros(0,5),'pressureRootResidualMs',zeros(0,1));
selected=false(3,2);
for axis=1:3,for side=1:2
    name=names{axis,side};selected(axis,side)=isfield(boundary,name)&&isfield(boundary.(name),'type')&&strcmp(boundary.(name).type,'riemannReservoir');
end,end
if ~any(selected,'all'),return;end
outside=find(f.area>0&xor(f.cellMinus>0,f.cellPlus>0));side=2*ones(numel(outside),1);side(f.cellMinus(outside)==0)=1;
pick=selected(sub2ind([3 2],double(f.axis(outside)),side));ids=outside(pick);side=side(pick);N=numel(ids);
out.faceID=ids;out.slot=zeros(numel(f.area),1);out.slot(ids)=(1:N)';out.primitive=zeros(N,5);out.outwardConvectiveFlux=zeros(N,5);
out.incoming=false(N,1);out.sampleRegion=zeros(N,1);out.reservoirPrimitive=zeros(N,5);out.pressureRootResidualMs=zeros(N,1);
owner=max(f.cellMinus(ids),f.cellPlus(ids));normal=f.normal(ids,:);normal(side==1,:)=-normal(side==1,:);
rho=Q(owner,1);u=Q(owner,2:4)./rho;T=(Q(owner,5)-.5*rho.*sum(u.^2,2))./(rho*gas.R/(gas.gamma-1));P=[rho u T];
for axis=1:3,for boundarySide=1:2
    if ~selected(axis,boundarySide),continue;end
    rows=find(f.axis(ids)==axis&side==boundarySide);if isempty(rows),continue;end
    b=boundary.(names{axis,boundarySide});[value,flux,tag]=resolveNSRiemannReservoir3D(P(rows,:),b,normal(rows,:),gas);
    out.primitive(rows,:)=value;out.outwardConvectiveFlux(rows,:)=flux;out.incoming(rows)=tag.isInflow;out.sampleRegion(rows)=tag.sampleRegion;
    out.reservoirPrimitive(rows,:)=tag.reservoirPrimitiveInput;out.pressureRootResidualMs(rows)=tag.pressureRootResidualMs;
end,end
end
