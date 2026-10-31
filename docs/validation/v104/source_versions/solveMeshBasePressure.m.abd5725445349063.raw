function result = solveMeshBasePressure(mesh,freestream,cfg)
%SOLVEMESHBASEPRESSURE 将有依据的底压关联映射到真实平圆底面。
% 模块用途：直接审查当前 V/F，而非信任几何名称或用户声明；对
% 孤立尖头直圆锥，内部求肩部锥流、转捩类别及层流物理厚度后，
% 将 SAND93-0280 的近似均匀底压赋给真实底面三角形。
% 本模块不修改 mesh.pressure，不把代表底热扩展成底面热分布。
% 坐标：主轴 X，鼻尖 +X，零攻角流速 -X；所有长度为 m。
%
% cfg.basePressure.hasRearSupport / hasGasInjection 须显式给出布尔值，
% operationSource 须说明无尾撑/无喷流声明的来源（STL不能证明工况）。
% 现版只接受显式 caloricallyPerfect、gamma=1.4 空气、零攻角及
% 9~10度尖头直圆锥层流；转捩使用项目 Reynolds 数经验阈值。
% 即使几何/类别合格，仍须通过同剖面物理 delta/rb<=0.05 筛选。
% 湍流及其他前体暂缺可信肩部薄层证据，保持未解。
if nargin<3, cfg=struct; end
if ~isstruct(mesh) || ~isfield(mesh,'V') || ~isfield(mesh,'F')
    error('RocketAeroHeat:MeshBase:Input','需要真实 mesh.V 和 mesh.F。');
end
validateattributes(mesh.V,{'numeric'},{'real','finite','2d','ncols',3,'nonempty'});
validateattributes(mesh.F,{'numeric'},{'real','finite','integer','positive', ...
    '2d','ncols',3,'nonempty','<=',size(mesh.V,1)});
V=double(mesh.V); F=double(mesh.F); nf=size(F,1);
result=blankResult(nf);

%% 1. 按真实顶点、连通性、法向和母线识别圆底尖锥
% 独立重算法向/面积，以免沿用已经过期的 localGeometry 或面积缓存。
[geometry,baseFace,area]=verifyGeometry(V,F);
result.geometry=geometry; result.baseFace=baseFace;
result.faceArea=area; result.summary.totalArea=sum(area);
result.summary.baseArea=sum(area(baseFace));
result.summary.baseFaces=nnz(baseFace);
result.faceStatus(baseFace)=uint8(2);
if ~geometry.verified
    result=reject(result,2,['真实几何不受支持：' geometry.status]); return;
end

%% 2. 审核外部操作声明和同一自由流物理范围
% 几何可以证明平底和没有附着尾撑面，不能证明发动机是否正在喷流。
if ~isfield(cfg,'basePressure') || ~isstruct(cfg.basePressure)
    result=reject(result,3,'缺少无尾撑/无喷流的显式操作声明及来源'); return;
end
bp=cfg.basePressure;
for key={'hasRearSupport','hasGasInjection'}
    if ~isfield(bp,key{1}) || ~isscalar(bp.(key{1})) || ...
            ~(islogical(bp.(key{1})) || isnumeric(bp.(key{1}))) || ...
            ~ismember(bp.(key{1}),[0 1])
        result=reject(result,3,['缺少有效 cfg.basePressure.' key{1}]); return;
    end
end
if ~isfield(bp,'operationSource') || ...
        ~(ischar(bp.operationSource) || (isstring(bp.operationSource) && isscalar(bp.operationSource))) || ...
        isempty(strtrim(char(bp.operationSource)))
    result=reject(result,3,'操作声明 operationSource 不能为空'); return;
end
result.operation=struct('hasRearSupport',logical(bp.hasRearSupport), ...
    'hasGasInjection',logical(bp.hasGasInjection), ...
    'valueKind','input','source',char(bp.operationSource));
if bp.hasRearSupport || bp.hasGasInjection
    result=reject(result,3,'当前有尾撑或底部喷流，不能采用无干扰底压关联'); return;
end
required={'pressure','temperature','density','speed','Mach','gamma','R', ...
    'alpha','beta','totalPressure','totalTemperature'};
for j=1:numel(required)
    key=required{j};
    if ~isfield(freestream,key) || ~isnumeric(freestream.(key)) || ...
            ~isscalar(freestream.(key)) || ~isfinite(freestream.(key))
        result=reject(result,4,['自由流缺少有效 ' key]); return;
    end
end
fs=freestream;
if abs(fs.alpha)>1e-8 || abs(fs.beta)>1e-8
    result=reject(result,4,'底压关联限零攻角和零侧滑角'); return;
end
[airMode,~]=selectAirThermalModel(fs,cfg);
if ~strcmp(airMode,'caloricallyPerfect') || abs(fs.gamma-1.4)>1e-8 || ...
        abs(fs.R/287.05287-1)>.005 || fs.Mach<2 || fs.Mach>20 || ...
        min([fs.pressure fs.temperature fs.density fs.speed])<=0
    result=reject(result,4,'现版底压映射仅适用显式定比热 gamma=1.4 空气及 2<=Mach<=20'); return;
end
cp=fs.gamma*fs.R/(fs.gamma-1);
stateResidual=max([abs(fs.density*fs.R*fs.temperature/fs.pressure-1), ...
    abs(fs.speed/(fs.Mach*sqrt(fs.gamma*fs.R*fs.temperature))-1), ...
    abs(fs.totalTemperature/(fs.temperature+fs.speed^2/(2*cp))-1), ...
    abs(fs.totalPressure/(fs.pressure*(fs.totalTemperature/fs.temperature)^(fs.gamma/(fs.gamma-1)))-1)]);
result.freestreamRelativeResidual=stateResidual;
if stateResidual>1e-5 || fs.temperature<80 || fs.totalTemperature>2500
    result=reject(result,4,'自由流气体/总焓不一致或超出80~2500K无化学空气范围'); return;
end
if geometry.halfAngleDeg<9-1e-8 || geometry.halfAngleDeg>10+1e-8
    result=reject(result,4,'现版真实底面接入限9~10度层流尖头直圆锥'); return;
end

%% 3. 程序内部求肩部状态与转捩类别；不借用用户指定底压
state=solveConicalShock(fs,geometry.halfAngleDeg);
result.internalShoulderState=state;
if ~state.attached
    result=reject(result,5,['内部肩部外流失效：' state.status]); return;
end
shoulder=struct;
for key={'pressure','temperature','density','speed','Mach'}
    shoulder.(key{1})=state.(key{1});
    shoulder.sources.(key{1})=['内部求解：' state.method '；真实尖头直锥肩部'];
end
shoulder.distance=geometry.generatorLength;
shoulder.dynamicViscosity=1.458e-6*state.temperature^1.5/(state.temperature+110.4);
shoulder.sources.distance='真实圆底/尖端的母线长度 hypot(轴向长度,底半径)';
shoulder.sources.dynamicViscosity='Sutherland 空气黏度，与生产输运实现一致';
result.shoulder=shoulder;
Re=shoulder.density*shoulder.speed*shoulder.distance/shoulder.dynamicViscosity;
transition=calcTransition(Re,true,cfg);
result.shoulderTransition=transition;
result.fieldMetadata.boundaryLayerRegime=fieldRecord('computed',transition.method, ...
    true,transition.scope);
result.fieldMetadata.boundaryLayerRegime.criticalReynolds=transition.criticalReynolds;
if transition.turbulentFace
    result=reject(result,6,['肩部按当前经验阈值已转捩；现版缺少受核验的' ...
        '湍流物理delta，不用动量厚度任意换算薄层']); return;
end
if Re<1e5
    result=reject(result,6,'肩部层流Re_e,s<1e5，低于报告底压关联范围'); return;
end

%% 4. 同一尖锥外流的轴对称层流物理剖面与薄层审核
% 直锥 Ue、pe 恒定，h=s*sin(theta)，xi=rho*mu*U*sin²(theta)*s³/3。
% 因而beta=0且满足沿程动量方程。该剖面只审核肩部薄层，不输出底热。
Tw=NaN;
if isfield(cfg,'flight') && isfield(cfg.flight,'Twall'), Tw=cfg.flight.Twall; end
if isfield(cfg,'heating') && isfield(cfg.heating,'wallTemperatureByFace')
    wall=cfg.heating.wallTemperatureByFace;
    ids=geometry.shoulderVertexIDs;
    shoulderFace=~baseFace & any(ismember(F,ids),2);
    wall=double(wall(:));
    if numel(wall)~=nf || any(~isfinite(wall(shoulderFace))) || ...
            max(wall(shoulderFace))-min(wall(shoulderFace))>1e-8
        result=reject(result,6,'非均匀/缺失肩部壁温不能用单一轴对称剖面审核'); return;
    end
    Tw=wall(find(shoulderFace,1));
    % 前体沿程温度不同会改变肩部剖面；不以局部值代替完整热边界。
    side=wall(~baseFace);
    if any(~isfinite(side)) || max(side)-min(side)>1e-8
        result=reject(result,6,'当前直锥剖面要求前体沿程等温壁边界'); return;
    end
end
if ~isnumeric(Tw) || ~isscalar(Tw) || ~isfinite(Tw) || Tw<80 || Tw>2500
    result=reject(result,6,'薄层审核需要80~2500K的给定等温壁边界'); return;
end
edge=shoulder;
edge.totalEnthalpy=cp*shoulder.temperature+shoulder.speed^2/2;
edge.sources.totalEnthalpy='内部Taylor-Maccoll肩部同物性静焓+动能';
s=shoulder.distance; rb=geometry.baseRadius;
station=struct('streamlineMetric',rb, ...
    'xi',shoulder.density*shoulder.dynamicViscosity*shoulder.speed*rb^2*s/3, ...
    'beta',0,'pressureGradient',0,'velocityGradient',0);
station.sources=struct('streamlineMetric','真实直锥肩部半径；h=s*sin(theta_c)', ...
    'xi','恒定锥面外流下解析积分int(rho_e mu_e U_e h² ds)', ...
    'beta','直锥外缘速度梯度为零的Levy-Lees局部相似参数', ...
    'pressureGradient','Taylor-Maccoll直锥面静压沿母线恒定', ...
    'velocityGradient','Taylor-Maccoll直锥面速度沿母线恒定', ...
    'wallTemperature','给定等温前体壁温，仅用于肩部物理厚度审核');
profile=solveAxisymmetricLaminarMassProfile(fs,edge,Tw,station,cfg);
result.shoulderBoundaryLayer=profile;
result.fieldMetadata.wallTemperatureInput=fieldRecord('prescribed', ...
    station.sources.wallTemperature,true,'前体等温壁边界，不是底面壁温求解值');
result.shoulderThicknessRatio=profile.delta/rb;
result.thinLayerMaximumRatio=.05;
result.thinLayerCriterionSource=['delta/rb<=0.05 为本实现显式保守工程门限；' ...
    '不是SAND93-0280报告给出的普适数值阈值'];
if ~profile.available || ~isfinite(profile.delta) || profile.delta<=0
    result=reject(result,6,['肩部物理剖面未通过审核：' profile.status]); return;
end
if result.shoulderThicknessRatio>.05
    result=reject(result,6,'肩部物理delta/rb>0.05，不采用薄边界层底压关联'); return;
end

%% 5. 调用已校核关联，并仅映射压力和标准Cp到真实底面
g=struct('family','cone','boundaryLayerRegime','laminar','baseRadius',rb, ...
    'halfAngleDeg',min(10,max(9,geometry.halfAngleDeg)), ...
    'noseBluntnessRatio',0,'flatCircularBase',true, ...
    'isolatedAxisymmetricBody',true,'hasRearSupport',false,'hasGasInjection',false, ...
    'source','当前V/F的圆底、直母线、同轴圆截面、尖端和闭合拓扑独立审核');
correlationCfg=struct;
if isfield(bp,'maxContinuumKnudsen'), correlationCfg.maxContinuumKnudsen=bp.maxContinuumKnudsen; end
correlation=solveAxisymmetricBaseCorrelation(fs,shoulder,g,Tw,correlationCfg);
result.correlation=correlation;
if ~correlation.pressureAvailable
    result=reject(result,7,['底压关联不适用：' correlation.status]); return;
end
result.pressure(baseFace)=correlation.basePressure;
result.Cp(baseFace)=correlation.basePressureCoefficient;
result.computedFace=baseFace; result.pressureAvailable=true; result.available=true;
result.faceStatus(baseFace)=uint8(0);
result.pressureMethod=correlation.pressureMethod;
result.CpMethod=[correlation.pressureMethod '；Cp=(pb-pinf)/qinf'];
result.pressureScope=['真实孤立平圆底尖锥；零攻角/侧滑、9~10度、gamma=1.4定比热空气；' ...
    '程序内部附体肩部、层流经验判据及物理薄层审核；' ...
    '报告p11支持近似均匀静压，不解析底唇细节或回流三维速度；无喷流/尾撑'];
result.faceMethod(baseFace)={result.pressureMethod};
result.faceScope(baseFace)={result.pressureScope};
result.fieldMetadata.pressure=fieldRecord('computed',result.pressureMethod,true,result.pressureScope);
result.fieldMetadata.Cp=fieldRecord('computed',result.CpMethod,true,result.pressureScope);
result.pressureEligibleForSurfaceCoverage=true;
result.summary.computedPressureArea=sum(area(baseFace));
result.summary.computedPressureFaces=nnz(baseFace);
result.summary.pressureAreaFraction=result.summary.computedPressureArea/result.summary.totalArea;
result.summary.basePressureAreaFraction=1;
result.status='真实底面压力已计算并映射；底面热流及相关温度仍未求解';
% 代表底热仍为地区标量，没有真实径向/面积平均证据，不纳入热覆盖。
result.regionRepresentative=struct('heatFlux',correlation.heatFluxRepresentative, ...
    'StantonNumber',correlation.baseStantonNumber,'method',correlation.heatingMethod, ...
    'scope',correlation.heatingScope,'eligibleForSurfaceCoverage',false);
end

function [g,base,area]=verifyGeometry(V,F)
% 模块用途：确认实际三角网格确为单一闭合、圆底、共轴尖头直圆锥。
nf=size(F,1); nv=size(V,1);
g=struct('verified',false,'flatBaseVerified',false,'status','尚未审核','baseRadius',NaN, ...
    'halfAngleDeg',NaN,'generatorLength',NaN,'axisCenterYZ',[NaN NaN], ...
    'shoulderVertexIDs',[],'circleRelativeResidual',NaN,'coneRelativeResidual',NaN, ...
    'source','从当前V/F重新计算；不信任cfg几何标签');
span=max(V)-min(V); scale=max(span); tol=max(scale*1e-8,64*eps(max(abs(V(:)))));
x=V(:,1); xmin=min(x); xmax=max(x); L=xmax-xmin;
base=all(abs(x(F)-xmin)<=tol,2);
raw=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
area=vecnorm(raw,2,2)/2;
if L<=tol || any(area<=64*eps*scale^2)
    g.status='轴长为零或存在退化三角形'; return;
end
N=raw./(2*area);
if ~any(base) || any(N(base,1)>-1+1e-8)
    g.status='没有具有正确尾向外法向的完整平底三角形'; return;
end
if numel(unique(F(:)))~=nv
    g.status='存在游离顶点，不能证明单一孤立物体'; return;
end
% 每条边必须使用恰好两次且方向相反，拒绝开放喷管、重叠面和非流形。
allEdges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
faceOfEdge=repmat((1:nf)',3,1);
[edges,~,edgeIndex]=unique(sort(allEdges,2),'rows');
count=accumarray(edgeIndex,1);
direction=sign(allEdges(:,2)-allEdges(:,1));
if any(count~=2) || any(accumarray(edgeIndex,direction)~=0)
    g.status='整体不是具有一致外向绕序的闭合二流形'; return;
end
owners=accumarray(edgeIndex,faceOfEdge,[],@(a){a});
adj=zeros(size(edges,1),2);
for j=1:numel(owners),adj(j,:)=owners{j};end
graphFaces=graph(adj(:,1),adj(:,2),[],nf);
if max(conncomp(graphFaces))~=1 || nv-size(edges,1)+nf~=2
    g.status='整体不是单连通闭合物体，可能含支撑、洞或额外组件'; return;
end
% 真正底面须为一个拓扑圆盘，边界为一个没有断点的简单环。
baseCount=accumarray(edgeIndex,double(base(faceOfEdge)));
boundary=edges(baseCount==1,:); rim=unique(boundary(:));
g.shoulderVertexIDs=rim;
rimDegree=accumarray(boundary(:),1,[nv 1]);
if numel(rim)<8 || any(rimDegree(rim)~=2)
    g.status='底面边界不是至少8点的单一闭合周向环'; return;
end
baseVertices=unique(F(base,:)); baseEdges=edges(baseCount>0,:);
if numel(baseVertices)-size(baseEdges,1)+nnz(base)~=1
    g.status='底面不是无孔圆盘'; return;
end
[~,localEdges]=ismember(boundary,rim);
if max(conncomp(graph(localEdges(:,1),localEdges(:,2),[],numel(rim))))~=1
    g.status='底面含多个独立边界圈'; return;
end
yz=V(rim,2:3); origin=mean(yz,1); yz0=yz-origin;
A=[2*yz0 ones(numel(rim),1)];
if rank(A)<3, g.status='底边界不能拟合非退化圆'; return; end
c=A\sum(yz0.^2,2); center=origin+c(1:2)';
rb=sqrt(c(3)+sum(c(1:2).^2));
radius=vecnorm(yz-center,2,2);
g.circleRelativeResidual=max(abs(radius/rb-1));
if ~isfinite(rb) || rb<=tol || g.circleRelativeResidual>1e-6
    g.status='真实底面边界不满足圆截面'; return;
end
theta=atan2(yz(:,2)-center(2),yz(:,1)-center(1));
[theta,order]=sort(theta); ring=rim(order); gap=diff([theta;theta(1)+2*pi]);
expectedEdges=sort([ring ring([2:end 1])],2);
if any(gap<=1e-8) || max(gap)>pi/4+1e-8 || ...
        ~isequal(sortrows(expectedEdges),sortrows(boundary))
    g.status='底边界不是按周向顺序连接的充分采样简单圆环'; return;
end
polygonArea=abs(sum(yz(order,1).*yz(order([2:end 1]),2)- ...
    yz(order,2).*yz(order([2:end 1]),1)))/2;
if abs(sum(area(base))/polygonArea-1)>1e-6
    g.status='底三角形未无重叠地覆盖真实边界圆盘'; return;
end
% 圆底核验与前体模型核验分开：圆底完整但前体不受支持时，可将
% 已确定的整片底面标为未解；圆底不完整时调用方必须拒绝应用，
% 不能只把剩余X极值平片置NaN而让倾斜残片暗中保留给定Cp。
g.flatBaseVerified=true;
tip=find(abs(x-xmax)<=tol);
if numel(tip)~=1 || norm(V(tip,2:3)-center)>rb*1e-6
    g.status='鼻端不是同轴孤立尖端，钝锥不能套用尖锥肩部解'; return;
end
% 底面内部顶点只检验圆盘包含性；全部前体顶点须位于同一直母线锥面。
sideVertices=setdiff(unique(F(~base,:)),baseVertices);
sideVertices=unique([sideVertices;rim;tip]);
rho=vecnorm(V(:,2:3)-center,2,2);
expected=rb*(xmax-x)/L;
g.coneRelativeResidual=max(abs(rho(sideVertices)-expected(sideVertices)))/rb;
if g.coneRelativeResidual>1e-6 || any(rho(baseVertices)>rb*(1+1e-6))
    g.status='前体不满足同轴尖头直母线；球鼻、锥柱、尾锥或支撑需独立模型'; return;
end
C=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
radial=[zeros(nf,1) C(:,2:3)-center];
if any(N(~base,1)<=0) || any(sum(N(~base,:).*radial(~base,:),2)<=0)
    g.status='侧壁法向不是向外/朝鼻方向，不能确认真实外表面'; return;
end
% 限制前体每片的周向跨度；拉链网格的中间圈点数变化，不应要求
% 分片体积恰等于底多边形直棱锥，须保留圆锥的弦离散误差。
maxFaceSpan=0;
for id=find(~base)'
    vertices=F(id,:); vertices=vertices(rho(vertices)>tol);
    a=sort(atan2(V(vertices,3)-center(2),V(vertices,2)-center(1)));
    d=diff([a;a(1)+2*pi]);
    faceSpan=2*pi-max(d); maxFaceSpan=max(maxFaceSpan,faceSpan);
end
if maxFaceSpan>pi/4+1e-8
    g.status='前体三角形跨越超过45度周向弧，无法核验圆锥面离散'; return;
end
% 与解析圆锥体积的弦误差界比较；平移坐标避免大坐标取消误差。
P=V-[xmin center];
volume=sum(dot(P(F(:,1),:),cross(P(F(:,2),:),P(F(:,3),:),2),2))/6;
idealVolume=pi*rb^2*L/3;
if volume/idealVolume>1+1e-6 || volume/idealVolume<cos(maxFaceSpan/2)^2-1e-6
    g.status='闭合体积不满足真实圆底尖锥，可能有内部或重叠面'; return;
end
g.verified=true;g.status='真实闭合圆底、共轴尖端及直母线均通过审核';
g.baseRadius=rb;g.halfAngleDeg=atan2d(rb,L);g.generatorLength=hypot(rb,L);
g.axisCenterYZ=center;g.baseX=xmin;g.tipX=xmax;g.basePolygonArea=polygonArea;
g.closedSignedVolume=volume;g.circleIdealArea=pi*rb^2;
g.maximumFaceCircumferentialSpanDeg=rad2deg(maxFaceSpan);
g.facetedVolumeFraction=volume/idealVolume;
end

function r=blankResult(nf)
% 模块用途：所有未求解量从开始就保持NaN，零值不充当无解标记。
r=struct('available',false,'pressureAvailable',false,'status','未计算', ...
    'baseFace',false(nf,1),'computedFace',false(nf,1), ...
    'pressure',nan(nf,1),'Cp',nan(nf,1),'heatFlux',nan(nf,1), ...
    'edgeTemperature',nan(nf,1),'recoveryTemperature',nan(nf,1), ...
    'solvedWallTemperature',nan(nf,1),'faceStatus',ones(nf,1,'uint8'), ...
    'faceMethod',{repmat({'未求解'},nf,1)},'faceScope',{repmat({'不属于本次已解区域'},nf,1)}, ...
    'pressureMethod','','CpMethod','','pressureScope','', ...
    'pressureEligibleForSurfaceCoverage',false,'heatEligibleForSurfaceCoverage',false, ...
    'radialHeatingDistributionSolved',false,'shoulderThicknessRatio',NaN);
r.faceStatusNames={'已计算底压','非底面区域','真实几何不受支持', ...
    '操作声明缺失或存在干扰','自由流/气体/角度范围外', ...
    '内部肩部外流失效','肩部边界层类别或物理薄层未通过','底压关联范围外'};
r.summary=struct('totalArea',NaN,'baseArea',0,'baseFaces',0, ...
    'computedPressureArea',0,'computedPressureFaces',0,'pressureAreaFraction',0, ...
    'basePressureAreaFraction',0,'computedHeatArea',0);
r.fieldMetadata=struct;
for key={'pressure','Cp','heatFlux','edgeTemperature','recoveryTemperature','solvedWallTemperature'}
    r.fieldMetadata.(key{1})=fieldRecord('unresolved','未求解',false,'没有有效区域计算结果');
end
r.reference=struct('report','SAND93-0280','doi','10.2172/10111104', ...
    'pressureEquation','Eq.(8b)','uniformPressureEvidence','报告打印页11，低速回流区静压近似均匀');
end
function r=reject(r,code,message)
% 模块用途：保留真实底面掩码和拒绝原因，绝不放宽判据或给定底压。
r.status=message;r.faceStatus(r.baseFace)=uint8(code);
r.faceScope(r.baseFace)={message};
end
function f=fieldRecord(kind,source,valid,scope)
f=struct('valueKind',kind,'source',source,'valid',logical(valid),'scope',scope);
end
