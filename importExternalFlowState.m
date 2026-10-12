function mesh = importExternalFlowState(mesh,freestream,data)
%IMPORTEXTERNALFLOWSTATE 导入与最终面元逐一对应的三维外流边界层外缘状态。
% data 每行对应 mesh.F 的同一面元，必需字段：
%   faceCenter [N,3] m、faceNormal [N,3]、pressure [N,1] Pa、
%   temperature [N,1] K、density [N,1] kg/m^3、velocity [N,3] m/s、
%   method（CFD 软件/物理模型名称）和 caseID（流场工况标识）。
% 中心、法向、理想气体状态方程及当前飞行工况都经过校验。
% 本程序不从压力“推测”速度，也不把不匹配的流场插到当前网格上。

%% 1. 检查面元几何和数据来源
if ~isstruct(mesh) || ~all(isfield(mesh,{'V','F','localGeometry'}))
    error('RocketAeroHeat:ExternalFlow:Mesh','需要最终 mesh 和 localGeometry。');
end
if ~isstruct(data) || ~isscalar(data)
    error('RocketAeroHeat:ExternalFlow:Data','data 必须是标量结构体。');
end
nf=size(mesh.F,1); g=mesh.localGeometry;
for name={'faceCenter','faceNormal','pressure','temperature','density','velocity', ...
        'method','caseID'}
    if ~isfield(data,name{1})
        error('RocketAeroHeat:ExternalFlow:Data','缺少 data.%s。',name{1});
    end
end
for name={'faceCenter','faceNormal','velocity'}
    key=name{1};
    validateattributes(data.(key),{'numeric'}, ...
        {'real','finite','size',[nf 3]},mfilename,['data.' key]);
end
for name={'pressure','temperature','density'}
    key=name{1};
    validateattributes(data.(key),{'numeric'}, ...
        {'real','finite','positive','size',[nf 1]},mfilename,['data.' key]);
end
if ~(ischar(data.method) || (isstring(data.method) && isscalar(data.method))) || ...
        strlength(string(data.method))==0
    error('RocketAeroHeat:ExternalFlow:Source','必须说明外流计算方法。');
end
if ~(ischar(data.caseID) || (isstring(data.caseID) && isscalar(data.caseID))) || ...
        strlength(string(data.caseID))==0
    error('RocketAeroHeat:ExternalFlow:Source','必须提供流场工况标识 caseID。');
end
validateattributes(g.center,{'numeric'},{'real','finite','size',[nf 3]});
validateattributes(g.normal,{'numeric'},{'real','finite','size',[nf 3]});
span=max(max(mesh.V)-min(mesh.V));
positionTolerance=max(1e-8*span,1e-10);
if any(vecnorm(data.faceCenter-g.center,2,2)>positionTolerance) || ...
        any(vecnorm(data.faceNormal-g.normal,2,2)>1e-6)
    error('RocketAeroHeat:ExternalFlow:FaceMapping', ...
        '外流面元中心或法向与当前网格不匹配；请重新映射流场。');
end

%% 2. 校验工况和完全气体状态，防止旧流场混入新飞行条件
for name={'Mach','alpha','beta','altitude','gamma','R', ...
        'pressure','temperature','density','speed','totalPressure','totalTemperature'}
    key=name{1};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:ExternalFlow:Freestream','缺少 freestream.%s。',key);
    end
    if ~isfield(data,'freestream') || ~isfield(data.freestream,key) || ...
            ~isequal(data.freestream.(key),freestream.(key))
        error('RocketAeroHeat:ExternalFlow:CaseMismatch', ...
            '流场与当前工况不一致：%s。',key);
    end
end
p=double(data.pressure(:)); T=double(data.temperature(:));
rho=double(data.density(:)); Uvec=double(data.velocity);
eosError=abs(p-rho.*freestream.R.*T)./p;
if any(eosError>0.01)
    error('RocketAeroHeat:ExternalFlow:EquationOfState', ...
        '导入数据与当前完全气体状态方程偏差超过 1%%。');
end

%% 3. 计算可供后续边界层模块读取的外缘派生量
U=vecnorm(Uvec,2,2); a=sqrt(freestream.gamma*freestream.R*T);
M=U./a; cp=freestream.gamma*freestream.R/(freestream.gamma-1);
T0=T+U.^2/(2*cp);
p0=p.*(T0./T).^(freestream.gamma/(freestream.gamma-1));
mu=1.458e-6*T.^(3/2)./(T+110.4);
edge=struct;
edge.method=['逐面元导入三维外流：' char(string(data.method))];
edge.pressureSource=edge.method;
edge.assumption='使用给定外流数据；派生量采用定比热完全气体和 Sutherland 黏度';
edge.caseID=char(string(data.caseID));
edge.pressure=p; edge.temperature=T; edge.density=rho;
edge.velocityVector=Uvec; edge.speed=U; edge.Mach=M;
edge.speedOfSound=a; edge.dynamicViscosity=mu;
edge.totalPressure=p0; edge.totalTemperature=T0;
edge.algebraicallyValid=true(nf,1);
edge.physicallyUnresolved=false(nf,1);
edge.freestream=freestream;
edge.status='已导入与当前网格及工况匹配的外流状态；仍需核查取样是否位于边界层外缘';

%% 4. 保存几何和工况快照，供每次使用前复核
mesh.externalFlow=struct('edge',edge,'faceCenter',g.center, ...
    'faceNormal',g.normal,'V',mesh.V,'F',mesh.F,'freestream',freestream, ...
    'method',edge.method,'caseID',edge.caseID);
mesh.boundaryLayerEdge=edge;
if isfield(mesh,'heating'), mesh=rmfield(mesh,'heating'); end
if isfield(mesh,'heatingCoverage'), mesh=rmfield(mesh,'heatingCoverage'); end
if isfield(mesh,'boundaryLayer'), mesh=rmfield(mesh,'boundaryLayer'); end
if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
if isfield(mesh,'surfaceStreamlines'), mesh=rmfield(mesh,'surfaceStreamlines'); end
end
