function boundary=resolveSU2ThermalBoundary(surface,spec)
%RESOLVESU2THERMALBOUNDARY 统一等温、绝热和空间壁温的输入身份，不冒充壁体温度解。
% 模块一：校验全部原节点的温度定义；轴向表不允许静默外推或漏掉鼻/底面。
% 模块二：逐面温度为原三角节点线性平均，保留输入来源；绝热不产生给定Tw。
% 模块三：官方Windows可执行文件尚无逐节点热边界调用器，空间模式显式标为不能执行。
assert(isstruct(spec)&&isscalar(spec)&&isfield(spec,'mode'),'缺少热边界模式。');
assert(all(ismember(fieldnames(spec),{'mode','temperature','x','source'})),'不支持的热边界字段。');
mode=validatestring(spec.mode,{'isothermal','adiabatic','axialTemperature'});
nv=size(surface.V,1);nf=size(surface.F,1);node=nan(nv,1);givenQ=nan(nf,1);
source='用户输入';if isfield(spec,'source'),source=char(spec.source);end
switch mode
    case 'isothermal'
        assert(isfield(spec,'temperature') && ~isfield(spec,'x'),'等温壁必须且只能提供统一温度。');
        validateattributes(spec.temperature,{'numeric'},{'real','scalar','finite','positive'});node(:)=spec.temperature;
        label='给定等温壁面边界';executable=true;
    case 'adiabatic'
        assert(~isfield(spec,'temperature') && ~isfield(spec,'x'),'绝热壁不能同时指定壁温。');
        givenQ(:)=0;label='绝热边界：给定q=0；壁面流体温度属于流场解，非材料热响应';executable=true;
    case 'axialTemperature'
        assert(isfield(spec,'temperature') && isfield(spec,'x'),'空间壁温需要项目X和配对温度表。');
        validateattributes(spec.x,{'numeric'},{'vector','real','finite','nonempty'});
        validateattributes(spec.temperature,{'numeric'},{'vector','real','finite','positive','numel',numel(spec.x)});
        xx=spec.x(:);tt=spec.temperature(:);
        assert(numel(xx)>=2 && all(diff(xx)>0),'X表必须严格递增且至少两点。');
        assert(xx(1)<=min(surface.V(:,1)) && xx(end)>=max(surface.V(:,1)), ...
            '温度表未覆盖全模型；须明确未测鼻端/底部延拓，不能自动外推。');
        node=interp1(xx,tt,surface.V(:,1),'linear');label='给定空间壁温：节点线性插值；引擎逐节点调用尚未接入';executable=false;
end
face=mean(node(surface.F),2);
boundary=struct('mode',mode,'source',source,'method',label,'specification',spec, ...
    'givenTemperatureByNode',node,'givenTemperatureByFace',face,'givenHeatFluxByFace',givenQ, ...
    'heatFluxIsBoundaryInput',strcmp(mode,'adiabatic'),'wallMaterialTemperatureSolved',false, ...
    'stockExecutableSupported',executable,'physicsQualified',false);
end
