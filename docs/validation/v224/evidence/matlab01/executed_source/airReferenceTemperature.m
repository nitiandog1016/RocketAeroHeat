function [Tref,info] = airReferenceTemperature(Te,Tw,Taw,freestream,cfg)
%AIRREFERENCETEMPERATURE 用相容焓模型计算参考温度及壁温导数。
% 模块作用：Riley NASA-TM-107838 附录H.6：
% h*=0.5*(he+hw)+0.22*(haw-he)。默认变cp反解h(T*)，不混用线性T。
% cfg.heating.referencePropertyModel='linearTemperature' 显式复现旧式；
% 定cp模式两式等价。范围外保持NaN，不借参考平均值掩盖无效端点。
if nargin<5, cfg=struct; end
[mode,~]=selectAirThermalModel(freestream,cfg);
model='referenceEnthalpy';
if isfield(cfg,'heating') && isfield(cfg.heating,'referencePropertyModel')
    model=char(cfg.heating.referencePropertyModel);
end
if ~ismember(model,{'referenceEnthalpy','linearTemperature'})
    error('RocketAeroHeat:ReferenceState:Model','未知参考性质模型。');
end
validateattributes(Te,{'numeric'},{'real'});
validateattributes(Tw,{'numeric'},{'real'});
validateattributes(Taw,{'numeric'},{'real'});
Te=double(Te); Tw=double(Tw); Taw=double(Taw);
if isscalar(Tw), Tw=repmat(Tw,size(Te)); end
if ~isequal(size(Te),size(Tw),size(Taw))
    error('RocketAeroHeat:ReferenceState:Size','外缘、壁温和恢复温度尺寸必须一致。');
end
Tref=nan(size(Te)); derivative=Tref; href=Tref;
valid=isfinite(Te)&isfinite(Tw)&isfinite(Taw)&Te>0&Tw>0&Taw>0;
info=struct('model',model,'thermalPropertyMode',mode, ...
    'source','Riley NASA-TM-107838 Appendix H.6/H.7', ...
    'method','Eckert 参考焓：h*=0.5(he+hw)+0.22(haw-he)');

%% 1. 定cp及显式旧公式；旧式对照仍拒绝变cp端点越域
if strcmp(mode,'thermallyPerfect')
    [~,he]=airThermalProperties(Te);
    [cpw,hw]=airThermalProperties(Tw);
    [~,haw]=airThermalProperties(Taw);
    valid=valid&isfinite(he)&isfinite(hw)&isfinite(haw);
    href=0.5*(he+hw)+0.22*(haw-he);
else
    validateattributes(freestream.gamma,{'numeric'},{'scalar','finite','>',1});
    validateattributes(freestream.R,{'numeric'},{'scalar','finite','positive'});
    cp=freestream.gamma*freestream.R/(freestream.gamma-1);
    href=cp*(Te+0.5*(Tw-Te)+0.22*(Taw-Te));
end
if strcmp(model,'linearTemperature') || strcmp(mode,'caloricallyPerfect')
    Tref(valid)=Te(valid)+0.5*(Tw(valid)-Te(valid))+0.22*(Taw(valid)-Te(valid));
    derivative(valid)=0.5;
    if strcmp(model,'linearTemperature')
        info.method='Eckert 线性参考温度（显式旧版对照）';
    else
        info.method='Eckert 参考焓（定cp极限等价于线性参考温度）';
    end
else
    %% 2. 单调焓反解及残差审核；低温/分段系数原样保持
    ids=find(valid); low=80*ones(size(ids)); high=2500*ones(size(ids));
    target=href(ids);
    for k=1:46
        middle=(low+high)/2;
        [~,hmid]=airThermalProperties(middle);
        below=hmid<target;
        low(below)=middle(below); high(~below)=middle(~below);
    end
    temperature=(low+high)/2;
    [cpr,checked]=airThermalProperties(temperature);
    good=isfinite(cpr)&cpr>0&abs(checked-target)<=2e-11*max(1,abs(target));
    Tref(ids(good))=temperature(good);
    % 固定外缘/恢复状态时 dh*/dTw=0.5 cp(Tw)。
    derivative(ids(good))=0.5*cpw(ids(good))./cpr(good);
end
%% 3. 返回独立方法身份，调用方不能把无效参考状态计为求解
info.targetReferenceEnthalpy=href;
if strcmp(mode,'thermallyPerfect')
    [~,info.referenceEnthalpy]=airThermalProperties(Tref);
else
    info.referenceEnthalpy=cp*Tref;
end
info.temperatureDerivativeByWallTemperature=derivative;
info.validFace=valid&isfinite(Tref)&isfinite(derivative);
info.referenceEnthalpy(~info.validFace)=NaN;
end
