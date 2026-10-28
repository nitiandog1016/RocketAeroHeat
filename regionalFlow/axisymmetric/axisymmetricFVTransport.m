function properties=axisymmetricFVTransport(T,gas)
%AXISYMMETRICFVTRANSPORT 按显式参考合同求局部层流mu(T)、k(T)，不裁温度。
% 模块一：常物性输入保持原mu/k；Sutherland参考参数由工况明确给定。
% 模块二：定比热、常Pr下k=cp*mu/Pr，EOS和声速仍由同一gamma/R确定。
% 模块三：面/源项温度外推超出定义域必须拒绝，不以参考温度或端点替代。
gas=axisymmetricFVGas(gas);
assert(isnumeric(T)&&isreal(T)&&~isempty(T)&&all(isfinite(T),'all')&& ...
    all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2),'all'), ...
    'AxisFV:TransportTemperature','局部输运温度超出显式接受域，未裁值。');
T=double(T);
if isfield(gas,'transport')
    tr=gas.transport;ratio=T/tr.referenceTemperature;
    mu=tr.referenceViscosity*ratio.^1.5.* ...
        ((tr.referenceTemperature+tr.sutherlandTemperature)./(T+tr.sutherlandTemperature));
    k=mu*(gas.cp/gas.Pr);
else
    mu=repmat(gas.mu,size(T));k=repmat(gas.k,size(T));
end
assert(all(isfinite(mu)&mu>0,'all')&&all(isfinite(k)&k>0,'all'), ...
    'AxisFV:TransportState','局部输运系数非有限或非正，未裁mu/k。');
properties=struct('mu',mu,'k',k);
end
