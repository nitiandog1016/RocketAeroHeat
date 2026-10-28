function gas=axisymmetricFVGas(gas)
%AXISYMMETRICFVGAS 明确定比热EOS、常Pr及常数/Sutherland输运合同。
% 模块一：旧gas.mu/k保持常数含义和完整输出，旧输入不增加默认字段。
% 模块二：显式transport时mu/k表示参考温度值，实际系数由温度独立求得。
% 不自动切换NASA7、补空气常数或推断湍流；每个研究工况须冻结参数来源。
required={'gamma','R','mu','Pr','temperatureRange'};
assert(isstruct(gas)&&isscalar(gas)&&all(isfield(gas,required)),'AxisFV:Gas','缺少完整气体合同。');
for j=1:4
    validateattributes(gas.(required{j}),{'numeric'},{'real','finite','scalar','positive'});
    gas.(required{j})=double(gas.(required{j}));
end
validateattributes(gas.temperatureRange,{'numeric'},{'real','finite','positive','numel',2});
gas.temperatureRange=double(gas.temperatureRange);
assert(gas.gamma>1&&gas.temperatureRange(2)>gas.temperatureRange(1),'AxisFV:Gas','gamma/温域无效。');
gas.cp=gas.gamma*gas.R/(gas.gamma-1);gas.cv=gas.R/(gas.gamma-1);k=gas.mu*gas.cp/gas.Pr;
if isfield(gas,'k'),assert(isreal(gas.k)&&isscalar(gas.k)&&isfinite(gas.k)&&abs(gas.k/k-1)<1e-12,'AxisFV:Gas','k与mu*cp/Pr不一致。');end
gas.k=k;
if isfield(gas,'transport')
    tr=gas.transport;fields={'model','referenceViscosity','referenceTemperature','sutherlandTemperature'};
    assert(isstruct(tr)&&isscalar(tr)&&all(isfield(tr,fields))&& ...
        all(ismember(fieldnames(tr),fields))&& ...
        ((ischar(tr.model)&&isrow(tr.model))||(isstring(tr.model)&&isscalar(tr.model)))&&strcmp(tr.model,'sutherland'), ...
        'AxisFV:GasTransport','输运须显式给出Sutherland模型和三个参考参数。');
    for key=fields(2:end)
        value=tr.(key{1});
        assert(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value>0, ...
            'AxisFV:GasTransport','Sutherland参数必须为有限正标量。');
        tr.(key{1})=double(value);
    end
    assert(abs(gas.mu/tr.referenceViscosity-1)<1e-12,'AxisFV:GasTransport', ...
        '显式变输运时gas.mu必须是相同参考温度的muRef，不能混作入口黏度。');
    gas.transport=tr;
end
end
