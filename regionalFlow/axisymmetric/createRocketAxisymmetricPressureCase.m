function item=createRocketAxisymmetricPressureCase(options)
%CREATEROCKETAXISYMMETRICPRESSURECASE 完整箭体压力源积分的独立新工况工厂。
% 模块一：继承完整母线/法向层、远场、气体、给定Tw及均匀五量初猜。
% 模块二：仅增加明确压力源积分选择；旧工况工厂及旧输入合同不改字节。
% 模块三：本工厂的新源/新链不覆盖旧链；内部暖初猜需在独立研究中另行完整绑定。
if nargin<1,options=struct;end
assert(isstruct(options)&&isscalar(options),'RocketPressureCase:Options','工况选项须为单一结构。');
base=options;policy='cell-constant';
if isfield(base,'pressureSourceQuadrature')
    policy=base.pressureSourceQuadrature;base=rmfield(base,'pressureSourceQuadrature');
end
assert(((ischar(policy)&&isrow(policy))||(isstring(policy)&&isscalar(policy)&&~ismissing(policy)))&& ...
    ismember(char(policy),{'linear-planar','cell-constant'}), ...
    'RocketPressureCase:Policy','压力源积分只接受linear-planar或cell-constant。');policy=char(policy);
item=createRocketAxisymmetricCollarCase(base);
item.boundary.pressureSourceQuadrature=policy;
item.caseFactory='createRocketAxisymmetricPressureCase';item.options=options;
item.inputContract.pressureSourceQuadraturePolicy=policy;
item.inputContract.pressureSourceQuadratureMeaning='explicit planar-area pressure integral from selected reconstruction; not a solved new boundary';
item.inputContract.changedFromV71=['same full collar case and uniform initialization; explicit pressure source=' policy];
item.inputContract.sourceQuadratureIsPhysicalAccuracyValidation=false;
end
