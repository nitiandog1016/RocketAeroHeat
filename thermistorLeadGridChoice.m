function [archive,options,label]=thermistorLeadGridChoice(level)
% 按侧栏选择返回带线试件的体网、复算选项及显示名称，旧逻辑值对应前两档。
switch double(level)
    case 0
        archive='leads3';options={};label='原混合网';
    case 1
        archive='leadsMidbody3';options={'--refine-leads'};label='中段加密';
    case 2
        archive='leadsThird3';options={'--third-midbody'};label='第三档中段细化';
    otherwise
        error('RAH:Lead:Grid','请选择已登记的带线试件体网。');
end
end
