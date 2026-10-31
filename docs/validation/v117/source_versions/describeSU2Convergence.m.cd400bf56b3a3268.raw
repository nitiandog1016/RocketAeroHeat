function lines=describeSU2Convergence(result)
%DESCRIBESU2CONVERGENCE 将真实进程、末次监测与显示检查点分开呈现。
% 本模块只格式化审核，不把退出0、残差通过或完整云图等同于物理验证。
r=result.receipt;process=sprintf('进程 exit %g',r.exitCode);
if r.cancelled,process=[process ' · 已取消'];elseif r.timedOut,process=[process ' · 预算超时'];else,process=[process ' · 已结束'];end
if ~isfield(result,'convergenceAudit')
    lines={process;'显示检查点迭代未知 · 旧快照缺少完整收敛审核';'研究值 · 未获网格及物理资格'};return;
end
a=result.convergenceAudit;
if isfinite(a.terminalIteration)
    terminal=sprintf('末次监测 %d 步 · %s',a.terminalIteration,pass(a.terminalResidualPassed));
else,terminal='末次监测未知';end
if a.checkpointIdentityConfirmed
    checkpoint=sprintf('显示检查点 %d 步 · %s',a.checkpointIteration,pass(a.checkpointResidualPassed));
else,checkpoint='显示检查点迭代未确认';end
detail=a.reason;
if isfinite(a.threshold)
    values=strtrim(sprintf(' %.5g',a.checkpointResiduals));
    detail=sprintf('log10残差 [%s] / 门槛 %.4g · %s',values,a.threshold,strjoin(a.fields,', '));
    if ~a.checkpointIdentityConfirmed,detail=a.reason;end
end
lines={[process ' · ' terminal];[checkpoint ' · 研究值，未获物理资格'];detail};
end
function text=pass(yes)
if yes,text='原残差门槛通过';else,text='原残差门槛未通过/未确认';end
end
