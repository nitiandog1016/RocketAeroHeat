function validateShortTimeChain(p,target)
% 模块用途：按实际 BDF1 步长逐帧核对两条链，固定种子、真实前一帧、物理时钟和终点。
assert(p.elapsedDU==target && p.fullEnergy.elapsedS==p.elapsedS && p.halfEnergy.elapsedS==p.elapsedS && ...
    p.fullEnergy.dt==.25 && p.halfEnergy.dt==.125,'RAH:Short:Job','时间或步长与请求不符。');
assert(strcmp(p.fullEnergy.frame,p.fullFrame) && strcmp(p.halfEnergy.frame,p.halfFrame) && ...
    strcmp(p.firstHalfEnergy.frame,p.firstHalfFrame) && strcmp(p.firstHalfEnergy.old,p.seed), ...
    'RAH:Short:Job','时间层身份或初帧不符。');
if target==.25
    assert(strcmp(p.schema,'short65-supported-time-pair-v1') && ...
        strcmp(p.fullEnergy.old,p.seed) && strcmp(p.halfEnergy.old,p.firstHalfFrame),'RAH:Short:Job','早期时间层的来源不符。');
    return;
end
assert(target==1 && strcmp(p.schema,'short65-supported-time-continuation-v1'),'RAH:Short:Job','接续时刻的格式不符。');
for name={'full','half'}
    branch=name{1};chain=p.history.(branch);endPoint=p.([branch 'Energy']);dt=endPoint.dt;old=p.seed;
    assert(numel(chain)==target/dt && endPoint.step==numel(chain),'RAH:Short:Job','接续链步数不符。');
    for j=1:numel(chain)
        row=chain(j);
        assert(row.step==j && row.dt==dt && row.elapsedDU==j*dt && strcmp(row.old,old) && ...
            abs(row.elapsedS-j*dt*p.condition.diameterM/p.condition.speedMps)<1e-14 && ...
            abs(row.relativeBalance)<.005,'RAH:Short:Job','接续链存在错帧、时间错误或能量账不符。');
        old=row.frame;
    end
    assert(strcmp(chain(end-1).frame,endPoint.old) && strcmp(old,endPoint.frame) && ...
        strcmp(chain(round(.25/dt)).frame,p.timeComparison.(branch).earlyFrame),'RAH:Short:Job','早期端点或最终前一层不符。');
end
end
