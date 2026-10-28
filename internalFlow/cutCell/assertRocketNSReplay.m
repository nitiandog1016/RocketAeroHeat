function report=assertRocketNSReplay(rate,info,expectedRate,expectedInfo,volume)
%ASSERTROCKETNSREPLAY 逐值核对状态残差与局部通量，仅全域归约使用前向舍入界。
% MATLAB多线程sum的加法次序可能改变最后数位；不放宽物理守恒或稳态门槛。
assert(isequaln(rate,expectedRate)&&isequal(fieldnames(info),fieldnames(expectedInfo)), ...
    'PeakSystem:Replay','当前残差或诊断字段未逐值复现。');
reduced={'volumeIntegratedResidual','balanceResidual'};
for key=fieldnames(info)'
    name=key{1};if ismember(name,reduced),continue;end
    assert(isequaln(info.(name),expectedInfo.(name)),'PeakSystem:Replay','当前NS字段%s未逐值复现。',name);
end
n=size(rate,1);gamma=n*eps/(1-n*eps);sumAbsolute=sum(abs(rate.*volume(:)),1);
bound=2*gamma*sumAbsolute+8*eps*abs(info.netBoundaryFlux);
differences=zeros(2,5);
for k=1:2
    differences(k,:)=abs(info.(reduced{k})-expectedInfo.(reduced{k}));
    assert(all(differences(k,:)<=bound),'PeakSystem:Replay','全域归约差超出两个加法序列的前向舍入界。');
end
report=struct('rateAndLocalInfoBitwiseIdentical',true,'globalReductionAbsoluteDifference',differences, ...
    'globalReductionForwardRoundoffBound',bound,'physicalAcceptanceToleranceChanged',false);
end
