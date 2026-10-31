function boundary=shockFittedCharacteristicSpeed3D(downstreamQ,normalOut,upstream,gas)
%SHOCKFITTEDCHARACTERISTICSPEED3D 用下游声学不变量与上游RH求局部法向激波速度。
% 模块一：外法向指向上游，采用J+=u_n+2a/(gamma-1)，不是直接指定激波后压力。
% 模块二：在M_rel>1分支有界求根；最后以同cp移动RH核对五量和不变量。
% 模块三：这是NASA2020方法1的常gamma局部法向形式，未实现完整三维特征导数方法。
assert(strcmp(gas.thermalPropertyMode,'caloricallyPerfect'),'ShockFit:CharacteristicGas', ...
    '局部声学不变量相容当前只实现常gamma，不回退NASA7。');
N=size(downstreamQ,1);
assert(isequal(size(normalOut),[N 3])&&all(isfinite(normalOut),'all')&& ...
    all(abs(sum(normalOut.^2,2)-1)<1e-12),'ShockFit:CharacteristicNormal','需要实际单位法向。');
down=shockFittedState3D(downstreamQ,gas,'conservative');
up=shockFittedState3D([upstream.density upstream.velocityVector' upstream.temperature],gas,'primitive');
assert(abs(up.pressure-upstream.pressure)<1e-10*upstream.pressure,'ShockFit:EOS','明确上游状态方程不一致。');
gamma=gas.gamma;un=normalOut*upstream.velocityVector;
invariant=sum(down.velocity.*normalOut,2)+2*down.soundSpeed/(gamma-1);
target=(invariant-un)/up.soundSpeed;
assert(all(target>2/(gamma-1)),'ShockFit:CharacteristicWave','下游不变量没有对应的压缩激波分支。');
lo=ones(N,1);hi=2*lo;
for iteration=1:40
    bad=curve(hi)<target;if ~any(bad),break;end
    hi(bad)=2*hi(bad);
end
assert(all(curve(hi)>=target)&all(isfinite(hi)),'ShockFit:CharacteristicBracket','压缩激波速度未括住，不提交状态。');
for iteration=1:80
    middle=(lo+hi)/2;below=curve(middle)<target;
    lo(below)=middle(below);hi(~below)=middle(~below);
end
mach=(lo+hi)/2;speed=un+up.soundSpeed*mach;
jump=movingShockJump3D(upstream,gas,-normalOut,-speed);
recovered=sum(jump.state.velocity.*normalOut,2)+2*jump.state.soundSpeed/(gamma-1);
residual=(recovered-invariant)/up.soundSpeed;
assert(max(abs(residual))<2e-11,'ShockFit:CharacteristicCompatibility','移动RH与原下游不变量不相容。');
boundary=struct('normalOut',normalOut,'normalSpeedMS',speed,'relativeShockMach',mach, ...
    'downstreamInvariantMS',invariant,'recoveredInvariantMS',recovered, ...
    'invariantResidualByUpstreamSoundSpeed',residual,'jump',jump, ...
    'method','常gamma局部法向J+与压缩RH相容/有界激波Mach根；不含完整三维特征导数', ...
    'pressureGivenAsBoundaryCondition',false,'productionQualified',false,'newProductionSolvedAreaM2',0);

    function value=curve(M)
        ratio=(gamma+1)*M.^2./(2+(gamma-1)*M.^2);
        pressure=1+2*gamma/(gamma+1)*(M.^2-1);
        value=M.*(1-1./ratio)+2/(gamma-1)*sqrt(pressure./ratio);
    end
end
