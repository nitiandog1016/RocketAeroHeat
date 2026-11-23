function p=readShortSupportTimePair(folder)
% 模块用途：绑定稳态基线、共同扰动初帧、两种 BDF1 步长以及两个独立终点。
file=fullfile(folder,'supportedTime65.json');
assert(strcmp(cm011FileSHA256(file),'42c79caa7d36d7e142407945919252bc826c526d8481f9a6b53dc47a71fdb499'));
p=jsondecode(fileread(file));baseline=readShortSupportAngularPair(folder);
assert(strcmp(p.schema,'short65-supported-time-pair-v1') && strcmp(p.steadyFrame,baseline.fineFrame) && ...
    strcmp(p.meshSHA256,baseline.fineMesh) && strcmp(p.bodySHA256,baseline.bodySHA256) && ...
    strcmp(p.transferSHA256,baseline.transferSHA256) && p.faces==baseline.faces && p.elapsedDU==.25);
assert(strcmp(p.fullEnergy.old,p.seed) && strcmp(p.firstHalfEnergy.old,p.seed) && strcmp(p.halfEnergy.old,p.firstHalfFrame) && ...
    p.fullEnergy.elapsedS==p.halfEnergy.elapsedS && p.fullEnergy.dt==2*p.halfEnergy.dt);
end
