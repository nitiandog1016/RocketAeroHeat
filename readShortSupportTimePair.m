function p=readShortSupportTimePair(folder,radiusD)
% 模块用途：绑定稳态基线、共同扰动初帧、两种 BDF1 步长以及两个独立终点。
if nargin<2,radiusD=14;end
baseline=readShortSupportAngularPair(folder);
if isequal(radiusD,14)
    file=fullfile(folder,'supportedTime65.json');
    assert(strcmp(cm011FileSHA256(file),'42c79caa7d36d7e142407945919252bc826c526d8481f9a6b53dc47a71fdb499'));
    frame=baseline.fineFrame;mesh=baseline.fineMesh;
else
    assert(isequal(radiusD,20),'RAH:Short:Domain','可选远场半径为 14D 或 20D。');
    file=fullfile(folder,'outer20','supportedTime65.json');
    assert(strcmp(cm011FileSHA256(file),'07df1aa77303c11e81d4ebd2cd92861501cd9d24e58df9fd999e4ab45ab576d3'), ...
        'RAH:Short:Domain','20D 的归档时间层与登记结果不一致。');
    frame='8a6efe752ffd5d913196e884e51880deb686de78b20a54159571466f9e36bac0';
    mesh='da1ed66545237069b472c09c33f35d0c72b86d6d6e83a7f2c2de2b05c63227ec';
end
p=jsondecode(fileread(file));
assert(strcmp(p.schema,'short65-supported-time-pair-v1') && strcmp(p.steadyFrame,frame) && ...
    strcmp(p.meshSHA256,mesh) && strcmp(p.bodySHA256,baseline.bodySHA256) && ...
    strcmp(p.transferSHA256,baseline.transferSHA256) && p.faces==baseline.faces && p.elapsedDU==.25);
assert(strcmp(p.fullEnergy.old,p.seed) && strcmp(p.firstHalfEnergy.old,p.seed) && strcmp(p.halfEnergy.old,p.firstHalfFrame) && ...
    p.fullEnergy.elapsedS==p.halfEnergy.elapsedS && p.fullEnergy.dt==2*p.halfEnergy.dt);
end
