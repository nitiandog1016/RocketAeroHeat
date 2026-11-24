function p=readShortSupportTimePair(folder,radiusD,target)
% 模块用途：绑定稳态基线、共同扰动初帧、两种 BDF1 步长以及两个独立终点。
if nargin<2,radiusD=14;end
if nargin<3,target=.25;end
baseline=readShortSupportAngularPair(folder);
if target==1
    assert(radiusD==14,'RAH:Short:Domain','一个对流时间的工况采用 14D 外域。');
    file=fullfile(folder,'later1','supportedTime65.json');
    assert(strcmp(cm011FileSHA256(file),'da201cd731a21927cd875e1d9a85d5987bc8e678a422b81c06e497c40fba1617'),'RAH:Short:Job','接续归档与登记终点不符。');
    frame=baseline.fineFrame;mesh=baseline.fineMesh;
elseif isequal(radiusD,14)
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
assert(strcmp(p.steadyFrame,frame) && ...
    strcmp(p.meshSHA256,mesh) && strcmp(p.bodySHA256,baseline.bodySHA256) && ...
    strcmp(p.transferSHA256,baseline.transferSHA256) && p.faces==baseline.faces);
validateShortTimeChain(p,target);
end
