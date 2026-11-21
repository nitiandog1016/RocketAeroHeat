function mesh=attachThermistorPairedControl(mesh,kind)
% 模块一：按已验收的模型、冻结工况和源帧配对。其他已算结果仍可独立使用。
kind=validatestring(kind,{'wired','historical','midbody','wiredRefined'});
[folder,guard]=openThermistorArchive('paired3'); %#ok<NASGU>
c=jsondecode(fileread(fullfile(folder,'comparison.json')));
id=mesh.evidenceIdentity.sourceFrame;condition=mesh.hybridCondition;
if strcmp(kind,'wiredRefined')
    s=mesh.thermistorLead;
    sourceMatches=strcmp(id.bodySHA256,c.bodies.wired) && strcmp(id.sha256,s.sourceFrameSHA256) && strcmp(id.wallSHA256,s.wallSHA256);
elseif strcmp(kind,'midbody')
    s=mesh.thermistorMidbody;
    sourceMatches=strcmp(id.bodySHA256,c.bodies.smooth) && strcmp(id.sha256,s.sourceFrameSHA256) && strcmp(id.wallSHA256,s.wallSHA256);
else
    sourceMatches=strcmp(id.sha256,c.sourceFrames.(kind)) && strcmp(id.bodySHA256,c.bodies.(kind)) && strcmp(id.wallSHA256,c.primaryTables.(kind));
end
if ~sourceMatches || ~isequal(condition.velocityMps,[-3 0 0]) || ...
        condition.ambientTemperatureK~=293.15 || condition.ambientPressurePa~=101325 || ...
        ~isfield(condition,'wallTemperatureK') || condition.wallTemperatureK~=305.15
    return;
end
n=size(mesh.F,1);w=mesh.hybridWall;v=mesh.V;f=mesh.F;
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
q=nan(n,1);pressure=q;reverse=q;reverseChange=q;
if any(strcmp(kind,{'wired','wiredRefined'}))
    % 模块二：两种几何按共同三角面比较。引线、接头和线端的光滑对照保持未知。
    mapFile=fullfile(folder,'shared.csv');referenceFolder=folder;
    if strcmp(kind,'wiredRefined')
        % 两种加密网共用原面配对表，参考压力和热流取加密光滑解的真实终帧。
        [referenceFolder,refGuard]=openThermistorArchive('midbody3'); %#ok<NASGU>
        ref=jsondecode(fileread(fullfile(referenceFolder,'summary.json')));
        refFlow=jsondecode(fileread(fullfile(referenceFolder,'flow.json')));
        assert(strcmp(ref.bodySHA256,c.bodies.smooth) && strcmp(ref.sourceFrameSHA256,refFlow.sourceFrameSHA256));
        tableHash=ref.wallSHA256;flowHash=refFlow.tableSHA256;
    else
        tableHash=c.primaryTables.smooth;flowHash=c.flow.smoothNative.tableSHA256;
    end
    bodyFile=fullfile(referenceFolder,'body.json');wallFile=fullfile(referenceFolder,'wall.csv');
    assert(strcmp(cm011FileSHA256(mapFile),c.sharedMapSHA256) && strcmp(cm011FileSHA256(bodyFile),c.bodies.smooth) && ...
        strcmp(cm011FileSHA256(wallFile),tableHash),'RAH:Bead:Pair','公共面来源不符。');
    pair=readtable(mapFile);b=jsondecode(fileread(bodyFile));table=readtable(wallFile);
    flowFile=fullfile(referenceFolder,'flow.csv');assert(strcmp(cm011FileSHA256(flowFile),flowHash));
    flow=readtable(flowFile);i=pair.wiredID;j=pair.smoothID;
    assert(numel(i)==c.commonFaces && numel(unique(i))==numel(i) && numel(unique(j))==numel(j) && ...
        isequal(flow.face_id,table.faceID+1) && all(table.regionID(j)==2));
    for k=1:3
        assert(isequal(v(f(i,k),:),b.V(b.F(j,k),:)),'RAH:Bead:Pair','公共三角面的顶点或绕序不符。');
    end
    assert(max(abs(area(i)-table.areaM2(j)))<1e-18 && abs(sum(area(i))-c.commonAreaM2)<1e-16);
    q(i)=table.heatFluxWm2(j);pressure(i)=table.pressurePa(j);reverse(i)=flow.reverse_both(j);
    dq=w.heatFluxWm2-q;dp=w.pressurePa-pressure;
    reverseChange(i)=mesh.thermistorLeadFlow.reverse(i)-reverse(i);
    order='带线减光滑';scope='仅限两个原面三个顶点、面积与几何位置相同的珠面。';
    if strcmp(kind,'wiredRefined')
        order='带线减光滑 · 两网均中段加密';
        c=struct('powerW',struct('pairedSmooth',ref.totalOutwardW,'wiredBeadAndJoints',s.beadAndJointsW), ...
            'sourceFrames',struct('wired',id.sha256,'smooth',ref.sourceFrameSHA256),'commonFaces',numel(i), ...
            'commonAreaM2',sum(area(i)),'sharedMapSHA256',c.sharedMapSHA256,'referenceWallSHA256',tableHash);
    end
elseif strcmp(kind,'midbody')
    % 模块三：新旧混合网拥有相同的原面、顺序和壁温，可直接逐面比较。
    wallFile=fullfile(folder,'wall.csv');flowFile=fullfile(folder,'flow.csv');bodyFile=fullfile(folder,'body.json');
    assert(strcmp(cm011FileSHA256(bodyFile),id.bodySHA256) && strcmp(cm011FileSHA256(wallFile),c.primaryTables.smooth) && ...
        strcmp(cm011FileSHA256(flowFile),c.flow.smoothNative.tableSHA256));
    body=jsondecode(fileread(bodyFile));base=readtable(wallFile);flow=readtable(flowFile);
    assert(isequal(v,body.V) && isequal(f,body.F) && isequal(base.faceID,(0:n-1)') && ...
        isequal(flow.face_id,base.faceID+1) && all(abs(base.areaM2-area)<1e-18));
    q=base.heatFluxWm2;pressure=base.pressurePa;reverse=flow.reverse_both;
    assert(abs(area'*q+c.powerW.pairedSmooth)<1e-12 && abs(area'*w.heatFluxWm2+s.totalOutwardW)<1e-12);
    dq=w.heatFluxWm2-q;dp=w.pressurePa-pressure;reverseChange=mesh.thermistorLeadFlow.reverse-reverse;
    order='中段加密减原混合网';scope='光滑热敏珠 6,624 个相同三角面。两次独立冷态求解，共用壁面、首层和工况。';
else
    % 模块三：光滑旧网格上的两解都按原 5,216 面投影。用面积分复核两套壁面功率。
    assert(strcmp(kind,'historical'));
    path=fullfile(folder,'historical.csv');assert(strcmp(cm011FileSHA256(path),c.historicalTransferSHA256));
    t=readtable(path);
    assert(isequal(t.faceID,(1:n)') && max(abs(t.areaM2-area))<1e-18 && ...
        abs(area'*t.pairedHeatFluxWm2+c.powerW.pairedSmooth)<1e-12 && ...
        abs(area'*w.heatFluxWm2+c.powerW.oldSmooth)<1e-12,'RAH:Bead:Pair','光滑原面功率或顺序不符。');
    q=t.pairedHeatFluxWm2;pressure=t.pairedPressurePa;reverse=t.pairedReverse;
    dq=q-w.heatFluxWm2;dp=pressure-w.pressurePa;reverseChange=reverse-t.originalReverse;
    order='混合网减原结构网';scope='同一光滑椭球的两次计算按原 5,216 面比较；壁热按面积守恒投影。';
end
mesh.thermistorPaired=struct('kind',kind,'summary',c,'heatFluxWm2',q,'deltaHeatWm2',dq, ...
    'deltaPressurePa',dp,'reverse',reverse,'reverseChange',reverseChange,'order',order,'scope',scope);
end
