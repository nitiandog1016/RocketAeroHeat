function mesh=attachThermistorLayerBudget(mesh,area)
% 模块一：将给定壁温下的分带收支绑定到两份实际终场及本工况的原面顺序。
g=mesh.thermistorThreeGrids;
if ~strcmp(g.summary.mode,'givenWall'),return;end
path=fullfile(fileparts(mfilename('fullpath')),'cases','thermistor_b','layerBands.json');
digest=cm011FileSHA256(path);
assert(strcmp(digest,'764fd03400467414bde957d28b5294224c986a20d7dd8fa27731eadd91bb53df'),'热敏珠分带数据的来源指纹不符。');
b=jsondecode(fileread(path));frames=g.summary.grids(2:3);
for j=1:2
    src=b.sourceGrids(j);
    if ~strcmp(src.frameSHA256,frames(j).frameSHA256) || ~strcmp(src.meshSHA256,frames(j).meshSHA256) || ...
            ~strcmp(src.wallSHA256,frames(j).wallSHA256)
        return;
    end
end
s=mesh.thermistorLead;c=mesh.hybridCondition;ids=b.faceBand(:);
assert(strcmp(b.schema,'thermistor-b-band-heat-v1') && strcmp(b.mode,'givenWall') && ...
    strcmp(b.bodySHA256,s.bodySHA256) && strcmp(b.caseSHA256,s.caseSHA256) && ...
    c.wallTemperatureK==b.wallK && c.ambientTemperatureK==b.airK && isequal(c.velocityMps,[-b.speedMps 0 0]));
% 模块二：逐带重算两档壁面放热及其差，检查分带数量、面积和收支分项。
for k=1:numel(b.bands)
    band=b.bands(k);use=ids==band.id;
    assert(band.id==k && nnz(use)==band.faces && abs(sum(area(use))-band.areaM2)<1e-17);
    watts=-area(use)'*g.heat(use,2:3);delta=band.delta;
    assert(max(abs(watts-[band.medium.faceOutW band.fine.faceOutW]))<1e-11 && abs(diff(watts)-delta.faceOutW)<1e-11);
    assert(abs(delta.outerConvW+delta.outerDiffW+delta.sideConvW+delta.sideDiffW+delta.mapCorrectionW-delta.residualW-delta.faceOutW)<1e-11);
end
assert(numel(ids)==numel(area) && all(ismember(ids,[b.bands.id])));
b.packageSHA256=digest;mesh.thermistorLayerBudget=b;
mesh.evidenceIdentity.sourceFrame.layerBudget=struct('packageSHA256',digest,'sourceGrids',b.sourceGrids);
end
