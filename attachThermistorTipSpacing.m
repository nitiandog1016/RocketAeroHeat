function mesh=attachThermistorTipSpacing(mesh)
% 模块一：按几何、体网、原始场与热边界的指纹，读取前引线端面局部细化参照。
root=fileparts(mfilename('fullpath'));name='sidewindTipSpacing';
key=jsondecode(fileread(fullfile(root,'cases','thermistor_b',[name '.json'])));
[folder,guard]=openThermistorArchive(name); %#ok<NASGU>
s=jsondecode(fileread(fullfile(folder,'summary.json')));base=mesh.thermistorSidewind;
spec=jsondecode(fileread(fullfile(folder,'case.json')));body=jsondecode(fileread(fullfile(folder,'body.json')));
assert(strcmp(base.meshSHA256,key.baseMeshSHA256) && strcmp(base.caseSHA256,key.baseCaseSHA256) && ...
    strcmp(base.sourceFrameSHA256,mesh.evidenceIdentity.sourceFrame.sha256) && strcmp(base.bodySHA256,key.bodySHA256));
for field={'sourceFrameSHA256','sourceReactionSHA256','caseSHA256','meshSHA256','bodySHA256','wallSHA256'}
    k=field{1};assert(strcmp(s.(k),key.(k)));
end
files={'body.json','bodySHA256';'case.json','caseSHA256';'wall.csv','wallSHA256';'matched_temperature.csv','sampleSHA256'};
for j=1:size(files,1)
    assert(strcmp(cm011FileSHA256(fullfile(folder,files{j,1})),key.(files{j,2})));
end
c=mesh.hybridCondition;
assert(isequal(body.V,mesh.V) && isequal(body.F,mesh.F) && isequal(c.velocityMps,[0 -3 0]) && ...
    isequal(spec.velocityMps(:)',c.velocityMps) && spec.wallK==c.wallTemperatureK && spec.wallK==299.15 && ...
    spec.airK==c.ambientTemperatureK && spec.airPa==c.ambientPressurePa && spec.outerRadiusOverD==20 && spec.sectorDeg==360);
assert(spec.frontCapFirstLayerRatio==.5 && s.energy.domainCopies==1 && strcmp(s.energy.inputFrameSHA256,s.sourceFrameSHA256) && ...
    abs(s.energy.relativeDifference)<.005 && s.relativeChange20<.002 && all(structfun(@(v)v<-8,s.residualWindow)) && ...
    all(structfun(@(v)v<-8,s.plus20Window)));
rows=readtable(fullfile(folder,'wall.csv'));sample=readtable(fullfile(folder,'matched_temperature.csv'));
w=mesh.hybridWall;f=mesh.F;v=mesh.V;
a=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(isequal(rows.faceID,w.faceID,sample.faceID) && isequal(rows.regionID,body.regionID) && ...
    all(abs(a-rows.areaM2)<1e-18) && all(isfinite(rows.heatFluxWm2)) && all(rows.pressurePa>0) && ...
    all(abs(rows.pressurePa-(spec.airPa+rows.pressureGaugePa))<1e-9));
cap=rows.regionID==7;
assert(all(sample.fineDistanceM>0) && all(sample.fineDistanceM<=sample.baseDistanceM+1e-18) && ...
    all(abs(sample.fineDistanceM(cap)-sample.baseDistanceM(cap)/2)<1e-17));
temps=[rows.nearWallCellTemperatureK sample.matchedK];
assert(all(isfinite(temps) & temps>=spec.airK & temps<=spec.wallK,'all'));
% 模块二：热流按同面比较，温差使用两网共有的原首层节点，另列实际细网首层温度。
g=struct('baseFrame',base.sourceFrameSHA256,'fineFrame',s.sourceFrameSHA256,'seedFrame',key.seedFrameSHA256, ...
    'fineReaction',s.sourceReactionSHA256,'baseReaction',base.sourceReactionSHA256,'baseMeshSHA256',base.meshSHA256, ...
    'fineMeshSHA256',s.meshSHA256,'fineWallSHA256',s.wallSHA256,'sampleSHA256',key.sampleSHA256, ...
    'outwardBaseW',-a'*w.heatFluxWm2,'outwardFineW',-a'*rows.heatFluxWm2, ...
    'basePeakFace',base.peakFaceID,'finePeakFace',s.peakFaceID,'commonLocation','原 20D 首层物理节点', ...
    'scope','带线 B 型热敏珠，3 m/s 沿 −Y 侧吹，20D 外域，给定壁温 299.15 K。');
% 同时保留第二强面与两者的实际差值，镜像面的极小差别可能改变排序。
[q0,i0]=mink(w.heatFluxWm2,2);[q1,i1]=mink(rows.heatFluxWm2,2);
assert(i0(1)==g.basePeakFace && i1(1)==g.finePeakFace && abs(g.outwardFineW-s.totalOutwardW)<1e-12);
g.baseSecondFace=i0(2);g.fineSecondFace=i1(2);
g.basePeakGapWm2=diff(q0);g.finePeakGapWm2=diff(q1);
for j=1:numel(s.regions)
    mask=rows.regionID==j;
    g.regions(j)=struct('name',s.regions(j).name,'faces',nnz(mask), ...
        'baseOutW',-a(mask)'*w.heatFluxWm2(mask),'fineOutW',-a(mask)'*rows.heatFluxWm2(mask));
    assert(abs(g.regions(j).fineOutW-s.regions(j).outwardW)<1e-12);
end
mesh.thermistorTipSpacing=struct('summary',g,'fineHeat',rows.heatFluxWm2,'deltaHeat',rows.heatFluxWm2-w.heatFluxWm2, ...
    'fineGauge',rows.pressureGaugePa,'deltaPressure',rows.pressurePa-w.pressurePa,'fineFirstK',rows.nearWallCellTemperatureK, ...
    'matchedK',sample.matchedK,'deltaMatchedK',sample.matchedK-w.nearWallCellTemperatureK,'fineDistanceM',sample.fineDistanceM);
end
