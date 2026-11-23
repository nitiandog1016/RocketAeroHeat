function [rows,s,t,deltaFlux]=applyThermistorThermometry(rows,s,c,internalK,action)
% 模块一：绑定论文测温模型、原始长引线场和采用真实修正边界的独立原生解。
t=thermistorResistanceModel(c,internalK);
root=fileparts(mfilename('fullpath'));
key=jsondecode(fileread(fullfile(root,'cases','thermistor_b','resistance12.json')));
assert(abs(t.scale-key.factor)<1e-13 && abs(t.reynoldsL-key.figure12.reynoldsL)<1e-10);
t.figure12=key.figure12;t.matchedCheck=key.matchedCheck;
[original,guard]=openThermistorArchive('longWire12'); %#ok<NASGU>
ref=jsondecode(fileread(fullfile(original,'summary.json')));
old=readtable(fullfile(original,'wall.csv'));
assert(strcmp(ref.sourceFrameSHA256,key.referenceFrame) && isequal(rows.faceID,old.faceID) && ...
    isequal(rows.regionID,old.regionID) && isequal(rows.areaM2,old.areaM2));
t.referenceFrame=ref.sourceFrameSHA256;t.checkFrame=key.frameSHA256;t.referenceWallSHA256=ref.wallSHA256;
t.referenceBeadW=ref.beadAndJointsW;t.sourceBoundarySHA256=s.boundarySHA256;t.sourceFieldSHA256=s.wallSHA256;
t.sourceEnergy=s.energy;
% 模块二：对正比的全壁超温施加一个系数，保留真实静压及热流正负号。
if strcmp(action,'solve')
    assert(strcmp(s.sourceFrameSHA256,ref.sourceFrameSHA256) && strcmp(s.wallSHA256,ref.wallSHA256));
    f=t.scale;rows.heatFluxWm2=f*rows.heatFluxWm2;
    rows.nearWallCellTemperatureK=c.airK+f*(rows.nearWallCellTemperatureK-c.airK);
    rows.wallInputK=c.airK+f*(rows.wallInputK-c.airK);
    s.totalOutwardW=f*s.totalOutwardW;s.beadAndJointsW=f*s.beadAndJointsW;
    s.beadChange20W=f*s.beadChange20W;
    for j=1:numel(s.regions)
        s.regions(j).outwardW=f*s.regions(j).outwardW;
        s.regions(j).wallRangeK=c.airK+f*(s.regions(j).wallRangeK-c.airK);
    end
    s.profile.endK=c.airK+f*(s.profile.endK-c.airK);
    for k={'oneLeadRootW','oneLeadThroughCutW','oneLeadSideW'},s.profile.(k{1})=f*s.profile.(k{1});end
    % 参照能量账原样保存，结果账中列出按恒物性模型推得的两个功率。
    s.energy=struct('wallReleaseW',f*t.sourceEnergy.wallReleaseW,'farOutW',f*t.sourceEnergy.farOutW, ...
        'relativeDifference',t.sourceEnergy.relativeDifference,'sourceFrameSHA256',t.referenceFrame,'method','参照同态功率按全壁超温比换算');
    t.fieldFactor=f;t.fieldMethod='长引线原生参照场的全域超温线性换算';
else
    assert(strcmp(action,'archive') && strcmp(s.sourceFrameSHA256,key.frameSHA256) && strcmp(s.boundarySHA256,key.boundarySHA256));
    t.fieldFactor=1;t.fieldMethod='实际修正节点壁温的 SU2 原生解及两次壁反力';
end
% 模块三：同步换算用于老界面引线差值的均温控制，另列测温修正相对原长引线场的变化。
rows.uniformHeatFluxWm2=t.scale*rows.uniformHeatFluxWm2;
s.controlBeadW=t.scale*s.controlBeadW;s.controlTotalW=t.scale*s.controlTotalW;
for j=1:numel(s.regions),s.regions(j).controlW=t.scale*s.regions(j).controlW;end
s.beadChangeW=s.beadAndJointsW-s.controlBeadW;
s.boundarySHA256=key.boundarySHA256;
deltaFlux=rows.heatFluxWm2-old.heatFluxWm2;
t.apparentHWm2K=s.beadAndJointsW/(c.sourceAreaM2*(t.internalK-c.airK));
t.surfaceHWm2K=s.beadAndJointsW/(c.sourceAreaM2*(t.surfaceK-c.airK));
t.internalMinusSurfaceK=internalK-t.surfaceK;
s.beadHWm2K=t.surfaceHWm2K;
end
