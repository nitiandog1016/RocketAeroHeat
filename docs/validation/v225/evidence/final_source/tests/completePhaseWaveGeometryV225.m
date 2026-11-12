function report=completePhaseWaveGeometryV225(folder)
%COMPLETEPHASEWAVEGEOMETRYV225 只验证已定位的环/尖端读取修复，不重算已保存波系。
d=load('outputs/v104_validation/baseline02/current_case.mat','current','fs');
saved=jsondecode(fileread('outputs/v225_validation/matlab02/waves.json'));
oracle=jsondecode(fileread('outputs/v225_validation/audit01/rocket_inlet.json'));
r=inspectPhaseAwareRocketInlet(d.current,d.fs,saved(4).result.maximumTurnDeg);
for name={'faces','tipVertexID','tipAdjacentFaces'},assert(r.(name{1})==oracle.(name{1}));end
for name={'areaM2','lastRingX','lastRingRadiusM','tipAdjacentAreaM2', ...
        'meridianInclinationDeg','windwardFirstTurnDeg','leewardFirstTurnDeg','maximumWeakShockTurnDeg'}
    assert(abs(r.(name{1})-oracle.(name{1}))<1e-12*max(1,abs(oracle.(name{1}))));
end
assert(r.bothMeridiansOutsideAttachedWeakBranch&&r.lastRingIsNotTip);
report=struct('passed',true,'scope','corrected original tip/ring geometry only', ...
    'waveRecalculation',false,'wholeMatlabSuiteRerun',false,'rocket',r);
fid=fopen(fullfile(folder,'report.json'),'w','n','UTF-8');clean=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));disp(report);
end
