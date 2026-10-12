function testExternalFlowImport()
% 外流状态逐面元匹配及工况过期检测。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); fs=calcFreestream(cfg.flight);
m.V=[0 0 0;0 1 0;0 0 1]; m.F=[1 2 3];
m.localGeometry.center=mean(m.V,1);
m.localGeometry.normal=[1 0 0];
d=struct('faceCenter',m.localGeometry.center, ...
    'faceNormal',m.localGeometry.normal,'pressure',fs.pressure, ...
    'temperature',fs.temperature,'density',fs.pressure/(fs.R*fs.temperature), ...
    'velocity',[0 fs.speed 0],'method','解析均匀测试流','caseID','test-M5', ...
    'freestream',fs);
m=importExternalFlowState(m,fs,d);
m=calcBoundaryLayerEdgeState(m,fs,cfg);
assert(contains(m.boundaryLayerEdge.method,'逐面元导入'));
assert(abs(m.boundaryLayerEdge.Mach-fs.Mach)<1e-10);
assert(abs(m.boundaryLayerEdge.totalPressure-fs.totalPressure)<1e-8*fs.totalPressure);
assert(isequal(m.boundaryLayerEdge.velocityVector,d.velocity));
bad=d; bad.faceCenter=bad.faceCenter+[1 0 0];
mustFail(@()importExternalFlowState(m,fs,bad));
bad=d; bad.density=2*bad.density;
mustFail(@()importExternalFlowState(m,fs,bad));
changed=fs; changed.alpha=10;
mustFail(@()calcBoundaryLayerEdgeState(m,changed,cfg));
fprintf('External-flow import passed: mapping, gas state, provenance and stale-case checks.\n');
end
function mustFail(callback)
failed=false; try, callback(); catch, failed=true; end
assert(failed);
end
