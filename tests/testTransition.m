function testTransition()
% 阈值两侧、等值和无效面元都应给出确定类别。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();');
Re=[1e5;5e5;1e6;NaN]; mask=[true;true;true;false];
t=calcTransition(Re,mask,cfg);
assert(isequal(t.laminarFace,[true;false;false;false]));
assert(isequal(t.turbulentFace,[false;true;true;false]));
assert(isnan(t.intermittency(4)));
assert(contains(t.method,'Reynolds'));
% 已发生的湍流沿有向流线传递；局部 Re 下降不自动回到层流。
path=struct('predecessorFace',[0;1;2;3;4], ...
    'pathCost',(0:4)','connectedFace',true(5,1));
pathRe=[1e5;6e5;2e5;3e5;1e5];
tracked=calcTransition(pathRe,true(5,1),cfg,path);
assert(isequal(tracked.turbulentFace,[false;true;true;true;true]));
assert(isequal(tracked.onsetFace,[false;true;false;false;false]));
assert(tracked.summary.inheritedTurbulentFaces==3);
% 热流空白不能把外流状态路径切断；隐藏的首次转捩仍由下游继承。
heatMask=[true;true;false;true;true];
carried=calcTransition(pathRe,true(5,1),cfg,path,heatMask);
assert(isequal(carried.turbulentFace,logical([0;1;0;1;1])));
assert(isequal(carried.onsetFace,logical([0;1;0;0;0])));
assert(carried.turbulentCarrierFace(3));
assert(all(carried.inheritedThroughUnavailableFace(4:5)));
assert(isnan(carried.intermittency(3)));
hiddenOnset=[true;false;true;true;true];
carried=calcTransition(pathRe,true(5,1),cfg,path,hiddenOnset);
assert(isequal(carried.turbulentFace,logical([0;0;1;1;1])));
assert(~any(carried.onsetFace));
assert(carried.summary.onsetsOutsideHeatDomain==1);
% 中间无解面不能跨越传播湍流状态。
maskWithGap=[true;true;false;true;true];
broken=calcTransition(pathRe,maskWithGap,cfg,path);
assert(~broken.turbulentFace(4) && ~broken.turbulentFace(5));
failed=false;
try
    calcTransition(pathRe,maskWithGap,cfg,path,true(5,1));
catch ME
    failed=strcmp(ME.identifier,'RocketAeroHeat:Transition:Mask');
end
assert(failed);
bad=Re; bad(1)=NaN;
failed=false; try, calcTransition(bad,mask,cfg); catch, failed=true; end
assert(failed);
fprintf('Transition passed: Recrit=%.0f, threshold equality is turbulent.\n', ...
    t.criticalReynolds);
end
