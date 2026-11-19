function results=testComputedWallCoverage()
%TESTCOMPUTEDWALLCOVERAGE 验证原面交集、缺项及有符号热功率，从普通计算入口回读结果。
% 一、三个面积不等的原面。给定压力排除，热流为零也算出；用分开掩码防止错配。
a=[1;2;4];
r=summarizeComputedWall(a,[false;true;true],[true;false;true],[0;NaN;-3]);
assert(r.jointComputedFaceCount==1 && r.jointComputedAreaM2==4 && r.jointMissingAreaM2==3);
assert(r.calculatedHeatPowerW==-12 && isnan(r.wholeSurfaceHeatPowerW) && ~r.fullSurfaceComputed);
disjoint=summarizeComputedWall(a,[true;false;false],[false;true;true],[17;2;-1]);
assert(disjoint.jointComputedFaceCount==0 && disjoint.jointMissingAreaM2==7 && disjoint.calculatedHeatPowerW==0);
assert(contains(disjoint.coverageMessage,'已算面热功率 0 W') && isnan(disjoint.wholeSurfaceHeatPowerW));
% 全壁热流都有解时，即使整壁缺少压力，仍可报告总热功率。
heatOnly=summarizeComputedWall(a,false(1,3),true(1,3),[0 2 -1]);
assert(heatOnly.wholeSurfaceHeatPowerW==0 && ~heatOnly.fullSurfaceComputed);
absent=summarizeComputedWall(a,true(3,1),false(3,1),nan(3,1));
assert(isnan(absent.calculatedHeatPowerW) && isnan(absent.wholeSurfaceHeatPowerW));
try
    summarizeComputedWall(a,true(3,1),true(3,1),[0;NaN;0]);
    error('Test:MissingRejection','应拒绝已计算面上的非有限热流。');
catch ex
    assert(strcmp(ex.identifier,'RocketAeroHeat:Coverage:HeatStatus'),ex.message);
end

% 二、标准球经过现有工程计算入口。按原面字段独立求交，并对照工作台的显示文字。
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=2;cfg.flight.alpha=0;cfg.flight.beta=0;
cfg.heating.airThermalModel='caloricallyPerfect';
mesh=buildCanonicalBody('sphere',struct('resolution',16));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
[current,~,run]=runAeroHeatSingleCase(mesh,cfg);
data=buildAerothermalViewerData(current);
p=current.pressure.fieldStatus(:)==1; q=current.heating.eligibleFace(:);area=current.localGeometry.area(:);
assert(nnz(p&q)==run.jointComputedFaceCount && ~run.fullSurfaceComputed && strcmp(run.status,'completed'));
assert(abs(sum(area(p&q))-run.jointComputedAreaM2)<1e-12 && abs(sum(area(~(p&q)))-run.jointMissingAreaM2)<1e-12);
assert(strcmp(run.coverageMessage,data.wallCoverage.coverageMessage) && contains(data.caseNote,run.coverageMessage));
% 三、同一普通入口进入混合模式，用已知整壁字段检查汇总交付，不执行新流场。
n=size(mesh.F,1);cfg.hybridSurface=struct('condition',struct('name','数据回传测试'), ...
    'region',repmat("test-wall",n,1),'wallTemperatureK',300*ones(n,1), ...
    'engineeringSolver',@knownFields,'localSolver',@unexpectedLocalSolve);
cfg.analysisMode='hybridSurface';
[model,~,hybrid]=runAeroHeatSingleCase(mesh,cfg);
view=buildAerothermalViewerData(model);
assert(hybrid.fullSurfaceComputed && hybrid.jointComputedFaceCount==n && hybrid.jointMissingAreaM2==0);
assert(abs(hybrid.wholeSurfaceHeatPowerW+2*sum(area))<1e-12 && hybrid.regionalAddedAreaM2==0);
assert(strcmp(view.wallCoverage.coverageMessage,hybrid.coverageMessage));
results=struct('partialWithZeroAndNegative',r,'disjoint',disjoint,'completeHeatOnly',heatOnly, ...
    'noHeat',absent,'standardSphereEntry',data.wallCoverage,'hybridInterfaceFixture',view.wallCoverage);
end

function rows=knownFields(g,~,~)
% 提供已知字段，仅测试几何、状态和统计量在常用混合入口的传递。
n=size(g.F,1);
rows=table((0:n-1)',ones(n,1),-2*ones(n,1),repmat("测试输入",n,1), ...
    'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
end

function result=unexpectedLocalSolve(varargin) %#ok<STOUT,INUSD>
error('Test:UnexpectedLocalSolve','已有完整测试字段，不应调用区域求解器。');
end
