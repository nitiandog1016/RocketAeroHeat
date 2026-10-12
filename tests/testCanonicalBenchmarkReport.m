function testCanonicalBenchmarkReport()
% 模块用途：完整链路报告必须区分给定/计算/未解，底面不得伪造热流。
addpath(fileparts(fileparts(mfilename('fullpath'))));
folder=tempname;
r=runCanonicalBenchmarks(folder,struct('kinds',{{'hemisphereCylinder'}}, ...
    'resolutions',48,'alphas',0,'exportWeb',false));
assert(r.shockTable.passed && ~r.experimentalValidation);
e=r.cases(1);names=fieldnames(e.fields);
for k=1:numel(names)
    f=e.fields.(names{k});
    assert(abs(f.computedFraction+f.prescribedFraction+f.missingFraction-1)<1e-12);
    regionArea=0;computedArea=0;
    for j=1:numel(e.regions)
        region=e.regions(j);regionArea=regionArea+region.area;
        computedArea=computedArea+region.area*region.fields.(names{k}).computedFraction;
    end
    assert(abs(regionArea/e.area-1)<1e-12);
    assert(abs(computedArea/e.area-f.computedFraction)<1e-12);
end
base=e.regions(strcmp({e.regions.name},'base'));
assert(base.fields.heatFlux.missingFraction==1);
assert(base.fields.pressure.prescribedFraction==1);
nose=e.regions(strcmp({e.regions.name},'sphericalNose'));
assert(nose.fields.heatFlux.computedFraction>0.8,'已验证半球头应进入现有低阶球鼻模型');
assert(e.fields.wallInput.prescribedFraction==1);
saved=jsondecode(fileread(fullfile(folder,'report.json')));
assert(saved.cases.faces==e.faces && exist(fullfile(folder,'report.md'),'file')==2);
fprintf('Canonical full-chain coverage report passed; output %s\n',folder);
end
