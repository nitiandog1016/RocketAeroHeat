function summary=validateSphericalNoseWarmRegion(outputFolder,baselineFile,reuseDefault)
%VALIDATESPHERICALNOSEWARMREGION 默认2度/细1度完整球鼻有效区校核。
% 模块作用：保存各默认/细表的真实算例与性能，复核每个有效热根、
% 插值面元和来源，按相同网格面积记录v53后的增量。不是完整后体解。
if nargin<1,outputFolder='outputs/v54_validation';end
if nargin<2,baselineFile='outputs/v53_validation/sphere_laminar_results.mat';end
if nargin<3,reuseDefault=false;end
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<2 && ~isfile(baselineFile)
    % 完整工程包把旧版本证据放在validation目录，仍使用原始已测基线。
    baselineFile=fullfile(root,'validation','inherited_v53','v53','sphere_laminar_results.mat');
end
if ~exist(outputFolder,'dir'),mkdir(outputFolder);end
base=load(baselineFile,'cfg','fs','mesh','meta','nasa');
cfg=base.cfg;fs=base.fs;mesh=base.mesh;meta=base.meta;
cfg.heating=rmfield(cfg.heating,'sphericalNoseAngleStepDeg');
eligible=base.nasa.edgeState.algebraicallyValid;
summary=struct('tested',false,'experimentalValidation',false,'interfaceSolved',false, ...
    'baselineFile',baselineFile,'baselineComputedFaces',nnz(base.nasa.faceMask), ...
    'baselineComputedAreaM2',sum(mesh.area(base.nasa.faceMask)), ...
    'baselineEligibleFaces',nnz(eligible),'baselineEligibleAreaM2',sum(mesh.area(eligible)), ...
    'totalMeshAreaM2',sum(mesh.area),'totalMeshFaces',size(mesh.F,1));

%% 1. 实际默认2度角度表；沿用已核验几何，不先预设覆盖结果
if reuseDefault
    % 验证任务被中断后，只在几何/工况/配置完全相同的情况下复用已保存
    % 的默认表；细表必须实际计算。不能以中断的日志声称收敛已通过。
    cached=load(fullfile(outputFolder,'sphere_laminar_default2.mat'));
    assert(isequaln(cached.cfg,cfg) && isequaln(cached.fs,fs) && ...
        isequaln(cached.mesh,mesh),'默认表缓存与本次算例不一致，不能复用');
    result=cached.result;seconds2=cached.summary.default2.seconds;
else
    tic;result=solveSphericalNoseLaminarHeating(mesh,fs,300,cfg);seconds2=toc;
end
fprintf('Default observed: angle %.9g, eligible %d / computed %d / calls %d / stations %d / %.6g s.\n', ...
    result.angularStepDeg,nnz(result.edgeState.algebraicallyValid),nnz(result.faceMask), ...
    result.profileCalls,result.uniqueStations,seconds2);
assert(result.angularStepDeg==2,'必须实测未配置角度步长时的生产默认');
summary.default2=checkRegion(result,mesh,seconds2);
assert(isequal(result.edgeState.algebraicallyValid,eligible), ...
    '比较不得通过改变几何/外缘范围扩大有效区');
save(fullfile(outputFolder,'sphere_laminar_default2.mat'),'cfg','fs','mesh','meta','result','summary');

%% 2. 高角度和所有同一真面元上的1度细表；不能跨无解节点比较
fineCfg=cfg;fineCfg.heating.sphericalNoseAngleStepDeg=1;
tic;fine=solveSphericalNoseLaminarHeating(mesh,fs,300,fineCfg);seconds1=toc;
fprintf('Fine observed: angle %.9g, eligible %d / computed %d / calls %d / stations %d / %.6g s.\n', ...
    fine.angularStepDeg,nnz(fine.edgeState.algebraicallyValid),nnz(fine.faceMask), ...
    fine.profileCalls,fine.uniqueStations,seconds1);
summary.fine1=checkRegion(fine,mesh,seconds1);
assert(isequal(result.edgeState.algebraicallyValid,fine.edgeState.algebraicallyValid));
common=result.faceMask & fine.faceMask;
high=common & result.polarAngleDeg>=45;
assert(nnz(high)>0 && all(common(eligible)),'不能用较小交集冒充同一完整有效区收敛');
summary.allAreaConvergence=compare(result,fine,common);
summary.highAngleConvergence=compare(result,fine,high);
% 误差标准沿用5/2.5度既有工程表校核范围；2/1应有更严格实际结果。
assert(summary.allAreaConvergence.heatRelativeDifference<.003 && ...
    summary.allAreaConvergence.recoveryDifferenceK<.15, ...
    '默认2度/细1度NASA7整个有效区未达到既有角度表收敛标准');
assert(summary.highAngleConvergence.heatRelativeDifference<.003 && ...
    summary.highAngleConvergence.recoveryDifferenceK<.15);
summary.additionalComputedFaces=nnz(result.faceMask & ~base.nasa.faceMask);
summary.additionalComputedAreaM2=sum(mesh.area(result.faceMask & ~base.nasa.faceMask));
summary.additionalTotalAreaFraction=summary.additionalComputedAreaM2/sum(mesh.area);
summary.baselineTotalAreaFraction=summary.baselineComputedAreaM2/sum(mesh.area);
summary.newTotalAreaFraction=summary.default2.computedAreaM2/sum(mesh.area);
summary.tested=true;
fprintf('Area increment: %d faces, %.12g m2, model total %.9g -> %.9g%%; eligible %.9g%%.\n', ...
    summary.additionalComputedFaces,summary.additionalComputedAreaM2, ...
    100*summary.baselineTotalAreaFraction,100*summary.newTotalAreaFraction, ...
    100*summary.default2.eligibleAreaFraction);
fprintf('Default2/fine1 all: q %.9g relative / Taw %.9g K; high45+: q %.9g / Taw %.9g K.\n', ...
    summary.allAreaConvergence.heatRelativeDifference,summary.allAreaConvergence.recoveryDifferenceK, ...
    summary.highAngleConvergence.heatRelativeDifference,summary.highAngleConvergence.recoveryDifferenceK);
save(fullfile(outputFolder,'sphere_laminar_default2.mat'),'cfg','fs','mesh','meta','result','summary');
save(fullfile(outputFolder,'sphere_laminar_angle_convergence.mat'),'summary','fineCfg','fine');
fid=fopen(fullfile(outputFolder,'sphere_laminar_angle_convergence.json'),'w');
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true));clear guard;
end
function metrics=checkRegion(r,mesh,seconds)
% 模块作用：统计实得区域再断言，保持未解区及算法和温度种类真实。
eligible=r.edgeState.algebraicallyValid;valid=r.faceMask;
metrics=struct('seconds',seconds,'angularStepDeg',r.angularStepDeg, ...
    'computedFaces',nnz(valid),'eligibleFaces',nnz(eligible), ...
    'computedAreaM2',sum(mesh.area(valid)),'eligibleAreaM2',sum(mesh.area(eligible)), ...
    'eligibleAreaFraction',sum(mesh.area(valid))/sum(mesh.area(eligible)), ...
    'profileCalls',r.profileCalls,'uniqueStations',r.uniqueStations, ...
    'maximumAdiabaticResidualWm2',0,'maximumHighAngleDeg',max(r.polarAngleDeg(valid)));
assert(r.available && all(valid(eligible)) && all(r.angularProfileTable.statusCode==0), ...
    '当前受支持球鼻有效区仍存在未解热根/面元，不能宣称该区全部求解');
assert(~r.interfaceSolved && ~r.entropySwallowingSolved && ...
    ~r.nonSimilarTermsSolved && ~r.experimentalValidation && ~r.historyCompatible);
assert(strcmp(r.fieldMetadata.wallTemperatureInput.valueKind,'prescribed'));
for key={'heatFlux','recoveryTemperature','heatTransferCoefficient'}
    m=r.fieldMetadata.(key{1});assert(strcmp(m.valueKind,'computed') && ...
        ~isempty(m.method) && ~isempty(m.scope));
    value=r.(key{1});assert(all(isfinite(value(valid))) && all(isnan(value(~valid))));
end
assert(all(r.heatFlux(valid)>0) && all(r.heatTransferCoefficient(valid)>0));
for j=1:numel(r.stationResults)
    n=r.stationResults{j};assert(n.wallProfile.available && n.adiabaticProfile.available);
    ad=n.adiabaticProfile;
    assert(ad.wallShearStress>=0 && all(ad.profile.speed>=-1e-8));
    assert(abs(ad.heatFlux)<=1e-5*max(abs(n.wallProfile.heatFlux),1));
    metrics.maximumAdiabaticResidualWm2=max(metrics.maximumAdiabaticResidualWm2,abs(ad.heatFlux));
    if n.polarAngleDeg>0
        t=n.adiabaticBracketTrials;ends=n.adiabaticBracket;
        for e=1:2
            row=find(abs(t(:,1)-ends(e))<1e-8,1,'last');
            assert(~isempty(row) && t(row,2)==1);
        end
    end
end
for id=find(valid)'
    pair=r.profileNodeIDs(id,:);w=r.profileNodeWeight(id);
    assert(all(r.angularProfileTable.statusCode(pair)==0) && r.statusCode(id)==0);
    for key={'heatFlux','recoveryTemperature','wallShearStress','massThickness','delta','xi','beta'}
        col=find(strcmp(r.angularProfileTable.fieldNames,key{1}));
        expected=(1-w)*r.angularProfileTable.values(pair(1),col)+ ...
            w*r.angularProfileTable.values(pair(2),col);
        assert(abs(r.(key{1})(id)-expected)<=1e-10*max(abs(expected),1));
    end
end
end
function m=compare(a,b,mask)
m=struct('faces',nnz(mask), ...
    'heatRelativeDifference',max(abs(a.heatFlux(mask)./b.heatFlux(mask)-1)), ...
    'recoveryDifferenceK',max(abs(a.recoveryTemperature(mask)-b.recoveryTemperature(mask))), ...
    'heatTransferRelativeDifference',max(abs(a.heatTransferCoefficient(mask)./b.heatTransferCoefficient(mask)-1)), ...
    'massThicknessRelativeDifference',max(abs(a.massThickness(mask)./b.massThickness(mask)-1)), ...
    'deltaRelativeDifference',max(abs(a.delta(mask)./b.delta(mask)-1)));
end
