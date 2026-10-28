function path=exportNonsimilarSphericalStudy(mesh,mapped,folder)
%EXPORTNONSIMILARSPHERICALSTUDY 用既有WebGL/探针/截面查看条件研究映射。
% 本页与生产结果分开导出，不改main默认链。条件q明确标名，不计生产覆盖；
% 恢复温度和壁温解未求时保持NaN，外缘给定模型与求导热流分别说明。
if ~isstruct(mapped)||~isfield(mapped,'kind')|| ...
        ~strcmp(mapped.kind,'conditional-spherical-face-map')|| ...
        ~isequal(mapped.physicalSolution,false)||~isequal(mapped.eligibleForSurfaceCoverage,false)|| ...
        ~isequal(mesh.V,mapped.geometryVertices)||~isequal(mesh.F,mapped.geometryFaces)
    error('RocketAeroHeat:NonSimilarView:Binding','研究映射与当前真实V/F不匹配。');
end
nf=size(mesh.F,1);mask=mapped.faceMask;
if ~islogical(mask)||~isequal(size(mask),[nf 1])
    error('RocketAeroHeat:NonSimilarView:Mask','条件映射掩码必须逐面元对应。');
end
for key={'heatFlux','wallShearStress','edgePressure','edgeTemperature'}
    values=mapped.(key{1});
    if ~isnumeric(values)||~isreal(values)||~isequal(size(values),[nf 1])|| ...
            any(~isfinite(values(mask)))||any(~isnan(values(~mask)))
        error('RocketAeroHeat:NonSimilarView:Values','条件字段与真实数值掩码不一致。');
    end
end
if ~isnumeric(mapped.wallTemperatureInput)||~isreal(mapped.wallTemperatureInput)|| ...
        ~isequal(size(mapped.wallTemperatureInput),[nf 1])|| ...
        any(~isfinite(mapped.wallTemperatureInput))||any(mapped.wallTemperatureInput<=0)|| ...
        any(~isnan(mapped.wallTemperatureSolved))||any(~isnan(mapped.recoveryTemperature))
    error('RocketAeroHeat:NonSimilarView:TemperatureIdentity','给定壁温与未求温度字段身份不一致。');
end
data=buildAerothermalViewerData(mesh);template=data.fields(1);data.fields=struct([]);
data.researchIdentity=struct('physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'experimentalValidation',false,'newProductionSolvedAreaM2',0, ...
    'inletAngleDeg',mapped.inletAngleDeg,'acceptedIntervalDeg',mapped.acceptedIntervalDeg, ...
    'source',mapped.source,'currentSourceHashes',mapped.currentSourceHashes, ...
    'historicalPreRunSourceBinding',false);
data.caseNote=sprintf(['条件研究，非生产/实验解：M%.4g，高度%.4g km，攻角0°，' ...
    '整段Tw=%.6g K；有限角局部相似入口%.4g°；仅映射%.4g–%.4g°；新增生产面积0。'], ...
    mapped.condition.Mach,mapped.condition.altitude/1000,mapped.wallTemperatureInput(1), ...
    mapped.inletAngleDeg,mapped.acceptedIntervalDeg);
data.resultCondition=struct('caption',data.caseNote);data.flowDirection=[-1 0 0];
data.missingReason=mapped.missingReason;data.heatingGaps=struct([]);
keep=mapped.faceMask;nf=size(mesh.F,1);area=data.faceArea;
add('conditionalHeatFlux','条件研究热流（不计生产覆盖）','kW/m²',mapped.heatFlux/1000,double(keep), ...
    mapped.method,[mapped.scope '；' mapped.interpolation ...
    '；面元常值采样：同一X截面可能穿过质心极角不同的三角面，局部锯齿不是横流物理预测']);
add('conditionalWallShear','条件研究壁面剪切','Pa',mapped.wallShearStress,double(keep), ...
    '同一非相似u剖面壁面导数，Sutherland黏度',mapped.scope);
add('conditionalEdgePressure','条件外缘压力','kPa',mapped.edgePressure/1000,double(keep), ...
    '给定固定熵工程外缘：球面修正牛顿p及正激波H0/熵', ...
    '工程外缘模型；不是求得真实激波层或表面载荷压力');
add('conditionalEdgeTemperature','条件外缘静温','K',mapped.edgeTemperature,double(keep), ...
    '同一NASA7正激波驻点H0/熵恢复',mapped.scope);
add('wallInput','整段壁温边界输入','K',mapped.wallTemperatureInput,2*ones(nf,1), ...
    '整段给定均匀Tw历史，不是求解壁温','仅适用于与保存条件试验完全一致的Tw历史');
add('wallSolved','求解壁温','K',nan(nf,1),zeros(nf,1), ...
    '未执行独立壁温求解','未提供材料/热响应或独立壁温解');
add('conditionalRecoveryTemperature','条件研究恢复温度','K',nan(nf,1),zeros(nf,1), ...
    '未求解绝热边界根，不能借用总温或局部相似结果','本非相似试验仅求给定Tw下q与tau');
% 分区字段须对应本页实际研究字段，不能残留模板中的旧wallSolved汇总。
for j=1:numel(data.regions)
    mask=strcmp(data.region,data.regions(j).name);data.regions(j).fields=struct;
    for k=1:numel(data.fields)
        field=data.fields(k);state=field.status;
        data.regions(j).fields.(field.key)=struct( ...
            'computedArea',sum(area(mask & state==1)), ...
            'prescribedArea',sum(area(mask & state==2)), ...
            'missingArea',sum(area(mask & state==0)));
    end
end
if ~exist(folder,'dir'),mkdir(folder);end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
for name={'index.html','viewer.css','section.js','export.js','viewer.js'}
    copyfile(fullfile(root,'web',name{1}),fullfile(folder,name{1}));
end
fid=fopen(fullfile(folder,'case-data.js'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'window.ROCKET_CASE=%s;\n',jsonencode(data,'ConvertInfAndNaN',true));clear cleanup;
fid=fopen(fullfile(folder,'results.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(data,'ConvertInfAndNaN',true));
path=fullfile(folder,'index.html');

    function add(key,label,unit,values,status,method,scope)
        f=template;f.key=key;f.label=label;f.unit=unit;f.values=values;f.status=status;
        f.methodNames={method};f.methodIndex=zeros(nf,1);f.scope=scope;
        f.scopeNames={};f.scopeIndex=[];f.history=[];f.availability='conditionalResearch';
        f.stateNote='本页为条件研究；计算面积仅对应数值映射，不计生产覆盖，不代表实验验证。';
        f.missingReason=mapped.missingReason;f.missingReason(status>0)={''};
        if all(status==0),f.missingReason(:)={method};end
        f.computedAreaFraction=sum(area(status==1))/sum(area);
        f.prescribedAreaFraction=sum(area(status==2))/sum(area);
        f.missingAreaFraction=sum(area(status==0))/sum(area);
        if isempty(data.fields),data.fields=f;else,data.fields(end+1)=f;end
    end
end
