function data=buildPhaseAwareNoseHeatingViewerData(mesh)
%BUILDPHASEAWARENOSEHEATINGVIEWERDATA 同一原面压力/热流/剖面积分及方法显示。
% 模块用途：复用原工作台字段合同；给定Tw、未知固体Tw、恢复温度估计分开。
data=buildPhaseAwareNoseViewerData(mesh);h=mesh.phaseAwareNoseHeating;
assert(strcmp(h.edgeSourceSHA256,mesh.phaseAwareNoseEdge.source.sha256) ...
    &&strcmp(h.meshSHA256,mesh.phaseAwareNoseEdge.meshSHA256), ...
    'RAH:PhaseNoseHeatViewer:Source','热流与当前外缘来源不同。');
digest=phaseAwarePayloadSHA(h.sourceText,h.heatFlux,h.wallShearStress,h.delta,h.deltaStar, ...
    h.thetaPhysical,h.massThickness,h.skinFriction,h.validFace,h.wallTemperature, ...
    h.lowerStation,h.upperStation,h.upperWeight);
assert(strcmp(digest,h.sourceSHA256)&&isequal(mesh.phaseAwareWallTemperature,h.wallTemperature), ...
    'RAH:PhaseNoseHeatViewer:Payload','热流或壁温已改变。');
data.defaultField='heatFlux';data.caseNote=['单工况鼻部热候选：' h.scope];
data.missingReason=h.reason;data.region=repmat({'其他待求解区域'},size(mesh.F,1),1);
data.region(h.validFace)={'相态球鼻热工程候选'};
data.regions=struct('name',{'相态球鼻热工程候选','其他待求解区域'}, ...
    'faces',{nnz(h.validFace),nnz(~h.validFace)}, ...
    'area',{sum(data.faceArea(h.validFace)),sum(data.faceArea(~h.validFace))});
keys={'heatFlux','wallShearStress','delta','deltaStar','thetaPhysical','massThickness','skinFriction'};
labels={'入壁对流热流','壁面剪应力','有限边缘厚度','密度加权位移厚度','密度加权动量厚度','质量厚度','局部摩擦系数'};
units={'W/m²','Pa','m','m','m','m','1'};
for k=1:numel(keys)
    f=data.fields(7);f.key=keys{k};f.label=labels{k};f.unit=units{k};f.values=h.(keys{k});
    f.status=double(h.validFace&isfinite(f.values));f.methodNames={h.method};f.scope=h.scope;
    f.missingReason=h.reason;f.computedStateLabel='未验证鼻部热工程候选';
    if k==1,data.fields(7)=f;else,data.fields(end+1)=f;end
end
data.researchIdentity.heatSourceSHA256=h.sourceSHA256;data.researchIdentity.heatSummary=h.summary;
data.heatingGaps=struct('reason','无适用相态球鼻热剖面，其他区域仍待求解','faces',nnz(~h.validFace), ...
    'area',sum(data.faceArea(~h.validFace)),'areaFraction',sum(data.faceArea(~h.validFace))/data.totalArea);
data.profileMapping=struct('lowerStation',h.lowerStation,'upperStation',h.upperStation, ...
    'upperWeight',h.upperWeight,'anglesRad',h.anglesRad,'note','相邻节点工程插值；非三维沿程剖面推进');
end
