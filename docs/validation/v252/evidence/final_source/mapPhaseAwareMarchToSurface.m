function h=mapPhaseAwareMarchToSurface(mesh,e,chain)
%MAPPHASEAWAREMARCHTOSURFACE 已有真实球帽原面上的非相似剖面候选映射。
% 模块一核对原面及连续外缘来源；模块二仅在相邻成功站间插值。
% 全部原面继续计入面积分母，鼻外不补零，输入和物理资格不升级。
assert(strcmp(e.meshSHA256,phaseAwareMeshIdentity(mesh))&&strcmp(chain.domain.sourceSHA256,e.source.sha256), ...
    'RAH:PhaseMarch:Mapping','几何或外缘身份不符。');
digest=phaseAwarePayloadSHA(e.sourceText,e.pressurePa,e.temperatureK,e.speedMps,e.velocityMps, ...
    e.recoveryTemperatureK,e.validFace,e.reasonCode,e.arcLengthM,e.Mach,e.faceAreaM2);
assert(strcmp(digest,e.source.sha256)&&strcmp(chain.sourceSHA256,phaseAwarePayloadSHA(chain.sourceText,jsonencode(chain.profiles))), ...
    'RAH:PhaseMarch:Mapping','来源数据已被修改。');
% 来源文本还须与实际映射所用元数据相符，避免只验剖面而漏验角站/壁温/可用标记。
source=jsondecode(chain.sourceText);
assert(strcmp(jsonencode(source.domain),jsonencode(chain.domain))&& ...
    strcmp(jsonencode(source.station),jsonencode(chain.stations))&& ...
    strcmp(jsonencode(source.edge),jsonencode(chain.edges))&& ...
    strcmp(jsonencode(source.wallHistory),jsonencode(chain.wallTemperatureHistory))&& ...
    strcmp(chain.domain.geometryKind,'rotated-spherical-research'), ...
    'RAH:PhaseMarch:Mapping','实际区域/角站/外缘/壁温与来源文本不符。');
for k=1:numel(chain.available)
    p=chain.profiles{k};expected=isstruct(p)&&isfield(p,'available')&&p.available;
    assert(chain.available(k)==expected,'RAH:PhaseMarch:Mapping','可用标记与实际剖面不符。');
end
n=size(mesh.F,1);angles=[chain.stations.thetaRad];R=chain.domain.radiusM;theta=e.arcLengthM/R;
valid=logical(e.validFace(:));assert(~any(mesh.localGeometry.featureFace(valid)|mesh.localGeometry.topologyIssue(valid)));
h=struct('schema','rah-phase-nonsimilar-surface-v1','validFace',false(n,1), ...
    'reason',{repmat({'鼻外/强三维/无适用真实外缘，保持未解'},n,1)},'lowerStation',zeros(n,1), ...
    'upperStation',zeros(n,1),'upperWeight',nan(n,1),'meshSHA256',e.meshSHA256, ...
    'edgeSourceSHA256',e.source.sha256,'chainSourceSHA256',chain.sourceSHA256,'anglesRad',angles, ...
    'method',chain.method,'scope',[chain.scope '；此真实20度箭体仅使用已具名旋转球帽近似'], ...
    'physicalQualification',false,'newQualifiedAreaM2',0,'wallTemperature',nan(n,1));
keys={'heatFlux','wallShearStress','delta','deltaStar','thetaPhysical','massThickness','energyThickness'};
for key=keys,h.(key{1})=nan(n,1);end
h.skinFriction=nan(n,1);
for id=find(valid)'
    j=find(angles<=theta(id),1,'last');if isempty(j),continue;end
    if j==numel(angles),j=j-1;end
    w=(theta(id)-angles(j))/(angles(j+1)-angles(j));assert(w>=-1e-12&&w<=1+1e-12);
    h.lowerStation(id)=j;h.upperStation(id)=j+1;h.upperWeight(id)=w;
    if ~all(chain.available(j:j+1)),h.reason{id}='上游非相似剖面拒绝，整条下游保持未知';continue;end
    a=chain.profiles{j};b=chain.profiles{j+1};
    for key=keys,h.(key{1})(id)=(1-w)*a.(key{1})+w*b.(key{1});end
    h.wallTemperature(id)=(1-w)*chain.wallTemperatureHistory(j)+w*chain.wallTemperatureHistory(j+1);
    c=phaseAwareAirCaloric(e.temperatureK(id));rho=e.pressurePa(id)/(c.R*e.temperatureK(id));
    h.skinFriction(id)=h.wallShearStress(id)/(.5*rho*e.speedMps(id)^2);
    h.validFace(id)=true;h.reason{id}='同相态非相似剖面工程候选；源/三维/实验资格未闭合';
end
h.summary=struct('originalFaces',n,'totalAreaM2',sum(e.faceAreaM2),'candidateFaces',nnz(h.validFace), ...
    'candidateAreaM2',sum(e.faceAreaM2(h.validFace)),'signedPowerW',sum(e.faceAreaM2(h.validFace).*h.heatFlux(h.validFace)), ...
    'newQualifiedAreaM2',0);
end
