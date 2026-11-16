function chain=marchPhaseAwareLaminarProfiles(edges,stations,Tw,seed,domain)
%MARCHPHASEAWARELAMINARPROFILES 保留真实剖面与壁温历史的附着层流沿程链。
% 模块一：核对物理入口和适用身份；不把调用者声明转成物理资格。
% 模块二：从保存入口逐站隐式推进，失败阻断下游，不重置厚度或剖面。
% 模块三：返回可供区域映射/独立审核的完整剖面与历史指纹。
required={'geometryKind','crossflowPresent','separationPresent','layerState', ...
    'sourceRole','sourceSHA256','angleOfAttackDeg','radiusM'};
assert(isstruct(domain)&&all(isfield(domain,required)),'RAH:PhaseMarch:Domain','缺少具名区域合同。');
assert(islogical(domain.crossflowPresent)&&~domain.crossflowPresent&& ...
    islogical(domain.separationPresent)&&~domain.separationPresent&&strcmp(domain.layerState,'laminar'), ...
    'RAH:PhaseMarch:Domain','本版不解强横流、逆流/分离或湍流。');
assert(ismember(domain.geometryKind,{'planar','axisymmetric','rotated-spherical-research'}) ...
    &&ismember(domain.sourceRole,{'prescribed-control','engineering-candidate'}) ...
    &&strlength(string(domain.sourceSHA256))==64,'RAH:PhaseMarch:Domain','不支持此几何或输入资格。');
if ~strcmp(domain.geometryKind,'rotated-spherical-research')
    assert(domain.angleOfAttackDeg==0,'RAH:PhaseMarch:Domain','不能把轴对称/平板范围套到有攻角箭体。');
end
N=numel(edges);assert(N==numel(stations)&&N>=2&&numel(Tw)==N&&all(isfinite(Tw)));
assert(all(diff([stations.xi])>0)&&stations(1).xi==seed.station.xi, ...
    'RAH:PhaseMarch:Order','入口或严格上游序列不一致。');
for key={'pressurePa','temperatureK','speedMps','totalEnthalpyJkg'}
    assert(abs(edges(1).(key{1})-seed.edge.(key{1}))<1e-12*max(1,abs(edges(1).(key{1}))), ...
        'RAH:PhaseMarch:Seed','保存入口与首站物理状态不一致。');
end
assert(seed.available&&seed.wallTemperature==Tw(1),'RAH:PhaseMarch:Seed','入口壁温/剖面无效。');
c0=phaseAwareAirCaloric(edges(1).temperatureK);S0=c0.s0-c0.R*log(edges(1).pressurePa);
for j=1:N
    c=phaseAwareAirCaloric(edges(j).temperatureK);
    assert(abs(c.s0-c.R*log(edges(j).pressurePa)-S0)<1e-5, ...
        'RAH:PhaseMarch:Entropy','本版同线外缘须同熵；不能把变熵层输入改标签。');
end
chain=struct('schema','rah-phase-nonsimilar-chain-v1','domain',domain,'edges',edges,'stations',stations, ...
    'wallTemperatureHistory',Tw(:)','profiles',{cell(1,N)},'available',false(1,N), ...
    'physicalQualification',false,'profileReadyForCFD',false,'sourceRole',domain.sourceRole, ...
    'method','TP-3271层流Levy–Lees非相似方程；xi后退一阶＋MATLAB bvp4c耦合五状态', ...
    'scope','同线常总焓/同熵、附着薄层；不解实际三维横流/分离/转捩；源资格不因求解而升级');
seed.streamwiseOrder=0;seed.energyThickness=seed.thicknessScale*trapz(seed.profile.eta,seed.profile.states(2,:).*(1-seed.profile.states(2,:).^2));
seed.enthalpyDefectThickness=NaN;chain.profiles{1}=seed;chain.available(1)=true;
for j=2:N
    r=solvePhaseAwareNonsimilarStation(edges(j),stations(j),Tw(j),chain.profiles{j-1});
    if r.available
        a=phaseAwareAirState('PT',edges(j).pressurePa,edges(j).temperatureK,'linearTable');
        r.Re_s=a.densityKgm3*edges(j).speedMps*stations(j).arcM/a.dynamicViscosityPas;
        if r.Re_s>=5e5||r.delta/domain.radiusM>=.05
            r.available=false;r.reason='Re_s或薄层delta/R超过原工程筛选';
        end
    end
    chain.profiles{j}=r;chain.available(j)=r.available;
    if ~r.available,break;end
end
chain.sourceText=jsonencode(struct('domain',domain,'codeSHA256',fileSHA('solvePhaseAwareNonsimilarStation.m'), ...
    'marchSHA256',fileSHA('marchPhaseAwareLaminarProfiles.m'),'seed',seed.numericalProfile, ...
    'edge',edges,'station',stations,'wallHistory',Tw));
chain.sourceSHA256=phaseAwarePayloadSHA(chain.sourceText,jsonencode(chain.profiles));
end
function h=fileSHA(name)
fid=fopen(fullfile(fileparts(mfilename('fullpath')),name),'rb');assert(fid>=0);c=onCleanup(@()fclose(fid));
b=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(b,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
