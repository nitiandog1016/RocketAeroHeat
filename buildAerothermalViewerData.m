function data = buildAerothermalViewerData(mesh)
%BUILDAEROTHERMALVIEWERDATA 将真实计算结果转换成离线三维工作台数据。
if isfield(mesh,'aerothermalResultMode')&&strcmp(mesh.aerothermalResultMode,'phaseAwareNoseHeating')
    data=buildPhaseAwareNoseHeatingViewerData(mesh);return;
end
if isfield(mesh,'aerothermalResultMode')&&strcmp(mesh.aerothermalResultMode,'phaseAwareNoseEdge')
    data=buildPhaseAwareNoseViewerData(mesh);return;
end
% 新相态区域优先使用自身来源，避免混入mesh中保留的旧NASA7/压力快照。
if isfield(mesh,'heating') && isfield(mesh.heating,'sidewallModel') && ...
        strcmp(mesh.heating.sidewallModel,'phaseAwareLaminarIntegral')
    data=buildPhaseAwareLaminarViewerData(mesh);
    return;
end
% 模块一：验证几何，统一面元编号；不会插值或修改计算结果。
V=double(mesh.V); F=double(mesh.F); n=size(F,1);
validateattributes(V,{'double'},{'finite','ncols',3,'nonempty'});
validateattributes(F,{'double'},{'finite','integer','positive','ncols',3,'nonempty','<=',size(V,1)});
a=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
area=sqrt(sum(a.^2,2))/2;
assert(all(area>0),'RocketAeroHeat:Viewer:Area','模型不能包含面积为零的退化面元。');
data=struct('kind','aerothermal-model','schemaVersion',1, ...
    'vertices',V,'faces',F-1,'faceArea',area,'fields',struct([]), ...
    'coordinateNote','SI · 鼻尖 +X · 零攻角来流 -X · 真实几何比例', ...
    'caseNote','单点工程结果与壁体热响应历史分别列示；缺少逐时刻压力时不推断瞬态压力');
data.region=repmat({'未分类'},n,1);
if isfield(mesh,'localGeometry') && isfield(mesh.localGeometry,'region')
    data.region=cellstr(string(mesh.localGeometry.region(:)));
end
data.missingReason=repmat({'尚无热流计算结果'},n,1);
data.heatingGaps=struct([]);
if isfield(mesh,'heatingCoverage')
    c=mesh.heatingCoverage;
    data.missingReason=reshape(c.reasonNames(double(c.reasonCode)+1),[],1);
    codes=unique(double(c.reasonCode(:)));codes=codes(codes>0);
    gaps=struct([]);
    for j=1:numel(codes)
        mask=double(c.reasonCode(:))==codes(j);
        gaps(j).reason=c.reasonNames{codes(j)+1};
        gaps(j).faces=nnz(mask);gaps(j).area=sum(area(mask));
        gaps(j).areaFraction=sum(area(mask))/sum(area);
    end
    data.heatingGaps=gaps;
end
% 模块二：逐变量记录方法、状态和适用范围。状态 0 未解、1 工程计算、2 给定。
if isfield(mesh,'pressure')
    p=mesh.pressure; state=ones(n,1); methods=repmat({p.model},n,1);
    if isfield(p,'freestream')
        fs=p.freestream;
        data.flowDirection=fs.unitFlowDirection;
        data.caseNote=sprintf(['单点工况：M = %.4g，攻角 %.4g°，侧滑角 %.4g°。' ...
            '输入壁温与热响应历史分别列示；历史热流不替代单点热流。'], ...
            fs.Mach,fs.alpha,fs.beta);
    end
    for name={'isPrescribedBase','isPrescribedLeeward'}
        if isfield(p,name{1})
            mask=logical(p.(name{1})(:)); state(mask)=2;
            if strcmp(name{1},'isPrescribedBase'),label='尾底';else,label='背风';end
            methods(mask)={['给定 ' label ' Cp 边界条件']};
        end
    end
    pressureScopes={};
    if isfield(p,'fieldStatus'),state=double(p.fieldStatus(:));end
    if isfield(p,'faceMethod'),methods=p.faceMethod(:);end
    if isfield(p,'faceScope'),pressureScopes=p.faceScope(:);end
    add('pressure','表面压力','kPa',p.staticPressure/1000,state,methods, ...
        '逐面元来源：修正牛顿、给定边界或经审核的区域关联；缺失保留未解',pressureScopes);
    add('cp','压力系数 Cp','1',p.Cp,state,methods,'与表面压力采用相同来源；Cp=(p-pinf)/qinf',pressureScopes);
end
expansionOutside=[];
% 混合入口已按原面完成求解，显示同一份有符号场和方法来源。
if isfield(mesh,'hybridWall')
    w=mesh.hybridWall;
    hybridScope='当前混合计算原面结果；方法见逐面来源，实验精度单独验证';
    if isfield(mesh,'thermistorSidewind'),hybridScope=mesh.evidenceIdentity.scope;end
    if isfield(mesh,'thermistorLead') && isfield(mesh.thermistorLead,'thermal')
        t=mesh.thermistorLead.thermal;
        hybridScope=sprintf('%s 本次壁温 %.12g K，参照终帧 %s。',t.scope,t.wallTemperatureK,t.nativeFrameSHA256);
    end
    data.missingReason(isfinite(w.heatFluxWm2))={''};
    add('pressure','壁面压力','kPa',w.pressurePa/1000,double(isfinite(w.pressurePa)),hybridMethod(w,'pressureMethod'),hybridScope);
    add('heatFlux','入壁热流','kW/m²',w.heatFluxWm2/1000,double(isfinite(w.heatFluxWm2)),hybridMethod(w,'heatFluxMethod'),hybridScope);
    if ismember('pressureGaugePa',w.Properties.VariableNames)
        delta=w.pressureGaugePa;
        sourceKey='pressureMethod';
        if ismember('pressureGaugeMethod',w.Properties.VariableNames),sourceKey='pressureGaugeMethod';end
        add('pressureGauge','相对环境壁压','Pa',delta,double(isfinite(delta)),hybridMethod(w,sourceKey), ...
            '同工况壁面绝对压力与环境参考压力之差');
    end
    if ismember('wallFluidTemperatureK',w.Properties.VariableNames)
        fluidT=w.wallFluidTemperatureK;
        add('wallFluidTemperature','壁节点流体温度','K',fluidT,double(isfinite(fluidT)), ...
            hybridMethod(w,'wallFluidTemperatureMethod'),'求解器实际输出的壁节点流体温度');
    end
    if ismember('nearWallCellTemperatureK',w.Properties.VariableNames)
        tCell=w.nearWallCellTemperatureK;
        firstScope='与原面对应的首层流体单元温度';
        if isfield(mesh,'thermistorLead'),firstScope=[firstScope '。' hybridScope];end
        add('nearWallCellTemperature','临壁流体单元温度','K',tCell,double(isfinite(tCell)), ...
            hybridMethod(w,'nearWallTemperatureMethod'),firstScope);
    end
    if isfield(mesh,'hybridSolvedWallTemperatureK')
        reason='此工况输入净对流功率，壁温由稳态热平衡求出。';
        add('wallInput','给定壁温','K',nan(n,1),zeros(n,1),repmat({'未给定壁温'},n,1),reason,{},repmat({reason},n,1));
        data.fields(end).stateNote=reason;
    else
        wallMethod='输入壁温';
        if isfield(mesh,'thermistorWire'),wallMethod='珠部和接头等温，直引线按一维长线导热解逐节点给定，交叠积分至原面';end
        if isfield(mesh,'thermometryFluxChange'),wallMethod='由电阻等效内部温度反推珠表面温度，长引线超温按同一比例作为空气侧边界';end
        add('wallInput','给定壁温','K',mesh.hybridWallTemperatureK,2*ones(n,1),repmat({wallMethod},n,1),hybridScope);
    end
    % 混合结果使用适配器实际返回的算法，不沿用默认工程压力、热流的解释文字。
    physical={'pressure','pressureGauge','heatFlux','wallFluidTemperature','nearWallCellTemperature'};
    for k=find(ismember({data.fields.key},physical))
        data.fields(k).stateNote=[data.fields(k).label '取自当前工况的逐面数值结果，算法见面元探针。'];
        if isfield(mesh,'evidenceIdentity') && isfield(mesh.evidenceIdentity,'origin') && strcmp(mesh.evidenceIdentity.origin,'research')
            data.fields(k).computedStateLabel='研究数值解';
        end
    end
    % 热流视图同时给出当前终场的原生边界账，保留全部原面字段及其数值。
    if isfield(mesh,'hybridEnergyBalance')
        b=mesh.hybridEnergyBalance;
        data.nativeEnergyBalance=b;
        k=find(strcmp({data.fields.key},'heatFlux'),1);
        if isfield(b,'schema') && strcmp(b.schema,'cm011-inclined-half-energy-v1')
            data.fields(k).stateNote=sprintf(['%s 三维半域折合整壁放热 %.3f µW，远场流出 %.3f µW。' ...
                '壁面减远场为 %.3g W，占壁面放热 %.3g%%。该工况的实验精度待验证。'],data.fields(k).stateNote, ...
                1e6*b.nativeWallReleaseW,1e6*b.nativeFarExcessOutW,b.wallMinusFarW,100*b.wallMinusFarFraction);
        else
            data.fields(k).stateNote=sprintf(['%s 壁面放热 %.3f µW，远场显热流出 %.3f µW。' ...
            '两者相差 %.3f µW（壁面放热的 %.2f%%）。' ...
            '离散方程中，内部面与轴对称源项合计 %.3f µW。该试件的实验精度待验证。'], ...
            data.fields(k).stateNote,1e6*b.nativeWallReleaseW,1e6*b.nativeFarExcessOutW, ...
            1e6*b.nativeSurfaceImbalanceW,100*b.nativeSurfaceImbalanceW/b.nativeWallReleaseW, ...
            1e6*b.innerFacePlusAxisSourceW);
        end
    end
end
% 热敏珠的差值场与主场共用原面，探针、截面与导出均能追溯到两种周向网格。
if isfield(mesh,'beadAngularDelta')
    angularDelta=mesh.beadAngularDelta;
    method=repmat({'同一冻结工况逐面相减，7.5° 网减 15° 网'},n,1);
    scope='周向离散敏感性，正值表示细网计算值较大。';
    add('angularHeatDifference','入壁热流差 · 7.5° − 15°','W/m²',angularDelta.heatFluxWm2,ones(n,1),method,scope);
    data.fields(end).stateNote=scope;data.fields(end).computedStateLabel='网格对照';
    add('angularPressureDifference','壁压差 · 7.5° − 15°','Pa',angularDelta.pressurePa,ones(n,1),method,scope);
    data.fields(end).stateNote=scope;data.fields(end).computedStateLabel='网格对照';
end
% 同一低速试件的三维体网格对照。研究帧沿用自己的残差、资格和功率，单列对照面积。
if isfield(mesh,'cm011Sector')
    sector=mesh.cm011Sector.summary; w=mesh.cm011Sector.wall;
    assert(isequal(w.faceID,(0:n-1)'),'RAH:CM011:SectorIDs','三维结果的原面编号不一致。');
    data.cm011SectorSummary=sector;
    data.researchIdentity=struct('productionQualified',false,'scope','CM-011 同工况的数值对照');
    label='三维研究帧（未达稳态验收）';
    if sector.numericalWallAccepted,label='三维数值稳态（实验精度待验）';end
    scope=sprintf(['低速等温试件，6° 三维对称扇区，原生保存步 %d。%s。' ...
        '同一原三角面入壁热功率 %.3f µW，二维对照 %.3f µW。'], ...
        sector.savedIndex,label,1e6*sector.fullTrianglePowerW,1e6*sector.axisymmetricTrianglePowerW);
    if isfield(sector,'nativeEnergyBalance')
        balance=sector.nativeEnergyBalance;
        scope=[scope sprintf(' 三维原生壁放热 %.3f µW，远场流出 %.3f µW，差额占壁放热 %.4f%%。', ...
            1e6*balance.nativeWallReleaseW,1e6*balance.nativeFarExcessOutW,100*balance.wallMinusFarFraction)];
    end
    source='SU2 三维不可压层流，按对称扇区镜像后映射到原三角面';
    entries={
        'sectorGaugePressure','三维对照 · 相对环境壁压','Pa',w.gaugePressurePa,'数值壁压';
        'sectorHeatFlux','三维对照 · 入壁热流','kW/m²',w.heatFluxWm2/1000,'等温壁 Fourier 热流';
        'sectorWallNode','三维对照 · 壁节点流体温度','K',w.wallNodeK,'原生壁节点温度';
        'sectorFirstNode','三维对照 · 首层流体温度','K',w.firstNodeK,'首层节点焓按定比热还原温度'};
    for i=1:size(entries,1)
        values=entries{i,4};methods=repmat({[source '，' entries{i,5}]},n,1);
        add(entries{i,1},entries{i,2},entries{i,3},values,double(isfinite(values)),methods,scope);
        data.fields(end).stateNote=scope;
        data.fields(end).computedStateLabel=label;
    end
    for pair={{'heatFlux','sectorHeatFlux'},{'pressureGauge','sectorGaugePressure'}}
        keys=pair{1};
        k1=find(strcmp({data.fields.key},keys{1}),1);k2=find(strcmp({data.fields.key},keys{2}),1);
        both=[data.fields(k1).values(:);data.fields(k2).values(:)];
        limits=[min(both),max(both)];
        data.comparisonColorRanges.(keys{1})=limits;
        data.comparisonColorRanges.(keys{2})=limits;
    end
end
if isfield(mesh,'boundaryLayerEdge')
    e=mesh.boundaryLayerEdge; methods=methodNames(e,n);
    state=double(e.algebraicallyValid(:));
    edgePressure=e.pressure;edgeTemperature=e.temperature;
    % 数值占位可用于诊断流向，但不是压力/温度场的有效计算结果。
    placeholder=false(n,1);
    if isfield(e,'thermalExpansionOutOfRangeFace') && ~isempty(e.thermalExpansionOutOfRangeFace)
        rawMask=e.thermalExpansionOutOfRangeFace(:);
        assert(numel(rawMask)==n && all(ismember(rawMask,[0 1])), ...
            'RocketAeroHeat:Viewer:Size','膨胀超域标志须与原面一一对应。');
        expansionOutside=logical(rawMask);placeholder=expansionOutside;
    end
    if isfield(e,'diagnosticPlaceholderFace') && ~isempty(e.diagnosticPlaceholderFace)
        rawMask=e.diagnosticPlaceholderFace(:);
        assert(numel(rawMask)==n && all(ismember(rawMask,[0 1])), ...
            'RocketAeroHeat:Viewer:Size','诊断占位标志须与原面一一对应。');
        placeholder=placeholder | logical(rawMask);
    end
    if isfield(e,'thermalStateValidFace')
        placeholder=placeholder | (logical(e.algebraicallyValid(:)) & ...
            ~logical(e.thermalStateValidFace(:)));
    end
    state(placeholder)=0;edgePressure(placeholder)=NaN;edgeTemperature(placeholder)=NaN;
    edgeReasons=repmat({'内部外流没有接受本面元外缘热力状态；没有有效外缘解'},n,1);
    edgeReasons(placeholder)={'本面元外缘热力状态未知；原记录只提供诊断状态'};
    edgeReasons(expansionOutside)={'本面元外缘热力状态未知；原膨胀计算超出热物性适用范围，留存状态只用于诊断流向'};
    pressureScopes=repmat({'内部工程外流外缘静压；代数有效不代表实验验证或区域接口已匹配'},n,1);
    [temperatureMethods,temperatureScopes]=edgeTemperatureProvenance(e,n);
    for i=find(placeholder)'
        methods{i}=['外缘压力未解：原热力状态未被接受；尝试模型：' methods{i}];
    end
    pressureScopes(placeholder)={'原结果未给出有效外缘压力，流向占位不计入已求面积'};
    % 局部相似鼻区使用自己的同压力/熵外缘。这里仅组合区域结果用于
    % 展示，不改写后体外流或流线，不把分区展示称作接口匹配已完成。
    if isfield(mesh,'heating') && isfield(mesh.heating,'noseModel') && ...
            strcmp(mesh.heating.noseModel,'sphericalLocalSimilarity')
        noseEdge=mesh.heating.noseDistribution.edgeState;
        sphere=noseEdge.geometryFace;keep=sphere & noseEdge.algebraicallyValid;
        if ~isempty(expansionOutside),expansionOutside(sphere)=false;end
        state(sphere)=0;edgePressure(sphere)=NaN;edgeTemperature(sphere)=NaN;
        state(keep)=1;edgePressure(keep)=noseEdge.pressure(keep);
        edgeTemperature(keep)=noseEdge.temperature(keep);
        methods(sphere)={noseEdge.method};temperatureMethods(sphere)={noseEdge.stateMethod};
        pressureScopes(sphere)={[noseEdge.scope '；球鼻/后体接口尚未匹配；只组合显示区域结果']};
        temperatureScopes(sphere)={[noseEdge.scope '；球鼻/后体接口尚未匹配；只组合显示区域结果']};
        edgeReasons(sphere & ~keep)={'球鼻区域外缘状态没有被接受；不以后体或占位状态补齐'};
    end
    add('edgePressure','边界层外缘压力','kPa',edgePressure/1000,state,methods, ...
        '内部工程外流/局部区域模型；代数有效不代表经过实验验证或区域接口已经匹配',pressureScopes,edgeReasons);
    add('edgeTemperature','边界层外缘温度','K',edgeTemperature,state,temperatureMethods, ...
        '内部工程外流的外缘静温；各面元的热物性与适用范围单独记录',temperatureScopes,edgeReasons);
end
if isfield(mesh,'heating')
    h=mesh.heating; state=double(h.eligibleFace(:)); methods=heatingMethods(h,n);
    detail={};
    for key={'boundaryLayerMethod','turbulentMethod','noseMethod','stagnationMethod','thermalPropertyMethod'}
        if isfield(h,key{1}),detail{end+1}=[key{1} ': ' char(h.(key{1}))];end %#ok<AGROW>
    end
    scope=strjoin([{h.validity},detail],'；');
    heatReasons=thermalMissingReasons(mesh,data.missingReason,state,n);
    add('heatFlux','对流热流密度','kW/m²',h.heatFlux/1000,state,methods,scope,{},heatReasons);
    if isfield(h,'kinneyAcreage')
        data.fields(end).computedStateLabel='工程方法对照';
        data.fields(end).stateNote=[data.fields(end).stateNote ...
            '；Kinney侧壁为离散入口假设下的公式对照，尚未独立验证入口/实验精度，不扩充正式求解面积'];
    end
    [recoveryMethods,recoveryScopes]=recoveryTemperatureProvenance(mesh,h,n);
    add('recoveryTemperature','恢复温度','K',h.recoveryTemperature,state,recoveryMethods, ...
        '绝热恢复温度；层流、湍流及球鼻使用各自温度模型，不能将热流关联式作为温度来源',recoveryScopes,heatReasons);
    add('wallInput','热流计算所用壁温','K',h.wallTemperatureByFace,2*ones(n,1), ...
        repmat({'热流边界输入（给定或耦合传入的面元温度）'},n,1), ...
        '此字段仅表示当前热流计算所用输入；即使由耦合程序传入，也不据此标为本次单点求解壁温');
    % 新侧壁方法保留参考流长和原积分法对照，可在现有云图/探针/截面直接查看。
    if isfield(h,'kinneyAcreage')
        k=h.kinneyAcreage;length=k.length;lengthState=double(k.acceptedFace);
        lengthMethods=repmat({[length.method '；' k.sourceURL '；' k.sourcePages]},n,1);
        lengthScope=[length.scope '；仅对当前侧壁热流已接受面元显示；' ...
            '默认离散起点零积分或局部常状态前缀是工程假设，未核验真实驻点/附着线入口'];
        add('heatingReferenceLength','热流所用参考流长','m',length.selectedLength, ...
            lengthState,lengthMethods,lengthScope);
        data.fields(end).computedStateLabel='工程参考量';
        [baselineMethods,baselineScopes]=baselineHeatingProvenance(h,n);
        add('baselineHeatFlux','原工程方法热流对照','kW/m²',k.baselineHeatFlux/1000, ...
            double(isfinite(k.baselineHeatFlux)),baselineMethods, ...
            '同一当前外流、路径和壁温下的原工程方法对照；鼻区与侧壁来源分别保留；不是独立实验或求解壁温',baselineScopes);
        % 两种方法共享同一范围；切换变量时可保持“对照统一色标”。
        both=[h.heatFlux(:);k.baselineHeatFlux(:)]/1000;both=both(isfinite(both));
        if ~isempty(both)
            limits=[min(0,min(both)),max(0,max(both))];
            if limits(1)==limits(2),limits=limits+[-1 1];end
            data.comparisonColorRanges=struct('heatFlux',limits,'baselineHeatFlux',limits);
        end
    end
end
% 风洞圆柱只显示可辨识的无量纲预测；实验点保存在测带表中。
if isfield(mesh,'lewisCylinder')
    nu=mesh.lewisCylinder.nusselt;
    scope='TM 1432 第20页空气四次式；Re_D=177000；中央66 mm、迎流角绝对值不超过55°。';
    reasons=cellstr(data.region);
    reasons(strcmp(reasons,'靠近风洞侧壁'))={'靠近风洞侧壁，未建模跨向变化'};
    reasons(strcmp(reasons,'中央跨度 · 模型角域外'))={'超出迎流层流级数的角度范围'};
    add('nusselt','局部努塞尔数 Nu_D','1',nu,double(isfinite(nu)), ...
        repmat({'Frössling 圆柱空气四次式，周向面带平均'},n,1),scope,{},reasons);
    data.fields(end).computedStateLabel='工程预测（研究对照）';
    if isfield(mesh,'lewisCondition')
        c=mesh.lewisCondition;w=mesh.lewisWall;
        pressureMethod=repmat({'Hiemenz 圆柱前部外缘速度，Bernoulli 静压，按周向面带平均'},n,1);
        pressureScope='中央附着区；Hiemenz 外缘取自不同雷诺数的试验，壁压尚无同工况独立验证。';
        add('pressure','绝对壁压','kPa',w.pressurePa/1000,double(isfinite(w.pressurePa)),pressureMethod,pressureScope,{},reasons);
        add('pressureGauge','相对上游壁压','Pa',w.gaugePa,double(isfinite(w.gaugePa)),pressureMethod,pressureScope,{},reasons);
        add('heatFlux','入壁热流','kW/m²',w.heatFluxWm2/1000,double(isfinite(w.heatFluxWm2)), ...
            repmat({'Frössling 带平均 Nu × k/D × (Tad − Tw)'},n,1),scope,{},reasons);
        add('wallInput','给定壁温','K',repmat(c.wallTemperatureK,n,1),2*ones(n,1), ...
            repmat({'研究工况输入 · 等温圆柱'},n,1),'设定值覆盖完整圆柱侧壁');
        % Tad 是试验换热系数统一采用的温度基准，不充当背面回流的局部恢复温度。
        add('heatReferenceTemperature','换热定义所用绝热壁温','K',repmat(c.adiabaticWallK,n,1),ones(n,1), ...
            repmat({'NASA TM-83557 式 (9)–(10)，来流总静温差与空气恢复因子'},n,1), ...
            '整组换热系数共用的定义基准；背面局部恢复温度未知');
        data.missingReason=reasons;data.missingReason(isfinite(w.heatFluxWm2))={''};
        if isfield(mesh,'lewisTrial')
            % 试算使用独立字段，绝不写入下方按 pressure / heatFlux 统计的覆盖和整壁功率。
            trial=mesh.lewisTrial;
            % 模块：探针和导出均引用本次载入的终态。旧版结构只有固定 CFL 的 120 次记录。
            if isfield(trial,'innerIterations') && isfield(trial,'schemeName')
                scheme=trial.schemeName;iterations=trial.innerIterations;
            else
                scheme='固定 CFL';iterations=120;
            end
            scopeTrial=sprintf('%s，首个物理步 Δt=%.6g s，%d 次内迭代，九分量残差未通过。源帧 %s', ...
                scheme,trial.physicalDtS,iterations,trial.sourceFrameSHA256);
            add('trialPressure','未收敛试算 · 壁压','kPa',trial.pressurePa/1000,ones(n,1), ...
                repmat({'SU2 INC_RANS + SST + LM；上游中心点定压，原面带积分'},n,1),scopeTrial,repmat({scopeTrial},n,1));
            data.fields(end).computedStateLabel='未收敛试算';
            add('trialHeatFlux','未收敛试算 · 入壁热流','kW/m²',trial.heatFluxWm2/1000,ones(n,1), ...
                repmat({'保存态等温圆柱法向 Fourier 通量；圆柱标记单独重建，原面带积分'},n,1),scopeTrial,repmat({scopeTrial},n,1));
            data.fields(end).computedStateLabel='未收敛试算';
        end
    else
        missing={'pressure','绝对壁压','kPa','本工况没有可用的绝对壁压解'; ...
             'heatFlux','入壁热流','kW/m²','缺少该次试验的导热系数与壁面至绝热壁温差'; ...
             'wallInput','输入壁温','K','原文说明恒温控制，图9没有给出该次设定温度'};
        for j=1:size(missing,1)
            reason=missing{j,4};
            add(missing{j,1},missing{j,2},missing{j,3},nan(n,1),zeros(n,1), ...
                repmat({reason},n,1),reason,{},repmat({reason},n,1));
        end
        data.missingReason=repmat(missing(2,4),n,1);
    end
end
% 模块三：显式区分输入壁温与壁温解；单一外流/热流工况本身不确定壁温。
% 尚无材料热响应解时仍允许探针查看“未求解”，不能以输入壁温补齐。
if isfield(mesh,'hybridWall')
    for term={'edgeTemperature','recoveryTemperature'}
        if any(strcmp({data.fields.key},term{1})),continue;end
        if strcmp(term{1},'edgeTemperature'),label='外缘静温';else,label='恢复温度';end
        reason=['本次混合计算未提供具有独立定义的' label];
        add(term{1},label,'K',nan(n,1),zeros(n,1),repmat({reason},n,1),reason,{},repmat({reason},n,1));
    end
end
if isfield(mesh,'hybridSolvedWallTemperatureK')
    t=solvedBeadThermal(mesh);
    add('wallSolved','稳态求解壁温','K',mesh.hybridSolvedWallTemperatureK,ones(n,1), ...
        repmat({t.wallMethod},n,1),t.scope);
    data.fields(end).stateNote=t.wallMethod;
    data.fields(end).computedStateLabel='研究数值解';
else
    add('wallSolved','求解壁温','K',nan(n,1),zeros(n,1), ...
    repmat({'未求解：当前未提供壁体热响应或其他独立壁温解'},n,1), ...
    '单点压力和对流热流计算使用壁温边界输入；求解壁温还需独立热平衡或材料热响应及边界条件');
    data.fields(end).availability='notProvidedInSnapshot';
end
% 模块四：采用求解器实际接受的时间轴，不插值、不恢复已失效的历程。
histories=buildAerothermalViewerHistory(mesh);
flowHistories=buildAerothermalFlowHistory(mesh);
if isempty(histories),histories=flowHistories;
elseif ~isempty(flowHistories),histories=[histories flowHistories];end
for j=1:numel(histories)
    h=histories(j); nt=numel(h.history.time); last=(nt-1)*n+(1:n);
    methods=reshape(h.history.methodNames(h.history.methodIndex(last)+1),[],1);
    existing=[];
    if ~isempty(data.fields),existing=find(strcmp({data.fields.key},h.key));end
    add(h.key,h.label,h.unit,h.history.values(last),h.history.status(last),methods,h.scope);
    data.fields(end).history=h.history;
    if ~isempty(existing)
        data.fields(existing)=data.fields(end);data.fields(end)=[];
    end
end
% 区域台账只接受当前场量的直接超域记录。历史场和独立球鼻外缘具有各自的来源。
if any(expansionOutside)
    for key={'edgePressure','edgeTemperature'}
        k=find(strcmp({data.fields.key},key{1}));
        if isempty(data.fields(k).history),data.fields(k).outOfDomainFace=expansionOutside;end
    end
end
% 模块五：按真实三角形面积分区汇总，给定值单独统计。
regions=unique(data.region,'stable'); rows=struct([]);
for j=1:numel(regions)
    mask=strcmp(data.region,regions{j});
    rows(j).name=regions{j}; rows(j).faces=nnz(mask); rows(j).area=sum(area(mask));
    rows(j).fields=struct;
    for k=1:numel(data.fields)
        field=data.fields(k); s=field.status(mask); ar=area(mask);
        rows(j).fields.(field.key)=struct('computedArea',sum(ar(s==1)), ...
            'prescribedArea',sum(ar(s==2)),'missingArea',sum(ar(s==0)));
    end
end
data.regions=rows; data.totalArea=sum(area);
% 模块六：冻结工况只取结果自己的自由流与实际热流壁温输入，绝不读待编辑配置。
data.resultCondition=snapshotCondition(mesh);
data.caseNote=[data.resultCondition.caption ...
    '。单点基准快照与热响应实际时间站分别列示；给定壁温不是求解壁温。'];
if isfield(mesh,'hybridWall')
    data.resultCondition.physics=mesh.hybridCondition;
    data.caseNote='混合求解单工况；逐面标注计算方法，正热流入壁。';
    hc=mesh.hybridCondition;
    if all(isfield(hc,{'ambientPressurePa','ambientTemperatureK','velocityMps'})) && ...
            isnumeric(hc.velocityMps) && numel(hc.velocityMps)==3 && all(isfinite(hc.velocityMps(:))) && ...
            isnumeric(hc.ambientPressurePa) && isscalar(hc.ambientPressurePa) && isfinite(hc.ambientPressurePa) && ...
            isnumeric(hc.ambientTemperatureK) && isscalar(hc.ambientTemperatureK) && isfinite(hc.ambientTemperatureK)
        u=reshape(hc.velocityMps,1,3);
        if norm(u)>0,data.flowDirection=u/norm(u);end
        data.caseNote=sprintf('混合求解 · 来流速度 %.6g m/s · 环境压力 %.8g Pa · 环境温度 %.6g K · 正热流入壁', ...
            norm(u),hc.ambientPressurePa,hc.ambientTemperatureK);
    end
    if isfield(hc,'caseLabel') && ischar(hc.caseLabel) && isrow(hc.caseLabel)
        data.caseNote=[hc.caseLabel ' · ' data.caseNote];
    end
    if all(isfield(hc,{'alphaDeg','betaDeg','ambientPressurePa','ambientTemperatureK'}))
        data.resultCondition.alphaDeg=hc.alphaDeg;data.resultCondition.betaDeg=hc.betaDeg;
        data.resultCondition.freestreamPressurePa=hc.ambientPressurePa;
        data.resultCondition.freestreamTemperatureK=hc.ambientTemperatureK;
    end
    data.resultCondition.caption=data.caseNote;
    data.defaultField='heatFlux';
end
% 带线试件按原面保留珠体、引线和接头分区。实验只对照扣除引线散热后的珠部平均量。
if isfield(mesh,'thermistorSidewind')
    s=mesh.thermistorSidewind;data.thermistorSidewind=s;
    data.coordinateNote='SI · 长轴 X · 来流 −Y · 正热流入壁';
    if isfield(mesh,'thermistorWire')
        w=mesh.thermistorWire;data.thermistorWire=w;
        data.caseNote=sprintf('B 珠 · 3 m/s 侧吹 · 线根 %.2f K · 名义线端 %.3f K · 长引线温度边界',mesh.hybridCondition.beadWallK,w.profile.endK);
        add('wireHeatFluxChange','引线边界热流差','kW/m²',mesh.thermistorWireFluxChange/1000,ones(n,1), ...
            repmat({'长引线三维入壁热流减同面均温控制解'},n,1),'同一原面和来流，均温参照与长引线终态分别记录');
        if isfield(w,'thermometry')
            m=w.thermometry;
            data.caseNote=sprintf('B 珠 · 3 m/s · 电阻等效温度 %.2f K · 推算珠表面 %.4f K · %s',m.internalK,m.surfaceK,m.fieldMethod);
            add('thermometryHeatFluxChange','测温修正热流差','kW/m²',mesh.thermometryFluxChange/1000,ones(n,1), ...
                repmat({'修正场减 v317 同面长引线热流'},n,1),'原长引线参照为 12 K 珠表面超温。原生核对另外推进了 20 步，同步比较见测温记录');
        end
    else
        data.caseNote=sprintf('带线 B 型热敏珠 · 全周侧吹 3 m/s · 给定壁温 %.6g K · 稳态',mesh.hybridCondition.wallTemperatureK);
    end
    data.resultCondition.caption=data.caseNote;
    k=find(strcmp({data.fields.key},'heatFlux'),1);
    data.fields(k).stateNote=sprintf('珠部含接头 %.5f mW，引线及端面 %.5f mW，七区合计向空气放热 %.5f mW。侧吹局部实验精度未知。', ...
        1000*s.beadAndJointsW,1000*(s.totalOutwardW-s.beadAndJointsW),1000*s.totalOutwardW);
    if isfield(mesh,'thermistorMeasuredPoint')
        m=mesh.thermistorMeasuredPoint;data.measuredPoint=m;
        data.caseNote='B 型热敏珠 · 图 7 · 3 m/s 名义侧吹 · 壁温 305.15 K';
        if isfield(m,'coldSource'),data.caseNote=sprintf('%s · 独立冷启动 %d 步',data.caseNote,m.coldSource.iterations);end
        data.resultCondition.caption=data.caseNote;
        data.fields(k).stateNote=sprintf('珠部平均 h %.2f W/(m²·K)，图 7 读数 %.2f（%.2f–%.2f），差 %+.2f%%。%s', ...
            m.beadHWm2K,m.experiment.hWm2K,m.experiment.graphReadRangeWm2K,100*m.relativeDifference,m.experiment.mounting);
        if isfield(m,'coldComparisonNote'),data.fields(k).stateNote=[data.fields(k).stateNote ' ' m.coldComparisonNote];end
    end
end
if isfield(mesh,'thermistorLead') || isfield(mesh,'thermistorMidbody')
    data.coordinateNote='SI · 长轴 X · 来流 −X · 正热流入壁';
    if isfield(mesh,'hybridSolvedWallTemperatureK')
        t=solvedBeadThermal(mesh);
        data.caseNote=sprintf('%s · 3 m/s · 整件净对流功率 %.9g mW · 求解壁温 %.9g K', ...
            mesh.hybridCondition.caseLabel,1000*t.netHeatingPowerW,t.wallTemperatureK);
    else
        data.caseNote=sprintf('%s · 3 m/s · 给定壁温 %.12g K · 30° 对称稳态', ...
            mesh.hybridCondition.caseLabel,mesh.hybridCondition.wallTemperatureK);
    end
    data.resultCondition.caption=data.caseNote;
    k=find(strcmp({data.fields.key},'heatFlux'),1);
    if isfield(mesh,'thermistorMidbody')
        % 中段加密有独立终帧；原混合网只作为同一批原三角面的对照字段。
        s=mesh.thermistorMidbody;data.thermistorMidbody=s;
        data.fields(k).stateNote=sprintf('整珠放热 %.5f mW，名义面积平均换热 %.2f W/(m²·K)，实验图读 %.2f W/(m²·K)。', ...
            1000*s.totalOutwardW,s.beadMeanHWm2K,s.experimentMeanHWm2K);
        scope='光滑珠体轴向对称稳态，壁剪切和首层速度投影到沿流向的母线切向。';
    else
        s=mesh.thermistorLead;data.thermistorLead=s;
        if isfield(s,'thermal')
            data.fields(k).stateNote=sprintf('珠体含接头 %.6f mW，引线含端面 %.6f mW。七区合计 %.6f mW，稳态热导 %.6f mW/K。', ...
                1000*s.beadAndJointsW,1000*(s.wireSidesW+s.wireEndsW),1000*s.totalOutwardW,1000*s.thermal.conductanceWK);
            if s.thermal.deltaTemperatureK==0
                data.fields(k).stateNote='壁面与环境同温，七区各面热流均为零。没有唯一热峰，按功率除以温差计算的平均换热系数未定义。';
            elseif isfield(s,'beadMeanComparison')
                obs=s.beadMeanComparison;
                data.fields(k).stateNote=[data.fields(k).stateNote sprintf( ...
                    ' 珠部计算平均换热系数 %.2f W/(m²·K)。图 7、8 在温差 %.9g K 处无匹配测点，实验偏差未知。', ...
                    obs.predictedHWm2K,obs.deltaTK)];
            end
        else
            data.fields(k).stateNote=sprintf('珠体含接头 %.4f mW，引线含端面 %.4f mW。珠部平均换热 %.2f W/(m²·K)，实验图读 %.2f W/(m²·K)。', ...
                1000*s.beadAndJointsW,1000*(s.wireSidesW+s.wireEndsW),s.beadMeanHWm2K,s.experimentMeanHWm2K);
        end
        % 珠部平均量随热流页展示，面元云图仍读取本域的逐面壁热。
        if isfield(s,'meanConvection')
            m=s.meanConvection;
            if m.available
                note=sprintf(['珠部平均换热系数 / 向空气放热：式（28）%.2f W/(m²·K) / %.4f mW，' ...
                    '三维积分 %.2f W/(m²·K) / %.4f mW。用表 1 极限 Nu 替换形状项，工程放热变化 %+.2f µW。'], ...
                    m.meanHWm2K,1000*m.outwardW,m.surfaceMeanHWm2K,1000*m.surfaceOutwardW,1e6*m.tableMinusFormulaW);
            else
                note=['珠部式（28）对照不可用。' m.reason];
            end
            data.fields(k).stateNote=[data.fields(k).stateNote ' ' note];
        end
        scope='带线试件对称稳态，母线正向取沿流向。端面及连接处按各自面元的切向计算。';
    end
    if isfield(mesh,'thermistorLeadFlow')
        flow=mesh.thermistorLeadFlow;
        add('leadDownstreamCf','沿流向壁面摩擦系数','1',flow.downstreamCf,ones(n,1),repmat({'原生壁面摩擦向量投影到沿流向母线切向，三角面平均'},n,1),scope);
        add('leadReverse','近壁反向流动','0 / 1',flow.reverse,ones(n,1),repmat({'壁剪切和首层流速的母线切向分量均为负时取 1'},n,1),scope);
        add('leadFirstYPlus','首层 y+ 估算','1',flow.estimatedYPlus,ones(n,1),repmat({'按原生壁面摩擦系数、实际首层法向距离和运动黏度估算'},n,1),scope);
        data.fields(end-1).stateNote=sprintf('1 表示两项均反向。%d 面，%.6f mm²。',flow.summary.reverseFaces,1e6*flow.summary.reverseAreaM2);
    end
end
% 热敏珠只具备整珠平均换热的实验读数；探针上的局部压力和热流均保留研究身份。
if isfield(mesh,'thermistorBead')
    bead=mesh.thermistorBead;
    data.thermistorBead=bead;
    data.coordinateNote='SI · 长轴 X · 来流 −X · 正热流入壁';
    powerMode=isfield(mesh.hybridCondition,'netHeatingPowerW');
    if powerMode
        thermal=bead.thermal;
        data.caseNote=sprintf('%s · 净加热功率 %.12g mW · 稳态求解壁温 %.9g K · 无引线名义珠体', ...
            data.caseNote,1000*thermal.netHeatingPowerW,thermal.wallTemperatureK);
    else
        data.caseNote=sprintf('%s · 给定壁温 %.12g K · 无引线名义珠体',data.caseNote,mesh.hybridCondition.wallTemperatureK);
    end
    if isfield(bead,'cellAngleDeg'),data.caseNote=sprintf('%s · 周向单元 %g°',data.caseNote,bead.cellAngleDeg);end
    data.resultCondition.caption=data.caseNote;
    k=find(strcmp({data.fields.key},'heatFlux'),1);
    if isfield(bead,'experiment')
        observation=bead.experiment;
        data.fields(k).stateNote=sprintf(['名义面积上的整珠平均换热系数 %.1f W/(m²·K)，论文图 7 约 %.1f，' ...
            '相差 %+.1f%%。图读范围 %.1f–%.1f W/(m²·K)。局部热流缺少逐面实验数据。'], ...
            bead.meanH,observation.hWm2K,100*bead.relativeDifference,observation.graphReadRangeWm2K);
    elseif bead.thermal.deltaTemperatureK==0
        data.fields(k).stateNote='珠面与来流同温，全部原面对流热流为 0。零温差下按功率除以温差得到的平均换热系数未定义。';
    else
        data.fields(k).stateNote=sprintf('恒物性对流扩散场按 %.12g K 壁面温差换算。名义面积平均换热系数 %.2f W/(m²·K)。', ...
            bead.thermal.deltaTemperatureK,bead.meanH);
    end
    if powerMode
        data.fields(k).stateNote=sprintf('%s 净功率 %.9g mW，原面放热积分 %.9g mW，空气侧热导 %.9g mW/K。', ...
            data.fields(k).stateNote,1000*thermal.netHeatingPowerW,1000*thermal.surfaceReleaseW,1000*thermal.conductanceWK);
    end
    if isfield(bead,'comparisonUnavailable')
        data.fields(k).stateNote=[data.fields(k).stateNote ' ' bead.comparisonUnavailable];
    end
    if isfield(bead,'baselineH')
        data.fields(k).stateNote=[data.fields(k).stateNote sprintf(' 15° 基准网为 %.2f W/(m²·K)，周向细化后变化 %+.2f%%。', ...
            bead.baselineH,100*(bead.meanH/bead.baselineH-1))];
    end
    k=find(strcmp({data.fields.key},'nearWallCellTemperature'),1);
    data.fields(k).label='首层流体节点温度';
    data.fields(k).stateNote='首层流体节点温度由三维温差标量求得，壁面温度由边界条件给定。';
    if powerMode,data.fields(k).stateNote='首层流体温度由参照温差标量按稳态集总壁温换算。';end
end
if isfield(mesh,'thermistorBead') || isfield(mesh,'thermistorLead')
    % 声明静态试件场的真实单值。空场不产生色标，固定的比较／手动色标另行保留。
    for k=1:numel(data.fields)
        f=data.fields(k);valid=f.values(f.status>0 & isfinite(f.values));
        if ~isempty(valid) && all(valid==valid(1))
            data.uniformFieldValues.(f.key)=struct('value',valid(1),'unit',f.unit);
        end
    end
end
% 热敏珠的配对控制随原面云图、探针和 CSV 共用同一组已核对字段。
if isfield(mesh,'thermistorPaired')
    pair=mesh.thermistorPaired;pc=pair.summary;data.thermistorPair=pc;
    if strcmp(pair.kind,'wiredRefined')
        note=sprintf('中段加密下，光滑珠体放热 %.5f mW，带线珠体含接头 %.5f mW。', ...
            1000*pc.powerW.pairedSmooth,1000*pc.powerW.wiredBeadAndJoints);
    else
        note=sprintf('3 m/s 珠部放热：原结构网光滑 %.5f mW，配对混合网光滑 %.5f mW，带线含接头 %.5f mW。', ...
            1000*pc.powerW.oldSmooth,1000*pc.powerW.pairedSmooth,1000*pc.powerW.wiredBeadAndJoints);
    end
    tab={'pairedSmoothHeat','配对光滑入壁热流','kW/m²',pair.heatFluxWm2/1000,'光滑混合网强等温标量的壁面离散反力';
         'pairedDeltaHeat',['入壁热流差 · ' pair.order],'W/m²',pair.deltaHeatWm2,'同位置两解的入壁热流相减';
         'pairedDeltaPressure',['壁压差 · ' pair.order],'Pa',pair.deltaPressurePa,'同位置两解的绝对壁压相减';
         'pairedSmoothReverse','配对光滑近壁反向','0 / 1',pair.reverse,'光滑壁面剪切与首层速度沿流向的分量均为负时取 1';
         'pairedReverseChange',['反向区变化 · ' pair.order],'1',pair.reverseChange,'两解的近壁反向标记相减，取 −1、0 或 1'};
    if strcmp(pair.kind,'midbody')
        note=sprintf('原混合网放热 %.5f mW，中段加密后 %.5f mW。',1000*pc.powerW.pairedSmooth,1000*mesh.thermistorMidbody.totalOutwardW);
        tab{1,2}='入壁热流 · 加密前';tab{4,2}='近壁反向 · 加密前';
    end
    status=double(isfinite(pair.heatFluxWm2));
    primary=data.fields(strcmp({data.fields.key},'heatFlux'));
    valid=[primary.values(isfinite(primary.values));pair.heatFluxWm2(status>0)/1000];
    if isfield(data,'comparisonColorRanges') && isfield(data.comparisonColorRanges,'heatFlux')
        valid=[valid(:);data.comparisonColorRanges.heatFlux(:)];
    end
    limits=[min(valid) max(valid)];data.comparisonColorRanges.heatFlux=limits;
    data.comparisonColorRanges.pairedSmoothHeat=limits;
    if isfield(data.comparisonColorRanges,'baselineHeatFlux'),data.comparisonColorRanges.baselineHeatFlux=limits;end
    data.comparisonColorRanges.pairedSmoothReverse=[0 1];data.comparisonColorRanges.pairedReverseChange=[-1 1];
    reasons=repmat({'光滑控制没有与这个接头、引线或端面重合的原三角面'},n,1);
    for j=1:size(tab,1)
        add(tab{j,1},tab{j,2},tab{j,3},tab{j,4},status,repmat(tab(j,5),n,1),pair.scope,{},reasons);
        data.fields(end).stateNote=[note ' ' pair.scope];
        if ismember(j,[2 3])
            bound=max(abs(tab{j,4}(status>0)));
            if bound>0,data.comparisonColorRanges.(tab{j,1})=[-bound bound];end
        end
    end
end
% 带线旧网与加密网在七区均可对照，和只在共同珠面定义的几何差别分开显示。
if isfield(mesh,'thermistorLeadGrid')
    g=mesh.thermistorLeadGrid;
    data.leadGridReference=rmfield(g,{'baseHeatFluxWm2','deltaHeatWm2','deltaPressurePa','deltaNearWallK'});
    scope=sprintf('两档带线网共用 10,272 个原面和首层位置，给定壁温 %.12g K。',mesh.hybridCondition.wallTemperatureK);
    tab={'leadGridBaseHeat','入壁热流 · 带线原网','kW/m²',g.baseHeatFluxWm2/1000,'原网强温差标量的壁节点反力按邻接面积守恒分配'; ...
         'leadGridDeltaHeat','入壁热流差 · 带线加密减原网','W/m²',g.deltaHeatWm2,'两套带线体网的原面入壁热流相减'; ...
         'leadGridDeltaPressure','壁压差 · 带线加密减原网','Pa',g.deltaPressurePa,'两套带线体网的同面绝对壁压相减'; ...
         'leadGridDeltaTemperature','首层温差 · 带线加密减原网','K',g.deltaNearWallK,'同一首层节点位置的面平均流体温度相减'};
    if mesh.hybridCondition.wallTemperatureK~=305.15,tab{1,5}='原网 305.15 K 壁反力按实际温差比换算';end
    for j=1:size(tab,1)
        add(tab{j,1},tab{j,2},tab{j,3},tab{j,4},ones(n,1),repmat(tab(j,5),n,1),scope);
        if all(tab{j,4}==tab{j,4}(1))
            data.uniformFieldValues.(tab{j,1})=struct('value',tab{j,4}(1),'unit',tab{j,3});
        elseif j>1
            bound=max(abs(tab{j,4}));if bound>0,data.comparisonColorRanges.(tab{j,1})=[-bound bound];end
        end
    end
    values=[mesh.hybridWall.heatFluxWm2;g.baseHeatFluxWm2]/1000;
    if isfield(data,'comparisonColorRanges') && isfield(data.comparisonColorRanges,'heatFlux')
        values=[values;data.comparisonColorRanges.heatFlux(:)];
    end
    limits=[min(values) max(values)];
    if limits(2)>limits(1)
        data.comparisonColorRanges.leadGridBaseHeat=limits;data.comparisonColorRanges.heatFlux=limits;
        if isfield(mesh,'thermistorPaired'),data.comparisonColorRanges.pairedSmoothHeat=limits;end
    end
end
% 外边界敏感性取两次冷态的同面结果，值、单位和相减方向随导出保留。
if isfield(mesh,'thermistorOuterDomain')
    domain=mesh.thermistorOuterDomain;
    if domain.available
        faceFields={'heat30Wm2','deltaHeatWm2','deltaPressurePa','deltaFirstK'};
        if strcmp(domain.mode,'sidewindGivenWall'),faceFields=[faceFields {'gauge30Pa','first30K'}];end
        data.outerDomainComparison=rmfield(domain,faceFields);
        note=sprintf('%s 20D / 30D 壁温 %.7f / %.7f K，整件放热 %.7f / %.7f mW。', ...
            domain.scope,domain.wall20K,domain.wall30K,1000*domain.outward20W,1000*domain.outward30W);
        tab={'outerHeat30','入壁热流 · 30D 外域','kW/m²',domain.heat30Wm2/1000,'30D 域强等温标量的原面壁反力'; ...
            'outerDeltaHeat','入壁热流差 · 30D 减 20D','W/m²',domain.deltaHeatWm2,'两域同一原面和热边界下的入壁热流相减'; ...
            'outerDeltaPressure','壁压差 · 30D 减 20D','Pa',domain.deltaPressurePa,'两次层流解的同面绝对壁压相减'; ...
            'outerDeltaFirst','首层温差 · 30D 减 20D','K',domain.deltaFirstK,'两个域相同首层节点上的面平均流体温度相减'};
        if strcmp(domain.mode,'sidewindGivenWall')
            tab{1,5}='30D 全周侧吹 6 K 强温差原生壁反力，按实际方位的原面面积守恒分配';
        elseif strcmp(domain.mode,'givenPower'),tab{1,5}='由 30D 空气侧热导反求壁温，按温差比换算该域壁反力';
        elseif domain.wall30K~=305.15,tab{1,5}='30D 域 305.15 K 壁反力按实际给定温差比换算';end
        for j=1:size(tab,1)
            add(tab{j,1},tab{j,2},tab{j,3},tab{j,4},ones(n,1),repmat(tab(j,5),n,1),domain.scope);
            data.fields(end).stateNote=note;
            if all(tab{j,4}==tab{j,4}(1))
                data.uniformFieldValues.(tab{j,1})=struct('value',tab{j,4}(1),'unit',tab{j,3});
            elseif j>1
                bound=max(abs(tab{j,4}));data.comparisonColorRanges.(tab{j,1})=[-bound bound];
            end
        end
        % 侧吹对照另外提供原生表压与首层温度，同一量的两域云图使用共同色标。
        if strcmp(domain.mode,'sidewindGivenWall')
            views={'outerGauge30','相对环境壁压 · 30D 外域','Pa','pressureGauge', ...
                domain.gauge30Pa,'30D 侧吹壁节点相对环境的表压按原面平均'; ...
                'outerFirst30','首层流体温度 · 30D 外域','K','nearWallCellTemperature', ...
                domain.first30K,'30D 对应首层节点的强温差标量按 6 K 换算并取面平均'};
            for j=1:size(views,1)
                baseField=data.fields(strcmp({data.fields.key},views{j,4}));
                values=views{j,5};range=[min([baseField.values(:);values]) max([baseField.values(:);values])];
                add(views{j,1},views{j,2},views{j,3},values,ones(n,1),repmat(views(j,6),n,1),domain.scope);
                data.fields(end).stateNote=note;
                if range(2)>range(1)
                    data.comparisonColorRanges.(views{j,1})=range;data.comparisonColorRanges.(views{j,4})=range;
                end
            end
        end
        limits=[mesh.hybridWall.heatFluxWm2;domain.heat30Wm2]/1000;
        if isfield(data,'comparisonColorRanges') && isfield(data.comparisonColorRanges,'heatFlux')
            limits=[limits;data.comparisonColorRanges.heatFlux(:)];
        end
        limits=[min(limits) max(limits)];
        if limits(2)>limits(1)
            data.comparisonColorRanges.heatFlux=limits;data.comparisonColorRanges.outerHeat30=limits;
            if isfield(data.comparisonColorRanges,'leadGridBaseHeat'),data.comparisonColorRanges.leadGridBaseHeat=limits;end
            if isfield(data.comparisonColorRanges,'pairedSmoothHeat'),data.comparisonColorRanges.pairedSmoothHeat=limits;end
        end
    else
        data.outerDomainComparison=domain;k=find(strcmp({data.fields.key},'heatFlux'),1);
        data.fields(k).stateNote=[data.fields(k).stateNote ' ' domain.note];
    end
end
% 热敏珠前引线端面近壁细化：分开细网首层与共同物理位置，仅在共同位置定义温差。
if isfield(mesh,'thermistorTipSpacing')
    t=mesh.thermistorTipSpacing;g=t.summary;data.tipSpacingComparison=g;
    note=sprintf('原网／细网前端面放热 %.5f／%.5f µW，整件 %.6f／%.6f mW。细网最强两面 %d、%d，放热密度相差 %.3g W/m²。', ...
        1e6*g.regions(7).baseOutW,1e6*g.regions(7).fineOutW,1e3*g.outwardBaseW,1e3*g.outwardFineW, ...
        g.finePeakFace,g.fineSecondFace,g.finePeakGapWm2);
    tab={'tipFineHeat','入壁热流 · 端面细网','kW/m²',t.fineHeat/1000,'细网强等温壁节点反力按原面面积守恒分配'; ...
        'tipDeltaHeat','入壁热流差 · 端面细网减原网','W/m²',t.deltaHeat,'两套 20D 全周解在同一原面的入壁热流相减'; ...
        'tipFineGauge','相对环境壁压 · 端面细网','Pa',t.fineGauge,'细网不可压流动的原面壁节点表压平均'; ...
        'tipDeltaPressure','壁压差 · 端面细网减原网','Pa',t.deltaPressure,'同一原面两网绝对壁压相减'; ...
        'tipFineFirst','实际首层流体温度 · 端面细网','K',t.fineFirstK,'细网实际第一层的三个对应节点按 6 K 温差换算并平均'; ...
        'tipMatched','原首层位置流体温度 · 端面细网','K',t.matchedK,'由细网终场直接读取保留的原首层三个物理节点并平均'; ...
        'tipDeltaMatched','同位置流体温差 · 细网减原网','K',t.deltaMatchedK,'在相同三个原首层物理节点读取两网温度，各自面平均后相减'; ...
        'tipHeight','实际首层法向距离 · 端面细网','µm',1e6*t.fineDistanceM,'三个壁顶点至对应细网首层节点的位移投影到原面单位法向后取平均'};
    for j=1:size(tab,1)
        key=tab{j,1};vals=tab{j,4};add(key,tab{j,2},tab{j,3},vals,ones(n,1),repmat(tab(j,5),n,1),g.scope);
        data.fields(end).stateNote=note;
        if startsWith(key,'tipDelta')
            limit=max(abs(vals));if limit>0,data.comparisonColorRanges.(key)=[-limit limit];end
        end
    end
    % 96 档返回原三角形平均值；两档温差使用前引线近壁细化后的同一组流体节点。
    if isfield(mesh,'thermistorAngularSpacing')
        az=mesh.thermistorAngularSpacing;ag=az.summary;data.angularSpacingComparison=ag;
        note=sprintf('48／96 档整件放热 %.6f／%.6f mW。96 档最强两面 %d、%d，相差 %.3g W/m²。', ...
            1e3*ag.outward48W,1e3*ag.outward96W,ag.finePeakFace,ag.fineSecondFace,ag.finePeakGapWm2);
        tab={'angular96Heat','入壁热流 · 周向 96 档','kW/m²',az.heat96/1000,'强等温反力按细面节点面积分摊，再按细三角线性形函数积分到原面'; ...
            'angular96Gauge','相对环境壁压 · 周向 96 档','Pa',az.gauge96,'细网原生壁节点表压在每个交叠多边形上积分，除以原面面积'; ...
            'angular96First','同位置首层温度 · 周向 96 档','K',az.common96,'96 档终场在 48 档端面细网的三个首层节点读取温度并取面平均'; ...
            'angularDeltaHeat','入壁热流差 · 96 减 48','W/m²',az.deltaHeat,'两档各自的壁反力守恒分配到相同原面后相减'; ...
            'angularDeltaPressure','壁压差 · 96 减 48','Pa',az.deltaPressure,'相同环境压力下，两网在同一原面上的平均表压相减'; ...
            'angularDeltaFirst','同位置首层温差 · 96 减 48','K',az.deltaCommon,'两网在完全相同的三个首层物理位置取温度平均后相减'};
        for j=1:size(tab,1)
            key=tab{j,1};vals=tab{j,4};add(key,tab{j,2},tab{j,3},vals,ones(n,1),repmat(tab(j,5),n,1),ag.scope);
            data.fields(end).stateNote=note;
            if startsWith(key,'angularDelta')
                limit=max(abs(vals));if limit>0,data.comparisonColorRanges.(key)=[-limit limit];end
            end
        end
    end
    % 同量绝对场共用色标，实际存在的对照字段参与范围计算。
    sets={{'heatFlux','outerHeat30','tipFineHeat','angular96Heat'},{'pressureGauge','outerGauge30','tipFineGauge','angular96Gauge'}, ...
        {'nearWallCellTemperature','outerFirst30','tipFineFirst','tipMatched','angular96First'}};
    for j=1:numel(sets)
        keys=sets{j};found=ismember({data.fields.key},keys);
        vals=vertcat(data.fields(found).values);bounds=[min(vals) max(vals)];
        for k=keys
            if ismember(k{1},{data.fields.key}),data.comparisonColorRanges.(k{1})=bounds;end
        end
    end
end
% 三档中段网距在同一热边界下比较，前两档取真实归档，第三档为当前原面主解。
if isfield(mesh,'thermistorThreeGrids')
    g=mesh.thermistorThreeGrids;info=g.summary;data.midbodyGridComparison=info;
    if info.available
        q=g.heat;scope=sprintf('3 m/s 带线试件，三档均采用 20D 外域，热边界 %.7g %s。',info.boundaryValue,info.boundaryUnit);
        method='同工况强等温标量的壁反力按原面积守恒分配';
        if strcmp(info.mode,'givenPower'),method='以每档实算热导反求等温壁温，按温差比换算该网原面壁反力';
        elseif info.boundaryValue~=305.15,method='各档 305.15 K 参照场的原面壁反力按实际给定温差比换算';end
        note=sprintf('%s 三档壁温 %.6f / %.6f / %.6f K。',scope,info.grids.wallK);
        tab={'gridHeat1','入壁热流 · 第一档','kW/m²',q(:,1)/1000,method; ...
             'gridHeat2','入壁热流 · 第二档','kW/m²',q(:,2)/1000,method; ...
             'gridHeat21','入壁热流差 · 第二档减第一档','W/m²',q(:,2)-q(:,1),[method '，第二档减第一档']; ...
             'gridHeat32','入壁热流差 · 第三档减第二档','W/m²',q(:,3)-q(:,2),[method '，第三档减第二档']};
        absolute=[min(q(:)) max(q(:))]/1000;bound=max(abs(diff(q,1,2)),[],'all');
        for j=1:size(tab,1)
            key=tab{j,1};values=tab{j,4};
            add(key,tab{j,2},tab{j,3},values,ones(n,1),repmat(tab(j,5),n,1),scope);
            data.fields(end).stateNote=note;
            if all(values==values(1))
                data.uniformFieldValues.(key)=struct('value',values(1),'unit',tab{j,3});
            elseif j<=2,data.comparisonColorRanges.(key)=absolute;
            else,data.comparisonColorRanges.(key)=[-bound bound];end
        end
        if absolute(2)>absolute(1),data.comparisonColorRanges.heatFlux=absolute;end
    else
        k=find(strcmp({data.fields.key},'heatFlux'),1);
        data.fields(k).stateNote=[data.fields(k).stateNote ' 当前功率超过部分参照网的支持范围。'];
    end
end
% 分带收支的云图取功率差除以带原面面积，积分可以恢复该带功率。正号表示向外。
if isfield(mesh,'thermistorLayerBudget')
    b=mesh.thermistorLayerBudget;data.layerBudget=rmfield(b,'faceBand');
    bandIndex=b.faceBand;bands=b.bands;deltas=[bands.delta];areas=[bands.areaM2];
    scopes=cell(numel(bands),1);
    for k=1:numel(bands)
        d=deltas(k);
        scopes{k}=sprintf(['%s · 第三档减第二档，305.15 K。全带放热变化 %+.3f µW，' ...
            '外侧对流 %+.3f、扩散 %+.3f、带间净流出 %+.3f、映射修正 %+.3f、残差 %+.5f µW。'], ...
            bands(k).name,1e6*d.faceOutW,1e6*d.outerConvW,1e6*d.outerDiffW,1e6*(d.sideConvW+d.sideDiffW), ...
            1e6*d.mapCorrectionW,1e6*d.residualW);
    end
    scope='3 m/s 带线热敏珠；同一带内显示面积平均，正值表示向外传热增加。';
    rows={'bandRelease','分带平均 · 壁面放热变化',[deltas.faceOutW],'原面放热按带积分'; ...
          'bandConvection','分带平均 · 外侧对流变化',[deltas.outerConvW],'共同近壁控制体外侧的原生对流通量'; ...
          'bandDiffusion','分带平均 · 外侧扩散变化',[deltas.outerDiffW],'共同近壁控制体外侧的原生扩散通量'; ...
          'bandSide','分带平均 · 带间净流出变化',[deltas.sideConvW]+[deltas.sideDiffW],'共同近壁控制体沿相邻带界的净通量'};
    bound=0;
    for j=1:size(rows,1)
        density=rows{j,3}./areas;bound=max(bound,max(abs(density)));
        method=[rows{j,4} '，第三档减第二档后除以本带原面面积'];
        add(rows{j,1},rows{j,2},'W/m²',density(bandIndex),ones(n,1),repmat({method},n,1),scope,scopes(bandIndex));
        data.fields(end).stateNote=scope;
    end
    for j=1:size(rows,1),data.comparisonColorRanges.(rows{j,1})=[-bound bound];end
end
% 全周启动步使用独立的时间、算法和储热账，供云图、截面、探针与离线导出共用。
if isfield(mesh,'thermistorStartup')
    s=mesh.thermistorStartup;e=s.energy;c=mesh.hybridCondition;
    data.thermistorStartup=s;
    data.coordinateNote='SI · 长轴 X · 来流 −X · 正热流入壁';
    data.caseNote=sprintf('%s · 距构造初场 %.8g µs · 给定壁温 %.6g K · 360° 全周', ...
        data.caseNote,1e6*c.elapsedPhysicalTimeS,c.wallTemperatureK);
    data.resultCondition.caption=data.caseNote;
    data.resultCondition.elapsedPhysicalTimeS=c.elapsedPhysicalTimeS;
    data.resultCondition.timeCaption=sprintf('距构造初场 %.8g µs',1e6*c.elapsedPhysicalTimeS);
    k=find(strcmp({data.fields.key},'heatFlux'),1);
    data.fields(k).stateNote=sprintf(['启动瞬态原面结果。珠面放热 %.8g mW，远场流出 %.8g mW，' ...
        '空气储热率 %.6g nW，收支差 %.4g pW。'],1e3*e.wallReleaseW,1e3*e.farOutW,1e9*e.airStorageW,1e12*e.imbalanceW);
    k=find(strcmp({data.fields.key},'nearWallCellTemperature'),1);
    data.fields(k).label='首层流体节点温度';
    data.fields(k).stateNote='终场的流体温度，取原面所对应的全周首层节点。';
    data.uniformFieldValues.wallInput=struct('value',c.wallTemperatureK,'unit','K');
end
% 模块用途：按对照类型指定同刻两份实际场、差值方向及温度取样位置。
if isfield(mesh,'beadStepComparison')
    trial=mesh.beadStepComparison;s=trial.summary;w=trial.wall;
    data.defaultField='timeHeatDifference';
    data.caseNote=[data.caseNote ' · 主值为一个整步，对照为两个半步'];
    data.resultCondition.caption=data.caseNote;
    data.timeStepComparison=s;
    scope='同一全周网格与物理时刻。比较时间步长对启动瞬态的影响，逐面实验精度待验证。';
    fields={
        'timeHalfPressure','绝对壁压 · 两个半步','kPa',w.pressurePa/1000,'全周层流 BDF1，终帧壁压面内积分并加环境压力','pressure';
        'timeHalfGauge','相对环境壁压 · 两个半步','Pa',w.pressureGaugePa,'全周层流 BDF1，终帧表压按原三角面内积分平均','pressureGauge';
        'timeHalfHeat','入壁热流 · 两个半步','kW/m²',w.heatFluxWm2/1000,'全周 BDF1 强温差标量固定端点壁反力，按原面面积保守分配','heatFlux';
        'timeHalfNearWall','首层流体温度 · 两个半步','K',w.nearWallCellTemperatureK,'全周 BDF1 终帧首层节点标量换算为 K 后按原面平均','nearWallCellTemperature'};
    differenceMethod='同一原面、同一物理时刻，两个 BDF1 半步减一个 BDF1 整步';
    differenceLabel='时间步对照';
    delta={'timeHeatDifference','入壁热流差','W/m²',trial.heatDifferenceWm2;
           'timePressureDifference','壁压差','Pa',trial.pressureDifferencePa};
elseif isfield(mesh,'beadNearWallComparison')
    trial=mesh.beadNearWallComparison;s=trial.summary;r=trial.reference;w=mesh.hybridWall;
    data.nearWallComparison=s;data.defaultField='gridHeatDifference';
    data.caseNote=[data.caseNote ' · 主值为 56 层，对照为 48 层 · 两网各推进两个半步'];
    data.resultCondition.caption=data.caseNote;
    scope='同一启动时刻和原面。近壁最内八层等分，温度差在原第一层的相同物理位置计算。';
    fields={
        'gridRefPressure','绝对壁压 · 48 层','kPa',r.pressurePa/1000,'48 层全周层流 BDF1，壁节点压力按原面积分并加环境压力','pressure';
        'gridRefGauge','相对环境壁压 · 48 层','Pa',r.pressureGaugePa,'48 层全周层流 BDF1，表压按原面积分','pressureGauge';
        'gridRefHeat','入壁热流 · 48 层','kW/m²',r.heatFluxWm2/1000,'48 层 BDF1 强温差标量壁反力，按原面保守分配','heatFlux';
        'gridMatchedTemperature','共同位置温度 · 56 层','K',trial.matchedTemperatureK,'56 层第二层节点的终场标量换算为温度，按原面权重平均','';
        'gridRefTemperature','共同位置温度 · 48 层','K',r.nearWallCellTemperatureK,'48 层第一层节点的终场标量换算为温度，按相同原面权重平均','gridMatchedTemperature'};
    differenceMethod='同一原面、同一物理时刻，56 层减 48 层';
    differenceLabel='近壁网格对照';
    delta={'gridHeatDifference','入壁热流差','W/m²',w.heatFluxWm2-r.heatFluxWm2;
           'gridPressureDifference','壁压差','Pa',w.pressureGaugePa-r.pressureGaugePa;
           'gridTemperatureDifference','共同位置温度差','K',trial.matchedTemperatureK-r.nearWallCellTemperatureK};
    k=find(strcmp({data.fields.key},'nearWallCellTemperature'),1);
    data.fields(k).label='首层流体温度 · 56 层';
    data.fields(k).stateNote='加密网实际首层。共同位置对照使用加密网第二层及原网第一层。';
end
% 主场和对照场使用共同比较色标，各差场使用关于零对称的比较色标。
if isfield(mesh,'beadStepComparison') || isfield(mesh,'beadNearWallComparison')
    for j=1:size(fields,1)
        values=fields{j,4};key=fields{j,1};ref=fields{j,6};
        add(key,fields{j,2},fields{j,3},values,ones(n,1),repmat(fields(j,5),n,1),scope);
        data.fields(end).stateNote=scope;data.fields(end).computedStateLabel='研究数值解';
        if isempty(ref),continue;end
        pair=[values;data.fields(strcmp({data.fields.key},ref)).values];
        data.comparisonColorRanges.(key)=[min(pair) max(pair)];
        data.comparisonColorRanges.(ref)=data.comparisonColorRanges.(key);
    end
    methods=repmat({differenceMethod},n,1);
    for j=1:size(delta,1)
        key=delta{j,1};values=delta{j,4};
        add(key,delta{j,2},delta{j,3},values,ones(n,1),methods,scope);
        data.fields(end).stateNote=[differenceMethod '。' scope];data.fields(end).computedStateLabel=differenceLabel;
        limit=max(abs(values));if limit>0,data.comparisonColorRanges.(key)=[-limit limit];end
    end
end
if isfield(mesh,'beadNearWallComparison')
    % 在对应差场的说明中给出量纲、归一化参照及整珠功率，便于读者核对局部图。
    k=find(strcmp({data.fields.key},'gridHeatDifference'),1);
    data.fields(k).stateNote=sprintf(['56 层减 48 层。面积加权 RMS %.5g W/m²，局部最大绝对差 %.5g W/m²' ...
        '（占 48 层峰值 %.3g%%）。整珠放热 %.7g → %.7g mW，变化 %+.4g%%。'], ...
        s.heatFluxWm2.weightedRMS,s.heatFluxWm2.maximumAbsolute,100*s.heatFluxWm2.maximumFraction, ...
        1000*s.referenceWallReleaseW,1000*s.fineSurface.energy.wallReleaseW,100*s.relativePowerDifference);
    k=find(strcmp({data.fields.key},'gridPressureDifference'),1);
    data.fields(k).stateNote=sprintf(['56 层减 48 层。面积加权 RMS %.5g Pa，局部最大绝对差 %.5g Pa。' ...
        '比较用动压 %.5g Pa。'],s.pressurePa.weightedRMS,s.pressurePa.maximumAbsolute,s.pressurePa.scale);
end
if isfield(mesh,'lewisCylinder')
    data.defaultField='nusselt';data.flowDirection=[0 -1 0];
    data.coordinateNote='SI · X 为圆柱跨向 · +Y 为迎流驻点 · +Z 为正圆周角';
    data.caseNote='NASA Lewis 光滑圆柱 · Re_D = 177000 · 低湍流空气 · 直径 66 mm';
    data.resultCondition.caption=data.caseNote;
    data.resultCondition.physics=rmfield(mesh.lewisCylinder,'nusselt');
    if isfield(mesh,'lewisCondition')
        c=mesh.lewisCondition;data.defaultField='heatFlux';
        data.caseNote=sprintf('给定边界的 66 mm 圆柱 · Re_D = 177000 · 来流 %.5g m/s · 上游 %.7g Pa、%.6g K · 给定壁温 %.6g K', ...
            c.velocityMps,c.referencePressurePa,c.airTemperatureK,c.wallTemperatureK);
        data.resultCondition.physics=c;data.resultCondition.caption=data.caseNote;
        data.resultCondition.freestreamPressurePa=c.referencePressurePa;
        data.resultCondition.freestreamTemperatureK=c.airTemperatureK;
        data.resultCondition.wallInputType='uniform';data.resultCondition.wallUniformValueK=c.wallTemperatureK;
        data.resultCondition.wallMinK=c.wallTemperatureK;data.resultCondition.wallMaxK=c.wallTemperatureK;
        data.resultCondition.wallKnownFaces=n;data.resultCondition.wallUnknownFaces=0;
    end
end
% 静态三角面快照按当前两列场重新统计。旧存档可直接导出，历史帧单独显示。
if isempty(histories)
    pf=data.fields(strcmp({data.fields.key},'pressure'));
    qf=data.fields(strcmp({data.fields.key},'heatFlux'));
    if isscalar(pf) && isscalar(qf)
        % 上述标准字段固定使用 kPa 和 kW/m²，功率汇总的输入单位为 W/m²。
        data.wallCoverage=summarizeComputedWall(area,pf.status==1,qf.status==1,1000*qf.values);
        data.caseNote=[data.wallCoverage.coverageMessage ' ' data.caseNote];
    end
end
% 仅传递明确保存的来源/证据；旧快照不凭方法名自动获得物理资格。
if isfield(mesh,'evidenceIdentity')
    if isfield(mesh,'thermistorImport')
        % 逐面导出保留文件身份和对应标准面，避免把本机文件路径写入离线页。
        b=mesh.thermistorImport;
        data.thermistorImport=struct('referenceFaceID',b.referenceFaceID,'referenceAreaM2',b.referenceAreaM2, ...
            'stlSHA256',b.source.sha256,'scaleToM',b.source.scale,'meshSHA256',mesh.evidenceIdentity.meshSHA256, ...
            'rotationMatrix',b.rotationMatrix,'rotationCenterM',b.rotationCenterM,'translationM',b.translationM);
    end
    data.evidenceIdentity=mesh.evidenceIdentity;
    % 对已明确标注研究资格的快照，曲线和 CSV 共用同一来源。
    if isfield(mesh.evidenceIdentity,'productionQualified') && isequal(mesh.evidenceIdentity.productionQualified,false)
        data.researchIdentity=mesh.evidenceIdentity;
    end
end
if isfield(mesh,'evidenceRecords'),data.evidenceRecords=mesh.evidenceRecords;end
    function add(key,label,unit,values,state,methods,scope,faceScopes,faceReasons)
        if nargin<8,faceScopes={};end
        if nargin<9,faceReasons={};end
        values=double(values(:)); state=double(state(:));sourceState=state;
        assert(numel(values)==n && numel(state)==n && numel(methods)==n, ...
            'RocketAeroHeat:Viewer:Size','字段 %s 与面元数量不一致。',key);
        state(~isfinite(values))=0; values(state==0)=NaN;
        if isempty(faceReasons)
            faceReasons=repmat({'原结果未记录该字段的未解原因；不以其他字段原因代替'},n,1);
            if strcmp(key,'wallSolved')
                faceReasons(:)={'当前单点快照未提供独立壁温解；未执行壁温求解，不表示热流全部未解'};
            elseif strcmp(key,'wallInput')
                faceReasons(:)={'该面元热流边界输入壁温未知；不由初始壁温补齐，也不是壁温求解失败'};
            elseif any(strcmp(key,{'pressure','cp'}))
                faceReasons(:)={'本面元表面压力结果未被接受；给定边界与区域计算分别审核'};
            end
        end
        assert(numel(faceReasons)==n,'RocketAeroHeat:Viewer:Size','字段未解原因与面元数不一致。');
        faceReasons=faceReasons(:);faceReasons(state>0)={''};
        for i=find(state==0)'
            if isempty(faceReasons{i}) || strcmp(faceReasons{i},'已计算')
                if sourceState(i)>0
                    faceReasons{i}='源状态标为有效但本字段没有有限值；不能由其他已求场量补齐';
                else
                    faceReasons{i}='本字段未解；原结果没有记录具体拒绝原因';
                end
            end
        end
        % 缺失原因同时进入逐面元范围；既有PNG/CSV来源列即可保留同一拒绝解释。
        if any(state==0)
            if isempty(faceScopes),faceScopes=repmat({scope},n,1);else,faceScopes=faceScopes(:);end
            for i=find(state==0)'
                faceScopes{i}=[faceScopes{i} '；未解原因：' faceReasons{i}];
            end
        end
        % 模块七：缺失字段保留原尝试模型，但不能把其名称当成已接受求解。
        methods=methods(:);
        for i=find(state==0)'
            if ~startsWith(methods{i},{'未解','未求解','未计算','来源未记录'})
                methods{i}=['未接受的尝试模型：' methods{i}];
            end
        end
        [names,~,ids]=unique(methods,'stable');
        scopeNames={};scopeIDs=[];
        if ~isempty(faceScopes)
            assert(numel(faceScopes)==n,'RocketAeroHeat:Viewer:Size', ...
                '字段 %s 的逐面元适用范围与面元数量不一致。',key);
            [scopeNames,~,scopeIDs]=unique(faceScopes(:),'stable');
            scopeIDs=scopeIDs-1;
        end
        entry=struct('key',key,'label',label,'unit',unit,'values',values, ...
            'status',state,'methodNames',{names},'methodIndex',ids-1,'scope',scope, ...
            'scopeNames',{scopeNames},'scopeIndex',scopeIDs, ...
            'missingReason',{faceReasons},'stateNote',fieldStateNote(key),'availability','', ...
            'computedStateLabel','', ...
            'computedAreaFraction',sum(area(state==1))/sum(area), ...
            'prescribedAreaFraction',sum(area(state==2))/sum(area), ...
            'missingAreaFraction',sum(area(state==0))/sum(area),'history',[]);
        if isempty(data.fields),data.fields=entry;else,data.fields(end+1)=entry;end
    end
end
function names=hybridMethod(w,key)
% 模块用途：从本次混合求解的原面表读取当前变量的算法。
if ismember(key,w.Properties.VariableNames),v=w.(key);else,v=w.method;end
names=cellstr(v);names(strlength(v)==0)={'未求解'};
end

function reasons=thermalMissingReasons(mesh,original,state,n)
% 模块八：区分本面外缘未知与未知上游传播；只解释已有拒绝，不改变有效掩码。
reasons=original(:);missing=state==0;
if ~isfield(mesh,'heatingCoverage')
    reasons(:)={'热边界没有有效解；原快照未记录具体拒绝原因'};
end
if isfield(mesh,'continuumRegime')
    c=mesh.continuumRegime;
    names={'downstreamStateUnresolvedFace','directStateUnresolvedFace'};
    text={'上游外缘热力状态未知/域外占位影响；即使本面外缘有效，也不能跨过未知上游恢复热流', ...
        '本面元外缘热力状态未知/域外占位；此项不是由有效外缘Kn超限推得，壁温输入状态另列'};
    for k=1:numel(names)
        if isfield(c,names{k}) && numel(c.(names{k}))==n
            mask=logical(c.(names{k})(:)) & missing;reasons(mask)={text{k}};
        end
    end
    if all(isfield(c,names)) && numel(c.(names{1}))==n && numel(c.(names{2}))==n
        inherited=missing & logical(c.downstreamStateUnresolvedFace(:)) & ~logical(c.directStateUnresolvedFace(:));
        if any(inherited)
            reasons=addSavedUpstreamReason(mesh,reasons,inherited,c);
        end
    end
end
end
function reasons=addSavedUpstreamReason(mesh,reasons,inherited,regime)
% 模块用途：沿已保存的同源前驱链说明缺失位置，供探针和导出共同使用。
if ~isfield(mesh,'surfaceStreamlines') || ~isfield(mesh,'boundaryLayerEdge'),return;end
path=mesh.surfaceStreamlines;edge=mesh.boundaryLayerEdge;
if ~all(isfield(path,{'predecessorFace','pathCost','connectedFace','velocitySource'})) || ...
        ~isfield(edge,'method') || ~strcmp(path.velocitySource,edge.method) || ...
        ~isscalar(string(edge.method)) || strlength(strtrim(string(edge.method)))==0,return;end
n=numel(reasons);
parent=double(path.predecessorFace(:));cost=path.pathCost(:);connected=logical(path.connectedFace(:));
if numel(parent)~=n || numel(cost)~=n || numel(connected)~=n,return;end
direct=logical(regime.directStateUnresolvedFace(:));
downstream=logical(regime.downstreamStateUnresolvedFace(:));
% 前驱成本必须严格递减。断链、环和来源不明的分支自然保留原说明。
source=zeros(n,1);links=zeros(n,1);
known=find(direct & connected & isfinite(cost));source(known)=known;
ids=find(downstream & ~direct & connected & isfinite(cost));
[~,order]=sort(cost(ids));
for face=reshape(ids(order),1,[])
    p=parent(face);
    if ~isfinite(p) || p<1 || p>n || p~=fix(p) || cost(p)>=cost(face) || source(p)==0,continue;end
    source(face)=source(p);links(face)=links(p)+1;
    if inherited(face)
        reasons{face}=sprintf('%s。沿保存的准流线回溯 %d 步，最近的外缘状态直接缺失原面为 %d。', ...
            reasons{face},links(face),source(face));
    end
end
end
function note=fieldStateNote(key)
% 模块九：当前字段的输入/求解身份，避免把热流覆盖率当作所有场量覆盖率。
switch key
    case {'pressure','cp'}
        note='表面载荷压力：修正牛顿、给定Cp或区域关联；与边界层外缘压力分开。';
    case 'edgePressure'
        note='内部外流外缘静压；不是修正牛顿/给定Cp表面压力。局部外缘有效不保证沿程热流可解。';
    case 'edgeTemperature'
        note='内部外流外缘静温；不是给定壁温、恢复温度或求解壁温。';
    case 'heatFlux'
        note='当前单点对流热流；本面外缘未知与未知上游传播分别阻断，不以占位或零值补齐。';
    case 'recoveryTemperature'
        note='绝热恢复温度；当前热边界有效面元的温度模型结果，不是给定或求解壁温。';
    case 'wallInput'
        note='热流边界输入；给定面积100%也不表示求解了壁体温度场。';
    case 'wallSolved'
        note='需要独立壁温求解及边界条件；未执行壁温求解不表示对流热流全部未解。';
    otherwise
        note='状态与覆盖率仅对应当前变量，不代表其他场量均已求解。';
end
end
function t=solvedBeadThermal(mesh)
% 从实际试件取得集总热平衡记录，供冻结工况和壁温探针共用。
if isfield(mesh,'thermistorLead'),t=mesh.thermistorLead.thermal;
else,t=mesh.thermistorBead.thermal;end
end

function result=snapshotCondition(mesh)
% 模块十：可复核工况元数据；缺失记录明确未知，均匀/逐面元身份取实际输入数组。
result=struct('source','结果快照；未读取待编辑cfg','Mach',NaN,'alphaDeg',NaN, ...
    'betaDeg',NaN,'altitudeM',NaN,'freestreamPressurePa',NaN,'freestreamTemperatureK',NaN, ...
    'wallInputType','notRecorded','wallUniformValueK',NaN,'wallMinK',NaN,'wallMaxK',NaN, ...
    'wallKnownFaces',0,'wallUnknownFaces',size(mesh.F,1),'wallIdentity','热流边界输入；不是本次求解壁温');
if isfield(mesh,'pressure') && isfield(mesh.pressure,'freestream')
    fs=mesh.pressure.freestream;
    from={'Mach','alpha','beta','altitude','pressure','temperature'};
    to={'Mach','alphaDeg','betaDeg','altitudeM','freestreamPressurePa','freestreamTemperatureK'};
    for k=1:numel(from)
        if isfield(fs,from{k}) && isnumeric(fs.(from{k})) && isscalar(fs.(from{k})) && isfinite(fs.(from{k}))
            result.(to{k})=double(fs.(from{k}));
        end
    end
end
if isfield(mesh,'hybridSolvedWallTemperatureK')
    % 只在功率模式登记求解温度；给定壁温的计数及范围维持无输入状态。
    t=solvedBeadThermal(mesh);
    result.wallInputType='netHeatingPower';result.wallIdentity=t.wallMethod;
    result.netHeatingPowerW=t.netHeatingPowerW;result.solvedWallTemperatureK=t.wallTemperatureK;
    result.caption=sprintf('净加热功率 %.12g W，稳态求解壁温 %.12g K',t.netHeatingPowerW,t.wallTemperatureK);
    return;
end
wallText='壁温输入类型/数值未记录';
if isfield(mesh,'hybridWallTemperatureK') || (isfield(mesh,'heating') && isfield(mesh.heating,'wallTemperatureByFace'))
    if isfield(mesh,'hybridWallTemperatureK'),values=double(mesh.hybridWallTemperatureK(:));
    else,values=double(mesh.heating.wallTemperatureByFace(:));end
    known=isfinite(values);
    result.wallKnownFaces=nnz(known);result.wallUnknownFaces=nnz(~known);
    if any(known)
        result.wallMinK=min(values(known));result.wallMaxK=max(values(known));
        if all(known) && all(values==values(1))
            result.wallInputType='uniform';result.wallUniformValueK=values(1);
            wallText=sprintf('统一壁温边界输入 %.6g K',values(1));
        else
            result.wallInputType='perFace';
            wallText=sprintf('逐面元壁温边界输入 %.6g–%.6g K；未知 %d/%d 面', ...
                result.wallMinK,result.wallMaxK,result.wallUnknownFaces,numel(values));
        end
    else
        result.wallInputType='unknown';wallText=sprintf('壁温边界输入全部未知（%d 面）',numel(values));
    end
end
result.caption=sprintf(['单点基准快照：M %s，α %s°，β %s°，高度 %s km；' ...
    'p∞ %s Pa，T∞ %s K；%s（输入，不是壁温解）'], ...
    recorded(result.Mach),recorded(result.alphaDeg),recorded(result.betaDeg), ...
    recorded(result.altitudeM/1000),recorded(result.freestreamPressurePa), ...
    recorded(result.freestreamTemperatureK),wallText);
end
function text=recorded(value)
if isfinite(value),text=sprintf('%.8g',value);else,text='未记录';end
end
function [names,scopes]=edgeTemperatureProvenance(e,n)
% 模块八：温度独立追溯到外流静态状态和逐面元热物性；不重算静温。
names=methodNames(e,n);
summary='内部工程外流静温；不是壁温；代数有效不等于已通过实验验证';
if isfield(e,'scope'),summary=[summary '；' char(e.scope)];end
if isfield(e,'assumption'),summary=[summary '；' char(e.assumption)];end
scopes=repmat({summary},n,1);
if isfield(e,'thermalPropertyModeByFace') && isfield(e,'thermalPropertyModeNames')
    ids=double(e.thermalPropertyModeByFace(:))+1;
    assert(numel(ids)==n,'RocketAeroHeat:Viewer:Size','外缘热物性编号与面元数量不一致。');
    for i=1:n
        if ids(i)>=1 && ids(i)<=numel(e.thermalPropertyModeNames) && ids(i)==fix(ids(i))
            property=char(e.thermalPropertyModeNames{ids(i)});
            names{i}=[names{i} '；外缘静温热物性：' property];
        else
            names{i}=[names{i} '；外缘静温热物性编号无有效名称'];
        end
    end
end
if isfield(e,'thermalExpansionOutOfRangeFace') && ~isempty(e.thermalExpansionOutOfRangeFace)
    mask=logical(e.thermalExpansionOutOfRangeFace(:));
    assert(numel(mask)==n,'RocketAeroHeat:Viewer:Size','外缘低温适用标志与面元数量不一致。');
    for i=find(mask)'
        scopes{i}=[scopes{i} '；本面元低于焓模型适用域，定比热外缘只用于流线占位，不能据此宣称热流已解'];
    end
end
end
function [names,scopes]=recoveryTemperatureProvenance(mesh,h,n)
% 模块九：恢复温度与热流具有不同来源；只解码实际结果，不作新物理推断。
% 层流/湍流由已保存的热边界编号选择恢复因子；球鼻沿用驻点总温，
% 此处明确它是球鼻低阶延拓所用统一温度，不能说是局部球面温度解。
names=repmat({'恢复温度原始结果；恢复因子/热物性细节未记录'},n,1);
scopes=repmat({'只保留原计算值；原结果缺少恢复温度来源细节，不能据热流方法反推'},n,1);
if ~isfield(h,'faceMethodID'),return;end
ids=double(h.faceMethodID(:));
assert(numel(ids)==n,'RocketAeroHeat:Viewer:Size','热边界方法编号与面元数量不一致。');
mode='未记录';
if isfield(h,'thermalPropertyMethod'),mode=char(h.thermalPropertyMethod);end
propertyScope='恢复温度热物性适用范围未记录';
if isfield(h,'thermalPropertyScope'),propertyScope=char(h.thermalPropertyScope);end
if strcmp(mode,'thermallyPerfect')
    formula='干空气焓恢复：h(Taw)=h(Te)+r*Ue²/2；NASA7/低温振动冻结反解温度';
elseif strcmp(mode,'caloricallyPerfect')
    formula='定比热恢复：Taw=Te+r*Ue²/(2*cp)，cp=gamma*R/(gamma-1)';
else
    formula='恢复温度结果；热物性模式未记录（不推断焓恢复或定比热）';
end
for i=1:n
    if ismember(ids(i),[1 2 4 5 7 8 9 10])
        if ismember(ids(i),[1 4 7 9]),factor='层流 r=√Pr';else,factor='湍流 r=Pr^(1/3)';end
        if isfield(h,'Prandtl')
            factor=[factor sprintf('，Pr=%.6g',h.Prandtl)];
        else
            factor=[factor '（Pr 数值未记录）'];
        end
        names{i}=[formula '；' factor];
        scopes{i}=['附着边界层绝热恢复温度的工程假设；仅输出当前热边界有效面元；' propertyScope];
        if ismember(ids(i),[4 5 9 10])
            scopes{i}=[scopes{i} '；外缘来自受限前体脱体工程闭合，其压力/速度及恢复温度均为低阶估算'];
        end
    elseif ids(i)==3
        names{i}='球鼻热边界恢复温度：沿用驻点绝热总温 T0（r=1）';
        if isfield(mesh,'stagnationHeating')
            s=mesh.stagnationHeating;
            if isfield(s,'thermalPropertyMode') && strcmp(s.thermalPropertyMode,'thermallyPerfect')
                names{i}=[names{i} '；干空气总焓反解温度'];
            elseif isfield(s,'thermalPropertyMode') && strcmp(s.thermalPropertyMode,'caloricallyPerfect')
                names{i}=[names{i} '；定比热总温'];
            end
            if isfield(s,'shockStateMethod'),names{i}=[names{i} '；' char(s.shockStateMethod)];end
        end
        scopes{i}=['已验证球鼻热边界的统一驻点温度假设；不是各面元局部球面恢复温度求解；' ...
            '热流几何延拓不改变该温度；' propertyScope];
        if isfield(h,'noseScope'),scopes{i}=[scopes{i} '；' char(h.noseScope)];end
    elseif ids(i)==6 && isfield(h,'noseDistribution')
        nose=h.noseDistribution;
        names{i}=nose.recoveryTemperatureMethod;
        scopes{i}=[nose.scope '；同一局部相似剖面的绝热热边界；不是输入壁温或统一驻点总温'];
    elseif ids(i)==0
        names{i}='未解：当前无有效恢复温度';
        scopes{i}='热边界未计算或超出适用范围；保留 NaN，不以输入壁温或总温补齐';
    end
end
end
function [names,scopes]=baselineHeatingProvenance(h,n)
% 原方法对照逐面恢复替换前来源。旧草稿缺少完整来源时拒绝，不能猜测。
k=h.kinneyAcreage;
assert(isfield(k,'baselineProvenance') && isfield(k,'baselineFaceMethodID'), ...
    'RocketAeroHeat:Viewer:BaselineProvenance', ...
    '原方法对照缺少替换前来源，请用原基线重新生成方法对照，不能从新方法推断。');
original=h;
for key=fieldnames(k.baselineProvenance)'
    original.(key{1})=k.baselineProvenance.(key{1});
end
original.faceMethodID=k.baselineFaceMethodID;
ids=double(original.faceMethodID(:));
assert(numel(ids)==n && all(isfinite(ids) & ids==fix(ids) & ids>=0 & ids<=6), ...
    'RocketAeroHeat:Viewer:BaselineProvenance','原方法编号必须逐面保存且属于原工程方法。');
assert(isequal(ids>0,isfinite(k.baselineHeatFlux(:))), ...
    'RocketAeroHeat:Viewer:BaselineProvenance','原方法编号与热流有值状态不一致。');
names=heatingMethods(original,n);
scopes=repmat({'原工程方法未解；不由新的方法或其他面元补齐'},n,1);
for i=find(ids>0)'
    if ismember(ids(i),[3 6])
        scopes{i}=original.noseScope;
    elseif ismember(ids(i),[1 4])
        scopes{i}=original.boundaryLayerValidity;
    else
        scopes{i}=original.turbulentScope;
    end
    if ismember(ids(i),[4 5])
        scopes{i}=[scopes{i} '；外缘来自受限前体脱体工程闭合，未经完整三维激波层验证'];
    end
    scopes{i}=[scopes{i} '；同一工况原算法对照，不是独立实验验证'];
end
end
function names=heatingMethods(h,n)
% 模块六：展开真实的局部边界层来源，区分 Head 积分和显式平板回退。
names=methodNames(h,n);
if ~isfield(h,'boundaryLayerMethod'),return;end
ids=double(h.faceMethodID(:));
for i=1:n
    if ismember(ids(i),[1 4])
        names{i}=[names{i} '；' h.boundaryLayerMethod '；Eckert 参考性质 + Reynolds 比拟'];
    elseif ismember(ids(i),[2 5])
        local=h.turbulentMethod;
        if isfield(h,'turbulentBoundaryLayer') && isstruct(h.turbulentBoundaryLayer) && ...
                isfield(h.turbulentBoundaryLayer,'sourceID')
            t=h.turbulentBoundaryLayer;
            local=[local '；本面元：' t.sourceNames{double(t.sourceID(i))+1}];
        end
        names{i}=[names{i} '；' local '；Eckert 参考性质 + Reynolds 比拟'];
    elseif ids(i)==3
        names{i}=[names{i} '；' h.stagnationMethod];
    elseif ids(i)==6
        names{i}=[names{i} '；' h.noseDistribution.heatFluxMethod];
    elseif ismember(ids(i),[7 8 9 10]) && isfield(h,'kinneyAcreage')
        k=h.kinneyAcreage;
        names{i}=[names{i} '；' k.sourceURL '；' k.sourcePages ...
            '；Eckert参考性质 + Reynolds比拟'];
        if k.length.prefixAssumedFace(i)
            names{i}=[names{i} '；入口前缀使用局部常状态工程假设'];
        elseif k.length.originAssumedFace(i)
            names{i}=[names{i} '；离散起点I=0假设，真实驻点/附着线入口尚未核验'];
        end
    end
end
end
function names=methodNames(result,n)
% 模块七：将原始逐面元算法编号解码；保留原名称与整体算法回退记录。
names=repmat({'来源未记录'},n,1);
if isfield(result,'method'),names(:)={char(result.method)};end
if isfield(result,'faceMethodID') && isfield(result,'faceMethodNames')
    ids=double(result.faceMethodID(:))+1;
    valid=ids>=1 & ids<=numel(result.faceMethodNames);
    names(valid)=reshape(result.faceMethodNames(ids(valid)),[],1);
end
% 将外流转角的几何来源也带入逐面元拾取；轴向推进单独注明。
if isfield(result,'normalReconstruction')
    frame=result.normalReconstruction;
    labels=repmat({'三角面法向入射角'},n,1);
    labels(frame.appliedFace)={'校验圆截面/直母线重构入射角'};
    if isfield(result,'faceMethodID')
        labels(result.faceMethodID==6)={'原生截面轴向波系推进'};
    end
    names=cellstr(string(names)+"；"+string(labels));
end
end
