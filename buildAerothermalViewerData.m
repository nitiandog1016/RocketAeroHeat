function data = buildAerothermalViewerData(mesh)
%BUILDAEROTHERMALVIEWERDATA 将真实计算结果转换成离线三维工作台数据。
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
if isfield(mesh,'boundaryLayerEdge')
    e=mesh.boundaryLayerEdge; methods=methodNames(e,n);
    state=double(e.algebraicallyValid(:));
    edgePressure=e.pressure;edgeTemperature=e.temperature;
    % 数值占位可用于诊断流向，但不是压力/温度场的有效计算结果。
    placeholder=false(n,1);
    if isfield(e,'thermalExpansionOutOfRangeFace')
        placeholder=logical(e.thermalExpansionOutOfRangeFace(:));
    end
    if isfield(e,'thermalStateValidFace')
        placeholder=placeholder | (logical(e.algebraicallyValid(:)) & ...
            ~logical(e.thermalStateValidFace(:)));
    end
    state(placeholder)=0;edgePressure(placeholder)=NaN;edgeTemperature(placeholder)=NaN;
    edgeReasons=repmat({'内部外流没有接受本面元外缘热力状态；没有有效外缘解'},n,1);
    edgeReasons(placeholder)={'本面元外缘热力状态未知/域外占位；定比热状态仅供流向诊断，不是有效外缘解'};
    pressureScopes=repmat({'内部工程外流外缘静压；代数有效不代表实验验证或区域接口已匹配'},n,1);
    [temperatureMethods,temperatureScopes]=edgeTemperatureProvenance(e,n);
    for i=find(placeholder)'
        methods{i}=['外缘压力未解：域外定比热状态仅供流向诊断；尝试模型：' methods{i}];
    end
    pressureScopes(placeholder)={'未取得所选热物性范围内的外缘压力；占位不计入已求压力面积'};
    % 局部相似鼻区使用自己的同压力/熵外缘。这里仅组合区域结果用于
    % 展示，不改写后体外流或流线，不把分区展示称作接口匹配已完成。
    if isfield(mesh,'heating') && isfield(mesh.heating,'noseModel') && ...
            strcmp(mesh.heating.noseModel,'sphericalLocalSimilarity')
        noseEdge=mesh.heating.noseDistribution.edgeState;
        sphere=noseEdge.geometryFace;keep=sphere & noseEdge.algebraicallyValid;
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
% 模块三：显式区分输入壁温与壁温解；单一外流/热流工况本身不确定壁温。
% 尚无材料热响应解时仍允许探针查看“未求解”，不能以输入壁温补齐。
add('wallSolved','求解壁温','K',nan(n,1),zeros(n,1), ...
    repmat({'未求解：当前未提供壁体热响应或其他独立壁温解'},n,1), ...
    '单点压力和对流热流计算使用壁温边界输入；求解壁温还需独立热平衡或材料热响应及边界条件');
data.fields(end).availability='notProvidedInSnapshot';
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
wallText='壁温输入类型/数值未记录';
if isfield(mesh,'heating') && isfield(mesh.heating,'wallTemperatureByFace')
    values=double(mesh.heating.wallTemperatureByFace(:));known=isfinite(values);
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
if isfield(e,'thermalExpansionOutOfRangeFace')
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
