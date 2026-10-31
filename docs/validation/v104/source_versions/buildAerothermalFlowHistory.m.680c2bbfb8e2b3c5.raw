function entries=buildAerothermalFlowHistory(mesh)
%BUILDAEROTHERMALFLOWHISTORY 将实际保存的外流快照对齐到接受的时间站。
% 模块一：耦合历程只读接受步元数据；局部反馈读原始计算站，均不插值。
% 外流是否有效独立于材料温度历程，不能因热壁失效抹掉有效压力。
entries=struct([]);n=size(mesh.F,1);stations={};t=[];
if isfield(mesh,'wallThermal') && isfield(mesh.wallThermal,'coupling')
    w=mesh.wallThermal;
    if ~isfield(w.coupling,'boundaryMetadata'),return;end
    meta=w.coupling.boundaryMetadata;t=double(w.time(:)');stations=cell(size(t));
    assert(numel(meta)==numel(t),'RocketAeroHeat:Viewer:FlowTime','接受步元数据数量错误。');
    for j=1:numel(t)
        if isempty(meta{j}) || ~isfield(meta{j},'flowSnapshot'),continue;end
        assert(meta{j}.time==t(j),'RocketAeroHeat:Viewer:FlowTime','外流快照不属于当前接受时刻。');
        stations{j}=meta{j}.flowSnapshot;
        checkCondition(stations{j},meta{j}.flight);
    end
elseif isfield(mesh,'flightHeatingHistory') && isfield(mesh.flightHeatingHistory,'flowStations')
    h=mesh.flightHeatingHistory;t=double(h.time(:)');stations=h.flowStations;
    assert(numel(stations)==numel(t),'RocketAeroHeat:Viewer:FlowTime','外流快照数量错误。');
    for j=1:numel(t)
        if isempty(stations{j}),continue;end
        condition=struct;
        for key={'Mach','altitude','alpha','beta'},condition.(key{1})=h.flight.(key{1})(j);end
        checkCondition(stations{j},condition);
    end
    if isfield(mesh,'wallThermal')
        assert(isequal(t,double(mesh.wallThermal.time(:)')), ...
            'RocketAeroHeat:Viewer:FlowTime','外流与壁温时间轴不同，不自动插值。');
    end
end
if isempty(t) || all(cellfun(@isempty,stations)),return;end
assert(all(isfinite(t)) && t(1)==0 && all(diff(t)>0), ...
    'RocketAeroHeat:Viewer:FlowTime','外流历史时间轴无效。');
% 保存快照时的精确网格用于拒绝重排/变形后误用旧历程。
if isfield(mesh,'flightHeatingHistory') && isfield(mesh.flightHeatingHistory,'couplingInput')
    input=mesh.flightHeatingHistory.couplingInput;
    assert(isequal(input.vertices,mesh.V) && isequal(input.faces,mesh.F), ...
        'RocketAeroHeat:Viewer:FlowMesh','外流历史与当前网格不一致。');
end

%% 2. 合并每站方法字典；缺失站单独标记，不能借用前/后站数值
nt=numel(t);keys={'pressure','cp','edgePressure','edgeTemperature'};
for k=1:numel(keys)
    values=nan(n,nt);status=zeros(n,nt);methods=repmat({'该时间站未保存此场量'},n,nt);
    notes=repmat({'该时间站外流快照缺失'},1,nt);scopes=notes;
    directions=nan(nt,3);sample=[];
    for j=1:nt
        s=stations{j};if isempty(s),continue;end
        assert(s.faceCount==n,'RocketAeroHeat:Viewer:FlowSize','快照与面元数量不一致。');
        index=find(strcmp({s.fields.key},keys{k}));
        if isempty(index),continue;end
        assert(isscalar(index),'RocketAeroHeat:Viewer:FlowSize','快照包含重复场量。');
        f=s.fields(index);sample=f;
        assert(numel(f.values)==n && numel(f.status)==n && numel(f.methodIndex)==n, ...
            'RocketAeroHeat:Viewer:FlowSize','快照字段尺寸错误。');
        ids=double(f.methodIndex(:));
        assert(all(isfinite(ids) & ids==fix(ids) & ids>=0 & ids<numel(f.methodNames)) && ...
            all(ismember(f.status,[0 1 2])),'RocketAeroHeat:Viewer:FlowSize','快照状态或方法编号错误。');
        values(:,j)=f.values;status(:,j)=f.status;
        methods(:,j)=reshape(f.methodNames(ids+1),[],1);
        notes{j}=s.timeNote;scopes{j}=f.scope;directions(j,:)=s.flowDirection;
    end
    if isempty(sample),continue;end
    status(~isfinite(values))=0;values(status==0)=NaN;
    [names,~,ids]=unique(methods(:),'stable');
    history=struct('time',t,'layout','time-major','values',values(:), ...
        'status',status(:),'methodNames',{names},'methodIndex',ids-1, ...
        'faceNote',{repmat({'外流独立于热壁有效掩码；给定底压/背风压力仍为给定条件'},n,1)}, ...
        'timeNote',{notes},'timeScope',{scopes},'flowDirections',directions(:));
    entry=struct('key',sample.key,'label',sample.label,'unit',sample.unit, ...
        'scope','逐站内部准定常外流；实际保存时刻，不是非定常 CFD 或时间插值结果', ...
        'history',history);
    if isempty(entries),entries=entry;else,entries(end+1)=entry;end %#ok<AGROW>
end
end
function checkCondition(snapshot,condition)
% 同一时刻标签之外，数值快照还必须属于同一组飞行条件。
expected=[condition.Mach condition.altitude condition.alpha condition.beta];
assert(isequal(snapshot.condition,expected), ...
    'RocketAeroHeat:Viewer:FlowCondition','外流快照工况与该时间站不一致。');
end
