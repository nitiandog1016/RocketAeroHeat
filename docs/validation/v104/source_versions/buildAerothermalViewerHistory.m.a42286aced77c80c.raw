function entries = buildAerothermalViewerHistory(mesh)
%BUILDAEROTHERMALVIEWERHISTORY 导出完整热历史，保持时间、状态和方法的对应。
% 模块一：统一时间轴；扁平数组按“每个时间站的全部面元”顺序保存，
% 避免 MATLAB JSON 将单面元/单时刻矩阵压缩后产生维度歧义。
entries=struct([]);
if ~isfield(mesh,'wallThermal'),return;end
w=mesh.wallThermal; n=size(mesh.F,1); t=double(w.time(:)'); nt=numel(t);
assert(nt>0 && all(isfinite(t)) && t(1)==0 && all(diff(t)>0), ...
    'RocketAeroHeat:Viewer:HistoryTime','热历史时间必须有限且严格递增。');
assert(numel(w.validFace)==n,'RocketAeroHeat:Viewer:HistorySize','热历史有效面元数量错误。');
valid=repmat(logical(w.validFace(:)),1,nt);
scope='实际接受的离散时间站；无时间插值；未标注时间的单点字段保留其原工况。';
if isfield(w,'scope'),scope=[scope '；' w.scope];end
timeNote=repmat({'固定工况热响应；工况见单点结果说明'},1,nt);
if isfield(w,'flightHistory')
    flight=w.flightHistory;
    assert(isequal(double(flight.time(:)'),t),'RocketAeroHeat:Viewer:HistoryTime','飞行条件与热响应时间轴不一致。');
    for name={'Mach','altitude','alpha','beta'}
        assert(numel(flight.(name{1}))==nt && all(isfinite(flight.(name{1}))), ...
            'RocketAeroHeat:Viewer:HistorySize','飞行条件记录尺寸或数值无效。');
    end
    for j=1:nt
        timeNote{j}=sprintf('热响应工况：M=%.4g，H=%.6g m，攻角 %.4g°，侧滑 %.4g°', ...
            flight.Mach(j),flight.altitude(j),flight.alpha(j),flight.beta(j));
    end
end
thermal=repmat({w.method},n,nt);
for key={'lumpedFace','oneDimensionalFace'}
    if isfield(w,key{1})
        mask=logical(w.(key{1})(:));
        if strcmp(key{1},'lumpedFace'),label='本面元：集中热容';else,label='本面元：一维厚度导热';end
        thermal(mask,:)=cellstr(string(thermal(mask,:))+"；"+label);
    end
end
% 模块二：对流历史使用逐站记录的热边界编号；不把末站方法当成所有站。
convection=thermal;
if isfield(w,'coupling') && isfield(w.coupling,'faceMethodIDHistory')
    c=w.coupling; ids=double(c.faceMethodIDHistory);
    checkSize(ids);
    labels={};
    if isfield(c,'finalHeating') && isfield(c.finalHeating,'faceMethodNames')
        labels=c.finalHeating.faceMethodNames;
    end
    convection=appendIDs(convection,ids,labels,'热边界');
    if isfield(c,'turbulentSourceIDHistory')
        labels={};
        if isfield(c.finalHeating,'turbulentBoundaryLayer') && ...
                isfield(c.finalHeating.turbulentBoundaryLayer,'sourceNames')
            labels=c.finalHeating.turbulentBoundaryLayer.sourceNames;
        end
        sources=double(c.turbulentSourceIDHistory);
        sources(~ismember(ids,[2 5]))=-1; % 层流面不显示无关的湍流“未计算”状态。
        convection=appendIDs(convection,sources,labels,'湍流来源');
    end
elseif isfield(w,'heatingHistoryMethod')
    convection=cellstr(string(convection)+"；"+w.heatingHistoryMethod);
    if isfield(mesh,'flightHeatingHistory') && isfield(mesh.flightHeatingHistory,'faceMethodID')
        fh=mesh.flightHeatingHistory;
        assert(isequal(double(fh.time(:)'),t),'RocketAeroHeat:Viewer:HistoryTime','逐站热边界时间轴不一致。');
        labels={};
        if isfield(mesh,'heating') && isfield(mesh.heating,'faceMethodNames'),labels=mesh.heating.faceMethodNames;end
        convection=appendIDs(convection,double(fh.faceMethodID),labels,'热边界');
    end
else
    convection=cellstr(string(convection)+"；局部换热边界反馈（详见热响应适用范围）");
end
reason=repmat({'热响应无有效历程；未以给定初温补齐'},n,1);
reason(w.validFace)={'按当前时间站的有限数值和求解器有效标志显示'};
if isfield(w,'coupling') && isfield(w.coupling,'lostFace')
    for i=find(w.coupling.lostFace(:))'
        reason{i}=sprintf('耦合途中失效，整条历程留空；首次失效 t=%.6g s',w.coupling.firstFailureTime(i));
    end
end
% 模块三：逐变量记录；初始壁温是给定条件，其后的壁温才是求解结果。
put('wallSolved','热响应壁温','K',w.temperature,1,thermal,true);
put('wallFinalFlux','热响应对流热流','kW/m²',w.convectiveHeatFlux,1000,convection,false);
optional={ 'radiativeHeatFlux','wallRadiation','外表面辐射散热','kW/m²',1000; ...
    'netHeatFlux','wallNetFlux','材料净输入热流','kW/m²',1000; ...
    'backTemperature','wallBackTemperature','内表面温度','K',1; ...
    'backHeatFlux','wallBackFlux','内侧净散热','kW/m²',1000; ...
    'recoveryTemperatureHistory','wallRecovery','热响应恢复温度','K',1};
for j=1:size(optional,1)
    if isfield(w,optional{j,1})
        methods=thermal;
        if strcmp(optional{j,1},'recoveryTemperatureHistory'),methods=convection;end
        if strcmp(optional{j,1},'radiativeHeatFlux')
            methods=cellstr(string(thermal)+"；Stefan–Boltzmann：εσ(Tw⁴−Trad⁴)，向外为正");
        elseif strcmp(optional{j,1},'netHeatFlux')
            methods=cellstr(string(convection)+"；材料净输入=外对流−外辐射−内侧散热");
        elseif strcmp(optional{j,1},'backHeatFlux') && isfield(w,'backBoundary')
            methods=cellstr(string(thermal)+"；"+w.backBoundary.method+"；内侧散出为正");
        end
        put(optional{j,2},optional{j,3},optional{j,4},w.(optional{j,1}),optional{j,5},methods,strcmp(optional{j,1},'backTemperature'));
    end
end
% 累计热量只导出已经过能量核对的结果，不在查看器里重新积分。
if isfield(mesh,'aerothermalExposure')
    e=mesh.aerothermalExposure;
    assert(isequal(double(e.time(:)'),t),'RocketAeroHeat:Viewer:ExposureTime','累计热量与壁温时间轴不一致。');
    assert(numel(e.validFace)==n,'RocketAeroHeat:Viewer:HistorySize','累计热量的有效面元数量错误。');
    savedValid=valid; valid=valid & repmat(logical(e.validFace(:)),1,nt);
    if isfield(e,'scope'),scope=[scope '；' e.scope];end
    if isfield(e,'energyConvention'),scope=[scope '；' e.energyConvention];end
    for item={'cumulativeConvectiveHeat','cumulativeNetHeat','cumulativeRadiativeHeat','cumulativeBackHeat'}
        key=item{1};if ~isfield(e,key),continue;end
        labels=struct('cumulativeConvectiveHeat','累计对流热量','cumulativeNetHeat','材料累计净热量', ...
            'cumulativeRadiativeHeat','累计外辐射散热','cumulativeBackHeat','累计内侧散热');
        put(key,labels.(key),'MJ/m²',e.(key),1e6,repmat({[e.method '；' w.method]},n,nt),false);
    end
    valid=savedValid; %#ok<NASGU>
end
    function put(key,label,unit,values,divisor,methods,initialGiven)
        checkSize(values); values=double(values)/divisor;
        status=double(valid & isfinite(values)); values(status==0)=NaN;
        if initialGiven
            status(status(:,1)>0,1)=2;methods(status(:,1)==2,1)={'给定材料初始温度'};
        end
        if strcmp(key,'wallBackFlux') && isfield(w,'backBoundary') && ~w.backBoundary.active
            status(status>0)=2;methods(status>0)={'给定绝热内侧边界：热流为零'};
        end
        [names,~,indices]=unique(methods(:),'stable');
        history=struct('time',t,'layout','time-major','values',values(:), ...
            'status',status(:),'methodNames',{names},'methodIndex',indices-1,'faceNote',{reason},'timeNote',{timeNote});
        entry=struct('key',key,'label',label,'unit',unit,'scope',scope,'history',history);
        if isempty(entries),entries=entry;else,entries(end+1)=entry;end
    end
    function checkSize(values)
        assert(isequal(size(values),[n nt]),'RocketAeroHeat:Viewer:HistorySize','历史字段与网格/时间尺寸不匹配。');
    end
    function out=appendIDs(out,ids,labels,prefix)
        checkSize(ids);
        for id=unique(ids(:))'
            if id<0,continue;end
            label=sprintf('%s记录编号 %g（名称未记录）',prefix,id);
            if id>=0 && id==fix(id) && id<numel(labels),label=[prefix '：' labels{id+1}];end
            mask=ids==id;out(mask)=cellstr(string(out(mask))+"；"+label);
        end
    end
end
