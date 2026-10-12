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
    add('pressure','表面压力（牛顿/给定）','kPa',p.staticPressure/1000,state,methods, ...
        '修正牛顿压力与给定背风/底压；不是尾迹求解结果');
    add('cp','压力系数 Cp','1',p.Cp,state,methods,'与表面压力采用相同来源');
end
if isfield(mesh,'boundaryLayerEdge')
    e=mesh.boundaryLayerEdge; methods=methodNames(e,n);
    state=double(e.algebraicallyValid(:));
    add('edgePressure','边界层外缘压力','kPa',e.pressure/1000,state,methods, ...
        '内部工程外流，代数有效不代表经过实验验证');
    add('edgeTemperature','边界层外缘温度','K',e.temperature,state,methods, ...
        '内部工程外流；不是壁面温度');
end
if isfield(mesh,'heating')
    h=mesh.heating; state=double(h.eligibleFace(:)); methods=heatingMethods(h,n);
    detail={};
    for key={'boundaryLayerMethod','turbulentMethod','noseMethod','stagnationMethod','thermalPropertyMethod'}
        if isfield(h,key{1}),detail{end+1}=[key{1} ': ' char(h.(key{1}))];end %#ok<AGROW>
    end
    scope=strjoin([{h.validity},detail],'；');
    add('heatFlux','对流热流密度','kW/m²',h.heatFlux/1000,state,methods,scope);
    add('recoveryTemperature','恢复温度','K',h.recoveryTemperature,state,methods, ...
        '恢复焓/恢复温度模型；不是求解壁温');
    add('wallInput','热流计算所用壁温','K',h.wallTemperatureByFace,2*ones(n,1), ...
        repmat({'给定/耦合传入的热流边界温度'},n,1), ...
        '此字段是当前热流计算的输入，不等于末时刻热响应');
end
% 模块三：采用求解器实际接受的时间轴，不插值、不恢复已失效的历程。
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
% 模块三：按真实三角形面积分区汇总，给定值单独统计。
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
    function add(key,label,unit,values,state,methods,scope)
        values=double(values(:)); state=double(state(:));
        assert(numel(values)==n && numel(state)==n && numel(methods)==n, ...
            'RocketAeroHeat:Viewer:Size','字段 %s 与面元数量不一致。',key);
        state(~isfinite(values))=0; values(state==0)=NaN;
        [names,~,ids]=unique(methods,'stable');
        entry=struct('key',key,'label',label,'unit',unit,'values',values, ...
            'status',state,'methodNames',{names},'methodIndex',ids-1,'scope',scope, ...
            'computedAreaFraction',sum(area(state==1))/sum(area), ...
            'prescribedAreaFraction',sum(area(state==2))/sum(area), ...
            'missingAreaFraction',sum(area(state==0))/sum(area),'history',[]);
        if isempty(data.fields),data.fields=entry;else,data.fields(end+1)=entry;end
    end
end
function names=heatingMethods(h,n)
% 模块五：展开真实的局部边界层来源，区分 Head 积分和显式平板回退。
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
    end
end
end
function names=methodNames(result,n)
% 模块四：将原始逐面元算法编号解码；保留原名称。
names=repmat({'来源未记录'},n,1);
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
