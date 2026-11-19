function [current,freestream,report]=runHybridSurfaceCase(mesh,cfg,progressFcn)
%RUNHYBRIDSURFACECASE 同一工况先算工程场，再对缺项区域调用黏性求解器。
% 模块一：三角形顺序定义零起点原面编号。调用方提供与几何同序的物理区域、
% 壁温和工况，以及两个求解函数。工程函数 (mesh,condition,Tw) 返回计算表；
% 局部函数 (mesh,condition,Tw,regionIDs,targetIDs,engineering) 返回
% accepted、reason 和 wall。regionIDs/targetIDs 均为原面编号，前者含全区。
% 各方法自行核对适用性、边界相容与数值接受条件。wall 只收计算值，
% 列为 faceID、pressurePa、heatFluxWm2、method。缺项用 NaN 或省略该面。
% 工程表的 Properties.UserData 可携带同一工况计算的区域交界状态，
% 完整表原样传入局部方法，由该方法核对界面物理与区域身份。
if nargin<3,progressFcn=[];end
setup=cfg.hybridSurface;
validateattributes(mesh.V,{'numeric'},{'real','finite','ncols',3,'nonempty'});
validateattributes(mesh.F,{'numeric'},{'integer','positive','ncols',3,'nonempty','<=',size(mesh.V,1)});
f=mesh.F; v=mesh.V; n=size(f,1);
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
validateattributes(area,{'numeric'},{'finite','positive'});
region=string(setup.region(:)); tw=setup.wallTemperatureK(:);
assert(numel(region)==n && all(~ismissing(region) & strlength(strtrim(region))>0),'RAH:Hybrid:Region','每个原面须有区域。');
validateattributes(tw,{'numeric'},{'real','positive','finite','numel',n});
assert(isstruct(setup.condition)&&isscalar(setup.condition)&&~isempty(fieldnames(setup.condition)), ...
    'RAH:Hybrid:Condition','需要明确的物理工况。');
% 提供给两个求解器的是同一次运行的输入几何，清除旧快照内可能遗留的结果。
geometry=struct('V',v,'F',f,'region',region);
ids=(0:n-1)'; wall=table(ids,nan(n,1),nan(n,1),strings(n,1), ...
    'VariableNames',{'faceID','pressurePa','heatFluxWm2','method'});
started=tic;
engineering=setup.engineeringSolver(geometry,setup.condition,tw);
idx=checkWall(engineering,ids);
wall(idx,:)=engineering(:,wall.Properties.VariableNames);
pMethod=wall.method;qMethod=wall.method;
pMethod(idx)=fieldMethod(engineering,'pressureMethod');qMethod(idx)=fieldMethod(engineering,'heatFluxMethod');
nearT=nan(n,1);nearMethod=strings(n,1);
[temperatures,methods]=adjacentFluid(engineering);nearT(idx)=temperatures;nearMethod(idx)=methods;
eng=isfinite(wall.pressurePa)&isfinite(wall.heatFluxWm2);
pending=unique(region(~eng),'stable'); attempts=struct('region',{},'targetFaceIDs',{},'accepted',{},'reason',{},'diagnostics',{});
% 模块二：每个有缺项的区域只调用一次。已有完整工程解优先保留；缺项原面
% 只有在同一次被接受的区域结果给出成对 p/q 时才替换，保留合法零和负热流。
for k=1:numel(pending)
    inRegion=region==pending(k); target=ids(inRegion&~eng);
    local=setup.localSolver(geometry,setup.condition,tw,ids(inRegion),target,engineering);
    validateattributes(local.accepted,{'logical'},{'scalar'});
    reason=string(local.reason);
    assert(isscalar(reason)&&~ismissing(reason)&&strlength(strtrim(reason))>0,'RAH:Hybrid:Reason','区域须记录停止或接受原因。');
    % 保留适配器实际计算的接口与守恒量，拒绝的区域也能回查。
    diagnostics=struct;
    if isfield(local,'diagnostics'),diagnostics=local.diagnostics;end
    attempts(end+1)=struct('region',pending(k),'targetFaceIDs',target,'accepted',local.accepted, ...
        'reason',reason,'diagnostics',diagnostics); %#ok<AGROW>
    if ~local.accepted,continue;end
    idx=checkWall(local.wall,ids(inRegion));
    take=ismember(local.wall.faceID,target)&isfinite(local.wall.pressurePa)&isfinite(local.wall.heatFluxWm2);
    wall(idx(take),:)=local.wall(take,wall.Properties.VariableNames);
    sources=fieldMethod(local.wall,'pressureMethod');pMethod(idx(take))=sources(take);
    sources=fieldMethod(local.wall,'heatFluxMethod');qMethod(idx(take))=sources(take);
    [temperatures,methods]=adjacentFluid(local.wall);nearT(idx(take))=temperatures(take);nearMethod(idx(take))=methods(take);
end
if ~isequal(pMethod,wall.method) || ~isequal(qMethod,wall.method)
    wall.pressureMethod=pMethod;wall.heatFluxMethod=qMethod;
end
if any(isfinite(nearT))
    wall.nearWallCellTemperatureK=nearT;wall.nearWallTemperatureMethod=nearMethod;
end
% 模块三：完整性以每个原面均有双变量计算值为准，独立验证由证据台账处理。
p=isfinite(wall.pressurePa); q=isfinite(wall.heatFluxWm2); joint=p&q;
total=sum(area); power=NaN;
if all(q),power=sum(area.*wall.heatFluxWm2);end
report=struct('status','completed','faceCount',n,'totalArea',total,'engineeringFirstAreaM2',sum(area(eng)), ...
    'regionalAddedAreaM2',sum(area(joint&~eng)),'jointComputedAreaM2',sum(area(joint)), ...
    'fullSurfaceComputed',all(joint),'wholeSurfaceHeatPowerW',power,'attempts',attempts, ...
    'pressureComputedAreaFraction',sum(area(p))/total,'heatComputedAreaFraction',sum(area(q))/total, ...
    'elapsedSeconds',toc(started));
current=struct('V',v,'F',f,'localGeometry',struct('region',{cellstr(region)}), ...
    'hybridWall',wall,'hybridCondition',setup.condition,'hybridWallTemperatureK',tw, ...
    'evidenceIdentity',struct('origin','research','physicalCondition',setup.condition));
freestream=[];
if ~isempty(progressFcn)
    progressFcn(struct('phase','coverage','status','completed','fraction',1,'message','混合壁面计算完成', ...
        'elapsedSeconds',report.elapsedSeconds,'errorIdentifier','','errorMessage',''));
end
end

function index=checkWall(rows,allowed)
% 模块四：原面编号显式绑定求解区，支持乱序和部分结果。
keys={'faceID','pressurePa','heatFluxWm2','method'};
assert(istable(rows)&&all(ismember(keys,rows.Properties.VariableNames)),'RAH:Hybrid:Table','结果表缺少原面或物理量列。');
id=rows.faceID;
assert(isnumeric(id)&&isreal(id)&&size(id,2)==1&&all(ismember(id,allowed))&&numel(unique(id))==numel(id), ...
    'RAH:Hybrid:FaceID','存在重复或区域外原面。');
for key={'pressurePa','heatFluxWm2'}
    x=rows.(key{1});
    assert(isnumeric(x)&&isreal(x)&&size(x,2)==1&&~any(isinf(x)),'RAH:Hybrid:Field','场量只接受有限实数或 NaN。');
end
assert(all(isnan(rows.pressurePa)|rows.pressurePa>0),'RAH:Hybrid:Pressure','计算压力须为正。');
assert(isstring(rows.method)&&size(rows.method,2)==1&&all(~ismissing(rows.method)&strlength(strtrim(rows.method))>0), ...
    'RAH:Hybrid:Method','每行须注明求解方法。');
index=double(id)+1;
end

function names=fieldMethod(rows,key)
% 允许压力和热流使用不同方法；历史适配器继续使用共同的 method 列。
names=rows.method;
if ismember(key,rows.Properties.VariableNames)
    names=rows.(key);
    assert(isstring(names)&&isequal(size(names),size(rows.method))&&all(~ismissing(names)&strlength(strtrim(names))>0), ...
        'RAH:Hybrid:Method','逐变量方法须与原面对应。');
end
end

function [values,methods]=adjacentFluid(rows)
% 可选临壁流体结果按同一次接受的原面记录传递，未提供时维持未知。
values=nan(height(rows),1);methods=strings(height(rows),1);
if ~ismember('nearWallCellTemperatureK',rows.Properties.VariableNames),return;end
values=rows.nearWallCellTemperatureK;
assert(isnumeric(values) && isreal(values) && isequal(size(values),size(rows.faceID)) && ...
    all(isnan(values) | (isfinite(values)&values>0)),'RAH:Hybrid:Temperature','临壁流体温度须为正 K 或 NaN。');
methods=fieldMethod(rows,'nearWallTemperatureMethod');
end
