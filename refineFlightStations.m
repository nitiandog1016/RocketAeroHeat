function [refined,info] = refineFlightStations(flightHistory,maxTimeStep)
%REFINEFLIGHTSTATIONS 在用户时间站之间生成较密的内置计算站。
% 模块作用：保留所有用户给出的时刻，并对高度、Mach、攻角和侧滑角
% 作分段线性插值，使相邻计算站的时间间隔不超过 maxTimeStep。
% 这只规定飞行条件的输入轨迹；每个新站仍需独立求定常外流和热流。
% 不推断用户站点之外的轨迹，也不计算非定常激波或流动记忆。

%% 1. 检查时间、字段及最大计算间隔
fields={'time','altitude','Mach','alpha','beta'};
for j=1:numel(fields)
    name=fields{j};
    if ~isfield(flightHistory,name)
        error('RocketAeroHeat:FlightRefine:Field','缺少 %s。',name);
    end
    validateattributes(flightHistory.(name),{'numeric'}, ...
        {'real','finite','vector','nonempty'},mfilename,name);
end
validateattributes(maxTimeStep,{'numeric'}, ...
    {'real','finite','scalar','positive'});
inputTime=double(flightHistory.time(:)');
if numel(inputTime)<2 || inputTime(1)~=0 || any(diff(inputTime)<=0)
    error('RocketAeroHeat:FlightRefine:Time', ...
        '用户时间站须从 0 开始，至少两站且严格递增。');
end
for j=2:numel(fields)
    if numel(flightHistory.(fields{j}))~=numel(inputTime)
        error('RocketAeroHeat:FlightRefine:Length', ...
            '%s 与时间站数量不一致。',fields{j});
    end
end

%% 2. 每段均匀细分；明确保留原始站点的精确时间
substeps=max(1,ceil(diff(inputTime)/maxTimeStep));
stationCount=1+sum(substeps);
if ~isfinite(stationCount) || stationCount>1e7
    error('RocketAeroHeat:FlightRefine:Size', ...
        '时间加密后站点过多，请增大最大计算间隔。');
end
time=zeros(1,stationCount);
time(1)=inputTime(1);
inputIndex=zeros(1,numel(inputTime));
inputIndex(1)=1;
cursor=1;
for j=1:numel(substeps)
    count=substeps(j);
    fraction=(1:count)/count;
    time(cursor+(1:count))=inputTime(j)+ ...
        fraction*(inputTime(j+1)-inputTime(j));
    cursor=cursor+count;
    time(cursor)=inputTime(j+1);
    inputIndex(j+1)=cursor;
end
if any(diff(time)<=0)
    error('RocketAeroHeat:FlightRefine:Precision', ...
        '计算时间站在浮点精度下不再严格递增，请增大最大间隔。');
end

%% 3. 对每个飞行输入分段线性插值，并记录来源与步长
refined=struct('time',time);
for j=2:numel(fields)
    name=fields{j};
    values=double(flightHistory.(name));
    values=values(:)';
    interpolated=interp1(inputTime,values,time,'linear');
    interpolated(inputIndex)=values;
    refined.(name)=interpolated;
end
info=struct;
info.method='用户飞行时间站分段线性插值 + 最大时间间隔均匀细分';
info.scope=['仅在相邻用户站之间插值飞行条件；' ...
    '不模拟非定常外流、激波运动或未给出的轨迹'];
info.maxTimeStep=maxTimeStep;
info.inputStationCount=numel(inputTime);
info.evaluationStationCount=stationCount;
info.inputStationIndex=inputIndex;
info.segmentSubsteps=substeps;
info.largestActualStep=max(diff(time));
end
