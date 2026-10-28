function [current,freestream,report] = runAeroHeatSingleCase(mesh,cfg,progressFcn)
%RUNAEROHEATSINGLECASE 显式执行一个当前工况的工程气动热计算。
% 模块一：核对已完成前处理的网格和当前边界，避免浏览页面自动求解。
% 输入 mesh 须已完成外法向、表面坐标和局部几何；cfg.flight 使用 m、K、deg。
% 模块二：依次调用既有压力、载荷、驻点、内部外流、波系、流线和热流。
% 与飞行历史入口不同，本模块保留当前给定的 wallTemperatureByFace。
% 模块三：向可选 progressFcn(event) 报告真实阶段，不估计阶段内部进度。
% event 含 phase/status/fraction/message/elapsedSeconds/errorIdentifier/errorMessage。
% fraction 为已完成阶段占比，不代表物理解面积或预计耗时。阶段失败原样抛错；
% 只有全链完成才返回 current，不能以中间网格冒充完成快照。
% 模块四：report 保存接受的设置与面积审计；执行完成不等于全面积求解完成。
% 未解值、给定值和场量方法均由既有后端保留，不作补零或模型回退。
if nargin<3,progressFcn=[];end
if ~isempty(progressFcn) && ~isa(progressFcn,'function_handle')
    error('RocketAeroHeat:SingleCase:ProgressCallback','progressFcn 须为空或函数句柄。');
end
started=tic;
phases={'inputs','freestream','regions','pressure','loads','stagnation', ...
    'externalFlow','axisymmetricMarch','streamlines','heating','coverage'};
messages={'检查当前网格和边界条件','计算自由来流','更新当前来流区域', ...
    '计算表面压力与所选底压','积分当前压力载荷','定位驻点区域', ...
    '计算内部工程外流','推进轴向激波/膨胀','追踪当前表面准流线', ...
    '计算当前壁温下的工程热边界','核对当前场量面积和完成快照'};
stageCount=numel(phases);events=struct([]);stage=1;
current=mesh;freestream=[];report=struct;
try
    for stage=1:stageCount
        announce('started',(stage-1)/stageCount,'','');
        switch phases{stage}
            case 'inputs'
                if ~isstruct(cfg) || ~isscalar(cfg) || ~isfield(cfg,'flight')
                    error('RocketAeroHeat:SingleCase:Config','需要当前 cfg.flight。');
                end
                if ~isstruct(mesh) || ~isscalar(mesh) || ...
                        ~all(isfield(mesh,{'V','F','surfaceCoord','localGeometry'}))
                    error('RocketAeroHeat:SingleCase:PreparedMesh', ...
                        '请先完成网格、外法向、表面坐标和局部几何。');
                end
                % 当前非均匀壁温须与实际面元逐一对应；NaN 明确保留未知。
                surfaceWallTemperature(size(mesh.F,1),cfg);
            case 'freestream'
                freestream=calcFreestream(cfg.flight);
                if freestream.Mach<=1
                    error('RocketAeroHeat:SingleCase:Mach', ...
                        '当前内部激波/膨胀工程链要求 M>1；不能以亚声速边界假装完成。');
                end
            case 'regions'
                % 原分类器独立重算 V/F 几何，拒绝已改变网格的旧局部数据。
                current=classifyAerothermalRegion(current,freestream,cfg);
            case 'pressure'
                current=modifiedNewtonianPressure(current,freestream,cfg);
            case 'loads'
                current=integrateAerodynamicLoads(current,cfg);
            case 'stagnation'
                current=findStagnationRegion(current,freestream,cfg);
            case 'externalFlow'
                current=solveExternalFlow(current,freestream,cfg);
            case 'axisymmetricMarch'
                current=marchAxisymmetricShockExpansion(current,freestream,cfg);
            case 'streamlines'
                current=traceSurfaceStreamlines(current,cfg);
            case 'heating'
                current=calcSurfaceHeating(current,freestream,cfg);
            case 'coverage'
                % 面积由实际三角形取得；预设底压或背风Cp不计入求解压力面积。
                area=current.localGeometry.area(:);totalArea=sum(area);
                ps=current.pressure.fieldStatus(:);
                if numel(ps)~=numel(area) || ...
                        any(~ismember(ps,[0 1 2])) || ...
                        ~isequal(ps==0,~isfinite(current.pressure.staticPressure(:)))
                    error('RocketAeroHeat:SingleCase:PressureStatus', ...
                        '当前压力值与逐面元未解/计算/给定状态不一致。');
                end
                report.status='completed';
                report.cfg=cfg;
                report.condition=cfg.flight;
                report.faceCount=size(current.F,1);
                report.totalArea=totalArea;
                report.pressureComputedAreaFraction=sum(area(ps==1))/totalArea;
                report.pressurePrescribedAreaFraction=sum(area(ps==2))/totalArea;
                report.pressureMissingAreaFraction=sum(area(ps==0))/totalArea;
                report.heatComputedAreaFraction=current.heatingCoverage.summary.computedAreaFraction;
                report.loadsAvailable=current.aerodynamicLoads.available;
                report.scope=['当前工况工程链执行完成；面积审计不是区域求解，' ...
                    '仍需检查各变量的未解、给定与模型适用范围。'];
        end
        announce('completed',stage/stageCount,'','');
    end
catch ME
    % 失败通知不能遮蔽真实后端错误；通知回调自身出错也不返回部分结果。
    try
        announce('failed',(stage-1)/stageCount,ME.identifier,ME.message);
    catch
    end
    rethrow(ME);
end
report.elapsedSeconds=toc(started);
report.events=events;

    function announce(status,fraction,errorIdentifier,errorMessage)
        event=struct('phase',phases{stage},'status',status,'fraction',fraction, ...
            'message',messages{stage},'elapsedSeconds',toc(started), ...
            'errorIdentifier',errorIdentifier,'errorMessage',errorMessage);
        if strcmp(status,'failed'),event.message=[event.message '：失败'];end
        if isempty(events),events=event;else,events(end+1,1)=event;end %#ok<AGROW>
        if ~isempty(progressFcn),progressFcn(event);end
    end
end
