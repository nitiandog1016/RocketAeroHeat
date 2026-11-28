function check_flight_cases(out)
% 模块一：绑定简化试件原面。生成器提供环带拓扑，计算顶点取已发布的米制 STL。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
tunnel=createCaseConfig('simpleLauncher');def=jsondecode(fileread(fullfile(root,'cases','simple_launcher','geometry.json')));
assert(strcmp(cm011FileSHA256(tunnel.geometry.stlFile),def.fileSHA256));
raw=loadSTLMesh(tunnel,true);mesh=buildRocketTestBody(def.options);
mesh.V=double(single(mesh.V));
for k=1:3
    assert(isequal(mesh.V(mesh.F(:,k),:),raw.V(raw.F(:,k),:)),'原面几何或顺序不符。');
end
cfg=createCaseConfig();cfg.computationBackend='engineering';cfg.geometry=tunnel.geometry;
mesh=calcFaceGeometry(mesh);mesh=buildSurfaceCoordinates(mesh,cfg);mesh=analyzeLocalGeometry(mesh,cfg);
area=mesh.localGeometry.area(:);n=numel(area);assert(n==def.faces && abs(sum(area)-def.areaM2)<1e-12);
fixed=ready(tunnel,raw); generic=ready(cfg,mesh);
mach=[.30 .85 1 1.20 2.00];
% 模块二：使用相同的环境高度、姿态及给定壁温。分别试固定试件入口和通用工程入口。
runs=cell(1,numel(mach));
for i=1:numel(mach)
    af=cfg;af.flight.altitude=5000;af.flight.alpha=0;af.flight.beta=0;af.flight.Twall=300;af.flight.Mach=mach(i);
    fs=calcFreestream(af.flight);timer=tic;phase='inputs';
    entry=struct('mach',mach(i),'speedMps',fs.speed,'soundSpeedMps',fs.speedOfSound,'freeTemperatureK',fs.temperature, ...
        'freePressurePa',fs.pressure,'altitudeM',af.flight.altitude,'wallGivenK',af.flight.Twall, ...
        'fixedBoundaryError','','engineeringError','','lastStage','','solverReturned',false,'publicationError','', ...
        'resultPageAccepted',false,'wholeSurfacePowerKnown',false);
    requested=tunnel.launcherBoundary;
    requested.velocityMps=[-fs.speed 0 0];requested.wallTemperatureK=af.flight.Twall;
    requested.ambientTemperatureK=fs.temperature;requested.ambientPressurePa=fs.pressure;
    try
        aeroHeatWorkflowState('confirmBoundary',fixed,struct('freestream',requested));
    catch ex
        entry.fixedBoundaryError=ex.identifier;
        if ~strcmp(ex.identifier,'RAH:Launcher:Boundary'),rethrow(ex);end
    end
    state=aeroHeatWorkflowState('setBoundary',generic,struct('cfg',af));
    state=aeroHeatWorkflowState('confirmBoundary',state,struct('freestream',fs));
    state=aeroHeatWorkflowState('beginSolve',state);
    p=zeros(n,1,'uint8');q=false(n,1);lastError='';
    try
        [result,actual,report]=runAeroHeatSingleCase(mesh,af,@track);
        entry.solverReturned=true;
    catch ex
        lastError=ex.identifier;entry.engineeringError=lastError;entry.engineeringMessage=ex.message;
        if isempty(ex.identifier),rethrow(ex);end
    end
    if entry.solverReturned
        % 模块三：只计“计算”状态的有限壁压及合格热流，另列给定压力与交集。
        assert(actual.Mach==fs.Mach && isequal(result.F,mesh.F) && isequal(result.V,mesh.V));
        p=result.pressure.fieldStatus(:);q=logical(result.heating.eligibleFace(:));
        assert(all(ismember(p,[0 1 2])) && isequal(isfinite(result.pressure.staticPressure(:)),p~=0) && ...
            isequal(q,isfinite(result.heating.heatFlux(:))));
        assert(abs(sum(area(p==1 & q))-report.jointComputedAreaM2)<1e-12);
        entry.wholeSurfacePowerKnown=isfinite(report.wholeSurfaceHeatPowerW);
        entry.pressureMethod=result.pressure.model;entry.heatMethod=result.heating.method;
        entry.pressureScope=unique(result.pressure.faceScope);entry.heatScope=result.heating.validity;
        token=struct('result',result,'freestream',actual,'sourceRevision',state.solveInputRevision,'solveId',state.solveId);
        try
            finished=aeroHeatWorkflowState('completeSolve',state,token);
            entry.resultPageAccepted=finished.steps(5).available && ~isempty(finished.result);
        catch ex
            entry.publicationError=ex.identifier;
            if ~strcmp(ex.identifier,'RocketAeroHeat:WorkflowIncompleteSurface'),rethrow(ex);end
        end
    else
        assert(~isempty(lastError));entry.publicationError='noReturnedSnapshot';
    end
    entry.lastStage=phase;entry.elapsedSeconds=toc(timer);
    entry.coverage=coverage(area,p,q);
    for group=def.regions'
        ids=group.firstFace:group.lastFace;
        entry.regions.(group.name)=coverage(area(ids),p(ids),q(ids));
    end
    % 面号由行号给出；本地只保存状态掩码，供独立面积复算。
    save(fullfile(out,sprintf('masks_%02d.mat',i)),'p','q','-v7');
    runs{i}=entry;
    fprintf('M %.2f | %s | joint %d/%d | %.5f m2\n',mach(i),phase,entry.coverage.jointFaces,n,entry.coverage.jointAreaM2);
    writeReport();
end
    function track(event)
        phase=event.phase;
    end
    function writeReport()
        rec=struct('surfaceSHA256',def.fileSHA256,'faces',n,'areaM2',sum(area),'cases',{runs(1:i)});
        writelines(string(jsonencode(rec,'PrettyPrint',true)),fullfile(out,'cases.json'),'Encoding','UTF-8');
    end
end

function state=ready(cfg,m)
% 走五步状态机，把同一个已核对的表面准备到边界确认步骤。
state=aeroHeatWorkflowState('init',[],struct('cfg',cfg));
state=aeroHeatWorkflowState('setModel',state,struct('rawMesh',m));
state=aeroHeatWorkflowState('setOrientation',state,struct('orientedMesh',m));
state=aeroHeatWorkflowState('setMesh',state,struct('mesh',m));
end

function c=coverage(a,p,q)
% 拒绝整例时各项计数代表本次返回的可用面数；给定量始终单计。
c=struct('faces',numel(a),'areaM2',sum(a),'computedPressureFaces',nnz(p==1), ...
    'computedPressureAreaM2',sum(a(p==1)),'prescribedPressureFaces',nnz(p==2),'prescribedPressureAreaM2',sum(a(p==2)), ...
    'heatFaces',nnz(q),'heatAreaM2',sum(a(q)),'jointFaces',nnz(p==1&q),'jointAreaM2',sum(a(p==1&q)));
end
