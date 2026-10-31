function fig=rocketAeroHeatWorkbench(cfg)
%ROCKETAEROHEATWORKBENCH Fluent式五步骤科研工作台，计算仍由MATLAB执行。
% 模块一：一个侧边栏依次提供模型、网格、边界、显式计算及结果入口。
% 模块二：集中状态控制输入版本和下游失效；查看页面不触发物理重算。
% 模块三：几何与结果采用离线WebGL，未求解预览没有虚构物理字段。
% 模块四：保留既有SI坐标、朝向和网格算法；各区域求解状态来自后端。
if nargin<1,cfg=createCaseConfig();end
% 统一可见默认值，避免只查看再确认就因补默认字段而误判输入改变。
cfg.heating.sphericalNoseAngleStepDeg=getOption(cfg.heating,'sphericalNoseAngleStepDeg',2);
cfg.pressure.aftBaseCp=getOption(getOption(cfg,'pressure',struct),'aftBaseCp',0);
cfg.pressure.leewardCp=getOption(getOption(cfg,'pressure',struct),'leewardCp',0);
state=aeroHeatWorkflowState('init',[],struct('cfg',cfg));
preview=[];rotation=eye(3);rotationCenter=[0 0 0];controls=struct;
webUI=[];report=[];logLines={'等待输入模型。'};preparing=false;methodDetailsOpen=false;
fig=uifigure('Name','RocketAeroHeat · 科研工作台', ...
    'Position',[30 30 1580 940],'Color',[.96 .97 .98]);
frame=uigridlayout(fig,[3 2]);frame.RowHeight={40,'1x',32};
frame.ColumnWidth={340,'1x'};frame.Padding=[8 6 8 6];frame.RowSpacing=6;
brand=uilabel(frame,'Text','RocketAeroHeat','FontSize',18,'FontWeight','bold');
brand.Layout.Row=1;brand.Layout.Column=1;
title=uilabel(frame,'FontSize',15,'FontWeight','bold');title.Layout.Row=1;title.Layout.Column=2;
sidebar=uipanel(frame,'Title','当前工况工作流','BackgroundColor',[.96 .97 .98]);
sidebar.Layout.Row=2;sidebar.Layout.Column=1;
side=uigridlayout(sidebar,[7 1]);side.RowHeight={32,32,32,32,32,48,'1x'};
side.Padding=[8 8 8 8];side.RowSpacing=5;buttons=gobjects(5,1);
for k=1:5
    step=k;buttons(k)=uibutton(side,'Text',sprintf('%d  %s',k,state.steps(k).label), ...
        'HorizontalAlignment','left','ButtonPushedFcn',@(~,~)go(step));
    buttons(k).Layout.Row=k;
end
inputSummary=uilabel(side,'WordWrap','on');inputSummary.Layout.Row=6;
editor=uipanel(side,'BorderType','none');editor.Layout.Row=7;
center=uipanel(frame,'BorderType','none');center.Layout.Row=2;center.Layout.Column=2;
centerGrid=uigridlayout(center,[1 1]);centerGrid.Padding=[0 0 0 0];
geometryHost=uipanel(centerGrid,'BorderType','none');geometryHost.Layout.Row=1;geometryHost.Layout.Column=1;
resultHost=uipanel(centerGrid,'BorderType','none','Visible','off');resultHost.Layout.Row=1;resultHost.Layout.Column=1;
geometryView=createAerothermalGeometryView(geometryHost,[]);
footer=uilabel(frame,'WordWrap','on');footer.Layout.Row=3;footer.Layout.Column=[1 2];
fig.UserData=struct('workflowState',state,'workflowButtons',buttons, ...
    'workflowSidebar',sidebar,'geometryView',geometryView,'resultHost',resultHost, ...
    'controls',controls,'webUI',[],'caseReport',[],'title',title,'status',footer,'editorPanel',editor);
fig.UserData.actions=struct('selectStep',@go,'loadModel',@()loadInput([],[]), ...
    'pca',@()pcaOrientation([],[]),'flip',@()flipOrientation([],[]), ...
    'confirmOrientation',@()confirmOrientation([],[]),'generateMesh',@()generate([],[]), ...
    'confirmBoundary',@()confirmBoundary([],[]),'run',@()solveCase([],[]),'getState',@getState);
renderEditor();sync();showGeometry(true);

%% 模块五：切页只读取当前状态；所有结果入口依赖已完成的当前版本。
    function go(step)
        if state.busy || preparing || ~isgraphics(fig),return;end
        try
            state=aeroHeatWorkflowState('selectStep',state,struct('step',step));
            renderEditor();sync();
            if step==5,showResults();else,showGeometry(false);end
        catch ex,showError(ex);end
    end
    function s=getState(),s=state;end
    function sync()
        if ~isgraphics(fig),return;end
        % 父面板也必须恢复；仅重建子控件不会解除父面板的禁用状态。
        if isprop(editor,'Enable'),editor.Enable=onoff(~state.busy && ~preparing);end
        for j=1:5
            tag='待前置';
            if state.steps(j).available,tag='待处理';end
            if state.steps(j).completed,tag='已完成';end
            if strcmp(state.steps(j).status,'busy'),tag='计算中';end
            if strcmp(state.steps(j).status,'failed'),tag='失败';end
            buttons(j).Text=sprintf('%d  %s  · %s',j,state.steps(j).label,tag);
            buttons(j).Enable=onoff(state.steps(j).available && ~state.busy && ~preparing);
            buttons(j).BackgroundColor=[.96 .97 .98];
        end
        buttons(state.currentStep).BackgroundColor=[.84 .90 .94];
        title.Text=sprintf('%s · 输入版本 %d',state.steps(state.currentStep).label,state.revision);
        counts='尚无求解网格';
        if ~isempty(state.mesh),counts=sprintf('求解网格 %d 节点 / %d 面元',size(state.mesh.V,1),size(state.mesh.F,1));end
        inputSummary.Text=[counts newline 'SI · 主轴 X · 鼻尖 +X · 零攻角来流 −X'];
        footer.Text=logLines{end};
        fig.UserData.workflowState=state;fig.UserData.controls=controls;
        fig.UserData.webUI=webUI;fig.UserData.caseReport=report;
        % 工作区也跟随当前有效状态；上游改变后不留下旧场量冒充当前结果。
        assignin('base','cfg',state.cfg);assignin('base','rawMeshOriginal',state.rawMesh);
        assignin('base','rawMesh',state.orientedMesh);assignin('base','quality',state.quality);
        current=state.mesh;if ~isempty(state.result),current=state.result;end
        assignin('base','mesh',current);assignin('base','freestream',state.freestream);
        if isfield(state.cfg,'flight'),assignin('base','flight',state.cfg.flight);end
        if state.busy || preparing
            list=findall(editor);
            for q=1:numel(list),if isprop(list(q),'Enable'),list(q).Enable='off';end,end
        end
        drawnow limitrate;
    end
    function message(text)
        logLines=[logLines;{char(text)}];if numel(logLines)>80,logLines=logLines(end-79:end);end
        if isfield(controls,'log') && isgraphics(controls.log),controls.log.Value=logLines;end
        sync();
    end
    function showError(ex)
        message(['操作未完成：' ex.message]);
    end
    function showGeometry(reset)
        if ~isgraphics(fig),return;end
        resultHost.Visible='off';geometryHost.Visible='on';
        m=preview;
        if state.currentStep>=2
            m=state.orientedMesh;if ~isempty(state.mesh),m=state.mesh;end
        end
        caption='仅几何预览；尚无压力、热流或温度解';
        if state.currentStep==1 && ~isempty(state.rawMesh) && isempty(state.orientedMesh)
            caption='朝向预览待确认：请核对主轴X、鼻尖+X；鼠标旋转仅改变视角';
        elseif ~isempty(state.mesh),caption='已生成真实求解网格；显示旋转不改变坐标或边界';end
        updateAerothermalGeometryView(geometryView,m,caption,struct('showWire',state.currentStep==2,'resetCamera',reset));
    end
    function showResults()
        if isempty(state.result),return;end
        geometryHost.Visible='off';resultHost.Visible='on';
        if isempty(webUI)
            tabs=uitabgroup(resultHost,'Units','normalized','Position',[0 0 1 1]);
            webUI=createAerothermalWebTab(tabs,state.result);
        else
            updateAerothermalWebTab(webUI,state.result);
        end
        sync();
    end

%% 模块六：参数表按当前步骤集中在同一侧边栏；编辑与确认分开。
    function renderEditor()
        delete(editor.Children);controls=struct;
        switch state.currentStep
            case 1
                grid=form(15);
                controls.file=textField(grid,1,'STL 文件',state.cfg.geometry.stlFile);
                controls.scale=numberField(grid,2,'原坐标 → m',state.cfg.geometry.scale,[eps inf]);
                controls.file.ValueChangedFcn=@modelEdited;controls.scale.ValueChangedFcn=@modelEdited;
                controls.browse=button(grid,3,'选择 STL 文件',@browseFile);
                controls.load=button(grid,4,'加载当前模型',@loadInput);
                controls.pca=button(grid,5,'PCA 对齐主轴',@pcaOrientation);
                controls.flip=button(grid,6,'翻转鼻向（Y轴180°）',@flipOrientation);
                controls.rx=numberField(grid,7,'绕 X 增量 / °',0,[-180 180]);
                controls.ry=numberField(grid,8,'绕 Y 增量 / °',0,[-180 180]);
                controls.rz=numberField(grid,9,'绕 Z 增量 / °',0,[-180 180]);
                controls.rotate=button(grid,10,'应用增量旋转预览',@manualOrientation);
                controls.reset=button(grid,11,'恢复导入坐标',@resetOrientation);
                controls.confirm=button(grid,12,'确认主轴 X、鼻尖 +X',@confirmOrientation);
                tip(grid,13,['STL不含单位。默认mm模型转换系数0.001；m模型为1。' newline ...
                    'PCA只识别轴，不能确定鼻向。实际旋转会使网格/结果失效；鼠标操作只改变观察视角。']);
                set([controls.pca controls.flip controls.rotate controls.reset controls.confirm], ...
                    'Enable',onoff(~isempty(state.rawMesh)));
            case 2
                grid=form(12);c=state.cfg.remesh;
                controls.hMin=numberField(grid,1,'最小尺寸 / mm',c.hMin*1000,[.01 inf]);
                controls.hMax=numberField(grid,2,'最大尺寸 / mm',c.hMax*1000,[.01 inf]);
                controls.spacing=numberField(grid,3,'初始轴向间距 / mm',c.initialAxialSpacing*1000,[.01 inf]);
                controls.angle=numberField(grid,4,'曲率目标角 / °',c.curvatureAngleTarget,[.1 90]);
                controls.minCircum=numberField(grid,5,'最少周向点数',c.minCircumPoints,[3 inf]);
                controls.maxCircum=numberField(grid,6,'最多周向点数',c.maxCircumPoints,[3 inf]);
                keys={'hMin','hMax','spacing','angle','minCircum','maxCircum'};
                for j=1:numel(keys),controls.(keys{j}).ValueChangedFcn=@meshEdited;end
                controls.generate=button(grid,7,'生成并检查表面网格',@generate);
                if isempty(state.quality),description='待生成；沿用既有自适应截面与拉链三角化算法。';
                else
                    q=state.quality;description=sprintf('AR 平均 %.4g / 最大 %.4g\n最小角 %.4g°，退化面 %d', ...
                        q.meanAspectRatio,q.maxAspectRatio,q.minimumAngle,q.numberDegenerate);
                end
                tip(grid,8,description);
                tip(grid,11,'参数修改即使旧网格和下游结果失效；不改变现有网格方法。质量差面保留真实指标。');
            case 3
                grid=form(24+12*methodDetailsOpen);c=state.cfg;
                controls.altitude=numberField(grid,1,'高度 / km',c.flight.altitude/1000,[0 86]);
                controls.mach=numberField(grid,2,'Mach（计算需 >1）',c.flight.Mach,[0 30]);
                controls.alpha=numberField(grid,3,'攻角 / °',c.flight.alpha,[-90 90]);
                controls.beta=numberField(grid,4,'侧滑角 / °',c.flight.beta,[-90 90]);
                controls.wall=numberField(grid,5,'给定壁温 / K',c.flight.Twall,[1 10000]);
                controls.air=choice(grid,6,'外流/热流物性',{'NASA7变比热','定比热空气'}, ...
                    {'thermallyPerfect','caloricallyPerfect'},c.heating.airThermalModel);
                controls.nose=choice(grid,7,'鼻区方法',{'驻点余弦对照','球鼻局部相似（0°）'}, ...
                    {'cosineExtension','sphericalLocalSimilarity'},c.heating.noseModel);
                controls.noseStep=numberField(grid,8,'鼻区角度表 / °',getOption(c.heating,'sphericalNoseAngleStepDeg',2),[.25 5]);
                controls.turbulent=choice(grid,9,'湍流侧壁方法',{'Head卷吸积分','固定形状因子对照'}, ...
                    {'head','fixedShape'},c.heating.turbulentModel);
                controls.transition=numberField(grid,10,'经验转捩 Re',c.heating.transitionRe,[1 inf]);
                controls.base=choice(grid,11,'尾底压力',{'给定Cp边界','SAND93层流尖锥'}, ...
                    {'prescribed','sand93Laminar'},getOption(c.basePressure,'model','prescribed'));
                controls.baseCp=numberField(grid,12,'给定尾底 Cp',getOption(getOption(c,'pressure',struct),'aftBaseCp',0),[-2 2]);
                controls.leeCp=numberField(grid,13,'给定背风 Cp',getOption(getOption(c,'pressure',struct),'leewardCp',0),[-2 2]);
                controls.support=choice(grid,14,'尾撑声明',{'未声明','无尾撑','有尾撑'}, ...
                    {'unknown','no','yes'},operation(c.basePressure,'hasRearSupport'));
                controls.injection=choice(grid,15,'底部喷流声明',{'未声明','无喷流','有喷流'}, ...
                    {'unknown','no','yes'},operation(c.basePressure,'hasGasInjection'));
                controls.source=textField(grid,16,'操作声明来源',getOption(c.basePressure,'operationSource',''));
                numeric={'altitude','mach','alpha','beta','wall','noseStep','transition','baseCp','leeCp'};
                controls.sidewall=choice(grid,17,'附着侧壁热流',{'边界层积分法','Kinney参考流长对照'}, ...
                    {'boundaryLayerIntegral','kinneyAcreage'},getOption(c.heating,'sidewallModel','boundaryLayerIntegral'));
                menu={'air','nose','turbulent','base','support','injection','sidewall'};
                for j=1:numel(numeric),controls.(numeric{j}).ValueChangedFcn=@boundaryEdited;end
                for j=1:numel(menu),controls.(menu{j}).ValueChangedFcn=@boundaryEdited;end
                controls.source.ValueChangedFcn=@boundaryEdited;
                controls.confirm=button(grid,18,'确认当前边界条件',@confirmBoundary);
                addMethodScope(grid,19,22,23);
                if isfield(c.heating,'wallTemperatureByFace')
                    tip(grid,23+12*methodDetailsOpen,sprintf('当前保留 %d 个面元壁温输入；修改统一壁温将切换为统一边界。',numel(c.heating.wallTemperatureByFace)));
                end
            case 4
                failed=~isempty(state.error.message);offset=6*failed;
                grid=form(19+12*methodDetailsOpen+offset);controls.run=button(grid,1,'计算当前工况',@solveCase);
                if failed
                    controls.failure=uitextarea(grid,'Editable','off','Value', ...
                        {'本次计算未完成：后端原始原因';state.error.message;state.error.identifier});
                    controls.failure.Layout.Row=[2 6];controls.failure.Layout.Column=[1 2];
                    controls.adjustBoundary=button(grid,7,'返回边界条件调整',@(~,~)go(3));
                end
                tip(grid,2+offset,sprintf('M %.5g，α %.5g°，β %.5g°\n高度 %.5g km，输入壁温 %.5g K', ...
                    state.cfg.flight.Mach,state.cfg.flight.alpha,state.cfg.flight.beta,state.cfg.flight.altitude/1000,state.cfg.flight.Twall));
                addMethodScope(grid,5+offset,8+offset,9+offset);
                offset=offset+12*methodDetailsOpen;
                tip(grid,9+offset,'仅按按钮启动计算。进度为已完成的程序阶段比例，不是面积覆盖率或剩余耗时。');
                controls.log=uitextarea(grid,'Editable','off','Value',logLines);controls.log.Layout.Row=[12+offset 19+offset];controls.log.Layout.Column=[1 2];
            case 5
                grid=form(10);d=buildAerothermalViewerData(state.result);
                h=d.fields(strcmp({d.fields.key},'heatFlux'));p=d.fields(strcmp({d.fields.key},'pressure'));
                tip(grid,1,sprintf('当前版本 %d 的计算快照\n热流计算面积 %.4g%%\n压力计算面积 %.4g%% / 给定 %.4g%%', ...
                    state.resultSourceRevision,100*h.computedAreaFraction,100*p.computedAreaFraction,100*p.prescribedAreaFraction));
                tip(grid,4,'点选面元核对场量、算法和范围。灰色未解；未知外缘热力状态与有效状态的Kn超限分开，两者均可阻断下游热流。输入壁温与求解壁温分开。');
                controls.export=button(grid,7,'导出完整离线工作台',@exportCase);
                tip(grid,8,'任意X截面曲线和导出使用同一快照。完成计算不等于全面积已解；修改上游使结果失效，浏览与导出不重算。');
        end
        sync();
    end
    function grid=form(rows)
        grid=uigridlayout(editor,[rows 2]);grid.ColumnWidth={138,'1x'};
        grid.RowHeight=repmat({31},1,rows);grid.Padding=[0 2 0 2];grid.RowSpacing=5;
        if isprop(grid,'Scrollable'),grid.Scrollable='on';end
    end
    function x=numberField(grid,row,label,value,limits)
        lab=uilabel(grid,'Text',label);lab.Layout.Row=row;lab.Layout.Column=1;
        x=uieditfield(grid,'numeric','Value',value,'Limits',limits);x.Layout.Row=row;x.Layout.Column=2;
    end
    function x=textField(grid,row,label,value)
        lab=uilabel(grid,'Text',label);lab.Layout.Row=row;lab.Layout.Column=1;
        x=uieditfield(grid,'text','Value',char(value));x.Layout.Row=row;x.Layout.Column=2;
    end
    function x=choice(grid,row,label,items,values,value)
        lab=uilabel(grid,'Text',label);lab.Layout.Row=row;lab.Layout.Column=1;
        x=uidropdown(grid,'Items',items,'ItemsData',values,'Value',value);x.Layout.Row=row;x.Layout.Column=2;
    end
    function x=button(grid,row,text,callback)
        x=uibutton(grid,'Text',text,'ButtonPushedFcn',callback);x.Layout.Row=row;x.Layout.Column=[1 2];
    end
    function tip(grid,row,text)
        x=uilabel(grid,'Text',text,'WordWrap','on','VerticalAlignment','top');
        x.Layout.Row=[row min(row+2,size(grid.RowHeight,2))];x.Layout.Column=[1 2];
    end

%% 模块六补充：方法说明可折叠；只读实际配置与几何，不求物理场或审核通过值。
    function addMethodScope(grid,summaryRow,toggleRow,detailRow)
        controls.methodSummary=uilabel(grid,'WordWrap','on','VerticalAlignment','top');
        controls.methodSummary.Layout.Row=[summaryRow summaryRow+2];controls.methodSummary.Layout.Column=[1 2];
        toggleText='算法与适用范围 ▸';if methodDetailsOpen,toggleText='算法与适用范围 ▾';end
        controls.methodToggle=button(grid,toggleRow,toggleText,@toggleMethodScope);
        if methodDetailsOpen
            controls.methodScope=uitextarea(grid,'Editable','off');
            controls.methodScope.Layout.Row=[detailRow detailRow+11];controls.methodScope.Layout.Column=[1 2];
        end
        updateMethodScope();
    end
    function toggleMethodScope(~,~)
        if state.busy || preparing,return;end
        methodDetailsOpen=~methodDetailsOpen;renderEditor();
    end
    function updateMethodScope()
        [summary,lines]=methodScopeText(state.cfg,state.mesh);
        if isfield(controls,'methodSummary'),controls.methodSummary.Text=summary;end
        if isfield(controls,'methodScope'),controls.methodScope.Value=lines;end
        % 所选关联不会读取给定底Cp；保留输入配置，但不能暗示它是求解值。
        if isfield(controls,'baseCp')
            controls.baseCp.Enable=onoff(strcmp(getOption(state.cfg.basePressure,'model','prescribed'),'prescribed'));
        end
    end

%% 模块七：真实模型输入及朝向确认，不以观察相机变化代替实际坐标变换。
    function modelEdited(~,~)
        if state.busy || preparing,return;end
        next=state.cfg;next.geometry.stlFile=controls.file.Value;next.geometry.scale=controls.scale.Value;
        state=aeroHeatWorkflowState('invalidateModel',state,struct('cfg',next));
        preview=[];report=[];rotation=eye(3);message('文件或尺度已编辑，旧网格/边界确认/结果失效；请重新加载。');showGeometry(true);
    end
    function browseFile(~,~)
        [name,path]=uigetfile('*.stl','选择几何模型');if isequal(name,0),return;end
        controls.file.Value=fullfile(path,name);modelEdited([],[]);loadInput([],[]);
    end
    function loadInput(~,~)
        if state.busy || preparing,return;end
        try
            next=state.cfg;next.geometry.stlFile=controls.file.Value;next.geometry.scale=controls.scale.Value;
            state=aeroHeatWorkflowState('invalidateModel',state,struct('cfg',next));
            preview=[];report=[];preparing=true;message('正在加载模型并核对SI尺度…');
            raw=quietCall(@()loadSTLMesh(next));raw=quietCall(@()calcFaceGeometry(raw));
            state=aeroHeatWorkflowState('setModel',state,struct('rawMesh',raw,'cfg',next, ...
                'modelSource',struct('file',next.geometry.stlFile,'scale',next.geometry.scale,'unit','m')));
            preview=raw;rotation=eye(3);rotationCenter=mean(raw.V,1);report=[];
            preparing=false;renderEditor();message('模型已按所选尺度转换为m；朝向尚未确认。');showGeometry(true);
        catch ex,preparing=false;renderEditor();showError(ex);showGeometry(true);end
    end
    function previewRotation(candidate)
        state=aeroHeatWorkflowState('invalidateOrientation',state);
        [preview,~]=quietPair(@()alignGeometryAxis(state.rawMesh,'matrix',candidate,rotationCenter));
        rotation=candidate;report=[];message('实际朝向预览已改变；确认后才可生成网格。');showGeometry(true);
    end
    function pcaOrientation(~,~)
        if isempty(state.rawMesh) || state.busy || preparing,return;end
        try
            [~,info]=quietPair(@()alignGeometryAxis(state.rawMesh,'pca'));
            rotationCenter=info.rotationCenter;previewRotation(info.rotationMatrix);
        catch ex,showError(ex);end
    end
    function flipOrientation(~,~)
        if isempty(state.rawMesh) || state.busy || preparing,return;end
        try,previewRotation([-1 0 0;0 1 0;0 0 -1]*rotation);catch ex,showError(ex);end
    end
    function manualOrientation(~,~)
        if isempty(state.rawMesh) || state.busy || preparing,return;end
        try
            rx=controls.rx.Value;ry=controls.ry.Value;rz=controls.rz.Value;
            X=[1 0 0;0 cosd(rx) -sind(rx);0 sind(rx) cosd(rx)];
            Y=[cosd(ry) 0 sind(ry);0 1 0;-sind(ry) 0 cosd(ry)];
            Z=[cosd(rz) -sind(rz) 0;sind(rz) cosd(rz) 0;0 0 1];
            previewRotation(Z*Y*X*rotation);
            controls.rx.Value=0;controls.ry.Value=0;controls.rz.Value=0;
        catch ex,showError(ex);end
    end
    function resetOrientation(~,~)
        if isempty(state.rawMesh) || state.busy || preparing,return;end
        try,previewRotation(eye(3));catch ex,showError(ex);end
    end
    function confirmOrientation(~,~)
        if isempty(preview) || state.busy || preparing,return;end
        try
            next=state.cfg;next.geometry.rotationMatrix=rotation;next.geometry.rotationCenter=rotationCenter;
            next.geometry.bodyAxis=[1 0 0];next.geometry.noseDirection=[1 0 0];
            next.geometry.isAligned=true;next.geometry.alignmentMethod='workbenchConfirmed';
            state=aeroHeatWorkflowState('setOrientation',state,struct('orientedMesh',preview,'cfg',next));
            message('用户已确认主轴X、鼻尖+X；可进入生成网格。');go(2);
        catch ex,showError(ex);end
    end

%% 模块八：保留原网格方法；参数变化、生成和求解输入版本明确关联。
    function next=meshConfiguration()
        next=state.cfg;next.remesh.hMin=controls.hMin.Value*.001;next.remesh.hMax=controls.hMax.Value*.001;
        next.remesh.initialAxialSpacing=controls.spacing.Value*.001;next.remesh.curvatureAngleTarget=controls.angle.Value;
        next.remesh.minCircumPoints=controls.minCircum.Value;next.remesh.maxCircumPoints=controls.maxCircum.Value;
    end
    function meshEdited(~,~)
        if state.busy || preparing,return;end
        try
            state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',meshConfiguration()));report=[];
            message('网格设置已修改，旧网格与下游结果失效；请重新生成。');showGeometry(false);
        catch ex,showError(ex);end
    end
    function generate(~,~)
        if state.busy || preparing,return;end
        try
            next=meshConfiguration();validateMeshOptions(next.remesh);
            state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',next));
            preparing=true;meshRevision=state.revision;oriented=state.orientedMesh;
            controls.generate.Enable='off';message('正在按既有算法重构表面网格…');drawnow;
            candidate=quietCall(@()autoRemeshSurface(oriented,next));
            candidate=quietCall(@()calcFaceGeometry(candidate));
            candidate=quietCall(@()fixNormalDirection(candidate,next));
            quality=quietCall(@()checkMeshQuality(candidate,next));
            candidate=quietCall(@()buildSurfaceCoordinates(candidate,next));
            candidate=analyzeLocalGeometry(candidate,next);
            edges=[candidate.F(:,[1 2]);candidate.F(:,[2 3]);candidate.F(:,[3 1])];
            [~,~,ids]=unique(sort(edges,2),'rows');
            if any(accumarray(ids,1)~=2) || quality.numberDegenerate>0
                error('RocketAeroHeat:Workbench:MeshTopology','生成网格含开口/非流形边或退化面，不接受为求解网格。');
            end
            if state.revision~=meshRevision,error('RocketAeroHeat:Workbench:StaleMesh','生成期间输入改变，拒绝旧网格。');end
            state=aeroHeatWorkflowState('setMesh',state,struct('mesh',candidate,'cfg',next,'quality',quality));
            preparing=false;
            renderEditor();message(sprintf('网格已生成并核对：%d面元；请添加边界条件。',size(candidate.F,1)));showGeometry(true);
        catch ex
            preparing=false;
            renderEditor();
            showError(ex);
        end
    end

%% 模块九：边界输入与热物性/方法选择均入集中状态，不立即求物理场。
    function next=boundaryConfiguration()
        next=state.cfg;altitude=controls.altitude.Value*1000;
        if altitude~=next.flight.altitude && isfield(next.flight,'atmosphere'),next.flight=rmfield(next.flight,'atmosphere');end
        next.flight.altitude=altitude;next.flight.Mach=controls.mach.Value;
        next.flight.alpha=controls.alpha.Value;next.flight.beta=controls.beta.Value;
        if controls.wall.Value~=next.flight.Twall && isfield(next.heating,'wallTemperatureByFace')
            next.heating=rmfield(next.heating,'wallTemperatureByFace');
        end
        next.flight.Twall=controls.wall.Value;next.heating.airThermalModel=controls.air.Value;
        next.heating.noseModel=controls.nose.Value;next.heating.sphericalNoseAngleStepDeg=controls.noseStep.Value;
        next.heating.turbulentModel=controls.turbulent.Value;next.heating.transitionRe=controls.transition.Value;
        next.heating.sidewallModel=controls.sidewall.Value;
        next.basePressure.model=controls.base.Value;next.basePressure.operationSource=controls.source.Value;
        next.pressure.aftBaseCp=controls.baseCp.Value;next.pressure.leewardCp=controls.leeCp.Value;
        for key={'hasRearSupport','hasGasInjection'}
            control=controls.support;if strcmp(key{1},'hasGasInjection'),control=controls.injection;end
            if strcmp(control.Value,'unknown')
                if isfield(next.basePressure,key{1}),next.basePressure=rmfield(next.basePressure,key{1});end
            else,next.basePressure.(key{1})=strcmp(control.Value,'yes');end
        end
    end
    function boundaryEdited(~,~)
        if state.busy || preparing,return;end
        try
            state=aeroHeatWorkflowState('setBoundary',state,struct('cfg',boundaryConfiguration()));report=[];
            updateMethodScope();message('边界设置已编辑；请重新确认并显式计算，旧结果不会冒充新工况。');showGeometry(false);
        catch ex,showError(ex);end
    end
    function confirmBoundary(~,~)
        if state.busy || preparing,return;end
        try
            next=boundaryConfiguration();fs=calcFreestream(next.flight);
            if fs.Mach<=1,error('RocketAeroHeat:Workbench:Mach','当前完整工程链要求M>1。');end
            state=aeroHeatWorkflowState('setBoundary',state,struct('cfg',next));
            state=aeroHeatWorkflowState('confirmBoundary',state,struct('freestream',fs));
            if isempty(state.result)
                message(sprintf('边界已确认：p∞ %.6gPa、T∞ %.6gK；尚未计算压力或热流。',fs.pressure,fs.temperature));
            else
                message('边界与当前计算快照一致，结果仍有效；再次计算需明确按计算按钮。');
            end
            go(4);
        catch ex,showError(ex);end
    end

%% 模块十：唯一物理计算入口；输入冻结，阶段失败不接受中间快照。
    function solveCase(~,~)
        if state.busy || preparing,return;end
        try
            state=aeroHeatWorkflowState('beginSolve',state);token=state.solveId;revision=state.solveInputRevision;
            renderEditor();message('正在显式计算当前工况；各阶段完成不代表全面积已求解。');
            frozen=state.frozenInputs;
            [solved,fs,newReport]=quietTriple(@()runAeroHeatSingleCase(frozen.mesh,frozen.cfg,@progress));
            state=aeroHeatWorkflowState('completeSolve',state,struct('result',solved, ...
                'freestream',fs,'solveId',token,'sourceRevision',revision));report=newReport;
            message(sprintf('当前计算快照已接受；热流计算面积 %.4g%%，未解区仍保留。',100*newReport.heatComputedAreaFraction));
            go(5);
        catch ex
            if exist('token','var') && state.busy
                state=aeroHeatWorkflowState('failSolve',state,struct('error',ex,'solveId',token,'sourceRevision',revision));
            end
            if isgraphics(fig),renderEditor();showError(ex);showGeometry(false);end
        end
    end
    function progress(event)
        if ~isgraphics(fig),error('RocketAeroHeat:Workbench:Closed','工作台已关闭，停止当前计算。');end
        if strcmp(event.status,'started') || strcmp(event.status,'failed'),message(event.message);end
        drawnow;
    end
    function exportCase(~,~)
        if isempty(state.result),return;end
        folder=uigetdir(state.cfg.projectRoot,'导出当前完整离线工作台');if isequal(folder,0),return;end
        try
            exportAerothermalWeb(state.result,fullfile(folder,'RocketAeroHeat_case'));
            message('已导出当前计算快照和逐场量来源；未重算物理。');
        catch ex,showError(ex);end
    end
end

function [summary,lines]=methodScopeText(cfg,mesh)
% 模块用途：说明所选方法的实际范围和已知输入冲突，不预先宣告区域求解通过。
% 此处只读取配置/做已有几何初检；Kn、肩部Re、物理厚度和BVP仍由显式计算审核。
h=cfg.heating;f=cfg.flight;bp=cfg.basePressure;
base=getOption(bp,'model','prescribed');nose=getOption(h,'noseModel','cosineExtension');
thermal=getOption(h,'airThermalModel','thermallyPerfect');
summary='壁温和背风Cp为给定边界；全面积热流尚未完成。';
current={};geometry=[];
if ~isempty(mesh) && (strcmp(base,'sand93Laminar') || strcmp(nose,'sphericalLocalSimilarity'))
    geometry=estimateNoseRadius(mesh);
end
lines={'范围说明是输入提示，不是本工况已通过审核。只有按“计算当前工况”才求物理场。';''; ...
    '压力与输入身份'; ...
    '迎风压力：修正牛顿法，按正激波驻点压力上限计算Cp；属于高超声速工程近似，不解遮挡、激波干扰或黏性压力。'; ...
    '背风Cp和“给定Cp边界”尾底Cp由用户指定；Cp=0表示p=p∞，不是程序求得尾迹压力为零。'; ...
    '壁温是热流边界输入。本次定常链路不求解壁体温度；外缘温度、恢复温度和求解壁温在结果中分别列出。';''};
if strcmp(base,'sand93Laminar')
    lines=[lines;{'尾底压力：SAND93-0280层流尖锥关联'; ...
        '只映射孤立、闭合平圆底、同轴尖头直圆锥，半锥角9–10°。真实V/F须审核通过；球鼻、锥柱、尾锥不是这个几何族。'; ...
        '零攻角、零侧滑；显式定比热γ=1.4空气，2≤Mach≤20；自由流T≥80K、总温≤2500K、给定等温前体壁温80–2500K。温度门限不是恒比热空气的实验准确度保证。'; ...
        '肩部Taylor–Maccoll外流须附体，按经验阈值仍为层流且Re≥10⁵；内部物理剖面需δ/rb≤0.05，连续流筛选Kn<0.01。δ/rb门限为项目保守筛选。'; ...
        '须明确无尾撑、无底部喷流并填写操作来源；STL不能证明发动机是否在喷流。'; ...
        '关联拒绝时不改用给定Cp：完整识别底面的拒绝面元为NaN；无法完整识别尾底目标时停止本次计算。代表底热不映射到面元，底热和相关底区温度仍未解。'}];
    if abs(f.alpha)>1e-8 || abs(f.beta)>1e-8
        current{end+1,1}=sprintf('当前α %.5g°、β %.5g°不满足SAND93零攻角/侧滑条件。',f.alpha,f.beta);
    end
    if ~strcmp(thermal,'caloricallyPerfect')
        current{end+1,1}='当前变比热物性不满足SAND93显式定比热要求。';
    end
    if f.Mach<2 || f.Mach>20,current{end+1,1}='当前Mach超出SAND93的2–20范围。';end
    if ~isempty(geometry) && geometry.available
        current{end+1,1}='当前网格球鼻几何初检可用；SAND93尖直锥底压不支持这类球鼻箭体。';
    end
    if ~strcmp(operation(bp,'hasRearSupport'),'no') || ...
            ~strcmp(operation(bp,'hasGasInjection'),'no') || isempty(strtrim(getOption(bp,'operationSource','')))
        current{end+1,1}='SAND93还需明确无尾撑/无喷流及操作声明来源。';
    end
else
    lines=[lines;{'尾底压力：给定Cp边界'; ...
        '用于非迎风尾盖；与背风Cp一样，属于输入，不是底区回流求解。底热和相关底区温度仍未解。'}];
end
lines=[lines;{''}];
if strcmp(nose,'sphericalLocalSimilarity')
    limit=getOption(getOption(cfg,'externalFlow',struct),'sphericalNoseMaxPolarAngleDeg',60);
    lines=[lines;{'鼻区：球鼻驻点熵外缘 + Levy–Lees局部相似'; ...
        sprintf('只在重新核验的光滑球面前缀计算；实际来流沿−X（零合成攻角），当前迎流极角筛选≤%.5g°，不延拓到后体。',limit); ...
        sprintf('附着层流、Kn<%.5g、Re<%.5g、δ/R≤0.05；壁温与有效物性80–2500K。单时刻最多16种给定壁温，只在相邻有效角度节点间插值。',getOption(h,'maxContinuumKnudsen',.01),getOption(h,'transitionRe',5e5)); ...
        '热流由同一BVP总焓壁面导数求得，恢复温度由q(Twall)=0求根；不把给定壁温或总温标成求解壁温。'; ...
        '不解沿程非相似、熵层吞入、横向涡量、鼻部/后体接口或转捩。拒绝的球鼻面元保持NaN，不用余弦对照偷偷补值。当前没有实验精度验证。'}];
    if abs(f.alpha)>1e-8 || abs(f.beta)>1e-8
        current{end+1,1}='当前攻角/侧滑非零，球鼻轴对称局部相似热边界不适用。';
    end
    if ~isempty(geometry)
        if geometry.available
            lines{end+1,1}=sprintf('当前仅几何初检：球鼻半径 %.6g m；完整面元范围和物理筛选仍待显式计算。',geometry.radius);
        else
            current{end+1,1}=['当前球鼻几何初检未通过：' geometry.status];
        end
    end
else
    lines=[lines;{'鼻区：Sutton–Graves驻点值 + 球面余弦延拓对照'; ...
        '仅已核验球形钝鼻；hypot(α,β)≤30°且通过前半球位置/法向检查，轴向及相对偏置驻点极角均≤80°。仍须驻点热流、连续流和物性审核通过。'; ...
        '这是驻点值锚定的低阶空间延拓，不是局部边界层剖面解，没有球面分布实验验证。'}];
end
lines=[lines;{'';'附着侧壁与转捩'; ...
    sprintf('层流沿准流线Thwaites型积分，热流用Eckert参考温度与Reynolds比拟；经验转捩Re=%.5g，不能预测真实扰动/粗糙度控制的转捩位置。',getOption(h,'transitionRe',5e5))}];
if strcmp(getOption(h,'sidewallModel','boundaryLayerIntegral'),'kinneyAcreage')
    lines=[lines;{'当前侧壁热流选Kinney/CBAero型参考流长对照（2007，II.I，印刷页4–5）。'; ...
        '原Thwaites/Head仅作入口、转捩和分离审核，热流改用公开Cf/St/参考密度公式；迎风几何流长、背风rho*mu*U加权流长。'; ...
        'Mangler因子默认层流3/湍流2，背风5°内smoothstep切换；首段局部常状态积分是项目假设，并非完整CBAERO复现或实验验证。'; ...
        '外缘未知、分离、域外及底面保持原未解；本方法仅验收单时刻，不用于未核验的热历史/壁体耦合。'}];
end
if strcmp(getOption(h,'turbulentModel','head'),'head')
    lines=[lines;{'当前湍流方法：Head卷吸 + von Kármán动量积分；H≥2.4为疑似分离，尚无再附模型，下游保留未解。'}];
else
    lines=[lines;{'当前湍流方法：固定H=1.4动量积分对照；不预测湍流分离、再附或三维横流。'}];
end
lines=[lines;{'Ludwieg–Tillmann关联仅在300≤Reθ≤10⁶使用；范围外明确标记局部平板回退，Head在范围内侧有过渡带。外流/路径/连续流失效或分离不由平板补齐。';'';'空气热物性'}];
if strcmp(thermal,'thermallyPerfect')
    lines=[lines;{'当前：固定组分干空气NASA7，150–2500K；80–150K为振动冻结低温近似，150–200K含少量外推；不含解离、电离或化学非平衡。'; ...
        '局部/轴向域外行为由后端逐项标记，不能因选择变比热就认定所有压力和外流算法均使用变比热。'}];
else
    lines=[lines;{'当前：定比热cp=γR/(γ−1)，高温物性误差未修正。选择定比热不自动满足SAND93的几何、层流、薄层或连续流条件。'}];
end
lines=[lines;{sprintf('当前连续流热模型要求Kn<%.5g；这是尺度筛选，不是滑移、过渡流或稀薄流求解。',getOption(h,'maxContinuumKnudsen',.01)); ...
    '只有有效外缘热力状态才能用于Kn流态判断。未知/域外的定比热流向占位不算已求外缘压力或温度，不由占位产生的大Kn宣称真实稀薄流。'; ...
    '外缘热力状态未知及其下游影响与有效状态的Kn超限分别记录；两类限制均阻断附着热流，不能在下游局部状态重新有效时跨过未知上游补值。'; ...
    ''; '后端来源：modifiedNewtonianPressure、solveMeshBasePressure、solveSphericalNoseLaminarHeating、calcBluntNoseHeating、calcSurfaceHeating、calcTurbulentBoundaryLayerState、selectAirThermalModel。'; ...
    '依据与证据边界详见docs/WORKFLOW_METHOD_SCOPE.md；逐面元最终方法/状态以当前结果探针为准。'}];
if ~isempty(current)
    summary=current{1};
    if numel(current)>1,summary=[summary newline current{2}];end
    lines=[{'当前已知输入冲突（尚未计算）：'};current;{''};lines];
else
    summary=[summary newline '已知输入提示不替代后端几何、连续流及剖面审核。'];
end
end

function value=onoff(flag)
if flag,value='on';else,value='off';end
end
function value=getOption(s,key,fallback)
value=fallback;if isstruct(s) && isfield(s,key),value=s.(key);end
end
function value=operation(s,key)
value='unknown';if isfield(s,key) && isscalar(s.(key)) && ismember(s.(key),[0 1])
    value='no';if logical(s.(key)),value='yes';end
end
end
function validateMeshOptions(r)
if r.hMin<=0 || r.hMin>=r.hMax || r.initialAxialSpacing<=0 || ...
        r.minCircumPoints<3 || r.maxCircumPoints<r.minCircumPoints || ...
        any([r.minCircumPoints r.maxCircumPoints]~=fix([r.minCircumPoints r.maxCircumPoints]))
    error('RocketAeroHeat:Workbench:MeshOptions','网格需0<hMin<hMax、正轴向间距和合法整数周向节点数。');
end
end
function out=quietCall(fn)
evalc('out=fn();');
end
function [first,second]=quietPair(fn)
evalc('[first,second]=fn();');
end
function [first,second,third]=quietTriple(fn)
evalc('[first,second,third]=fn();');
end
