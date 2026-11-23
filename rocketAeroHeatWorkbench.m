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
embeddedVolume=[];engineJob=[];engineTimer=[];researchUI=[];closing=false;
measuredDirectory='';measuredPauseFile='';lastWireDirectory='';
creasePreview=false;qualityPreview=false;candidatePreview=[];
fig=uifigure('Name','RocketAeroHeat · 科研工作台', ...
    'Position',[30 30 1580 940],'Color',[.96 .97 .98]);
fig.CloseRequestFcn=@closeWorkbench;
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
geometryReceiver=geometryView.html.DataChangedFcn;
geometryView.html.DataChangedFcn=@geometryMessage;
footer=uilabel(frame,'WordWrap','on');footer.Layout.Row=3;footer.Layout.Column=[1 2];
fig.UserData=struct('workflowState',state,'workflowButtons',buttons, ...
    'workflowSidebar',sidebar,'geometryView',geometryView,'resultHost',resultHost, ...
    'controls',controls,'webUI',[],'caseReport',[],'title',title,'status',footer,'editorPanel',editor);
fig.UserData.actions=struct('selectStep',@go,'loadModel',@()loadInput([],[]), ...
    'pca',@()pcaOrientation([],[]),'flip',@()flipOrientation([],[]), ...
    'confirmOrientation',@()confirmOrientation([],[]),'generateMesh',@()generate([],[]), ...
    'confirmBoundary',@()confirmBoundary([],[]),'run',@()solveCase([],[]),'getState',@getState, ...
    'loadNative',@()loadNative([],[]),'cancel',@()cancelEngine([],[]),'readResearch',@readResearch, ...
    'readTimeoutCheckpoint',@readTimeoutCheckpoint,'useSurfaceForSU2',@useSurfaceForSU2, ...
    'loadNativeFeature',@()loadNative([],[],true),'loadClassic',@loadClassic, ...
    'setWallTemperatureProfile',@setWallTemperatureProfile,'readHybridMesh',@readHybridMesh, ...
    'readVerifiedHEOS',@readVerifiedHEOS,'readSSTStudy',@readSSTStudy,'probeMeshQuality',@probeMeshQuality, ...
    'inspectMeshCandidate',@inspectMeshCandidate,'readBenchmarkReview',@readBenchmarkReview, ...
    'readIdealAirReview',@readIdealAirReview,'readStaticAirProbe',@readStaticAirProbe);
renderEditor();sync();showGeometry(true);

%% 模块五：切页只读取当前状态；所有结果入口依赖已完成的当前版本。
    function go(step)
        if state.busy || preparing || ~isgraphics(fig),return;end
        try
            state=aeroHeatWorkflowState('selectStep',state,struct('step',step));
            delete(findall(resultHost,'Tag','static-air-probe-host'));
            renderEditor();sync();
            if step==5,showResults();else,showGeometry(false);end
        catch ex,showError(ex);end
    end
    function s=getState(),s=state;end
    function sync()
        if ~isgraphics(fig),return;end
        if isempty(state.result) && ~isempty(webUI)
            % 输入失效后释放旧结果页，防止隐藏页探针继续指向上一次计算。
            if isgraphics(webUI.tab),delete(webUI.tab.Parent);end
            webUI=[];
        end
        if ~isempty(candidatePreview) && candidatePreview.revision~=state.revision,candidatePreview=[];qualityPreview=false;end
        % 父面板也必须恢复；仅重建子控件不会解除父面板的禁用状态。
        if isprop(editor,'Enable'),editor.Enable=onoff(~state.busy && ~preparing);end
        for j=1:5
            tag='待前置';
            if state.steps(j).available,tag='待处理';end
            if state.steps(j).completed,tag='已完成';end
            if strcmp(state.steps(j).status,'busy'),tag='计算中';end
            if strcmp(state.steps(j).status,'failed'),tag='失败';end
            if strcmp(state.steps(j).status,'paused'),tag='已保存';end
            if strcmp(state.steps(j).status,'research'),tag='研究快照';end
            buttons(j).Text=sprintf('%d  %s  · %s',j,state.steps(j).label,tag);
            buttons(j).Enable=onoff(state.steps(j).available && ~state.busy && ~preparing);
            buttons(j).BackgroundColor=[.96 .97 .98];
        end
        buttons(state.currentStep).BackgroundColor=[.84 .90 .94];
        title.Text=sprintf('%s · 输入版本 %d',state.steps(state.currentStep).label,state.revision);
        counts='尚无求解网格';
        if ~isempty(state.mesh),counts=sprintf('求解网格 %d 节点 / %d 面元',size(state.mesh.V,1),size(state.mesh.F,1));end
        inputSummary.Text=[counts newline 'SI · 主轴 X · 鼻尖 +X · 零攻角来流 −X'];
        if isBead()
            counts=sprintf('B 型热敏珠 · %d 个原面',size(state.rawMesh.F,1));
            if ~isempty(state.mesh),counts=sprintf('%s · %d 体单元',counts,state.mesh.thermistorPrepared.cells);end
            inputSummary.Text=[counts newline 'SI · 长轴 X · 来流 −X'];
            if isSidewind(),inputSummary.Text=[counts newline 'SI · 长轴 X · 侧吹来流 −Y'];end
        end
        footer.Text=logLines{end};
        fig.UserData.workflowState=state;fig.UserData.controls=controls;
        fig.UserData.webUI=webUI;fig.UserData.caseReport=report;
        fig.UserData.engineJob=engineJob;fig.UserData.researchUI=researchUI;
        fig.UserData.volumeInfo=[];
        if ~isempty(embeddedVolume)
            fig.UserData.volumeInfo=struct('nodes',size(embeddedVolume.vertices,1), ...
                'cells',numel(embeddedVolume.cells),'construction',embeddedVolume.meshConstruction);
        end
        % 工作区也跟随当前有效状态；上游改变后不留下旧场量冒充当前结果。
        assignin('base','cfg',state.cfg);assignin('base','rawMeshOriginal',state.rawMesh);
        assignin('base','rawMesh',state.orientedMesh);assignin('base','quality',state.quality);
        current=state.mesh;if ~isempty(state.result),current=state.result;end
        assignin('base','mesh',current);assignin('base','freestream',state.freestream);
        if isBead()
            assignin('base','flight',state.cfg.thermistorBoundary);
        elseif isCM011()
            % 工程链的飞行设置留在 cfg.flight；独立工作区变量显示当前试件边界。
            assignin('base','flight',getOption(state.cfg,'cm011Boundary',[]));
        elseif isfield(state.cfg,'flight'),assignin('base','flight',state.cfg.flight);end
        if state.busy || preparing
            list=findall(editor);
            for q=1:numel(list),if isprop(list(q),'Enable'),list(q).Enable='off';end,end
        end
        % 异步计算期间保留取消按钮；父面板禁用会一并禁用子控件。
        if state.busy && ~isempty(engineJob) && isfield(controls,'cancel')
            editor.Enable='on';controls.cancel.Enable=onoff(~engineJob.cancelRequested);
        end
        if state.busy && ~isempty(measuredPauseFile) && isfield(controls,'pauseMeasured')
            editor.Enable='on';controls.pauseMeasured.Enable=onoff(~isfile(measuredPauseFile));
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
        delete(findall(resultHost,'Tag','thermistor-bead-review'));fig.UserData.thermistorReview=[];
        delete(findall(resultHost,'Tag','lewis-cylinder-review'));fig.UserData.lewisCylinderReview=[];
        resultHost.Visible='off';geometryHost.Visible='on';
        fig.UserData.benchmarkReview=[];
        fig.UserData.idealAirReview=[];
        fig.UserData.meshQuality=[];
        if qualityPreview && state.currentStep==2 && (~isempty(embeddedVolume) || ~isempty(candidatePreview))
            % 网格质量来自真实四面体；灰色原壁面只用于空间定位。
            if isempty(candidatePreview)
                [m,info]=buildSU2MeshQualityPreview(embeddedVolume);
                caption=sprintf('低质量体单元 %d 个 · 最小SICN %.5g · 仅网格诊断',info.lowQualityCells,info.minimumSICN);
            else
                m=candidatePreview.mesh;info=candidatePreview.info;
                caption=sprintf('失败候选 · %s · %d个低质四面体 · 最小SICN %.5g · 只读，未接入计算',info.exteriorAlgorithm,info.lowQualityCells,info.minimumSICN);
            end
            updateAerothermalGeometryView(geometryView,m,caption,struct('showWire',true,'resetCamera',reset));
            s=geometryView.html.UserData;s.revision=s.revision+1;s.payload.revision=s.revision;s.rendered=[];
            s.payload.faceField=struct('name','四面体SICN','unit','1', ...
                'values',info.displayValues,'range',[0 .1],'productionQualified',false,'context','mesh-quality');
            geometryView.html.UserData=s;geometryView.html.Data=s.payload;
            fig.UserData.meshQuality=info;fig.UserData.creaseSlice=[];return;
        end
        if creasePreview && state.currentStep==2 && ~isempty(embeddedVolume)
            % 从已审核的真实体连接显示一列棱柱扇；切页仍保留原壁面与工况。
            [section,info]=buildSU2CreaseFanSlice(embeddedVolume);
            caption=sprintf('尾缘棱柱扇截面 · 原壁点%d · 仅几何，尚无压力/热流/温度解',info.sourceVertexID);
            updateAerothermalGeometryView(geometryView,section,caption,struct('showWire',true,'resetCamera',reset));
            fig.UserData.creaseSlice=info;return;
        end
        fig.UserData.creaseSlice=[];
        m=preview;
        if state.currentStep>=2
            m=state.orientedMesh;if ~isempty(state.mesh),m=state.mesh;end
        end
        caption='仅几何预览；尚无压力、热流或温度解';
        if state.currentStep==1 && ~isempty(state.rawMesh) && isempty(state.orientedMesh)
            caption='朝向预览待确认：请核对主轴X、鼻尖+X；鼠标旋转仅改变视角';
        elseif ~isempty(state.mesh),caption='已生成真实求解网格；显示旋转不改变坐标或边界';end
        if isBead()
            caption='B 型热敏珠 · 米制原面 · 沿长轴来流';
            if isSidewind(),caption='带线 B 型热敏珠 · 米制原面 · 侧向来流 −Y';end
            if state.currentStep==2 && ~isempty(state.mesh)
                p=state.mesh.thermistorPrepared;m=p.preview;
                caption=sprintf('实际三维体网的 30° 壁面扇区 · 周向单元 %g° · %d 个原生壁面',30/p.angularCells,p.wallFaces);
                if isFullBead(),caption=sprintf('360° 全周体网 · 56 层 · %d 个原生壁面',p.wallFaces);end
                if isSidewind(),caption=sprintf('侧吹全周体网 · %d 个周向间隔 · %d 个原生壁面',p.angularCells,p.wallFaces);end
            end
        end
        updateAerothermalGeometryView(geometryView,m,caption,struct('showWire',state.currentStep==2,'resetCamera',reset));
    end
    function showResults()
        delete(findall(resultHost,'Tag','thermistor-bead-review'));fig.UserData.thermistorReview=[];
        delete(findall(resultHost,'Tag','lewis-cylinder-review'));fig.UserData.lewisCylinderReview=[];
        if ~isempty(state.researchResult)
            geometryHost.Visible='off';resultHost.Visible='on';
            delete(resultHost.Children);webUI=[];
            % 只携带当前同一V/F原面的已有几何标签，不把工程流场资格移给研究场。
            displayResult=state.researchResult;
            if isfield(state.mesh,'localGeometry') && isfield(state.mesh.localGeometry,'region') && ...
                    isequal(displayResult.V,state.mesh.V) && isequal(displayResult.F,state.mesh.F)
                displayResult.geometryRegion=state.mesh.localGeometry.region;
                displayResult.geometryRegionSource='同一原V/F表面的已有几何标签；不代表流动区域资格';
            end
            researchUI=createSU2ResearchView(resultHost,displayResult);sync();return;
        end
        if isempty(state.result),return;end
        if ~isempty(researchUI),delete(resultHost.Children);researchUI=[];webUI=[];end
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
                if isBead()
                    grid=form(6);
                    if hasBeadLeads(),tip(grid,1,'B 型热敏珠 · 两端引线直径 0.0540 mm，各含 0.300 mm 名义直段和暴露端面。');
                    else,tip(grid,1,'B 型热敏珠 · 光滑名义珠体，沿 X 轴放置。长度和三角面采用公开试件的米制重建几何。');end
                    controls.nextMesh=button(grid,4,'选择周向网格',@(~,~)go(2));
                    controls.changeModel=button(grid,5,'更换模型',@leaveBead);sync();return;
                end
                grid=form(23);
                controls.file=textField(grid,1,'STL 文件',state.cfg.geometry.stlFile);
                controls.scale=numberField(grid,2,'原坐标 → m',state.cfg.geometry.scale,[eps inf]);
                controls.file.ValueChangedFcn=@modelEdited;controls.scale.ValueChangedFcn=@modelEdited;
                controls.browse=button(grid,3,'选择 STL 文件',@browseFile);
                controls.browse.Layout.Column=1;
                controls.beadModel=button(grid,3,'选用 B 型热敏珠',@loadBeadModel);controls.beadModel.Layout.Column=2;
                controls.load=button(grid,4,'加载当前模型',@loadInput);
                controls.pca=button(grid,5,'PCA 对齐主轴',@pcaOrientation);
                controls.flip=button(grid,6,'翻转鼻向（Y轴180°）',@flipOrientation);
                controls.rx=numberField(grid,7,'绕 X 增量 / °',0,[-180 180]);
                controls.ry=numberField(grid,8,'绕 Y 增量 / °',0,[-180 180]);
                controls.rz=numberField(grid,9,'绕 Z 增量 / °',0,[-180 180]);
                controls.rotate=button(grid,10,'应用增量旋转预览',@manualOrientation);
                controls.reset=button(grid,11,'恢复导入坐标',@resetOrientation);
                controls.confirm=button(grid,12,'确认主轴 X、鼻尖 +X',@confirmOrientation);
                controls.wiredBeadModel=button(grid,14,'选用带引线 B 型热敏珠',@(~,~)loadBeadModel([],[],'leads'));
                controls.wiredBeadModel.Layout.Column=1;
                controls.midbodyBeadModel=button(grid,14,'光滑热敏珠 · 中段加密',@(~,~)loadBeadModel([],[],'midbody'));
                controls.midbodyBeadModel.Layout.Column=2;
                tip(grid,13,['STL不含单位。默认mm模型转换系数0.001；m模型为1。' newline ...
                    'PCA只识别轴，不能确定鼻向。实际旋转会使网格/结果失效；鼠标操作只改变观察视角。']);
                controls.native=button(grid,15,'建立原生箭体 · 内嵌SU2研究',@loadNative);
                controls.nativeFeature=button(grid,16,'建立带整流罩箭体 · 内嵌SU2研究',@(~,~)loadNative([],[],true));
                controls.classic=button(grid,17,'经典锥柱扩张体 · 绝热研究',@(~,~)loadClassic());
                controls.readBenchmarkReview=button(grid,18,'公开算例验证记录',@readBenchmarkDialog);
                controls.readIdealAirReview=button(grid,19,'变比热热流接口记录',@readIdealAirDialog);
                controls.staticAir=button(grid,20,'静止空气 · 节点探针',@readStaticAirProbe);
                controls.lewisCylinder=button(grid,21,'计算 NASA 圆柱换热对照',@showLewisCylinder);
                controls.thermistorMode=choice(grid,22,'B 型热敏珠', ...
                    {'查看 3 m/s · 15° 基准网','查看 3 m/s · 7.5° 网及差值','查看 6 m/s · 全周启动瞬态', ...
                     '查看 6 m/s · 同刻时间步对照','查看 6 m/s · 近壁网格对照','计算 3 m/s · 15° 基准网','计算 3 m/s · 7.5° 细网'}, ...
                    {'archive','archive7p5','startup6','timeSteps6','nearWall6','solve','solve7p5'},'archive');
                controls.thermistor=button(grid,23,'打开热敏珠全表面结果',@showThermistor);
                controls.thermistorMode.ValueChangedFcn=@updateThermistorButton;
                set([controls.pca controls.flip controls.rotate controls.reset controls.confirm], ...
                    'Enable',onoff(~isempty(state.rawMesh)));
                if isEmbedded(),set([controls.pca controls.flip controls.rotate controls.reset controls.confirm],'Enable','off');end
            case 2
                if isBead()
                    grid=form(8);
                    if isMixedBead()
                        if hasBeadLeads()
                            controls.leadGrid=choice(grid,1,'带线试件体网',{'原混合网','中段加密','第三档中段细化','360° 全周 · 3 m/s 侧吹','96 档侧吹 · 图 7 风洞测点'},[0 1 2 3 4],double(getOption(state.cfg,'thermistorRefineLeads',false)));
                            controls.leadGrid.ValueChangedFcn=@changeLeadGrid;
                        else
                            tip(grid,1,[mixedBeadLabel() ' · 30° 对称扇区 · 周向单元 7.5°']);
                        end
                        controls.generate=button(grid,2,'生成试件混合体网',@generate);
                        if ~isempty(state.mesh)
                            p=state.mesh.thermistorPrepared;
                            layer='首层 0.937 µm';if p.angularCells==96,layer='前端面首层 0.469 µm';end
                            tip(grid,3,sprintf('%d 个体单元 · %d 个全表面原面 · %s',p.cells,size(state.mesh.F,1),layer));
                            controls.nextBoundary=button(grid,6,'设置边界条件',@(~,~)go(3));
                        end
                        sync();return;
                    end
                    controls.beadCells=choice(grid,1,'工况与体网',{'3 m/s · 15° 扇区单元','3 m/s · 7.5° 扇区单元','6 m/s · 360° 全周 · 56 层'},[2 4 48],state.cfg.thermistorAngularCells);
                    controls.beadCells.ValueChangedFcn=@changeBeadGrid;
                    controls.generate=button(grid,2,'生成热敏珠体网',@generate);
                    description='外域 20 倍横径 · 子午 64 × 径向 48 · 30° 对称扇区';
                    if isFullBead()
                        controls.generate.Text='准备 56 层全周体网';
                        description='外域 20 倍横径 · 48 个方位 · 56 个径向区间 · 全周启动瞬态';
                    end
                    if ~isempty(state.mesh)
                        p=state.mesh.thermistorPrepared;
                        description=sprintf('%s\n%d 节点 · %d 体单元 · %d 个全珠原面',description,p.nodes,p.cells,size(state.mesh.F,1));
                        controls.nextBoundary=button(grid,6,'设置边界条件',@(~,~)go(3));
                    end
                    controls.volumeSummary=tip(grid,3,description);sync();return;
                end
                if isEmbedded()
                    [~,mode]=su2WorkbenchMeshOptions(state.cfg);
                    grid=form(22);
                    controls.volumeMode=choice(grid,1,'体网格模式',{'法向偏置层 · 研究','射线层 · 历史对照','原三角面混合网格 · 研究'}, ...
                        {'normalLayers','starRay','hybridTriangles'},mode);
                    controls.volumeMode.ValueChangedFcn=@nativeMeshEdited;
                    controls.generate=button(grid,2,'生成完整三维体网格',@generate);
                    controls.reuseMesh=button(grid,3,'复用本程序已生成的混合体网格',@readHybridMeshDialog);
                    controls.reuseMesh.Enable=onoff(strcmp(mode,'hybridTriangles'));
                    description='待生成；保留完整原三角面，网格参数改变后须重新确认边界。';
                    if ~isempty(embeddedVolume)
                        description=sprintf('已生成 %d 节点 / %d 体单元\n原表面 %d 面，近壁精度尚未合格。', ...
                            size(embeddedVolume.vertices,1),numel(embeddedVolume.cells),size(state.mesh.F,1));
                    end
                    controls.volumeSummary=tip(grid,4,description);
                    if strcmp(mode,'hybridTriangles')
                        opts=getOption(state.cfg.su2,'hybridMesh',struct);
                        heights=getOption(opts,'layerHeights',[0 cumsum(2e-5*1.4.^(0:5))]);
                        if (~isempty(embeddedVolume) && embeddedVolume.meshConstruction.tetsBelowSICN01>0) || ~isempty(candidatePreview)
                            label='定位低质量体单元';if qualityPreview,label='返回整箭表面';end
                            controls.meshQuality=button(grid,14,label,@toggleMeshQuality);
                            controls.meshQualityProbe=tip(grid,15,'点选有色单元：查看体编号、SICN及坐标。灰色为原壁面定位参考。');
                        end
                        controls.inspectCandidate=button(grid,18,'检查质量未通过的网格 · 只读',@inspectMeshCandidateDialog);
                        tip(grid,19,'只显示已完成且拓扑通过的失败候选，保持当前模型、计算网格和边界状态。失败候选不能用于求解。');
                        fan=~isempty(embeddedVolume) && isfield(embeddedVolume.meshConstruction,'constructionOperation') && ...
                            startsWith(embeddedVolume.meshConstruction.constructionOperation,'closed-convex-rim-');
                        transition=~isempty(embeddedVolume) && isfield(embeddedVolume.meshConstruction,'noseRetainedLayers');
                        if ~isempty(embeddedVolume) && isfield(embeddedVolume.meshConstruction,'hasTangentialOffsets')
                            construction=embeddedVolume.meshConstruction;controls.generate.Text='重新生成常规体网格';
                            controls.meshMethod=tip(grid,7,sprintf('一般混合网格：%d棱柱 / %d四面体；保留原壁面及尾缘扇，外层含切向偏移。几何审核通过，最小SICN %.5g；流场精度未验证。', ...
                                construction.prisms,construction.tetrahedra,construction.minimumTetSICN));
                        elseif ~isempty(embeddedVolume) && strcmp(embeddedVolume.meshConstruction.method,'analytic-octahedron-prisms-and-staircase-tetrahedra')
                            construction=embeddedVolume.meshConstruction;
                            controls.meshMethod=tip(grid,7,sprintf('解析八面体接口网格 · %d棱柱 / %d四面体 · 最小SICN %.5g · 热流精度待验', ...
                                construction.prisms,construction.tetrahedra,construction.minimumTetSICN));
                        elseif transition
                            % 局部终止网格不再全表面等层；直接展示实际计数，不能沿用统一层数说明。
                            construction=embeddedVolume.meshConstruction;
                            controls.generate.Text='重新生成常规体网格';
                            controls.meshMethod=tip(grid,7,sprintf('已复用局部棱柱终止网格：%d棱柱 / %d金字塔 / %d四面体。鼻部保留%d层（%.5g mm），其余%d层。SICN<0.1：%d；未获物理资格。', ...
                                construction.prisms,construction.pyramids,construction.tetrahedra,construction.noseRetainedLayers, ...
                                1e3*construction.noseRetainedHeightM,construction.normalLayers,construction.tetsBelowSICN01));
                        elseif fan
                            controls.generate.Text='重新生成常规体网格';
                            controls.meshMethod=tip(grid,7,sprintf('已复用尾缘双法向棱柱扇；%d层，参考首层距%.5g µm。外域仍有%d个SICN<0.1单元；未获物理资格。', ...
                                numel(heights)-1,1e6*(heights(2)-heights(1)),embeddedVolume.meshConstruction.tetsBelowSICN01));
                        else
                            tip(grid,7,sprintf('保留全部原三角面：%d层法向棱柱，首层顶点位移%.5g µm；Gmsh四面体外域。热流分辨率尚未验证。',numel(heights)-1,1e6*(heights(2)-heights(1))));
                        end
                        if fan
                            label='查看尾缘体网格截面';if creasePreview,label='返回整箭表面';end
                            controls.creasePreview=button(grid,10,label,@toggleCreasePreview);
                        else,creasePreview=false;end
                    elseif strcmp(mode,'normalLayers')
                        tip(grid,7,'24层近壁偏置 + 16层外域；首层母线法向距离2 µm，近壁增长率1.25。三维面法向距离略受周向多边形影响，不能据此宣称y+或热流精度通过。');
                    else
                        tip(grid,7,'8层射线体网格，保留旧快照参数，仅作历史集成对照；近壁方向和精度尚未合格。');
                    end
                    tip(grid,11,'混合模式接受通过闭合/定向/单元审核的原三角表面和光滑整流罩；复杂悬挑/尖翼未验证。旧两种模式仅支持无突起原生箭体。');
                    sync();return;
                end
                if isCM011()
                    grid=form(8);meshBackendChoice(grid,1);
                    controls.cm011HalfSectors=choice(grid,2,'半周角向单元',{'12（15°）','24（7.5°）'}, ...
                        [12 24],getOption(state.cfg,'cm011HalfSectors',12));
                    controls.cm011HalfSectors.ValueChangedFcn=@changeCM011Resolution;
                    controls.generate=button(grid,3,'核对原面并生成半域体网',@generate);
                    description='使用已确认的米制 CM-011 表面，保留导入面元编号。';
                    if ~isempty(state.mesh) && isfield(state.mesh,'cm011Prepared')
                        p=state.mesh.cm011Prepared;
                        description=sprintf('%d 个导入面 · %d 个体单元\n最大节点匹配偏差 %.3g m', ...
                            size(state.mesh.F,1),p.cells,p.maxVertexDistanceM);
                        if isfield(p,'maxFacetSampleDistanceM')
                            description=sprintf('%d 个导入面 · %d 个体单元\n面片抽样偏差 %.3g µm · 顶点至求解壁面最大距离 %.3g µm', ...
                                size(state.mesh.F,1),p.cells,1e6*p.maxFacetSampleDistanceM,1e6*p.maxVertexToSolverWallM);
                        end
                    end
                    controls.volumeSummary=tip(grid,4,description);
                    sync();return;
                end
                grid=form(13);c=state.cfg.remesh;
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
                meshBackendChoice(grid,11);
                controls.useSurfaceForSU2=button(grid,13,'当前表面 → 内嵌SU2混合网格',@(~,~)useSurfaceForSU2());
                controls.useSurfaceForSU2.Enable=onoff(~isempty(state.mesh));
            case 3
                if isBead()
                    if isSidewind()
                        if state.cfg.thermistorAngularCells==96
                            grid=form(7);note=tip(grid,1,'侧吹空气 3 m/s · 沿 −Y · 293.15 K · 101325 Pa');note.Layout.Row=[1 2];
                            mode='uniform';if isfield(state.cfg.thermistorBoundary,'leadThermalModel'),mode='longWire';end
                            if isfield(state.cfg.thermistorBoundary,'resistanceTemperatureK'),mode='resistance';end
                            controls.wireBoundary=choice(grid,3,'试件热边界',{'整件等温','长引线导热分布','电阻测温 · 长引线'}, ...
                                {'uniform','longWire','resistance'},mode);
                            controls.wireBoundary.ValueChangedFcn=@boundaryEdited;
                            if strcmp(mode,'uniform'),tip(grid,4,'珠面、接头与引线均给定 305.15 K。');
                            elseif strcmp(mode,'longWire'),tip(grid,4,'珠面和线根给定 305.15 K。引线按论文长线近似降温，截面温度由模型给出。');
                            else,tip(grid,4,'给定电阻等效内部温度 305.15 K。按论文附录反推珠表面温度，并计算长引线温度分布。');end
                            controls.confirm=button(grid,7,'确认侧吹边界',@confirmBoundary);
                        else
                            grid=form(6);tip(grid,1,'侧吹空气 3 m/s · 沿 −Y · 293.15 K · 101325 Pa');
                            tip(grid,3,'珠体、接头、引线及端面统一给定 299.15 K。');
                            controls.confirm=button(grid,5,'确认侧吹边界',@confirmBoundary);
                        end
                        sync();return;
                    end
                    wiredPower=hasBeadLeads() && getOption(state.cfg,'thermistorRefineLeads',false);
                    if isMixedBead() && ~wiredPower
                        grid=form(6);tip(grid,1,'空气 3 m/s · 沿 −X · 293.15 K · 101325 Pa');
                        tip(grid,2,[mixedBeadLabel() ' · 全部暴露面给定 305.15 K。']);
                        controls.confirm=button(grid,5,'确认试件热边界',@confirmBoundary);sync();return;
                    end
                    if isFullBead()
                        grid=form(6);tip(grid,1,'6 m/s 空气 · 沿 −X · 293.15 K · 101325 Pa');
                        tip(grid,2,'珠面等温 303.15 K · 两个物理半步 · 距构造初场 7.809360865 µs');
                        controls.confirm=button(grid,5,'确认此工况与给定壁温',@confirmBoundary);sync();return;
                    end
                    grid=form(9);selected='none';temp=305.15;watts=0.004;
                    if wiredPower,selected='published';end
                    boundary=state.cfg.thermistorBoundary;
                    if isstruct(boundary)
                        if isfield(boundary,'netHeatingPowerW'),selected='power';watts=boundary.netHeatingPowerW;
                        else,selected='published';temp=boundary.wallTemperatureK;end
                    end
                    if wiredPower
                        items={'给定均匀壁温','整件净对流功率'};values={'published','power'};powerLabel='整件净对流功率 / W';
                    else
                        items={'请选择热边界','给定均匀壁温','给定净加热功率'};values={'none','published','power'};powerLabel='整珠净功率 / W';
                    end
                    controls.beadCondition=choice(grid,1,'热边界',items,values,selected);
                    controls.beadCondition.ValueChangedFcn=@boundaryEdited;
                    controls.beadWall=numberField(grid,2,'给定壁温 / K',temp,[293.15 305.15]);
                    controls.beadWall.ValueChangedFcn=@boundaryEdited;controls.beadWall.Enable=onoff(strcmp(selected,'published'));
                    limit=thermistorNetPowerLimit(state.cfg.thermistorAngularCells,state.cfg);
                    controls.beadPower=numberField(grid,3,powerLabel,watts,[0 limit]);
                    controls.beadPower.ValueDisplayFormat='%.10g';
                    controls.beadPower.ValueChangedFcn=@boundaryEdited;controls.beadPower.Enable=onoff(strcmp(selected,'power'));
                    if wiredPower
                        tip(grid,4,sprintf('给定壁温 293.15–305.15 K。净功率覆盖珠面、接头和引线，上限 %.8g mW。',1000*limit));
                    else
                        tip(grid,4,sprintf('3 m/s 空气，293.15 K，101325 Pa。净功率已扣除引线、辐射等损失。本网上限约 %.6g mW。',1000*limit));
                    end
                    controls.confirm=button(grid,7,'确认热敏珠边界',@confirmBoundary);
                    controls.confirm.Enable=onoff(wiredPower || ~isempty(state.cfg.thermistorBoundary));sync();return;
                end
                if isCM011()
                    grid=form(8);selected='none';
                    if ~isempty(state.cfg.cm011Boundary),selected='inclined';end
                    controls.cm011Condition=choice(grid,1,'试件工况',{'请选择工况','CM-011 · +10° 等温'}, ...
                        {'none','inclined'},selected);
                    controls.cm011Condition.ValueChangedFcn=@boundaryEdited;
                    tip(grid,2,'空气 300 K · 101325 Pa · 0.1 m/s。攻角 +10°，侧滑 0°。无滑移等温壁 300.1 K。');
                    controls.confirm=button(grid,5,'确认试件边界条件',@confirmBoundary);
                    controls.confirm.Enable=onoff(~isempty(state.cfg.cm011Boundary));sync();return;
                end
                if isEmbedded()
                    grid=form(61);c=state.cfg;
                    controls.altitude=numberField(grid,2,'高度 / km',c.flight.altitude/1000,[0 86]);
                    controls.mach=numberField(grid,3,'参考Mach',c.flight.Mach,[1.01 10]);
                    controls.alpha=numberField(grid,4,'攻角 / °',c.flight.alpha,[-89 89]);
                    controls.beta=numberField(grid,5,'侧滑角 / °',c.flight.beta,[-89 89]);
                    controls.wall=numberField(grid,6,'给定壁温 / K',c.flight.Twall,[1 10000]);
                    controls.iterations=numberField(grid,7,'最大迭代数',getOption(c.su2,'iterations',1000),[1 30000]);
                    controls.seconds=numberField(grid,8,'时间预算 / s',getOption(c.su2,'maximumSeconds',120),[1 720]);
                    mode='overwrite';if getOption(c.su2,'namedCheckpoints',false),mode='indexed';end
                    controls.checkpointMode=choice(grid,9,'检查点保存',{'覆盖写 · 默认','按迭代保留'}, {'overwrite','indexed'},mode);
                    controls.outputEvery=numberField(grid,10,'保存间隔 / 步',getOption(c.su2,'outputEvery',1000),[1 30000]);
                    thermal=getOption(c.su2,'wallThermalMode','isothermal');
                    controls.thermalMode=choice(grid,11,'壁面热边界',{'给定统一壁温','绝热 · 给定热流0','沿X变化的壁温'}, ...
                        {'isothermal','adiabatic','axialTemperature'},thermal);
                    controls.wall.Enable=onoff(strcmp(thermal,'isothermal'));
                    for key={'altitude','mach','alpha','beta','wall','iterations','seconds','checkpointMode','outputEvery','thermalMode'},controls.(key{1}).ValueChangedFcn=@boundaryEdited;end
                    controls.wallProfile=button(grid,12,'导入X—壁温分布表',@importWallTemperatureProfile);
                    controls.thermodynamics=choice(grid,1,'热力学模型',{'定比热 · 默认','变比热理想空气 · SST研究','HEOS Air · 研究'}, ...
                        {'idealGas','idealAir','HEOS-Air'},getOption(c.su2,'thermodynamicModel','idealGas'));
                    controls.thermodynamics.ValueChangedFcn=@thermodynamicsEdited;
                    controls.confirm=button(grid,13,'确认内嵌引擎边界',@confirmBoundary);
                    controls.propertySummary=tip(grid,14,propertySummaryText(c));
                    controls.storageEstimate=tip(grid,18,checkpointBudgetText(c));
                    text='CSV两列：X_m、T_K；覆盖整个模型，温度为给定条件，节点间线性变化。';
                    if isfield(c.su2,'wallTemperatureProfile')
                        p=c.su2.wallTemperatureProfile;
                        text=sprintf('已导入%d站：X %.5g—%.5g m，给定壁温 %.5g—%.5g K。',numel(p.x),min(p.x),max(p.x),min(p.temperature),max(p.temperature));
                    end
                    controls.thermalSummary=tip(grid,21,text);
                    tip(grid,24,'绝热q=0属于给定边界；流体节点温度不等于材料壁温解。输入变化后须重新确认与计算。');
                    if isfield(c.flight,'atmosphere') && ~isfield(c.flight.atmosphere,'geometricAltitude')
                        controls.altitude.Enable='off';
                        tip(grid,45,sprintf('指定静态来流：T %.6g K，p %.6g Pa；高度栏不参与此工况。',c.flight.atmosphere.temperature,c.flight.atmosphere.pressure));
                    end
                    % 初猜不是流场解；显式研究选项只对已核对的混合域和等温壁开放。
                    controls.initialMode=choice(grid,28,'初始场',{'自由流 · 默认','壁面相容总焓 · 研究','同工况检查点 · 继续研究','定比热源场转HEOS','层流或已审核SST源场','实际壁距独立初场 · SST研究'}, ...
                        {'freestream','wallEnthalpy','checkpoint','propertyTransfer','sstTransfer','sstWallDistance'},getOption(c.su2,'initialMode','freestream'));
                    controls.initialThickness=numberField(grid,29,'初猜厚度 / mm',1000*getOption(c.su2,'initialThicknessM',.001),[.001 1000]);
                    controls.checkpointLimit=numberField(grid,30,'检查点上限 / MiB',getOption(c.su2,'checkpointLimitMiB',100),[1 1536]);
                    for key={'initialMode','initialThickness','checkpointLimit'},controls.(key{1}).ValueChangedFcn=@boundaryEdited;end
                    controls.checkpointSource=textField(grid,31,'源检查点目录',getOption(c.su2,'checkpointSource',''));
                    controls.checkpointSource.ValueChangedFcn=@boundaryEdited;
                    controls.chooseCheckpoint=button(grid,32,'选择已完成的同工况任务',@chooseCheckpointInitial);
                    tip(grid,33,'继续研究须保持原体网格、物理及数值配置并启用具名保存。初猜厚度须在完整棱柱层内，不代表已解边界层。保存上限只限制检查点估计。');
                    controls.heosSource=textField(grid,36,'定比热源目录',getOption(c.su2,'heosSource',''));
                    controls.heosSource.ValueChangedFcn=@boundaryEdited;
                    controls.chooseHEOS=button(grid,37,'选择HEOS转换的定比热源场',@chooseHEOSSource);
                    tip(grid,38,'HEOS固定原p/T/速度，Mach栏为定比热速度参考；新EOS Mach另列。统一Tw、具名保存；源步号不累加。');
                    tip(grid,41,'HEOS仅层流；初态转换≤180s，状态审核≤120s，与流场预算分别计时。没有材料壁温解；低温输运/网格/收敛未合格。SST设置见下方。');
                    controls.flowModel=choice(grid,47,'黏性流动模型',{'层流 · 默认','SST-2003m · 研究'}, ...
                        {'laminar','SST-2003m'},getOption(c.su2,'flowModel','laminar'));
                    controls.flowModel.ValueChangedFcn=@flowModelEdited;
                    controls.turbulenceIntensity=numberField(grid,48,'入口湍流强度 I',getOption(c.su2,'turbulenceIntensity',.01),[1e-6 .5]);
                    controls.turbulentViscosityRatio=numberField(grid,49,'入口 μt/μ',getOption(c.su2,'turbulentViscosityRatio',10),[1e-6 1e4]);
                    controls.turbulentPrandtl=numberField(grid,50,'湍流 Prt',getOption(c.su2,'turbulentPrandtl',.9),[.01 10]);
                    controls.turbulenceProvenance=textField(grid,51,'湍流参数来源',getOption(c.su2,'turbulenceProvenance','全湍流受控研究假设；非飞行测量，转捩未知'));
                    for key={'turbulenceIntensity','turbulentViscosityRatio','turbulentPrandtl','turbulenceProvenance'},controls.(key{1}).ValueChangedFcn=@boundaryEdited;end
                    tip(grid,53,'SST采用统一给定Tw。变比热空气使用实际壁距初场及独立阶段保存；定比热还支持已有源场。I为比例，p/T热流阶段与更新k/ω分别审核。');
                    tip(grid,57,'全湍流是计算假设；没有转捩或壁温解。SST和HEOS不可组合；研究值未通过网格/外域/公开热流验证，不计新增合格面积。');
                    controls.boundaryLayout=struct('grid',grid,'items',grid.Children, ...
                        'rows',arrayfun(@(x)x.Layout.Row(1),grid.Children));
                    layoutEmbeddedBoundary();
                    sync();return;
                end
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
                if isBead()
                    grid=form(11);p=state.mesh.thermistorPrepared;
                    if isMixedBead()
                        controls.readBead=button(grid,1,'读取同网格已算结果',@(~,~)solveCase([],[],'archive'));
                        controls.run=button(grid,2,'从冷态计算当前试件',@solveCase);
                        if isSidewind()
                            if isfield(state.cfg.thermistorBoundary,'resistanceTemperatureK')
                                controls.run.Text='按电阻测温计算全表面';
                                controls.readBead.Text='读取修正壁温的原生核对';
                                tip(grid,3,'两种来源共用同一电阻测点。结果页分别按内部超温和表面超温列出珠部平均换热系数。');
                            elseif isfield(state.cfg.thermistorBoundary,'leadThermalModel')
                                controls.readBead.Text='读取长引线已算全场';controls.run.Text='计算长引线全表面响应';
                                controls.wireSeed=textField(grid,3,'均温初猜目录（可空）','');
                                controls.beadDirectory=textField(grid,4,'引线任务目录',lastWireDirectory);
                                controls.browseRun=button(grid,5,'选择引线任务目录',@chooseBeadDirectory);
                                controls.resumeWire=button(grid,6,'接续引线任务',@(~,~)solveCase([],[],'resumeWire'));
                                controls.pauseMeasured=button(grid,7,'保存当前阶段后暂停',@pauseMeasured);controls.pauseMeasured.Enable='off';
                                note=tip(grid,8,'新算可复用本程序的 96 档均温冷启动。初猜目录留空时由程序生成。');note.Layout.Row=[8 9];
                                controls.log=uitextarea(grid,'Editable','off','Value',logLines);
                                controls.log.Layout.Row=[10 11];controls.log.Layout.Column=[1 2];
                            elseif p.angularCells==96
                                controls.run.Text='用 6 K 参照场换算 12 K 测点';controls.readBead.Text='读取已算 12 K 原生结果';
                                controls.coldMeasured=button(grid,3,'从冷态新算 12 K 测点',@(~,~)solveCase([],[],'cold'));
                                controls.beadDirectory=textField(grid,4,'冷启动任务目录',measuredDirectory);
                                controls.browseRun=button(grid,5,'选择任务目录',@chooseBeadDirectory);
                                controls.resumeBead=button(grid,6,'接续或查看该次冷启动',@(~,~)solveCase([],[],'resume'));
                                controls.pauseMeasured=button(grid,7,'保存当前阶段后暂停',@pauseMeasured);
                                controls.pauseMeasured.Enable='off';
                                note=tip(grid,8,'冷启动从均匀来流计算。暂停在当前阶段保存并核验后生效，目录可用于接续。');
                                note.Layout.Row=[8 9];
                                controls.log=uitextarea(grid,'Editable','off','Value',logLines);
                                controls.log.Layout.Row=[10 11];controls.log.Layout.Column=[1 2];
                            else,tip(grid,3,'完整 360° · 沿 −Y 稳态侧吹 · 给定壁温 299.15 K');end
                            sync();return;
                        end
                        if isfield(state.cfg.thermistorBoundary,'netHeatingPowerW')
                            controls.readBead.Text='用已算参照场求壁温';controls.run.Text='冷态重算参照场并求壁温';
                        elseif state.cfg.thermistorBoundary.wallTemperatureK~=305.15
                            controls.readBead.Text='按给定壁温换算参照场';controls.run.Text='冷态重算参照场并换算';
                        end
                        tip(grid,3,[mixedBeadLabel() ' · 3 m/s 轴向稳态 · ' beadBoundaryText(state.cfg.thermistorBoundary)]);
                        sync();return;
                    end
                    if isFullBead()
                        controls.run=button(grid,1,'从构造初场计算两个半步',@solveCase);
                        tip(grid,2,'6 m/s · 给定壁温 303.15 K · 目标 7.809360865 µs。已有初场约 49 分钟，首次计算另需准备初场。');
                        controls.beadDirectory=textField(grid,4,'已有计算目录',fullfile(state.cfg.projectRoot,'outputs','v294_validation'));
                        controls.browseRun=button(grid,5,'选择计算目录',@chooseBeadDirectory);
                        controls.readBead=button(grid,6,'读取已完成计算',@(~,~)solveCase([],[],'read'));
                        controls.resumeBead=button(grid,7,'从已审核阶段接续',@(~,~)solveCase([],[],'resume'));
                        controls.log=uitextarea(grid,'Editable','off','Value',logLines);
                        controls.log.Layout.Row=[8 11];controls.log.Layout.Column=[1 2];sync();return;
                    end
                    controls.readBead=button(grid,1,'用参照场计算当前工况',@(~,~)solveCase([],[],'archive'));
                    controls.run=button(grid,2,'从冷初场计算当前试件',@solveCase);
                    tip(grid,3,sprintf('B 型热敏珠 · 周向单元 %g° · %d 体单元\n3 m/s 空气 · %s',30/p.angularCells,p.cells,beadBoundaryText(state.cfg.thermistorBoundary)));
                    controls.log=uitextarea(grid,'Editable','off','Value',logLines);
                    controls.log.Layout.Row=[6 11];controls.log.Layout.Column=[1 2];sync();return;
                end
                if isCM011()
                    grid=form(12);
                    controls.run=button(grid,1,'计算 CM-011 全表面压力与热流',@solveCase);
                    tip(grid,2,'+10° 低速等温试件 · 内部三维不可压层流与能量方程，Y≥0 半域。');
                    controls.log=uitextarea(grid,'Editable','off','Value',logLines);
                    controls.log.Layout.Row=[5 12];controls.log.Layout.Column=[1 2];sync();return;
                end
                if isEmbedded()
                    grid=form(20);controls.run=button(grid,1,'启动内嵌SU2 · 研究计算',@solveCase);
                    controls.cancel=button(grid,2,'取消当前引擎任务',@cancelEngine);
                    controls.cancel.Enable=onoff(state.busy);
                    controls.readResearch=button(grid,3,'重开本程序同工况研究结果',@readResearchDialog);
                    controls.readTimeoutCheckpoint=button(grid,4,'查看同工况超时检查点 · 研究',@readTimeoutDialog);
                    controls.readVerifiedHEOS=button(grid,5,'查看同表面/工况的已有HEOS场',@readVerifiedHEOSDialog);
                    controls.readSSTStudy=button(grid,6,'查看同表面/工况的SST数值对照',@readSSTStudyDialog);
                    tip(grid,7,'已有场须核对原文件、物理条件、EOS与面元；沿用实际网格/迭代和数值方案。此操作只查看，不继续计算。');
                    controls.readBenchmarkReview=button(grid,10,'公开算例验证记录',@readBenchmarkDialog);
                    controls.readIdealAirReview=button(grid,11,'变比热热流接口记录',@readIdealAirDialog);
                    controls.log=uitextarea(grid,'Editable','off','Value',logLines);controls.log.Layout.Row=[12 20];controls.log.Layout.Column=[1 2];
                    sync();return;
                end
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
                if ~isempty(state.researchResult)
                    grid=form(12);r=state.researchResult;
                    tip(grid,1,sprintf('输入版本%d · 内嵌NS研究快照\n%d个原面元 · 新增合格面积0',state.resultSourceRevision,size(r.F,1)));
                    tip(grid,4,'压力、物面平均q、流体节点T可拾取；保守热负荷密度另列。给定Tw只作输入；真实壁温未求。');
                    tip(grid,8,'点击三维表面或拖动X滑块查看该截面周向分布；云图与曲线共用色标，可导出带来源CSV。研究值仍须收敛、网格与物理审核。');
                    controls.staticAir=button(grid,12,'静止空气 · 节点探针',@readStaticAirProbe);
                    sync();return;
                end
                resultRows=11+double(isfield(state.result,'thermistorOuterDomain') || isfield(state.result,'thermistorLayerBudget'));
                if isfield(state.result,'thermistorWire'),resultRows=18;end
                grid=form(resultRows+double(isfield(state.result,'thermistorAngularSpacing')));d=buildAerothermalViewerData(state.result);
                h=d.fields(strcmp({d.fields.key},'heatFlux'));p=d.fields(strcmp({d.fields.key},'pressure'));
                tip(grid,1,sprintf('当前版本 %d 的计算快照\n热流计算面积 %.4g%%\n压力计算面积 %.4g%% / 给定 %.4g%%', ...
                    state.resultSourceRevision,100*h.computedAreaFraction,100*p.computedAreaFraction,100*p.prescribedAreaFraction));
                if isFullBead()
                    origins=struct('run','本次计算','read','已保存计算','resume','接续计算');
                    tip(grid,4,sprintf('%s · 全周 56 层 · 7.809360865 µs · 珠面放热 %.6f mW', ...
                        origins.(report.startupAction),-1000*report.wholeSurfaceHeatPowerW));
                elseif isfield(state.result,'thermistorMidbody')
                    s=state.result.thermistorMidbody;
                    tip(grid,4,sprintf('光滑珠中段加密 · 整珠放热 %.5f mW · 平均换热 %.2f W/(m²·K)',1000*s.totalOutwardW,s.beadMeanHWm2K));
                elseif isSidewind()
                    s=state.result.thermistorSidewind;
                    origin='侧吹 −Y';
                    if isfield(state.result,'thermistorWire'),origin='长引线温度边界';end
                    if isfield(state.result,'thermistorWire') && isfield(state.result.thermistorWire,'thermometry')
                        m=state.result.thermistorWire.thermometry;origin=sprintf('电阻 %.2f K · 珠表面 %.3f K',m.internalK,m.surfaceK);
                    end
                    if isfield(state.result,'thermistorMeasuredPoint')
                        m=state.result.thermistorMeasuredPoint;
                        if isfield(m,'coldSource'),origin=sprintf('独立冷启动 · %d 步',m.coldSource.iterations);
                        elseif strcmp(m.action,'solve'),origin='6 K 参照场换算';else,origin='12 K 原生归档';end
                    end
                    detail='沿 X 选择截面查看周向差异。';
                    if isfield(state.result,'thermistorMeasuredPoint') && isfield(state.result.thermistorMeasuredPoint,'coldComparisonNote')
                        detail=state.result.thermistorMeasuredPoint.coldComparisonNote;
                    end
                    tip(grid,4,sprintf('%s · 七区放热 %.5f mW · 最强放热面 %d。%s',origin,1000*s.totalOutwardW,s.peakFaceID,detail));
                elseif hasBeadLeads()
                    s=state.result.thermistorLead;
                    msg=sprintf('珠体和接头 %.4f mW · 两条引线含端面 %.4f mW · 总放热 %.4f mW', ...
                        1000*s.beadAndJointsW,1000*(s.wireSidesW+s.wireEndsW),1000*s.totalOutwardW);
                    if isfield(state.result,'hybridSolvedWallTemperatureK')
                        msg=sprintf('求解壁温 %.6f K · %s',s.thermal.wallTemperatureK,msg);
                    else
                        msg=sprintf('给定壁温 %.9g K · %s',state.result.hybridCondition.wallTemperatureK,msg);
                    end
                    tip(grid,4,msg);
                elseif isBead()
                    origin='本次冷启动计算';
                    if state.result.evidenceIdentity.reusedArchive,origin='使用归档温差场';end
                    tip(grid,4,sprintf('%s · 周向单元 %g°。点选面元可查壁压、入壁热流、首层流体温度和各量的方法。', ...
                        origin,state.result.thermistorBead.cellAngleDeg));
                elseif isCM011()
                    tip(grid,4,'点选导入三角面查看压力、热流和流体温度。相对环境壁压适合比较同一圆周两侧。');
                else
                    tip(grid,4,'点选面元核对场量、算法和范围。灰色未解；未知外缘热力状态与有效状态的Kn超限分开，两者均可阻断下游热流。输入壁温与求解壁温分开。');
                end
                controls.export=button(grid,7,'导出完整离线工作台',@exportCase);
                if isfield(state.result,'thermistorPaired')
                    controls.export.Layout.Column=2;
                    controls.pairedHeat=button(grid,7,'查看配对热流差',@showPairedHeat);
                    if isfield(state.result,'thermistorMidbody'),controls.pairedHeat.Text='查看加密前后热流差';end
                    controls.pairedHeat.Layout.Column=1;
                end
                % 在同一结果侧栏并列两种珠部平均量，热流云图和截面继续显示逐面解。
                if isfield(state.result,'thermistorWire')
                    w=state.result.thermistorWire;
                    controls.wireComparison=uitable(grid,'RowName',{},'ColumnName',{'珠部对照','h / W/(m²·K)','放热 / mW'}, ...
                        'Data',{'图 7 测点',sprintf('%.2f',w.experiment.hWm2K),sprintf('%.4f',1000*w.experimentEquivalentW); ...
                                '整件等温',sprintf('%.2f',w.controlHWm2K),sprintf('%.4f',1000*w.controlBeadW); ...
                                '长引线',sprintf('%.2f',w.beadHWm2K),sprintf('%.4f',1000*w.beadAndJointsW)});
                    controls.wireComparison.Layout.Row=[8 11];controls.wireComparison.Layout.Column=[1 2];
                    labels={'整段','沿线分段'};ids={'whole','bands'};selected='whole';
                    if isfield(w,'thermometry')
                        m=w.thermometry;
                        controls.wireComparison.Data={'图 7 · 电阻温差',sprintf('%.2f',w.experiment.hWm2K),sprintf('%.4f',1000*w.experimentEquivalentW); ...
                            '原长线 · 表温差',sprintf('%.2f',m.referenceBeadW/(m.areaM2*12)),sprintf('%.4f',1000*m.referenceBeadW); ...
                            '修正 · 电阻温差',sprintf('%.2f',m.apparentHWm2K),sprintf('%.4f',1000*w.beadAndJointsW); ...
                            '修正 · 表温差',sprintf('%.2f',m.surfaceHWm2K),sprintf('%.4f',1000*w.beadAndJointsW)};
                        labels=[{'测温口径'} labels];ids=[{'thermometry'} ids];selected='thermometry';
                    end
                    controls.wireBudgetMode=choice(grid,12,'对照量',labels,ids,selected);
                    controls.wireBudgetMode.ValueChangedFcn=@updateWireBudget;
                    controls.wireBudget=uitable(grid,'RowName',{});controls.wireBudget.Layout.Row=[13 18];controls.wireBudget.Layout.Column=[1 2];
                    updateWireBudget([],[]);sync();return;
                elseif isfield(state.result,'thermistorMeasuredPoint')
                    m=state.result.thermistorMeasuredPoint;
                    controls.measuredPoint=uitable(grid,'RowName',{},'ColumnName',{'珠部平均量','h / W/(m²·K)','比图读 / %'}, ...
                        'Data',{'图 7 测点',sprintf('%.2f',m.experiment.hWm2K),'—';'侧吹 96 档',sprintf('%.2f',m.beadHWm2K),sprintf('%+.2f',100*m.relativeDifference); ...
                        '名义轴向',sprintf('%.2f',m.axial.hWm2K),sprintf('%+.2f',100*m.axial.relativeDifference); ...
                        '原文式（28）',sprintf('%.2f',m.correlation.meanHWm2K),sprintf('%+.2f',100*m.correlation.relativeDifference)});
                    controls.measuredPoint.Layout.Row=[8 11];controls.measuredPoint.Layout.Column=[1 2];sync();return;
                elseif isfield(state.result,'thermistorLead') && isfield(state.result.thermistorLead,'meanConvection') ...
                        && state.result.thermistorLead.meanConvection.available
                    m=state.result.thermistorLead.meanConvection;
                    controls.beadMeanMethods=uitable(grid,'RowName',{}, ...
                        'ColumnName',{'珠部方法','h / W/(m²·K)','放热 / mW'}, ...
                        'Data',{'坪内–佐藤（28）',sprintf('%.2f',m.meanHWm2K),sprintf('%.4f',1000*m.outwardW); ...
                                '三维层流积分',sprintf('%.2f',m.surfaceMeanHWm2K),sprintf('%.4f',1000*m.surfaceOutwardW)});
                    controls.beadMeanMethods.Layout.Row=[8 10];controls.beadMeanMethods.Layout.Column=[1 2];
                else
                    tip(grid,8,'任意X截面曲线和导出使用同一快照。完成计算不等于全面积已解；修改上游使结果失效，浏览与导出不重算。');
                end
                controls.staticAir=button(grid,11,'静止空气 · 节点探针',@readStaticAirProbe);
                if isfield(state.result,'thermistorLeadGrid')
                    if isfield(controls,'pairedHeat'),controls.pairedHeat.Text='查看引线影响';end
                    controls.staticAir.Layout.Column=2;
                    controls.gridHeat=button(grid,11,'查看网距影响',@(~,~)showPairedHeat([],[],true));
                    controls.gridHeat.Layout.Column=1;
                end
                if isfield(state.result,'thermistorOuterDomain') && state.result.thermistorOuterDomain.available
                    controls.outerHeat=button(grid,12,'查看外边界影响 · 30D 减 20D',@showOuterHeat);
                end
                if isfield(state.result,'thermistorTipSpacing')
                    controls.staticAir.Layout.Column=2;
                    controls.tipSpacing=button(grid,11,'前引线端面网距',@showTipSpacing);
                    controls.tipSpacing.Layout.Column=1;
                end
                if isfield(state.result,'thermistorAngularSpacing')
                    controls.angularSpacing=button(grid,13,'周向网距 · 48 与 96 档',@showAngularSpacing);
                end
                if isfield(state.result,'thermistorThreeGrids') && state.result.thermistorThreeGrids.summary.available
                    controls.staticAir.Layout.Column=2;
                    controls.compareGrids=button(grid,11,'对照三档网距',@showThirdGrid);
                    controls.compareGrids.Layout.Column=1;
                end
                if isfield(state.result,'thermistorLayerBudget')
                    controls.layerBudget=button(grid,12,'查看分带热量账',@showLayerBudget);
                end
        end
        sync();
    end
    function showPairedHeat(~,~,gridEffect)
        % 打开当前冻结结果的配对证据页，初次定位绝对热流差最大的可比面。
        if nargin==3 && gridEffect
            delta=state.result.thermistorLeadGrid.deltaHeatWm2;key='leadGridDeltaHeat';order='带线加密减原网';
        else
            pair=state.result.thermistorPaired;delta=pair.deltaHeatWm2;key='pairedDeltaHeat';order=pair.order;
        end
        inspectHeatDifference(delta,key,order);
    end
    function showOuterHeat(~,~)
        % 按当前热边界定位两次外域计算的最大绝对热流差。
        pair=state.result.thermistorOuterDomain;
        inspectHeatDifference(pair.deltaHeatWm2,'outerDeltaHeat','30D 减 20D');
    end
    function showTipSpacing(~,~)
        % 直接定位细网放热峰；云图显示该面相对原网的热流变化。
        t=state.result.thermistorTipSpacing;g=t.summary;id=g.finePeakFace;
        tab=webUI.evidenceTab;tab.Parent.SelectedTab=tab;
        webUI.evidenceView.selectField('tipDeltaHeat');webUI.evidenceView.inspectFace(id,[]);
        message(sprintf('细网面元 %d 入壁热流 %.3f W/m²，与第二强面 %d 相差 %.3g W/m²。前端面放热 %.5f µW。', ...
            id,t.fineHeat(id),g.fineSecondFace,g.finePeakGapWm2,1e6*g.regions(7).fineOutW));
    end
    function showAngularSpacing(~,~)
        % 用同一原面显示两档热流差，定位 96 档积分回原面后的最强放热区。
        a=state.result.thermistorAngularSpacing;g=a.summary;id=g.finePeakFace;
        tab=webUI.evidenceTab;tab.Parent.SelectedTab=tab;
        webUI.evidenceView.selectField('angularDeltaHeat');webUI.evidenceView.inspectFace(id,[]);
        message(sprintf('96 档最强面 %d，入壁热流 %.2f W/m²，比 48 档变化 %+.2f W/m²。整件放热 %.6f mW。', ...
            id,a.heat96(id),a.deltaHeat(id),1e3*g.outward96W));
    end
    function showThirdGrid(~,~)
        % 定位第三档相对第二档的最大热流变化，可在变量菜单切换到上一次细化。
        pair=state.result.thermistorThreeGrids;
        inspectHeatDifference(pair.heat(:,3)-pair.heat(:,2),'gridHeat32','第三档减第二档');
    end
    function showLayerBudget(~,~)
        % 打开分带结果并定位原面放热变化最大的带；每条带的实际功率可在探针中查看。
        b=state.result.thermistorLayerBudget;d=[b.bands.delta];[~,k]=max(abs([d.faceOutW]));
        face=find(b.faceBand==k,1);tab=webUI.evidenceTab;tab.Parent.SelectedTab=tab;
        webUI.evidenceView.selectField('bandRelease');webUI.evidenceView.inspectFace(face,[]);
        message(sprintf('%s，第三档相对第二档的全带放热变化 %+.3f µW。',b.bands(k).name,1e6*d(k).faceOutW));
    end
    function inspectHeatDifference(delta,key,order)
        ids=find(isfinite(delta));[~,k]=max(abs(delta(ids)));
        tab=webUI.evidenceTab;tab.Parent.SelectedTab=tab;
        webUI.evidenceView.selectField(key);webUI.evidenceView.inspectFace(ids(k),[]);
        message(sprintf('%s，最大绝对差位于面元 %d，入壁热流差 %+.3f W/m²。',order,ids(k),delta(ids(k))));
    end
    function grid=form(rows)
        grid=uigridlayout(editor,[rows 2]);grid.ColumnWidth={138,'1x'};
        grid.RowHeight=repmat({31},1,rows);grid.Padding=[0 2 0 2];grid.RowSpacing=5;
        if isprop(grid,'Scrollable'),grid.Scrollable='on';end
    end
    function updateWireBudget(~,~)
        % 模块用途：分别展示解析轴向导热和三维空气侧对流，右、左线各自按原面求和。
        w=state.result.thermistorWire;t=controls.wireBudget;
        if strcmp(controls.wireBudgetMode.Value,'thermometry')
            m=w.thermometry;t.ColumnName={'测温模型','数值'};
            t.Data={'内部等效温度 / K',sprintf('%.4f',m.internalK);'推算珠表面 / K',sprintf('%.4f',m.surfaceK); ...
                '长度雷诺数',sprintf('%.3f',m.reynoldsL);'内外温差比 / %',sprintf('%.3f',100*m.ratio); ...
                '图 12 读数 / %',sprintf('%.3f–%.3f',100*m.figure12.range);'剩余放热偏差 / mW',sprintf('%+.5f',1000*m.remainingGapW); ...
                '所需内外温差比 / %',sprintf('%.2f',100*m.originalBias);'论文 3% 情景剩余偏差',sprintf('%+.2f%%',100*m.upperScenarioBias)};
        elseif strcmp(controls.wireBudgetMode.Value,'whole')
            t.ColumnName={'热量 / mW','后引线','前引线'};
            p=w.profile;values=1000*[p.oneLeadRootW p.oneLeadRootW;p.oneLeadThroughCutW p.oneLeadThroughCutW; ...
                p.oneLeadSideW p.oneLeadSideW;w.regions(2).outwardW w.regions(6).outwardW;w.regions(1).outwardW w.regions(7).outwardW];
            labels={'根部轴向导热';'截点轴向导热';'一维段内放热';'三维线侧放热';'三维端面放热'};
            t.Data=[labels cellstr(compose('%.5f',values(:,1))) cellstr(compose('%.5f',values(:,2)))];
        else
            t.ColumnName={'线段 / µm','一维 / µW','后侧 / µW','前侧 / µW'};s=w.bands;
            t.Data=[cellstr(compose('%.0f–%.0f',1e6*[s.fromM]',1e6*[s.toM]')) cellstr(compose('%.3f',1e6*[s.theorySideW]')) ...
                cellstr(compose('%.3f',1e6*[s.backSideW]')) cellstr(compose('%.3f',1e6*[s.frontSideW]'))];
        end
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
    function x=tip(grid,row,text)
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
            next=state.cfg;next.computationBackend='engineering';embeddedVolume=[];
            next.geometry.stlFile=controls.file.Value;next.geometry.scale=controls.scale.Value;
            state=aeroHeatWorkflowState('invalidateModel',state,struct('cfg',next));
            preview=[];report=[];preparing=true;message('正在加载模型并核对SI尺度…');
            hash=cm011FileSHA256(next.geometry.stlFile);
            raw=quietCall(@()loadSTLMesh(next));raw=quietCall(@()calcFaceGeometry(raw));
            assert(strcmp(hash,cm011FileSHA256(next.geometry.stlFile)),'RAH:STL:Changed','STL 读取期间发生改变，请重新加载。');
            [ok,stlFile]=fileattrib(next.geometry.stlFile);assert(ok,'无法定位导入文件。');
            state=aeroHeatWorkflowState('setModel',state,struct('rawMesh',raw,'cfg',next, ...
                'modelSource',struct('file',stlFile.Name,'scale',next.geometry.scale,'unit','m','sha256',hash)));
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
        if isBead()
            try
                preparing=true;message('正在准备并核对所选热敏珠体网…');
                if isMixedBead(),candidate=prepareThermistorLeadMesh(state.orientedMesh,@()isgraphics(fig),state.cfg.thermistorVariant,getOption(state.cfg,'thermistorRefineLeads',false));
                else,candidate=prepareThermistorWorkbenchMesh(state.orientedMesh,state.cfg.thermistorAngularCells,@()isgraphics(fig));end
                state=aeroHeatWorkflowState('setMesh',state,struct('mesh',candidate));report=[];
                preparing=false;renderEditor();showGeometry(true);message('体网已就绪，可以设置边界条件。');
            catch ex,preparing=false;if isgraphics(fig),renderEditor();showError(ex);end,end
            return;
        end
        if isEmbedded(),generateNative();return;end
        if isCM011()
            try
                preparing=true;message('正在核对导入原面并生成 CM-011 半域体网…');drawnow;
                candidate=prepareCM011WorkbenchMesh(state.orientedMesh,state.modelSource,getOption(state.cfg,'cm011HalfSectors',12));
                state=aeroHeatWorkflowState('setMesh',state,struct('mesh',candidate));
                preparing=false;renderEditor();message('体网已生成，请选择并确认低速试件工况。');showGeometry(true);
            catch ex,preparing=false;renderEditor();showError(ex);end
            return;
        end
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
        next=state.cfg;
        if isBead()
            next.thermistorBoundary=[];
            if isSidewind()
                b=thermistorWorkbenchBoundary();b.wallTemperatureK=299.15;b.betaDeg=90;
                if next.thermistorAngularCells==96,b.wallTemperatureK=305.15;end
                if isfield(controls,'wireBoundary') && ~strcmp(controls.wireBoundary.Value,'uniform')
                    b.leadThermalModel='longWire';
                    if strcmp(controls.wireBoundary.Value,'resistance')
                        b=rmfield(b,'wallTemperatureK');b.resistanceTemperatureK=305.15;
                    end
                end
                next.thermistorBoundary=thermistorWorkbenchBoundary(b,next.thermistorAngularCells,next);return;
            end
            if isMixedBead() && ~isfield(controls,'beadCondition'),next.thermistorBoundary=thermistorWorkbenchBoundary();return;end
            if isFullBead()
                b=thermistorWorkbenchBoundary();b.speed=6;b.wallTemperatureK=303.15;
                next.thermistorBoundary=thermistorWorkbenchBoundary(b,48);return;
            end
            if ~strcmp(controls.beadCondition.Value,'none')
                b=thermistorWorkbenchBoundary();
                if strcmp(controls.beadCondition.Value,'power')
                    b=rmfield(b,'wallTemperatureK');b.netHeatingPowerW=controls.beadPower.Value;
                else
                    b.wallTemperatureK=controls.beadWall.Value;
                end
                next.thermistorBoundary=thermistorWorkbenchBoundary(b,next.thermistorAngularCells,next);
            end
            return;
        end
        if isCM011()
            if strcmp(controls.cm011Condition.Value,'inclined'),next.cm011Boundary=cm011WorkbenchBoundary();
            else,next.cm011Boundary=[];end
            return;
        end
        altitude=controls.altitude.Value*1000;
        if altitude~=next.flight.altitude && isfield(next.flight,'atmosphere'),next.flight=rmfield(next.flight,'atmosphere');end
        next.flight.altitude=altitude;next.flight.Mach=controls.mach.Value;
        next.flight.alpha=controls.alpha.Value;next.flight.beta=controls.beta.Value;
        if controls.wall.Value~=next.flight.Twall && isfield(next.heating,'wallTemperatureByFace')
            next.heating=rmfield(next.heating,'wallTemperatureByFace');
        end
        next.flight.Twall=controls.wall.Value;
        if isEmbedded()
            next.heating.airThermalModel='caloricallyPerfect';
            property=controls.thermodynamics.Value;
            if ~strcmp(property,'idealGas') || isfield(next.su2,'thermodynamicModel'),next.su2.thermodynamicModel=property;end
            thermal=controls.thermalMode.Value;
            if ~strcmp(thermal,'isothermal') || isfield(next.su2,'wallThermalMode'),next.su2.wallThermalMode=thermal;end
            next.su2.iterations=controls.iterations.Value;next.su2.maximumSeconds=controls.seconds.Value;
            named=strcmp(controls.checkpointMode.Value,'indexed');every=controls.outputEvery.Value;
            % 旧快照未含默认字段时保持结构不变，避免只确认就使同工况复开失效。
            if named || isfield(next.su2,'namedCheckpoints'),next.su2.namedCheckpoints=named;end
            if every~=1000 || isfield(next.su2,'outputEvery'),next.su2.outputEvery=every;end
            initial=controls.initialMode.Value;thickness=controls.initialThickness.Value/1000;
            if ~strcmp(initial,'freestream') || isfield(next.su2,'initialMode'),next.su2.initialMode=initial;end
            if ~strcmp(initial,'freestream') || isfield(next.su2,'initialThicknessM'),next.su2.initialThicknessM=thickness;end
            limit=controls.checkpointLimit.Value;
            if limit~=100 || isfield(next.su2,'checkpointLimitMiB'),next.su2.checkpointLimitMiB=limit;end
            source=controls.checkpointSource.Value;
            if ~isempty(source) || isfield(next.su2,'checkpointSource'),next.su2.checkpointSource=source;end
            heosSource=controls.heosSource.Value;
            if ~isempty(heosSource) || isfield(next.su2,'heosSource'),next.su2.heosSource=heosSource;end
            flow=controls.flowModel.Value;
            if ~strcmp(flow,'laminar') || isfield(next.su2,'flowModel'),next.su2.flowModel=flow;end
            if strcmp(flow,'SST-2003m')
                assert(ismember(property,{'idealGas','idealAir'}),'HEOS RANS尚未适配，请选择定比热或变比热理想空气SST。');
                assert(strcmp(thermal,'isothermal') && ismember(initial,{'sstTransfer','sstWallDistance'}) && named,'SST需要统一壁温、明确初场及具名保存。');
                if strcmp(initial,'sstTransfer')
                    assert(isfile(fullfile(source,'engine_process.json')),'请选择同网格层流或已审核SST源场；启动时核对原始身份。');
                else
                    assert(~isempty(embeddedVolume) && strcmp(embeddedVolume.sourceKind,'original-triangle-hybrid-volume'));
                    assert(embeddedVolume.meshConstruction.tetsBelowSICN01==0,'完整域质量未通过，不能使用新壁距初场。');
                    validateattributes(thickness,{'numeric'},{'scalar','finite','positive'});next.su2.initialThicknessM=thickness;
                end
                next.su2.turbulenceIntensity=controls.turbulenceIntensity.Value;
                next.su2.turbulentViscosityRatio=controls.turbulentViscosityRatio.Value;
                next.su2.turbulentPrandtl=controls.turbulentPrandtl.Value;
                next.su2.turbulenceProvenance=controls.turbulenceProvenance.Value;
                assert(~isempty(strtrim(next.su2.turbulenceProvenance)),'请说明湍流参数的来源或假设。');
            else
                assert(~ismember(initial,{'sstTransfer','sstWallDistance'}),'层流模型不能使用SST初态。');
            end
            if strcmp(property,'HEOS-Air')
                assert(strcmp(initial,'propertyTransfer'),'HEOS必须选择“定比热源场转HEOS”，不能直接使用旧能量或自由流。');
                assert(strcmp(thermal,'isothermal'),'HEOS工作台目前只支持给定统一壁温。');
                assert(named,'HEOS研究须启用具名检查点保存。');
                assert(isfile(fullfile(heosSource,'engine_process.json')),'请选择已完成的定比热源场；启动时核对全部网格与状态。');
            else
                assert(~strcmp(initial,'propertyTransfer'),'当前物性入口需要自身相容初态。');
            end
            if strcmp(property,'idealAir')
                assert(strcmp(flow,'SST-2003m') && strcmp(initial,'sstWallDistance') && strcmp(thermal,'isothermal') && named, ...
                    '变比热主入口需要SST、实际壁距初场、统一给定Tw与具名保存。');
                validateattributes(next.flight.Twall,{'numeric'},{'scalar','finite','>=',60,'<=',2000});
            end
            if strcmp(initial,'checkpoint')
                assert(isfile(fullfile(source,'engine_process.json')),'请选择具有真实进程记录的源检查点目录。');
                assert(named,'继续研究须启用具名检查点保存，以保留本段和累计身份。');
            end
            if strcmp(initial,'wallEnthalpy')
                assert(strcmp(thermal,'isothermal'),'壁面相容初猜目前只支持给定统一壁温。');
                assert(~isempty(embeddedVolume) && strcmp(embeddedVolume.sourceKind,'original-triangle-hybrid-volume'),'壁面相容初猜需要原三角混合体网格。');
                assert(~isfield(embeddedVolume.meshConstruction,'noseRetainedLayers') && ~isfield(embeddedVolume.meshConstruction,'sideNodeIDs'), ...
                    '局部终止或尾缘扇网格尚未适配该初猜；不能使用等层编号推算近壁厚度。');
                nv=size(embeddedVolume.sourceVertices,1);nl=embeddedVolume.meshConstruction.prisms/size(state.mesh.F,1);
                height=vecnorm(embeddedVolume.vertices(nl*nv+(1:nv),:)-embeddedVolume.sourceVertices,2,2);
                assert(all(height>=thickness),'初猜厚度超出完整棱柱区，请调整网格或厚度。');
            end
            validateattributes(every,{'numeric'},{'scalar','integer','positive'});
            validateattributes(next.su2.iterations,{'numeric'},{'scalar','integer','positive'});return;
        end
        next.heating.airThermalModel=controls.air.Value;
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
    function thermodynamicsEdited(~,~)
        % 用户显式选物性后设置唯一相容初态；不自动启动或改变已有源文件。
        if strcmp(controls.thermodynamics.Value,'HEOS-Air')
            controls.initialMode.Value='propertyTransfer';controls.checkpointMode.Value='indexed';
        elseif strcmp(controls.thermodynamics.Value,'idealAir')
            controls.flowModel.Value='SST-2003m';controls.initialMode.Value='sstWallDistance';controls.checkpointMode.Value='indexed';
        elseif strcmp(controls.initialMode.Value,'propertyTransfer')
            controls.initialMode.Value='freestream';
        end
        display=state.cfg;display.su2.thermodynamicModel=controls.thermodynamics.Value;
        controls.propertySummary.Text=propertySummaryText(display);boundaryEdited([],[]);
    end
    function flowModelEdited(~,~)
        % 显式切换模型才选择新物理初猜；HEOS冲突保留为拒绝，不能静默降级物性。
        if strcmp(controls.flowModel.Value,'SST-2003m')
            controls.initialMode.Value='sstTransfer';controls.checkpointMode.Value='indexed';
            if strcmp(controls.thermodynamics.Value,'idealAir'),controls.initialMode.Value='sstWallDistance';end
        elseif ismember(controls.initialMode.Value,{'sstTransfer','sstWallDistance'}),controls.initialMode.Value='freestream';end
        display=state.cfg;display.su2.flowModel=controls.flowModel.Value;
        controls.propertySummary.Text=propertySummaryText(display);boundaryEdited([],[]);
    end
    function text=propertySummaryText(c)
        if strcmp(getOption(c.su2,'thermodynamicModel','idealGas'),'idealAir')
            text='SU2 / 变比热理想空气 / SST-2003m；Lemmon 2000，OTH参考，60–2000 K。保持输入Mach、静温/静压与流向，使用新声速与实际速度；一般EOS Roe，同阶段审核。';
        elseif strcmp(getOption(c.su2,'thermodynamicModel','idealGas'),'HEOS-Air')
            text='SU2 8.6.0 / CoolProp 7.2 HEOS Air：层流、ROE/特征远场；输运仍为Sutherland、μcp/Pr0.72。M输入是原速度参考；EOS审核通过后才能查看。';
        elseif strcmp(getOption(c.su2,'flowModel','laminar'),'SST-2003m')
            text='SU2 8.6.0 / SST-2003m全湍流研究；定比热/Sutherland/Pr0.72。入口I、μt/μ与Prt显式给定；新增k/ω和能量约定独立审核，未通过物理资格。';
        else
            text='SU2 8.6.0：层流NS、定比热理想气体γ1.4、Sutherland、Pr0.72；无湍流/转捩模型。空间壁温使用原生输入适配版本。';
        end
    end
    function boundaryEdited(~,~)
        if state.busy || preparing,return;end
        if isEmbedded() && isfield(controls,'boundaryLayout'),layoutEmbeddedBoundary();end
        try
            state=aeroHeatWorkflowState('setBoundary',state,struct('cfg',boundaryConfiguration()));report=[];
            if isBead() && isfield(controls,'wireBoundary')
                renderEditor();message('引线热边界已更新，请确认当前工况。');
            elseif isBead()
                controls.beadWall.Enable=onoff(strcmp(controls.beadCondition.Value,'published'));
                controls.beadPower.Enable=onoff(strcmp(controls.beadCondition.Value,'power'));
                controls.confirm.Enable=onoff(~isempty(state.cfg.thermistorBoundary));
                message('热敏珠边界选择已更新，请确认工况。');showGeometry(false);return;
            end
            if isCM011()
                controls.confirm.Enable=onoff(~isempty(state.cfg.cm011Boundary));
                message('试件工况已更改，请重新确认。');showGeometry(false);return;
            end
            if isEmbedded() && isfield(controls,'storageEstimate'),controls.storageEstimate.Text=checkpointBudgetText(state.cfg);end
            if isEmbedded(),controls.wall.Enable=onoff(strcmp(controls.thermalMode.Value,'isothermal'));end
            updateMethodScope();message('边界设置已编辑；请重新确认并显式计算，旧结果不会冒充新工况。');showGeometry(false);
        catch ex
            state=aeroHeatWorkflowState('invalidateBoundary',state);report=[];showGeometry(false);showError(ex);
        end
    end
    function layoutEmbeddedBoundary()
        % 模块用途：把物性/流动模型放在顶部，仅显示当前模型需要的来源和参数，去掉空白行。
        b=controls.boundaryLayout;sst=strcmp(controls.flowModel.Value,'SST-2003m');
        heos=strcmp(controls.thermodynamics.Value,'HEOS-Air');axial=strcmp(controls.thermalMode.Value,'axialTemperature');
        initial=controls.initialMode.Value;source=ismember(initial,{'checkpoint','sstTransfer'});
        order=[1 47 2:6 11 12 28 29 31 32 36 37 48:51 7:10 30 13 14 18 21 24 33 38 41 45 53 57];
        hidden=[];
        if ~sst,hidden=[hidden 48:51 53 57];end
        if ~heos,hidden=[hidden 36:38 41];end
        if ~axial,hidden=[hidden 12 21];end
        if ~ismember(initial,{'wallEnthalpy','sstWallDistance'}),hidden=[hidden 29];end
        if ~source,hidden=[hidden 31:33];end
        order=order(~ismember(order,hidden) & ismember(order,b.rows));
        for j=1:numel(b.items)
            row=find(order==b.rows(j));
            if isempty(row),b.items(j).Visible='off';b.items(j).Layout.Row=1;
            else,b.items(j).Visible='on';b.items(j).Layout.Row=row;end
        end
        height=repmat({31},1,numel(order));
        for j=find(ismember(order,[14 18 21 24 33 38 41 45 53 57])),height{j}=90;end
        b.grid.RowHeight=height;
    end
    function confirmBoundary(~,~)
        if state.busy || preparing,return;end
        try
            next=boundaryConfiguration();
            if isBead()
                [~,fs]=thermistorWorkbenchBoundary(next.thermistorBoundary,next.thermistorAngularCells,next);
                state=aeroHeatWorkflowState('setBoundary',state,struct('cfg',next));
                state=aeroHeatWorkflowState('confirmBoundary',state,struct('freestream',fs));
                message(sprintf('已确认 %g m/s、293.15 K 空气，%s。',fs.speed,beadBoundaryText(next.thermistorBoundary)));go(4);return;
            end
            if isCM011()
                [~,fs]=cm011WorkbenchBoundary(next.cm011Boundary);
                state=aeroHeatWorkflowState('setBoundary',state,struct('cfg',next));
                state=aeroHeatWorkflowState('confirmBoundary',state,struct('freestream',fs));
                message('边界已确认：0.1 m/s，攻角 +10°，给定壁温 300.1 K。');go(4);return;
            end
            if isEmbedded(),fs=calcSU2Freestream(next);else,fs=calcFreestream(next.flight);end
            if isEmbedded() && strcmp(getOption(next.su2,'wallThermalMode','isothermal'),'axialTemperature')
                assert(isfield(next.su2,'wallTemperatureProfile'),'请先导入覆盖整个模型的X—壁温表。');
                resolveSU2ThermalBoundary(state.mesh,next.su2.wallTemperatureProfile);
            end
            if fs.Mach<=1,error('RocketAeroHeat:Workbench:Mach','当前完整工程链要求M>1。');end
            state=aeroHeatWorkflowState('setBoundary',state,struct('cfg',next));
            state=aeroHeatWorkflowState('confirmBoundary',state,struct('freestream',fs));
            if isempty(state.result) && isempty(state.researchResult)
                message(sprintf('边界已确认 · p∞ %.6g Pa · T∞ %.6g K · M %.6g · 速度 %.6g m/s',fs.pressure,fs.temperature,fs.Mach,fs.speed));
            else
                message('边界与当前计算快照一致，结果仍有效；再次计算需明确按计算按钮。');
            end
            go(4);
        catch ex
            state=aeroHeatWorkflowState('invalidateBoundary',state);report=[];showGeometry(false);showError(ex);
        end
    end

    function importWallTemperatureProfile(~,~)
        % 两列显式SI表导入后把数值冻结到当前边界，不依赖外部文件继续存在。
        if state.busy || preparing,return;end
        [file,folder]=uigetfile('*.csv','导入壁温表（列名X_m、T_K）');if isequal(file,0),return;end
        try
            table=readtable(fullfile(folder,file),'VariableNamingRule','preserve');
            assert(isequal(table.Properties.VariableNames,{'X_m','T_K'}),'CSV须严格包含X_m、T_K两列。');
            setWallTemperatureProfile(struct('mode','axialTemperature','x',table.X_m, ...
                'temperature',table.T_K,'source',['用户导入：' file]));
        catch ex,showError(ex);end
    end
    function setWallTemperatureProfile(profile)
        assert(~state.busy && ~preparing && isEmbedded() && ~isempty(state.mesh),'请先完成当前表面网格。');
        assert(strcmp(profile.mode,'axialTemperature'));resolveSU2ThermalBoundary(state.mesh,profile);
        next=state.cfg;next.su2.wallThermalMode='axialTemperature';next.su2.wallTemperatureProfile=profile;
        state=aeroHeatWorkflowState('setBoundary',state,struct('cfg',next));report=[];
        go(3);message('已导入全模型壁温分布；这是给定边界，旧结果已失效，请重新确认。');
    end

%% 模块十：唯一物理计算入口；输入冻结，阶段失败不接受中间快照。
    function solveCase(~,~,source)
        if state.busy || preparing,return;end
        if isEmbedded(),startEngine();return;end
        try
            if nargin<3,source='solve';end
            request=struct();
            if isFullBead()
                if strcmp(source,'solve'),source='run';end
                request.startupAction=source;
                if ~strcmp(source,'run'),request.runDirectory=strtrim(controls.beadDirectory.Value);end
            elseif isMixedBead()
                request.leadAction=source;
                if any(strcmp(source,{'resume','resumeWire'})),request.runDirectory=strtrim(controls.beadDirectory.Value);end
                if strcmp(source,'solve') && isfield(controls,'wireSeed'),request.wireSeedDirectory=strtrim(controls.wireSeed.Value);end
            end
            state=aeroHeatWorkflowState('beginSolve',state,request);token=state.solveId;revision=state.solveInputRevision;
            if isfield(state.frozenInputs.cfg,'thermistorInvocation')
                c=state.frozenInputs.cfg;
                if ~isfield(c.thermistorBoundary,'leadThermalModel'),measuredDirectory=c.thermistorRunDirectory;end
                measuredPauseFile=[c.thermistorRunDirectory '.' c.thermistorInvocation '.pause'];
            end
            if isBead() && isfield(state.frozenInputs.cfg,'thermistorInvocation') && isfield(state.cfg.thermistorBoundary,'leadThermalModel')
                lastWireDirectory=state.frozenInputs.cfg.thermistorRunDirectory;
            end
            renderEditor();
            frozen=state.frozenInputs;
            if isBead()
                frozen.cfg.thermistorResultSource=source;
                if isFullBead()
                    message('正在核对全周计算的网格、初场和物理时刻…');
                elseif isfield(state.cfg.thermistorBoundary,'resistanceTemperatureK'),message('正在核对测温模型、长引线壁温和原面热流…');
                elseif isfield(state.cfg.thermistorBoundary,'leadThermalModel') && ~strcmp(source,'archive'),message(['长引线计算目录 ' lastWireDirectory]);
                elseif ~isempty(measuredPauseFile),message(['当前冷启动目录 ' measuredDirectory]);
                elseif isSidewind() && state.cfg.thermistorAngularCells==96,message('正在核对 96 档全表面测点和珠部平均换热…');
                elseif isSidewind() && strcmp(source,'archive'),message('正在核对全周侧吹归档与七区原面结果…');
                elseif strcmp(source,'archive'),message('正在核对体网，并用归档温差场计算当前壁温结果…');
                else,message('正在使用已确认的体网从冷初场计算热敏珠…');end
            else,message('正在显式计算当前工况；各阶段完成不代表全面积已求解。');end
            [solved,fs,newReport]=quietTriple(@()runAeroHeatSingleCase(frozen.mesh,frozen.cfg,@progress));
            measuredPauseFile='';
            state=aeroHeatWorkflowState('completeSolve',state,struct('result',solved, ...
                'freestream',fs,'solveId',token,'sourceRevision',revision));report=newReport;
            if isfield(newReport,'coverageMessage')
                message(newReport.coverageMessage);
            else
                message(sprintf('当前区域快照已接受，热流候选面积 %.4g%%。',100*newReport.heatComputedAreaFraction));
            end
            go(5);
        catch ex
            measuredPauseFile='';
            if exist('token','var') && state.busy
                operation='failSolve';if strcmp(ex.identifier,'RAH:Bead:Paused'),operation='pauseSolve';end
                state=aeroHeatWorkflowState(operation,state,struct('error',ex,'solveId',token,'sourceRevision',revision));
            end
            if isgraphics(fig)
                renderEditor();
                if strcmp(ex.identifier,'RAH:Bead:Paused'),message(ex.message);else,showError(ex);end
                showGeometry(false);
            end
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

%% 模块十一：内嵌引擎在原五步骤中运行；研究结果独立于工程覆盖状态。
    function yes=isEmbedded()
        yes=strcmp(getOption(state.cfg,'computationBackend','engineering'),'embeddedSU2Research');
    end
    function yes=isCM011()
        yes=strcmp(getOption(state.cfg,'computationBackend','engineering'),'cm011Imported');
    end
    function label=beadBoundaryText(b)
        % 计算页和确认回执显示实际冻结的热边界与单位。
        if isfield(b,'netHeatingPowerW') && hasBeadLeads(),label=sprintf('整件净对流功率 %.12g mW',1000*b.netHeatingPowerW);
        elseif isfield(b,'netHeatingPowerW'),label=sprintf('整珠净功率 %.12g mW',1000*b.netHeatingPowerW);
        elseif isfield(b,'resistanceTemperatureK'),label=sprintf('电阻等效内部温度 %.6g K、长引线温度分布',b.resistanceTemperatureK);
        elseif isfield(b,'leadThermalModel'),label=sprintf('线根 %.6g K、长引线温度分布',b.wallTemperatureK);
        else,label=sprintf('给定壁温 %.12g K',b.wallTemperatureK);end
    end
    function yes=isBead()
        yes=~isempty(state.rawMesh) && strcmp(getOption(state.cfg,'computationBackend','engineering'),'thermistorB');
    end
    function yes=isFullBead()
        yes=isBead() && state.cfg.thermistorAngularCells==48 && ~isSidewind();
    end
    function yes=isSidewind()
        yes=isBead() && strcmp(getOption(state.cfg,'thermistorVariant','smooth'),'sidewind');
    end
    function yes=hasBeadLeads()
        yes=isBead() && any(strcmp(getOption(state.cfg,'thermistorVariant','smooth'),{'leads','sidewind'}));
    end
    function yes=isMixedBead()
        yes=isBead() && any(strcmp(getOption(state.cfg,'thermistorVariant','smooth'),{'leads','midbody','sidewind'}));
    end
    function label=mixedBeadLabel()
        if isSidewind(),label='带线珠全周侧吹';return;end
        if hasBeadLeads()
            label='两端带引线';
            level=getOption(state.cfg,'thermistorRefineLeads',false);
            if level,[~,~,name]=thermistorLeadGridChoice(level);label=['带线珠' name];end
        else,label='光滑珠中段加密';end
    end
    function chooseBeadDirectory(~,~)
        % 这里只选择现有运行目录，读取、接续与新算分别由计算页按钮触发。
        base=fullfile(state.cfg.projectRoot,'outputs');
        if isFullBead(),base=fullfile(base,'v294_validation');end
        folder=uigetdir(base,'选择热敏珠计算根目录');
        if ~isequal(folder,0),controls.beadDirectory.Value=folder;end
    end
    function pauseMeasured(~,~)
        % 暂停请求只属于本次侧栏操作。求解器先保存、核验正在运行的阶段。
        if isempty(measuredPauseFile),return;end
        flag=java.io.File(measuredPauseFile);flag.createNewFile();
        message(['等待当前阶段保存。任务目录 ' state.frozenInputs.cfg.thermistorRunDirectory]);
    end
    function loadBeadModel(~,~,variant)
        % 固定朝向的民用试件从模型步骤进入主流程，原面在网格步骤继续使用。
        if state.busy || preparing,return;end
        try
            if nargin<3,variant='smooth';end
            if strcmp(variant,'smooth'),cells=2;else,cells=4;end
            [raw,source]=thermistorWorkbenchModel(variant);
            next=state.cfg;next.computationBackend='thermistorB';next.thermistorAngularCells=cells;next.thermistorBoundary=[];
            next.thermistorVariant=variant;
            next.thermistorRefineLeads=false;
            state=aeroHeatWorkflowState('setModel',state,struct('rawMesh',raw,'orientedMesh',raw,'modelSource',source,'cfg',next));
            embeddedVolume=[];qualityPreview=false;creasePreview=false;candidatePreview=[];preview=raw;report=[];
            go(2);
            if isMixedBead(),message([mixedBeadLabel() '已加载，可以生成体网。']);
            else,message('已加载 B 型热敏珠，沿长轴来流。请选择周向分辨率。');end
        catch ex,showError(ex);end
    end
    function leaveBead(~,~)
        if state.busy || preparing,return;end
        next=state.cfg;next.computationBackend='engineering';next=rmfield(next,{'thermistorAngularCells','thermistorBoundary'});
        if isfield(next,'thermistorVariant'),next=rmfield(next,'thermistorVariant');end
        if isfield(next,'thermistorRefineLeads'),next=rmfield(next,'thermistorRefineLeads');end
        state=aeroHeatWorkflowState('invalidateModel',state,struct('cfg',next));
        report=[];preview=[];renderEditor();showGeometry(true);message('请选择新的输入模型。');
    end
    function changeBeadGrid(field,~)
        if state.busy || preparing || field.Value==state.cfg.thermistorAngularCells,return;end
        next=state.cfg;next.thermistorAngularCells=field.Value;next.thermistorBoundary=[];
        state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',next));
        report=[];renderEditor();showGeometry(true);message('周向分辨率已改变，请生成新网格并确认边界。');
    end
    function changeLeadGrid(field,~)
        % 带线几何共用同一原面表。更换体网后，清除原计算和热边界确认。
        refine=field.Value;
        if state.busy || preparing || refine==getOption(state.cfg,'thermistorRefineLeads',false),return;end
        next=state.cfg;next.thermistorRefineLeads=refine;next.thermistorBoundary=[];
        if refine==4,next.thermistorVariant='sidewind';next.thermistorAngularCells=96;
        elseif refine==3,next.thermistorVariant='sidewind';next.thermistorAngularCells=48;
        else,next.thermistorVariant='leads';next.thermistorAngularCells=4;end
        state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',next));
        report=[];renderEditor();showGeometry(true);message('带线体网选项已改变，请生成网格并确认边界。');
    end
    function meshBackendChoice(grid,row)
        controls.meshBackend=choice(grid,row,'网格方法',{'工程表面重构','CM-011 低速半域'}, ...
            {'engineering','cm011Imported'},getOption(state.cfg,'computationBackend','engineering'));
        controls.meshBackend.ValueChangedFcn=@changeMeshBackend;
    end
    function changeMeshBackend(field,~)
        if state.busy || preparing,return;end
        next=state.cfg;next.computationBackend=field.Value;next.cm011Boundary=[];
        state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',next));report=[];
        renderEditor();showGeometry(true);message('网格方法已选择，请生成当前网格。');
    end
    function changeCM011Resolution(field,~)
        % 周向离散变化会使网格、边界确认和结果页失效。
        if state.busy || preparing,return;end
        next=state.cfg;next.cm011HalfSectors=field.Value;next.cm011Boundary=[];
        state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',next));
        report=[];renderEditor();showGeometry(true);
        message('半周单元数已改变，请重新生成网格。');
    end
    function loadNative(~,~,featureEnabled)
        if state.busy || preparing,return;end
        if nargin<3,featureEnabled=false;end
        try
            [raw,truth]=buildRocketTestBody(struct('resolution',24,'protrusionEnabled',featureEnabled));
            next=state.cfg;next.computationBackend='embeddedSU2Research';
            next.su2=getOption(next,'su2',struct);next.heating.airThermalModel='caloricallyPerfect';
            if featureEnabled,next.su2.volumeMode='hybridTriangles';end
            next.geometry.isAligned=true;next.geometry.bodyAxis=[1 0 0];next.geometry.noseDirection=[1 0 0];
            state=aeroHeatWorkflowState('setModel',state,struct('rawMesh',raw,'orientedMesh',raw,'cfg',next, ...
                'modelSource',struct('kind','native-rocket','units','m','truth',truth)));
            embeddedVolume=[];preview=raw;report=[];
            message('原生完整箭体已建立，坐标由参数化生成器明确给定；下一步生成体网格。');go(2);
        catch ex,showError(ex);end
    end
    function useSurfaceForSU2(surface)
        % 将用户已确认的当前表面原样接入；测试可传入同一原几何，绝不换成原生示例。
        assert(~state.busy && ~preparing,'当前有任务运行。');
        if nargin<1,surface=state.mesh;end
        assert(~isempty(surface),'请先生成并确认当前表面网格。');
        next=state.cfg;next.computationBackend='embeddedSU2Research';next.su2=getOption(next,'su2',struct);
        next.su2.volumeMode='hybridTriangles';next.heating.airThermalModel='caloricallyPerfect';
        state=aeroHeatWorkflowState('setModel',state,struct('rawMesh',surface,'orientedMesh',surface,'cfg',next, ...
            'modelSource',struct('kind','retained-original-triangles','units','m')));
        embeddedVolume=[];preview=surface;report=[];go(2);
        message('已保留当前全部原三角面；下一步生成混合体网格，旧工程结果已失效。');
    end
    function loadClassic()
        % 经典实验尺寸和来流为显式模型输入；底盖保留但不称已复现实验支杆。
        if state.busy || preparing,return;end
        try
            [raw,truth]=buildSchaeferFergusonBody();condition=schaeferFergusonCondition();
            next=state.cfg;next.computationBackend='embeddedSU2Research';next.su2=getOption(next,'su2',struct);
            next.su2.volumeMode='hybridTriangles';next.su2.wallThermalMode='adiabatic';
            next.su2.hybridMesh=struct('layerHeights',[0 cumsum(2e-6*1.25.^(0:23))], ...
                'outerMinSize',.006,'farSize',.08,'maximumSeconds',120);
            next.flight.Mach=condition.Mach;next.flight.alpha=0;next.flight.beta=0;next.flight.atmosphere=condition.atmosphere;
            next.heating.airThermalModel='caloricallyPerfect';next.geometry.isAligned=true;
            next.geometry.bodyAxis=[1 0 0];next.geometry.noseDirection=[1 0 0];
            state=aeroHeatWorkflowState('setModel',state,struct('rawMesh',raw,'orientedMesh',raw,'cfg',next, ...
                'modelSource',struct('kind','schaefer-ferguson-1961','units','m','truth',truth,'condition',condition)));
            embeddedVolume=[];preview=raw;report=[];go(2);
            message('经典完整锥柱扩张体已建立；绝热研究边界、指定实验来流，底盖/支杆差异及低温物性尚未验证。');
        catch ex,showError(ex);end
    end
    function generateNative()
        try
            preparing=true;message('正在建立完整三维体网格并核对闭合/体积…');
            candidate=state.orientedMesh;
            [options,mode]=su2WorkbenchMeshOptions(state.cfg);
            if strcmp(mode,'hybridTriangles')
                folder=fullfile(fileparts(mfilename('fullpath')),'outputs','hybrid_mesh',char(java.util.UUID.randomUUID()));
                embeddedVolume=buildSU2HybridVolume(candidate,folder,options);
            else
                assert(isfield(state.modelSource,'truth'),'当前原表面须使用混合网格，不能替换成轴对称模型。');
                embeddedVolume=buildSU2RocketVolume(candidate,state.modelSource.truth,options);
            end
            creasePreview=false;qualityPreview=false;
            state=aeroHeatWorkflowState('setMesh',state,struct('mesh',candidate));
            preparing=false;renderEditor();showGeometry(true);
            message(sprintf('体网格%d单元；原表面%d面未改变。近壁精度未合格，请确认边界。',numel(embeddedVolume.cells),size(candidate.F,1)));
        catch ex,preparing=false;renderEditor();showError(ex);end
    end
    function nativeMeshEdited(~,~)
        % 模式属于求解输入；清除旧体网格、边界确认和研究结果，禁止混用。
        if state.busy || preparing,return;end
        try
            next=state.cfg;next.su2.volumeMode=controls.volumeMode.Value;
            if isequaln(next,state.cfg),return;end
            state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',next));
            embeddedVolume=[];creasePreview=false;qualityPreview=false;researchUI=[];webUI=[];report=[];delete(resultHost.Children);
            renderEditor();showGeometry(false);message('体网格模式已改变；请重新生成网格并确认边界。');
        catch ex,showError(ex);end
    end
    function readHybridMeshDialog(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs','hybrid_mesh'),'选择本程序已完成的混合体网格目录');
        if isequal(folder,0),return;end
        try,readHybridMesh(folder);catch ex,showError(ex);end
    end
    function readHybridMesh(folder)
        % 复用本程序体网格前重新审核SHA、原面、正体积与边界；不加载外部CFD场。
        assert(~state.busy && ~preparing && isEmbedded() && ~isempty(state.orientedMesh));
        [~,mode]=su2WorkbenchMeshOptions(state.cfg);assert(strcmp(mode,'hybridTriangles'),'请先选择原三角面混合模式。');
        candidate=readSU2HybridVolume(folder,state.orientedMesh);
        source=load(candidate.meshReceipt.input,'optionsJSON');
        next=state.cfg;next.su2.hybridMesh=jsondecode(source.optionsJSON);
        next.su2.hybridMesh.layerHeights=next.su2.hybridMesh.layerHeights(:)';
        state=aeroHeatWorkflowState('setMeshOptions',state,struct('cfg',next));
        embeddedVolume=candidate;state=aeroHeatWorkflowState('setMesh',state,struct('mesh',state.orientedMesh));
        creasePreview=false;qualityPreview=false;
        researchUI=[];webUI=[];report=[];delete(resultHost.Children);
        renderEditor();showGeometry(true);
        message(sprintf('已核对并复用%d体单元，原表面不变；加载原网格参数，请重新确认边界。',numel(candidate.cells)));
    end
    function toggleCreasePreview(~,~)
        assert(~state.busy && ~isempty(embeddedVolume) && state.currentStep==2);
        candidatePreview=[];
        creasePreview=~creasePreview;qualityPreview=false;renderEditor();sync();showGeometry(true);
    end
    function toggleMeshQuality(~,~)
        assert(~state.busy && (~isempty(embeddedVolume) || ~isempty(candidatePreview)) && state.currentStep==2);
        if ~isempty(candidatePreview),candidatePreview=[];qualityPreview=true;end
        qualityPreview=~qualityPreview;creasePreview=false;renderEditor();sync();showGeometry(true);
    end
    function inspectMeshCandidateDialog(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs'),'选择已终态且质量未过的网格');
        if isequal(folder,0),return;end
        try,inspectMeshCandidate(folder);catch ex,showError(ex);end
    end
    function inspectMeshCandidate(folder)
        % 只读支路不赋值embeddedVolume，不修改网格、边界确认或下游计算资格。
        assert(~state.busy && ~preparing && isEmbedded() && state.currentStep==2 && ~isempty(state.orientedMesh));
        [m,info]=inspectSU2MeshCandidate(folder,state.orientedMesh);
        candidatePreview=struct('mesh',m,'info',info,'revision',state.revision);qualityPreview=true;creasePreview=false;
        renderEditor();sync();showGeometry(true);message('已定位失败候选；当前计算输入保持，候选未被接收。');
    end
    function geometryMessage(src,event)
        message=src.Data;geometryReceiver(src,event);
        if isstruct(message) && isfield(message,'kind') && strcmp(message.kind,'geometry-face-probe')
            probeMeshQuality(message);
        end
    end
    function probeMeshQuality(message)
        % 只接受当前显示版本的真实面索引；切页、旧帧与背景面不能读作体单元。
        if ~qualityPreview || state.currentStep~=2 || isempty(fig.UserData.meshQuality),return;end
        current=geometryView.html.UserData;
        if ~isfield(message,'revision') || message.revision~=current.revision || ~isfield(message,'faceIndex'),return;end
        i=double(message.faceIndex)+1;info=fig.UserData.meshQuality;
        if ~isscalar(i) || ~isfinite(i) || i~=fix(i) || i<1 || i>numel(info.displaySourceCellIDs),return;end
        id=info.displaySourceCellIDs(i);
        if id==0,controls.meshQualityProbe.Text='原壁面定位参考；这里没有体单元SICN或物理求解值。';return;end
        row=find(info.sourceCellIDs==id,1);xyz=info.cellCenters(row,:);
        controls.meshQualityProbe.Text=sprintf('体单元 %d · SICN %.6g\n中心 X/Y/Z = %.6g / %.6g / %.6g m\n等边参考Jacobian；网格指标，不是物理解。', ...
            id,info.cellSICN(row),xyz);
        fig.UserData.meshQualityProbe=struct('cellID',id,'SICN',info.cellSICN(row),'center',xyz,'geometryOnly',true);
    end
    function startEngine()
        try
            assert(~isempty(embeddedVolume),'请先生成当前原生体网格。');
            budget=checkpointBudget(state.cfg);
            assert(budget.withinLimit,'预计检查点%.1f MiB超过%.0f MiB研究限制；请调整保存设置或显式上限。',budget.estimatedBytes/1024^2,budget.limitBytes/1024^2);
            state=aeroHeatWorkflowState('beginSolve',state);
            frozen=state.frozenInputs;fs=frozen.freestream;
            folder=fullfile(fileparts(mfilename('fullpath')),'outputs','embedded_workbench',char(java.util.UUID.randomUUID()));
            flight=struct('Mach',frozen.cfg.flight.Mach,'temperature',fs.temperature,'pressure',fs.pressure, ...
                'wallTemperature',frozen.cfg.flight.Twall,'alpha',frozen.cfg.flight.alpha,'beta',frozen.cfg.flight.beta, ...
                'iterations',getOption(frozen.cfg.su2,'iterations',1000), ...
                'outputEvery',getOption(frozen.cfg.su2,'outputEvery',1000), ...
                'namedCheckpoints',getOption(frozen.cfg.su2,'namedCheckpoints',false));
            if strcmp(getOption(frozen.cfg.su2,'wallThermalMode','isothermal'),'adiabatic')
                flight.thermalBoundary=struct('mode','adiabatic','source','工作台明确选择绝热边界');
            elseif strcmp(getOption(frozen.cfg.su2,'wallThermalMode','isothermal'),'axialTemperature')
                flight.thermalBoundary=frozen.cfg.su2.wallTemperatureProfile;
            end
            manifest=writeSU2RocketCase(embeddedVolume,folder,flight);
            if strcmp(getOption(frozen.cfg.su2,'initialMode','freestream'),'wallEnthalpy')
                writeSU2HybridInitialState(embeddedVolume,folder,frozen.cfg.su2.initialThicknessM);
                manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
            elseif strcmp(getOption(frozen.cfg.su2,'initialMode','freestream'),'checkpoint')
                writeSU2CheckpointInitialState(folder,frozen.cfg.su2.checkpointSource, ...
                    getOption(frozen.cfg.su2,'pythonExecutable','python'));
                manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));
            end
            projection=buildSU2WallProjection(embeddedVolume,frozen.mesh);
            save(fullfile(folder,'input_snapshot.mat'),'frozen','manifest','projection','budget','-v7');
            jobOptions=struct('threads',2, ...
                'maximumSeconds',getOption(frozen.cfg.su2,'maximumSeconds',120), ...
                'pythonExecutable',getOption(frozen.cfg.su2,'pythonExecutable','python'));
            if strcmp(getOption(frozen.cfg.su2,'thermodynamicModel','idealGas'),'idealAir')
                c=frozen.cfg.su2;
                jobOptions.idealAirInitialization=struct('model','ideal-air-SST','mode','wallDistance', ...
                    'thicknessM',c.initialThicknessM,'intensity',c.turbulenceIntensity, ...
                    'viscosityRatio',c.turbulentViscosityRatio,'Prt',c.turbulentPrandtl, ...
                    'provenance',c.turbulenceProvenance,'maximumOutputBytes',budget.limitBytes);
            elseif strcmp(getOption(frozen.cfg.su2,'thermodynamicModel','idealGas'),'HEOS-Air')
                jobOptions.propertyInitialization=struct('model','HEOS-Air','runtimeFlavor','rah-coolprop-recon-v3', ...
                    'sourceDirectory',char(java.io.File(frozen.cfg.su2.heosSource).getCanonicalPath()));
            elseif strcmp(getOption(frozen.cfg.su2,'flowModel','laminar'),'SST-2003m')
                c=frozen.cfg.su2;
                jobOptions.turbulenceInitialization=struct('model','SST-2003m', ...
                    'sourceDirectory',char(java.io.File(getOption(c,'checkpointSource','')).getCanonicalPath()), ...
                    'intensity',c.turbulenceIntensity,'viscosityRatio',c.turbulentViscosityRatio, ...
                    'Prt',c.turbulentPrandtl,'provenance',c.turbulenceProvenance);
                if strcmp(getOption(c,'initialMode','sstTransfer'),'sstWallDistance')
                    jobOptions.turbulenceInitialization=rmfield(jobOptions.turbulenceInitialization,'sourceDirectory');
                    jobOptions.turbulenceInitialization.mode='wallDistance';
                    jobOptions.turbulenceInitialization.thicknessM=c.initialThicknessM;
                end
            end
            engineJob=embeddedSU2Job('start',folder,jobOptions);
            engineTimer=timer('ExecutionMode','fixedSpacing','Period',.5,'BusyMode','drop','TimerFcn',@pollEngine);
            renderEditor();message('内嵌引擎已启动；等待实际进程回执。');start(engineTimer);
        catch ex
            if state.busy,state=aeroHeatWorkflowState('failSolve',state,struct('error',ex,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));end
            renderEditor();showError(ex);
        end
    end
    function pollEngine(~,~)
        if isempty(engineJob),return;end
        try
            [engineJob,receipt]=embeddedSU2Job('poll',engineJob);
            if ~engineJob.terminal
                if isfield(engineJob,'pipelineReceipt')
                    p=engineJob.pipelineReceipt;labels=struct('preparation','初态物性转换','solver','原生流场求解','eosAudit','保存场EOS审核');
                    label=p.stage;if isfield(labels,label),label=labels.(label);end
                    message(sprintf('%s · 本阶段%.1fs / %.0fs · 总%.1fs · 取消%d',label,p.phaseElapsedSeconds,p.phaseBudgetSeconds,p.elapsedSeconds,engineJob.cancelRequested));
                elseif ~isempty(receipt) && isfield(receipt,'elapsedSeconds')
                    suffix='';if isfield(receipt,'progress'),suffix=jsonencode(receipt.progress);end
                    if engineJob.cancelRequested,label='取消请求已提交，等待真实终态';else,label='内嵌引擎运行中';end
                    message(sprintf('%s · %.1fs · %s',label,receipt.elapsedSeconds,suffix));
                end
                return;
            end
            stopEngineTimer();
            if closing,delete(fig);return;end
            if isfield(engineJob,'pipelineReceipt')
                p=engineJob.pipelineReceipt;
                assert(p.exitCode==0 && ~p.cancelled && ~p.timedOut, ...
                    '变物性流程在%s未完成：%s；原生求解器已启动%d。目录：%s',p.stage,getOption(p,'error','详见阶段回执'),p.engineStarted,engineJob.folder);
            end
            assert(engineJob.supervisorExitCode==0 && receipt.exitCode==0 && ~receipt.timedOut && ~receipt.cancelled, ...
                'RocketAeroHeat:SU2:Terminal','引擎未成功：exit %g / 取消%d / 超时%d。目录：%s',receipt.exitCode,receipt.cancelled,receipt.timedOut,engineJob.folder);
            stable=struct2cell(receipt.inputsStable);assert(all(cellfun(@(v)isequal(v,true),stable)),'引擎输入身份改变。');
            data=readSU2SurfaceResult(fullfile(engineJob.folder,'surface_flow.vtk'));
            projection=buildSU2WallProjection(embeddedVolume,state.mesh);
            result=mapSU2WallResult(data,projection);result.V=state.mesh.V;result.F=state.mesh.F;
            result.receipt=receipt;result.flight=state.frozenInputs.cfg.flight;
            result.wallTemperatureInput=result.flight.Twall;result.sourceRevision=state.revision;result.solveId=state.solveId;
            % 后台转换会更新真实运行清单，冻结UI快照仍保留原参考Mach；不可用旧清单覆盖EOS身份。
            manifest=jsondecode(fileread(fullfile(engineJob.folder,'case_manifest.json')));
            if isfield(engineJob,'pipelineReceipt')
                result.flight=manifest.flight;result.workflowReceipt=engineJob.pipelineReceipt;
            end
            result=attachSU2ThermalBoundary(result,manifest.thermalBoundary);
            result.resultDirectory=engineJob.folder;
            result.convergenceAudit=auditSU2ResultConvergence(engineJob.folder,result);
            if isfield(receipt,'thermodynamicModel')
                [~,~,result.thermalStateAudit]=readSU2IdealAirStateAudit(engineJob.folder,receipt);
                result.actualFreestreamSpeed=receipt.thermodynamicModel.freestreamSpeed;
                assert(norm(manifest.expectedVelocity(:)-state.frozenInputs.freestream.flowVelocityVector(:))<2e-10, ...
                    '物性参考速度与已确认来流不一致。');
            elseif isfield(receipt,'physicsModel')
                result.thermalStateAudit=readSU2SSTStateAudit(engineJob.folder,receipt);
            end
            if isfield(engineJob,'pipelineReceipt')
                assert(result.convergenceAudit.checkpointIdentityConfirmed,result.convergenceAudit.reason);
            end
            save(fullfile(engineJob.folder,'mapped_research.mat'),'result','-v7');
            state=aeroHeatWorkflowState('completeResearch',state,struct('result',result,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            message('实际进程已结束；原面元研究快照已接收，未授予物理资格。');go(5);
        catch ex
            if ~isempty(engineJob) && ~engineJob.process.HasExited
                engineJob=embeddedSU2Job('cancel',engineJob);showError(ex);return;
            end
            stopEngineTimer();
            if closing,delete(fig);return;end
            if state.busy,state=aeroHeatWorkflowState('failSolve',state,struct('error',ex,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));end
            renderEditor();showError(ex);
        end
    end
    function cancelEngine(~,~)
        if isempty(engineJob) || ~state.busy,return;end
        engineJob=embeddedSU2Job('cancel',engineJob);message('取消请求已提交；保持窗口响应并等待引擎实际结束。');
    end
    function readResearchDialog(~,~)
        [file,folder]=uigetfile('mapped_research.mat','重开本程序同几何/边界研究结果');
        if isequal(file,0),return;end
        try,readResearch(fullfile(folder,file));catch ex,showError(ex);end
    end
    function budget=checkpointBudget(c)
        budget=estimateSU2CheckpointStorage(embeddedVolume,getOption(c.su2,'iterations',1000), ...
            getOption(c.su2,'outputEvery',1000),getOption(c.su2,'namedCheckpoints',false),getOption(c.su2,'checkpointLimitMiB',100), ...
            getOption(c.su2,'flowModel','laminar'),getOption(c.su2,'thermodynamicModel','idealGas'));
    end
    function text=checkpointBudgetText(c)
        if isempty(embeddedVolume),text='生成体网格后估计输出体量。';return;end
        b=checkpointBudget(c);text=sprintf('检查点预计最多 %.1f MiB / 限制%.0f MiB；含基础文件及具名副本。估计不含网格、快照、求解内存；2线程。',b.estimatedBytes/1024^2,b.limitBytes/1024^2);
    end
    function chooseCheckpointInitial(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs'),'选择同工况已完成任务');
        if isequal(folder,0),return;end
        controls.checkpointSource.Value=folder;boundaryEdited([],[]);
    end
    function chooseHEOSSource(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs'),'选择同原p/T/速度与统一Tw的定比热源场');
        if isequal(folder,0),return;end
        controls.heosSource.Value=folder;boundaryEdited([],[]);
    end
    function readVerifiedHEOSDialog(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs'),'查看本程序已完成并有EOS审核的HEOS场');
        if isequal(folder,0),return;end
        try,readVerifiedHEOS(folder);catch ex,showError(ex);end
    end
    function readVerifiedHEOS(folder)
        % 研究入口读取原运行证据；不伪造历史UI冻结输入，不要求重新求解。
        result=readSU2VerifiedHEOSCase(folder,state,embeddedVolume);
        state=aeroHeatWorkflowState('beginSolve',state);
        try
            result.sourceRevision=state.revision;result.solveId=state.solveId;
            state=aeroHeatWorkflowState('completeResearch',state,struct('result',result, ...
                'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            message('原HEOS场与当前原面/物理条件一致，已核对EOS及具名步号；只读查看，未重新计算或增加资格。');go(5);
        catch ex
            state=aeroHeatWorkflowState('failSolve',state,struct('error',ex,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            renderEditor();rethrow(ex);
        end
    end
    function readBenchmarkDialog(~,~)
        [name,folder]=uigetfile('*.json','查看公开算例审核',fullfile(fileparts(mfilename('fullpath')),'outputs'));
        if isequal(name,0),return;end
        try,readBenchmarkReview(fullfile(folder,name));catch ex,showError(ex);end
    end
    function updateThermistorButton(~,~)
        % 切换选择时只更新操作文字。只有按下按钮才可能启动本次计算。
        if startsWith(controls.thermistorMode.Value,'solve')
            controls.thermistor.Text='从冷初场计算热敏珠';
        else
            controls.thermistor.Text='打开热敏珠全表面结果';
        end
    end
    function showThermistor(~,~)
        % 独立民用试件在现有中央区显示，页内探针与导出共用同一份原面表。
        if state.busy || preparing,return;end
        mode=controls.thermistorMode.Value;preparing=true;
        try
            if startsWith(mode,'solve'),message('正在计算 B 型热敏珠的三维壁压和热流…');
            else,message('正在读取 B 型热敏珠原面结果…');end
            [bead,summary]=runThermistorBeadCase([],mode,@()isgraphics(fig));
            if ~isgraphics(fig),return;end
            prior=resultHost.Children;
            group=uitabgroup(resultHost,'Tag','thermistor-bead-review','Visible','off', ...
                'Units','normalized','Position',[0 0 1 1]);
            view=createAerothermalWebTab(group,bead);
            delete(prior);webUI=[];researchUI=[];
            geometryHost.Visible='off';resultHost.Visible='on';group.Visible='on';
            fig.UserData.thermistorReview=view;preparing=false;renderEditor();
            controls.thermistorMode.Value=mode;updateThermistorButton([],[]);
            message(sprintf('热敏珠共 %d 面，整珠放热 %.4f mW，正热流入壁。',height(bead.hybridWall),-1000*summary.wholeSurfaceHeatPowerW));
            if strcmp(mode,'nearWall6')
                title.Text=sprintf('B 型热敏珠 · 近壁网格对照 · %.6g µs',summary.elapsedPhysicalTimeS*1e6);
            elseif strcmp(mode,'timeSteps6')
                title.Text=sprintf('B 型热敏珠 · 同刻时间步对照 · %.6g µs',summary.elapsedPhysicalTimeS*1e6);
            elseif strcmp(mode,'startup6')
                title.Text=sprintf('B 型热敏珠 · 6 m/s 启动瞬态 · %.6g µs',summary.elapsedPhysicalTimeS*1e6);
            else
                title.Text='B 型热敏珠 · 三维场与实验对照';
            end
        catch ex
            preparing=false;if ~isgraphics(fig),return;end
            if exist('group','var') && isgraphics(group),delete(group);end
            renderEditor();controls.thermistorMode.Value=mode;updateThermistorButton([],[]);showError(ex);
        end
    end
    function showLewisCylinder(~,~)
        % 显式运行独立风洞工况，在当前中央区显示结果并保留用户的模型输入。
        if state.busy || preparing,return;end
        try
            [mesh,stations,summary]=runLewisCylinderBenchmark();
            prior=resultHost.Children;
            view=createLewisCylinderReview(resultHost,mesh,stations,summary);
            delete(prior);webUI=[];researchUI=[];
            geometryHost.Visible='off';resultHost.Visible='on';view.tabs.Visible='on';
            sync();title.Text='NASA Lewis 圆柱 · 风洞测站对照';
            fig.UserData.lewisCylinderReview=view;
        catch ex,showError(ex);end
    end
    function readBenchmarkReview(file)
        % 公开基准是独立验证记录，按原终态显示并保留当前五步骤输入。
        assert(~state.busy && ~preparing,'请等待当前任务结束后查看公开验证记录。');
        review=readSU2BenchmarkReview(file);
        delete(resultHost.Children);webUI=[];researchUI=[];
        geometryHost.Visible='off';resultHost.Visible='on';
        view=createSU2BenchmarkReview(resultHost,review);
        sync();title.Text='公开算例验证记录';fig.UserData.benchmarkReview=view;
    end
    function readIdealAirDialog(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs'),'查看已审核的变比热热流算例');
        if isequal(folder,0),return;end
        try,readIdealAirReview(folder);catch ex,showError(ex);end
    end
    function readStaticAirProbe(varargin)
        % 模块用途：先核对归档，再在现有中央区显示节点页；保留当前工况和输入版本。
        if state.busy || preparing,return;end
        try
            page=buildStaticAirProbePage();
            prior=resultHost.Children;
            layout=uigridlayout(resultHost,[1 1],'Tag','static-air-probe-host','Visible','off');
            layout.Padding=[0 0 0 0];
            uihtml(layout,'Tag','static-air-probe','DataChangedFcn',@staticAirMessage,'HTMLSource',page);
            delete(prior);webUI=[];researchUI=[];
            geometryHost.Visible='off';resultHost.Visible='on';layout.Visible='on';
            sync();title.Text='静止空气 · 节点探针';
        catch ex
            if exist('layout','var') && isgraphics(layout),delete(layout);end
            showError(ex);
        end
    end
    function staticAirMessage(src,~)
        % 内嵌页面通过主窗口保存文件，独立 HTML 使用浏览器下载相同的 CSV。
        event=src.Data;
        if ~isstruct(event) || ~isfield(event,'kind') || ~strcmp(event.kind,'static-air-export'),return;end
        [file,folder]=uiputfile('*.csv','导出所选节点',event.filename);if isequal(file,0),return;end
        try
            fid=fopen(fullfile(folder,file),'wb');assert(fid>=0,'无法创建所选 CSV。');
            guard=onCleanup(@()fclose(fid));bytes=unicode2native(event.csv,'UTF-8');
            assert(fwrite(fid,bytes,'uint8')==numel(bytes),'CSV 写入不完整。');
        catch ex,showError(ex);end
    end
    function readIdealAirReview(folder)
        % 独立热物性接口记录在主中央区显示，当前箭体输入继续由工作流管理。
        assert(~state.busy && ~preparing,'请等待当前任务结束。');
        receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
        readSU2IdealAirStateAudit(folder,receipt);
        delete(resultHost.Children);webUI=[];researchUI=[];
        geometryHost.Visible='off';resultHost.Visible='on';
        view=createSU2IdealAirReview(resultHost,folder);
        sync();title.Text='变比热热流接口记录';fig.UserData.idealAirReview=view;
    end
    function readSSTStudyDialog(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs'),'查看本程序已审核的SST数值对照');
        if isequal(folder,0),return;end
        try,readSSTStudy(folder);catch ex,showError(ex);end
    end
    function readSSTStudy(folder)
        % 保留真实数值方案而复用同一个主结果区；只读操作不启动求解器。
        result=readSU2VerifiedSSTStudy(folder,state,embeddedVolume);
        state=aeroHeatWorkflowState('beginSolve',state);
        try
            result.sourceRevision=state.revision;result.solveId=state.solveId;
            state=aeroHeatWorkflowState('completeResearch',state,struct('result',result, ...
                'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            message('已核对SST原面/物理条件与数值方案，只读查看；没有续算或增加资格。');go(5);
        catch ex
            state=aeroHeatWorkflowState('failSolve',state,struct('error',ex,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            renderEditor();rethrow(ex);
        end
    end
    function readTimeoutDialog(~,~)
        folder=uigetdir(fullfile(fileparts(mfilename('fullpath')),'outputs'),'选择本程序原超时运行目录');
        if isequal(folder,0),return;end
        try,readTimeoutCheckpoint(folder);catch ex,showError(ex);end
    end
    function readTimeoutCheckpoint(folder)
        assert(isEmbedded() && ~state.busy,'请在当前内嵌工况的计算步骤显式重开。');
        result=readSU2TimeoutCheckpoint(folder,state,embeddedVolume);
        state=aeroHeatWorkflowState('beginSolve',state);
        try
            result.sourceRevision=state.revision;result.solveId=state.solveId;
            state=aeroHeatWorkflowState('completeCheckpointResearch',state,struct('result',result, ...
                'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            message('已显式重开完整超时检查点；原任务仍为超时、结果仍未获物理资格。');go(5);
        catch ex
            state=aeroHeatWorkflowState('failSolve',state,struct('error',ex,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            renderEditor();rethrow(ex);
        end
    end
    function readResearch(file)
        % 只复用本程序保存的同几何/同边界结果；输入版本不匹配必须拒绝。
        assert(~state.busy && state.boundaryConfirmed && isEmbedded(),'请先确认当前内嵌引擎工况。');
        saved=load(file,'result');result=saved.result;
        assert(isequal(result.V,state.mesh.V) && isequal(result.F,state.mesh.F),'保存结果不是当前原面元几何。');
        folder=fileparts(file);input=load(fullfile(folder,'input_snapshot.mat'),'frozen');
        assert(isequaln(input.frozen.mesh,state.mesh) && isequaln(input.frozen.cfg,state.cfg) && ...
            isequaln(input.frozen.freestream,state.freestream),'保存结果的边界/物性/网格设置与当前工况不同。');
        receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
        assert(isequaln(receipt,result.receipt),'保存结果与实际进程回执不一致。');
        auditFolder=folder;
        % 便携小快照可引用原运行目录；只有同一真实回执才可用其完整历史补审核。
        if (~isfile(fullfile(auditFolder,'history.csv')) || ~isfile(fullfile(auditFolder,'restart_flow.csv'))) && isfield(result,'resultDirectory') && ...
                isfile(fullfile(result.resultDirectory,'engine_process.json'))
            original=jsondecode(fileread(fullfile(result.resultDirectory,'engine_process.json')));
            if isequaln(original,receipt),auditFolder=result.resultDirectory;end
        end
        result.convergenceAudit=auditSU2ResultConvergence(auditFolder,result);
        if isfield(receipt,'thermodynamicModel')
            assert(result.convergenceAudit.checkpointIdentityConfirmed,result.convergenceAudit.reason);
            [~,~,result.thermalStateAudit]=readSU2IdealAirStateAudit(auditFolder,receipt);
        elseif isfield(receipt,'physicsModel')
            assert(result.convergenceAudit.checkpointIdentityConfirmed,result.convergenceAudit.reason);
            result.thermalStateAudit=readSU2SSTStateAudit(auditFolder,receipt);
        end
        origin=struct('sourceRevision',result.sourceRevision,'solveId',result.solveId,'file',file);
        state=aeroHeatWorkflowState('beginSolve',state);
        try
            result.reopenedFrom=origin;result.sourceRevision=state.revision;result.solveId=state.solveId;
            state=aeroHeatWorkflowState('completeResearch',state,struct('result',result,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            message('已重开本程序同几何/边界研究快照，保留原进程与输入身份；未重新计算。');go(5);
        catch ex
            state=aeroHeatWorkflowState('failSolve',state,struct('error',ex,'solveId',state.solveId,'sourceRevision',state.solveInputRevision));
            renderEditor();rethrow(ex);
        end
    end
    function stopEngineTimer()
        if ~isempty(engineTimer) && isvalid(engineTimer),stop(engineTimer);delete(engineTimer);end
        engineTimer=[];
    end
    function closeWorkbench(~,~)
        if ~isempty(measuredPauseFile),pauseMeasured([],[]);end
        if ~isempty(engineJob) && ~engineJob.process.HasExited
            closing=true;cancelEngine([],[]);return;
        end
        stopEngineTimer();delete(fig);
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
