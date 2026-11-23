function page=addShortSupportGridTab(tabs,geometry)
% 模块一：绑定已验收的双支管几何、两个真实体网和原生计算回执。
root=fileparts(mfilename('fullpath'));folder=fullfile(root,'cases','short_sphere');
file=fullfile(folder,'supportedGrid65.json');
assert(strcmp(cm011FileSHA256(file),'049332d1171345da4d5d6d9fc1fe099d07c1c83065d210bec77f2f8d46607820'));
report=jsondecode(fileread(file));
assert(strcmp(cm011FileSHA256(fullfile(folder,'supported65.zip')),report.bodySHA256) && ...
    strcmp(geometry.summary.caseSHA256,report.caseSHA256) && geometry.summary.extensionD==1);
fluxFile=fullfile(folder,'supportedFlux65.json');
assert(strcmp(cm011FileSHA256(fluxFile),'1c4fb8cff2b83f23ac4343b602b9b28a1c9d0f3026a4443928f366cc0b3c3a05'));
flux=jsondecode(fileread(fluxFile));
assert(strcmp(flux.caseSHA256,report.caseSHA256) && strcmp(flux.meshSHA256,report.hybrid.meshSHA256) && ...
    isequal({flux.frames.baselineSHA256},{report.hybrid.stages(end).frames.sha256}));
tab=uitab(tabs,'Title','流体体网与计算状态');
box=uigridlayout(tab,[6 2]);box.RowHeight={34,30,'1x',50,55,32};box.ColumnWidth={'1x','1x'};
box.Padding=[14 8 14 8];heading=uilabel(box,'Text','Short Test 65 · 封闭有限支管 · Re = 882','FontWeight','bold');
heading.Layout.Column=[1 2];
method=uidropdown(box,'Items',{'混合外区','结构外区对照','混合外区 · 质量通量对照'}, ...
    'ItemsData',{'hybrid','structured','flux'},'Value','hybrid');
zoom=uidropdown(box,'Items',{'球体与两侧接头','+Y 管根近壁层','全外域'},'ItemsData',{'body','joint','outer'},'Value','body');
ax=uiaxes(box);ax.Layout.Column=[1 2];hold(ax,'on');plotLine=plot(ax,nan,nan,'Color',[.28 .40 .49],'LineWidth',.5);
point=plot(ax,nan,nan,'o','MarkerFaceColor',[.74 .30 .12],'MarkerEdgeColor','white', ...
    'MarkerSize',9,'Visible','off','DisplayName','首次越界的最低值节点（径向投影）');
xlabel(ax,'沿支管轴 Y/D');ylabel(ax,'子午径向距离 r/D');axis(ax,'equal');
stats=uilabel(box,'WordWrap','on');stats.Layout.Column=[1 2];
outcome=uilabel(box,'WordWrap','on');outcome.Layout.Column=[1 2];
build=uibutton(box,'Text','从几何生成所选体网','ButtonPushedFcn',@buildGrid);
message=uilabel(box,'Text','','WordWrap','on');
method.ValueChangedFcn=@(~,~)drawMesh();zoom.ValueChangedFcn=@(~,~)drawMesh();drawMesh();
page=struct('tab',tab,'method',method,'zoom',zoom,'line',plotLine,'outcome',outcome,'stats',stats, ...
    'build',build,'message',message,'report',report,'flux',flux,'point',point);
    function kind=gridKind()
        % 两种质量通量计算共用同一套混合网，重建时按实际网格类型分派。
        kind=method.Value;
        if strcmp(kind,'flux'),kind='hybrid';end
    end
    function drawMesh()
        % 模块二：绘制实际体网母平面，诊断页另标首次越界节点绕支管轴的投影。
        kind=gridKind();d=report.(kind);p=d.section.pointsD;e=d.section.edges;
        x=reshape(p(e(:),1),size(e));y=reshape(p(e(:),2),size(e));blank=nan(size(e,1),1);
        plotLine.XData=reshape([x blank]',1,[]);plotLine.YData=reshape([y blank]',1,[]);
        switch zoom.Value
            case 'joint',xlim(ax,[.42 .59]);ylim(ax,[.06 .21]);
            case 'outer',xlim(ax,[-14.5 14.5]);ylim(ax,[-.25 14.5]);
            otherwise,xlim(ax,[-1.7 1.7]);ylim(ax,[-.04 1.1]);
        end
        geo=d.geometry;last=d.stages(end);frame=last.frames(end);
        stats.Text=sprintf('%d 个体单元，%d 个节点。%d 个求解壁片对应 %d 个原面，法向首层 %.2f–%.2f µm。', ...
            geo.cells,geo.nodes,geo.wallFaces,geo.transfer.targetFaces,1e6*geo.minFirstM,1e6*geo.maxFirstM);
        point.Visible='off';outcome.Tooltip='';
        outcome.Text=sprintf(['本次稳态推进 %d 次，正常保存终帧。主温差标量下限为 −1e−12，' ...
            '末帧最小 %.4g，%d 个节点低于下限。已接受壁压与热流面积 0 m²。'], ...
            sum([d.stages.iterations]),frame.minScalar,frame.below);
        message.Text='';
        if strcmp(method.Value,'flux')
            q=flux.first;xyz=q.xyzD;f=flux.frames(end);
            point.XData=xyz(2);point.YData=hypot(xyz(1),xyz(3));point.Visible='on';
            stats.Text=sprintf('首次越界在第 %d 次更新，共 %d 个节点。橙点为最低值节点 %d 的径向投影，距最近壁面节点 %.3fD。', ...
                q.update,numel(q.belowIDs),q.worstID,q.nearestWallNodeD);
            outcome.Text=sprintf('%s · 终帧 %d 的最小温差标量 %.4g，%d 个节点低于 %.0e。残差未收敛，合格壁压与热流面积 %g m²。', ...
                flux.method,f.index,f.minScalar,f.below,flux.lowerBound,flux.acceptedAreaM2);
            outcome.Tooltip=sprintf('温差标量 = (T − T∞) / %.0f K。终帧最低温差 %.3g K；原迎风格式同帧最低标量 %.4g，越界 %d 个。', ...
                flux.scalarUnitK,flux.scalarUnitK*f.minScalar,f.baselineMinimum,f.baselineBelow);
        end
    end
    function buildGrid(~,~)
        % 模块三：由同一源码重新生成体网，计算任务和浏览存档各有独立身份。
        kind=gridKind();job=fullfile(root,'outputs',['short65_grid_' char(java.util.UUID.randomUUID())]);
        build.Enable='off';method.Enable='off';zoom.Enable='off';
        unlock=onCleanup(@enableControls); %#ok<NASGU>
        message.Text='正在生成体网并检查原面对应…';drawnow;
        try
            if ~isfolder(fileparts(job)),mkdir(fileparts(job));end
            runCM011Python({fullfile(root,'scripts','shortSupportMesh.py'),job,'--outer',kind},[job '.log'],@()drawnow('limitrate'),true);
            if ~isgraphics(message),return;end
            p=jsondecode(fileread(fullfile(job,'prepared.json')));
            assert(strcmp(cm011FileSHA256(fullfile(job,'mesh.su2')),p.files.mesh_su2) && ...
                strcmp(p.files.mesh_su2,report.(kind).meshSHA256) && strcmp(p.supportedSHA256,report.bodySHA256));
            message.Text=sprintf('重建通过 · %d 单元',p.geometry.cells);message.UserData=job;
        catch ex
            if isgraphics(message),message.Text=ex.message;end
        end
    end
    function enableControls()
        if isgraphics(build),build.Enable='on';method.Enable='on';zoom.Enable='on';end
    end
end
