function view=createShortSphereReview(parent,mesh,stations,report)
% 模块一：将完整球面云图、面元证据和本次试验的逐站误差放入同一工作区。
tabs=uitabgroup(parent,'Tag','short-sphere-review','Visible','off','Units','normalized','Position',[0 0 1 1]);
try
    web=createAerothermalWebTab(tabs,mesh);tab=uitab(tabs,'Title','七个风洞测站');
    layout=uigridlayout(tab,[5 1]);layout.RowHeight={42,'1x',200,62,32};
    layout.Padding=[16 10 16 10];layout.RowSpacing=6;
    uilabel(layout,'Text',sprintf('Short 银球 · Test 65    |    Re = 882    |    %.4f m/s    |    直径 12.7 mm',report.condition.speedMps),'FontWeight','bold');
    wake=isfield(report,'elapsedPhysicalTimeS');label='SU2 稳态基流';
    if wake,label=sprintf('瞬时场 · %.4f ms',1000*report.elapsedPhysicalTimeS);end
    ax=uiaxes(layout);hold(ax,'on');plot(ax,stations.angleDeg,stations.computedNu,'-o','Color',[.15 .32 .49],'DisplayName',label);
    plot(ax,stations.angleDeg,stations.reportedNu,'d','Color',[.62 .30 .12],'MarkerFaceColor',[.62 .30 .12],'DisplayName','Short 表 III');
    xlabel(ax,'距迎风驻点的极角 / °');ylabel(ax,'局部 Nu（当地表温物性）');xlim(ax,[-4 180]);ax.YGrid='on';legend(ax,'Location','best');
    t=uitable(layout,'Data',round(stations{:,{'angleDeg','reportedNu','computedNu','deltaNu','relativeError'}}.*[1 1 1 1 100],3), ...
        'ColumnName',{'极角 / °','实验 Nu','计算 Nu','计算 − 实验','相对误差 / %'},'RowName',[]);
    note=sprintf('放热功率：网格 %.4f W，折合完整球面积 %.4f W。试验量热 %.4f W，探针梯度积分 %.4f W。\n七站 RMSE %.2f Nu。局部 Nu 来自线温梯度分布和整体量热，论文未给逐站误差条。', ...
        report.wholeOutW,report.equivalentSphereW,report.observedW,report.gradientIntegratedW,report.stationRMSE);
    uilabel(layout,'Text',note,'WordWrap','on');
    export=uibutton(layout,'Text','导出七站数据与误差 CSV','ButtonPushedFcn',@saveStations);
    tabs.SelectedTab=web.tab;
    view=struct('tabs',tabs,'web',web,'table',t,'stations',stations,'report',report,'mesh',mesh,'export',export,'exportTo',@writeStations);
    if wake
        % 分开绘制每物理步的速度和实际读取壁反力时的换热，避免插补热采样。
        wt=uitab(tabs,'Title','尾迹响应');box=uigridlayout(wt,[4 1]);box.RowHeight={44,'1x','1x',55};
        box.Padding=[16 10 16 10];g=report.balance;pair=report.wakeComparison;
        uilabel(box,'WordWrap','on','Text',sprintf('距扰动初场 %.6f ms   |   Δt = %.6f ms   |   第 %d 个物理步', ...
            1000*report.elapsedPhysicalTimeS,1000*report.dtSeconds,report.physicalStep));
        vel=uiaxes(box);plot(vel,[report.wake.timeDU],[report.wake.nonAxisymmetricRmsOverU],'-','Color',[.15 .32 .49]);
        vel.YGrid='on';xlabel(vel,'tU/D');ylabel(vel,'尾迹非轴对称速度 RMS / U');
        h=report.history;nu=[h.localNu];thermal=uiaxes(box);hold(thermal,'on');
        plot(thermal,[h.timeDU],nu(6,:)-report.baselineNu(6),'-o','DisplayName','119.8°','Color',[.15 .32 .49]);
        plot(thermal,[h.timeDU],nu(7,:)-report.baselineNu(7),'-s','DisplayName','166.7°','Color',[.61 .30 .12]);
        xlabel(thermal,'tU/D');ylabel(thermal,'瞬时 Nu − 稳态 Nu');legend(thermal,'Location','best');thermal.YGrid='on';
        uilabel(box,'WordWrap','on','Text',sprintf('放热 %.8f W，远场 %.8f W，空气储热率 %.4g W。两步长同刻热流最大差 %.5g W/m²。', ...
            g.wallOutW,g.farOutW,g.airStorageW,pair.heatMaxWm2));
        view.wakeTab=wt;view.wakeVelocity=vel;view.wakeHeat=thermal;
    end
    if isfield(report,'resolution')
        view.resolution=createResolutionTab(tabs,report.resolution,report.condition.diameterM);
    end
catch ex
    delete(tabs);rethrow(ex);
end

    function saveStations(~,~)
        [name,folder]=uiputfile('*.csv','导出风洞测站','Short65-stations.csv');if isequal(name,0),return;end
        writeStations(fullfile(folder,name));
    end
    function writeStations(path)
        % 模块二：同时导出归一化定义和对应的稳态源帧，可将单行还原到本试验。
        data=stations;data.test=repmat(65,height(data),1);data.sourceDOI=repmat(string(mesh.shortSphere.doi),height(data),1);
        data.frameSHA256=repmat(string(report.frameSHA256),height(data),1);
        data.method=repmat("三维稳态基流的子午线节点，按实验角度线性取样；Nu 使用当地表温导热系数",height(data),1);
        if wake
            data.method(:)="三维 BDF1 端点同态壁行反力，按实验角度线性取样，Nu 使用当地表温导热系数";
            data.physical_time_s=repmat(report.elapsedPhysicalTimeS,height(data),1);data.dt_s=repmat(report.dtSeconds,height(data),1);
            data.baselineNu=report.baselineNu(:);
        end
        writetable(data,path,'Encoding','UTF-8');
    end
end

function page=createResolutionTab(tabs,data,diameter)
% 模块三：按真实射线显示六个速度剖面，区分插值零点、原网夹区与壁面测站。
tab=uitab(tabs,'Title','尾迹空间分辨率');
grid=uigridlayout(tab,[6 1]);grid.RowHeight={44,32,'1x',58,80,36};
grid.Padding=[16 10 16 10];grid.RowSpacing=6;
axisData=data.axis;
uilabel(grid,'WordWrap','on','Text',sprintf( ...
    't = %.2f ms   |   从球尾量起的轴线回流长度 %.3f mm   |   零速点的相邻原节点位于 %.3f、%.3f mm', ...
    1000*data.timeS,1000*diameter*axisData.lengthD,1000*diameter*axisData.bracketD));
sections=data.sections;locations=[sections.xD];
choice=uidropdown(grid,'Items',cellstr(compose('子午截面   x/D = %g',locations(:))), ...
    'ItemsData',1:numel(sections),'Value',1,'ValueChangedFcn',@(~,~)showSection());
ax=uiaxes(grid);hold(ax,'on');ax.YGrid='on';
curve=plot(ax,nan,nan,'-o','Color',[.15 .32 .49],'MarkerSize',4,'DisplayName','原网格线交点');
marks=plot(ax,nan,nan,'s','Color',[.66 .32 .12],'MarkerFaceColor',[.66 .32 .12],'DisplayName','10% 与 90% 恢复位置');
xlabel(ax,'距轴线 r/D');ylabel(ax,'下游速度 −u_x / U');legend(ax,'Location','southeast');
note=uilabel(grid,'WordWrap','on');
stations=data.stations;
tableData=[[stations.angleDeg]' [stations.firstTangentialSpeedOverU]' [stations.computedNu]' [stations.measuredNu]'];
uitable(grid,'Data',round(tableData,5),'ColumnName',{'测站极角 / °','首层切向速度 / U','瞬时 Nu','实验 Nu'},'RowName',[]);
uilabel(grid,'WordWrap','on','Text','速度取自 +Z 测线相邻原节点的线性插值。首层切向以迎风点沿壁朝向尾极为正。测站热流可从七站页导出。');
showSection();
page=struct('tab',tab,'choice',choice,'curve',curve,'marks',marks,'note',note);
    function showSection()
        % 模块四：切换截面时同步更新网距和恢复带；图中保留原始剖面取样点。
        profile=sections(choice.Value);band=profile.band;
        curve.XData=profile.radiusD;curve.YData=profile.downstreamSpeedOverU;
        marks.XData=[band.r10D band.r90D];
        marks.YData=band.minSpeedOverU+[.1 .9]*(1-band.minSpeedOverU);
        xlim(ax,[0 1.4*band.r90D]);ylim(ax,[-.1 1.2]);
        note.Text=sprintf(['速度恢复带宽 %.4f D，带内 %d 条独立射线。射线横距 %.4f–%.4f D，' ...
            '供值边轴向跨度 %.4f–%.4f D。恢复幅度取速度谷底至自由流之差。'], ...
            band.widthD,band.innerRays,band.minRadialGapD,band.maxRadialGapD,band.minAxialSpanD,band.maxAxialSpanD);
    end
end
