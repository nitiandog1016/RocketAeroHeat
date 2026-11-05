function view=createSU2IdealAirReview(host,folder)
%CREATESU2IDEALAIRREVIEW 在统一工作台查看原生变比热小例及真实方法来源。
% 模块一。通过正常收敛与热状态读取入口，保持给定壁温和流体温度身份。
% 模块二。显示原壁节点p/q曲线和定量表，导出同一组原数据与方法。
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
[~,~,summary]=readSU2IdealAirStateAudit(folder,receipt);
index=jsondecode(fileread(fullfile(folder,'ideal_air_review.json')));
wall=readtable(fullfile(folder,index.wallPath));
runtime=describeSU2Runtime(receipt);
result=struct('receipt',receipt,'reconstructedEngineHeatPowerRaw',summary.nativeHeatPower);
convergence=auditSU2ResultConvergence(folder,result);
assert(convergence.checkpointIdentityConfirmed,convergence.reason);
assert(convergence.checkpointIteration==summary.checkpointLocalIteration);
assert(summary.dimension==2,'此定量曲线入口用于二维接口验证。三维结果使用原面工作台。');
grid=uigridlayout(host,[4 1]);grid.RowHeight={68,62,'1x',36};grid.Padding=[8 8 8 8];
heading=uilabel(grid,'FontSize',14,'FontWeight','bold','WordWrap','on', ...
    'Text',sprintf('%s\n局部 %d · %d 次更新 · 接口通过 · 物理资格待验', ...
    runtime.physicsCaption,summary.checkpointLocalIteration,summary.successfulUpdates));
heading.Layout.Row=1;
method=uilabel(grid,'FontSize',11,'WordWrap','on','Text',runtime.method);method.Layout.Row=2;
tabs=uitabgroup(grid);tabs.Layout.Row=3;
plotTab=uitab(tabs,'Title','壁面分布');plotGrid=uigridlayout(plotTab,[2 1]);
[x,order]=sort(wall.x_m);axesList=gobjects(2,1);
fields={'p_Pa','q_W_m2'};labels={'压力 (Pa)','入壁热流 (MW/m²)'};scale=[1 1e-6];
for j=1:2
    ax=uiaxes(plotGrid);axesList(j)=ax;
    plot(ax,x,wall.(fields{j})(order)*scale(j),'Color',[.17 .39 .57],'LineWidth',1.4);
    xlabel(ax,'X (m)');ylabel(ax,labels{j});ax.XGrid='on';ax.YGrid='on';
    title(ax,'有限更新研究场','FontWeight','normal');
end
tableTab=uitab(tabs,'Title','节点数值与温度');tableGrid=uigridlayout(tableTab,[2 1]);tableGrid.RowHeight={55,'1x'};
uilabel(tableGrid,'WordWrap','on','Text',sprintf( ...
    '给定壁温 %.6g K · 壁面流体节点 %.6g–%.6g K · 求解固体壁温未知\n原生二维热负荷 %.7g W/m · 原−8残差门槛未过 · 新增合格面积0', ...
    summary.givenWallTemperatureK,summary.wallFluidTemperatureRangeK,summary.nativeHeatPower));
tableUI=uitable(tableGrid,'Data',wall,'RowName',{},'ColumnName', ...
    {'原节点','X / m','Y / m','Z / m','p / Pa','q / W·m⁻²','流体T / K','给定Tw / K','独立y+'});
bar=uigridlayout(grid,[1 2]);bar.Layout.Row=4;bar.ColumnWidth={'1x',145};bar.Padding=[0 0 0 0];
uilabel(bar,'Text','二维接口验证 · 60–2000 K · OTH能熵参考','FontSize',11);
button=uibutton(bar,'Text','导出节点CSV','ButtonPushedFcn',@exportDialog);
view=struct('summary',summary,'runtime',runtime,'convergence',convergence,'wall',wall, ...
    'heading',heading,'method',method,'tabs',tabs,'tableTab',tableTab,'table',tableUI,'axes',axesList, ...
    'export',button,'exportTo',@exportTo);
    function exportDialog(~,~)
        [name,path]=uiputfile('*.csv','导出变比热研究节点');if isequal(name,0),return;end
        exportTo(fullfile(path,name));
    end
    function exportTo(file)
        data=wall;data.solvedWallTemperature_K=nan(height(data),1);
        data.Method=repmat(string(runtime.method),height(data),1);
        data.SavedLocalIteration=repmat(summary.checkpointLocalIteration,height(data),1);
        data.PhysicsQualified=false(height(data),1);writetable(data,file);
        meta=struct('runtime',runtime,'summary',summary,'convergence',convergence,'sourceFolder',folder);
        fid=fopen([file '.json'],'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(meta));
    end
end
