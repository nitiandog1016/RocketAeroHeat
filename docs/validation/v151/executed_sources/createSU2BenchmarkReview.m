function view=createSU2BenchmarkReview(host,review)
%CREATESU2BENCHMARKREVIEW 在现有工作台呈现公开算例的实际终态和定量历史。
% 模块一：清楚区分源步号、最后记录步号和新检查点，压力/热流资格单列。
% 模块二：科研图保留原残差门槛与原生热量单位，导出原审核及图像。
grid=uigridlayout(host,[4 1]);grid.RowHeight={60,82,'1x',36};grid.Padding=[8 8 8 8];
if review.runSucceeded,status='原生任务正常结束';else,status='原生任务失败';end
if review.hasSavedNewField,frame=sprintf('新保存场局部 %d',review.savedLocalIteration);
else,frame='本次尚无新保存场';end
if isempty(review.sourceLocalIteration),source='独立相容初场';
else,source=sprintf('源场局部 %d',review.sourceLocalIteration);end
heading=uilabel(grid,'FontSize',14,'FontWeight','bold','WordWrap','on', ...
    'Text',sprintf('NASA ACC-SBLI · M %.4g · %s\n%s · 退出码 %d', ...
    review.conditions.mach,review.model,status,review.exitCode));heading.Layout.Row=1;
if isfield(review,'acceptance')
    passed=sum(arrayfun(@(g)isequal(g.passed,true),review.acceptance.gates));
    heading.Text=sprintf('%s\n预设验收 %d / %d 通过',heading.Text,passed,numel(review.acceptance.gates));
end
summary=uilabel(grid,'WordWrap','on','FontSize',11,'VerticalAlignment','top', ...
    'Text',sprintf(['%s · 最后历史局部 %d · %s\n' ...
    '给定壁温 %.4g K · 耗时 %.4g s · 新增合格面积 0\n' ...
    '原生 HF 使用 W/m 平面积分。整周旋转热功率采用 W，随具名流场单独审核。'], ...
    source,review.lastRecordedIteration,frame, ...
    review.conditions.wallTemperatureK,review.elapsedSeconds));summary.Layout.Row=2;
tabs=uitabgroup(grid);tabs.Layout.Row=3;
historyTab=uitab(tabs,'Title','收敛历史');historyGrid=uigridlayout(historyTab,[1 1]);
picture=uiimage(historyGrid,'ImageSource',review.imageFile,'ScaleMethod','fit');
comparisonTab=[];gateTable=[];
if ~isempty(review.currentFieldAssessment)
    % 模块三：用原弧长测站和有量纲数值展示公开数据对照。
    comparisonTab=uitab(tabs,'Title','压力与热流对照');comparisonGrid=uigridlayout(comparisonTab,[2 1]);
    fields={'pressure','heat'};labels={'压力 (Pa)','入壁热流 (kW/m²)'};factors=[1 .001];
    for i=1:2
        ax=uiaxes(comparisonGrid);data=review.currentFieldAssessment.stations.(fields{i});
        plot(ax,data.S,data.predicted*factors(i),'-','Color',[.17 .39 .57],'LineWidth',1.4);hold(ax,'on');
        plot(ax,data.S,data.measured*factors(i),'o','Color',[.25 .25 .25],'MarkerSize',4);
        xlabel(ax,'壁面弧长 S (m)');ylabel(ax,labels{i});ax.XGrid='on';ax.YGrid='on';
        legend(ax,{'原生保存场','公开实验'},'Location','best');
        title(ax,sprintf('SST-Vm · 局部 %d · 研究结果',review.savedLocalIteration),'FontWeight','normal');
    end
end
if isfield(review,'acceptance')
    % 模块四：每项门槛、实值和待验原因同时显示，资格仍由完整验证决定。
    gateTab=uitab(tabs,'Title','验收指标');gateGrid=uigridlayout(gateTab,[1 1]);
    gates=review.acceptance.gates;rows=cell(numel(gates),5);
    for i=1:numel(gates)
        g=gates(i);
        if isempty(g.value),valueText='';else,valueText=sprintf('%.6g',g.value);end
        rows(i,1:3)={g.name,valueText,sprintf('%s %.5g',g.comparison,g.limit)};
        if isempty(g.passed),rows{i,4}='待验';rows{i,5}=g.reason;
        elseif g.passed,rows{i,4}='通过';rows{i,5}='';
        else,rows{i,4}='未通过';rows{i,5}='';end
    end
    gateTable=uitable(gateGrid,'Data',rows,'ColumnName',{'验收项','当前值','门槛','结果','说明'}, ...
        'ColumnWidth',{245,120,110,70,'auto'},'RowName',{});
end
bar=uigridlayout(grid,[1 2]);bar.Layout.Row=4;bar.ColumnWidth={'1x',160};bar.Padding=[0 0 0 0];
uilabel(bar,'Text','公开验证记录 · 当前箭体输入由侧栏管理','FontSize',11);
button=uibutton(bar,'Text','导出审核及科研图','ButtonPushedFcn',@exportDialog);
view=struct('heading',heading,'summary',summary,'picture',picture,'export',button, ...
    'review',review,'exportTo',@exportTo,'tabs',tabs,'comparisonTab',comparisonTab,'gateTable',gateTable);
    function exportDialog(~,~)
        [name,folder]=uiputfile('*.json','导出公开基准审核');if isequal(name,0),return;end
        exportTo(fullfile(folder,name));
    end
    function exportTo(file)
        copyfile(review.reviewFile,file);copyfile(review.imageFile,[file '.png']);
        if ~isempty(review.currentFieldAssessment)
            field=review.currentFieldAssessment.stations;names={'pressure','heat'};units={'Pa','W/m2'};
            parts=cell(1,2);
            for j=1:2
                data=field.(names{j});n=numel(data.S);
                parts{j}=table(repmat(names(j),n,1),data.S,data.predicted,data.measured,repmat(units(j),n,1), ...
                    repmat({review.model},n,1),repmat(review.savedLocalIteration,n,1), ...
                    repmat(review.validationPassed,n,1),'VariableNames', ...
                    {'Quantity','S_m','Predicted','Measured','Unit','Method','SavedIteration','ValidationPassed'});
            end
            output=vertcat(parts{:});
            writetable(output,[file '.stations.csv']);
        end
    end
end
