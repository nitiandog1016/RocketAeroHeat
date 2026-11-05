function view=createSU2BenchmarkReview(host,review)
%CREATESU2BENCHMARKREVIEW 在现有工作台呈现公开算例的实际终态和定量历史。
% 模块一：清楚区分源步号、最后记录步号和新检查点，压力/热流资格单列。
% 模块二：科研图保留原残差门槛与原生热量单位，导出原审核及图像。
grid=uigridlayout(host,[4 1]);grid.RowHeight={60,74,'1x',36};grid.Padding=[8 8 8 8];
if review.runSucceeded,status='原生任务正常结束';else,status='原生任务失败';end
if review.hasSavedNewField,frame=sprintf('新保存场局部 %d',review.savedLocalIteration);
else,frame='本次尚无新保存场';end
heading=uilabel(grid,'FontSize',14,'FontWeight','bold','WordWrap','on', ...
    'Text',sprintf('NASA ACC-SBLI · M %.4g · %s\n%s · 退出码 %d', ...
    review.conditions.mach,review.model,status,review.exitCode));heading.Layout.Row=1;
summary=uilabel(grid,'WordWrap','on','FontSize',11,'VerticalAlignment','top', ...
    'Text',sprintf(['源场局部 %d · 最后历史局部 %d · %s\n' ...
    '给定壁温 %.4g K · 耗时 %.4g s · 新增合格面积 0\n' ...
    '原生 HF 使用 W/m 平面积分。整周旋转热功率采用 W，随具名流场单独审核。'], ...
    review.sourceLocalIteration,review.lastRecordedIteration,frame, ...
    review.conditions.wallTemperatureK,review.elapsedSeconds));summary.Layout.Row=2;
picture=uiimage(grid,'ImageSource',review.imageFile,'ScaleMethod','fit');picture.Layout.Row=3;
bar=uigridlayout(grid,[1 2]);bar.Layout.Row=4;bar.ColumnWidth={'1x',160};bar.Padding=[0 0 0 0];
uilabel(bar,'Text','公开验证记录 · 当前箭体输入由侧栏管理','FontSize',11);
button=uibutton(bar,'Text','导出审核及科研图','ButtonPushedFcn',@exportDialog);
view=struct('heading',heading,'summary',summary,'picture',picture,'export',button, ...
    'review',review,'exportTo',@exportTo);
    function exportDialog(~,~)
        [name,folder]=uiputfile('*.json','导出公开基准审核');if isequal(name,0),return;end
        exportTo(fullfile(folder,name));
    end
    function exportTo(file)
        copyfile(review.reviewFile,file);copyfile(review.imageFile,[file '.png']);
    end
end
