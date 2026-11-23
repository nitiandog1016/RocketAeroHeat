function check_balance_ui(out)
% 模块一：通过工作台实际选中双支管试件，在同一窗体读取热量账。
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supports';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
review=fig.UserData.shortSphereReview;p=review.fluidGrid.balance;review.tabs.SelectedTab=p.tab;drawnow;
assert(contains(p.note.Text,'0.6400%') && contains(p.note.Text,'0 m²'));
d=p.data;assert(isequal(p.ledger.Data(:,2),[d.early.wallW;d.late.wallW]));
assert(isequal(p.ledger.Data(:,4),1000*[d.early.physicalGapW;d.late.physicalGapW]));
% 模块二：逐条绑定原始读数及计算更新号，切换流场与温差曲线。
for mode={'flow','heat'}
    name=mode{1};p.mode.Value=name;p.mode.ValueChangedFcn(p.mode,[]);drawnow;
    if strcmp(name,'flow'),cols=2:5;else,cols=6:8;end
    lines=flipud(findobj(p.ax,'Type','line'));assert(numel(lines)==numel(cols));
    for k=1:numel(lines)
        assert(isequal(lines(k).XData(:),d.trend(:,1)) && isequal(lines(k).YData(:),d.trend(:,cols(k))));
    end
end
review.tabs.SelectedTab=review.fluidGrid.tab;drawnow;assert(isgraphics(review.fluidGrid.line));
review.tabs.SelectedTab=review.geometryTab;drawnow;assert(height(review.geometry.stations)==7);
proof=struct('mesh',d.meshSHA256,'end',d.early.sourceFrame,'trendRows',size(d.trend,1), ...
    'table',p.ledger.Data,'note',p.note.Text);
f=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(f>0);done=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
