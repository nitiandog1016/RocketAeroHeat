function page=addShortSupportBalanceTab(tabs,folder,reference)
% 模块一：由同一试件、混合网和冷态绑定有界标量计算与两份端点收支。
file=fullfile(folder,'supportedCapacity65.json');
assert(strcmp(cm011FileSHA256(file),'9560d6951503cbe91fb06902bbbe9692239694832d489acf87f7aa476a143599'));
d=jsondecode(fileread(file));
assert(strcmp(d.caseSHA256,reference.caseSHA256) && strcmp(d.meshSHA256,reference.meshSHA256) && ...
    strcmp(d.coldSHA256,reference.coldSHA256) && strcmp(d.stages(end).frame,d.early.sourceFrame));
tab=uitab(tabs,'Title','收敛与热量账');
layout=uigridlayout(tab,[5 1]);layout.RowHeight={30,32,'1x',82,57};layout.Padding=[14 8 14 8];
uilabel(layout,'Text',sprintf('Short Test 65 · 节点余量通量修正 · %d 次更新 · %d 个节点越界', ...
    d.updates,d.boundViolations),'FontWeight','bold');
mode=uidropdown(layout,'Items',{'压力与三分量速度','弱焓与温差标量'}, ...
    'ItemsData',{'flow','heat'},'Value','heat');
ax=uiaxes(layout);xlabel(ax,'稳态更新次数');ylabel(ax,'log₁₀ RMS 残差');grid(ax,'on');
% 模块二：从两份真实快照显示给定壁温下的壁放热与标量远场输出，正号表示向外。
values=[d.early d.late];
rows=[[d.updates;d.updates+20] [values.wallW]' [values.scalarFarW]' 1000*[values.physicalGapW]' 100*[values.physicalGapFraction]'];
ledger=uitable(layout,'Data',rows,'RowName',{}, ...
    'ColumnName',{'更新次数','壁放热 (W)','远场输出 (W)','差额 (mW)','差额 (%)'});
note=uilabel(layout,'WordWrap','on','Text',sprintf(['收支验收上限 0.5%%，本次计算差额 %.4f%%。' ...
    '按流动质量通量复算远场为 %.7f W。合格壁压与热流面积 %g m²。'], ...
    100*d.early.physicalGapFraction,d.early.flowFluxComparisonW,d.acceptedAreaM2));
note.Tooltip='差额 =（标量远场输出 − 壁放热）/ 壁放热。流动质量通量的复算值单列供定位原因。';
mode.ValueChangedFcn=@(~,~)draw();draw();
page=struct('tab',tab,'mode',mode,'ax',ax,'ledger',ledger,'note',note,'data',d);
    function draw()
        % 模块三：物理变量与未限幅标量残差分别画线，同一工况共用 −8 门槛。
        cla(ax);hold(ax,'on');
        if strcmp(mode.Value,'flow')
            cols=2:5;labels={'压力','X 速度','Y 速度','Z 速度'};
        else
            cols=6:8;labels={'弱焓','实际温差方程','原未限幅温差方程'};
        end
        lines=gobjects(numel(cols),1);
        for k=1:numel(cols)
            lines(k)=plot(ax,d.trend(:,1),d.trend(:,cols(k)),'LineWidth',1.4,'DisplayName',labels{k});
        end
        yline(ax,-8,'--','HandleVisibility','off','Color',[.32 .32 .32]);
        xlim(ax,[0 d.updates]);legend(ax,lines,labels,'Location','southwest');
    end
end
