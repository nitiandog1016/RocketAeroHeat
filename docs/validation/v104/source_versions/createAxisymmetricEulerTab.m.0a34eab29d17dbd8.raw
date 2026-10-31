function ui = createAxisymmetricEulerTab(tabs,mesh,freestream)
%CREATEAXISYMMETRICEULERTAB 在统一框架中计算/展示内部 Euler 流场。
% 独立验证页：该阶段只显示无粘场，不改写已采用工程外流的气动热。
% 未收敛结果可查看；近壁候选显式标记为诊断，不作为已验证外缘。
fig=ancestor(tabs,'figure'); tab=uitab(tabs,'Title','Euler 外流验证');
layout=uigridlayout(tab,[5 1]);layout.RowHeight={35,35,'1x',50,68};
bar=uigridlayout(layout,[1 13]);bar.Padding=[0 0 0 0];
bar.ColumnWidth={55,70,55,70,60,55,90,70,65,85,145,'1x',125};
uilabel(bar,'Text','基础格数'); nx=uieditfield(bar,'numeric','Value',120,'Limits',[16 500],'RoundFractionalValues','on');
uilabel(bar,'Text','径向格数'); nr=uieditfield(bar,'numeric','Value',64,'Limits',[16 300],'RoundFractionalValues','on');
uilabel(bar,'Text','母线细分'); refine=uieditfield(bar,'numeric','Value',1,'Limits',[1 8], ...
    'RoundFractionalValues','on','ValueDisplayFormat','%.0f');
uilabel(bar,'Text','外边界 R / L'); far=uieditfield(bar,'numeric','Value',0.65,'Limits',[0.2 5]);
uilabel(bar,'Text','迭代上限');maxIter=uieditfield(bar,'numeric','Value',50000,'Limits',[1 200000], ...
    'RoundFractionalValues','on','ValueDisplayFormat','%.0f');
run=uibutton(bar,'Text','计算内部 Euler 场','ButtonPushedFcn',@calculate);
metric=uidropdown(bar,'Items',{'子午面压力 p / kPa','子午面马赫数','三维壁面压力 p / kPa','壁面压力曲线','收敛残差', ...
    '三维近壁切向速度 / m/s','三维近壁温度 / K','近壁速度及法向分量','近壁速度梯度 / 1/s','近壁总压恢复系数','壁压估计对比'}, ...
    'ValueChangedFcn',@plotResult);
saveButton=uibutton(bar,'Text','导出此页图','ButtonPushedFcn',@savePlot);
secondBar=uigridlayout(layout,[1 7]);secondBar.Padding=[0 0 0 0];secondBar.ColumnWidth={195,205,195,'1x',140,120,120};
scheme=uidropdown(secondBar,'Items',{'HLL + 总焓输运（试验）','原始 HLL（对照）'}, ...
    'ItemsData',{'HLLenthalpy','HLL'},'Value','HLLenthalpy');
gridType=uidropdown(secondBar,'Items',{'X 截面网格（原有）','法向网格：光滑钝鼻（试验）','鼻部法向过渡（试验）'}, ...
    'ItemsData',{'axial','normalOffset','noseBlend'},'Value','axial');
spatialOrder=uidropdown(secondBar,'Items',{'一阶常量重构（默认）','线性＋单边梯度（试验）', ...
    '线性＋镜像梯度（试验）','线性＋镜像＋总焓（试验）'},'ItemsData',[1 2 3 4],'Value',1);
guide=uilabel(secondBar,'Text','零攻角/侧滑。凹角可能不适用法向网格；Euler 尚未接入热流。','WordWrap','on');
reuse=uicheckbox(secondBar,'Text','沿用同网格初场','Value',true);
logPressure=uicheckbox(secondBar,'Text','压力对数色标','Value',true,'ValueChangedFcn',@plotResult);
fullField=uicheckbox(secondBar,'Text','查看整个远场','Value',false,'ValueChangedFcn',@plotResult);
panel=uipanel(layout,'BorderType','none');ax=uiaxes(panel,'Units','normalized','Position',[.02 .035 .95 .94]);
status=uilabel(layout,'Text','尚未计算。Euler 场会同时给出收敛与远场边界检查。','WordWrap','on');
description=uilabel(layout,'Text', ...
    ['算法：可选 HLL 质量/动量 + 上风总焓能量通量（试验），或原始 HLL 对照。' ...
    '一阶有限体积、滑移壁、局部伪时间、径向 η^1.8 加密。' ...
    '定比热完全气体；本页不是 NASA7 热力模型。激波由守恒方程捕捉，厚度受网格耗散影响。' ...
    '尾端须超声速流出；无底部尾迹、黏性或高温化学。该阶段尚不将场结果用于热流。'],'WordWrap','on');
ui=struct('tab',tab,'axes',ax,'runButton',run,'metric',metric,'axialCells',nx, ...
    'radialCells',nr,'farRadius',far,'status',status,'description',description,'guide',guide, ...
    'logPressure',logPressure,'fullField',fullField,'fluxScheme',scheme, ...
    'maxIterations',maxIter,'reuseInitialState',reuse,'gridType',gridType,'wallRefinement',refine, ...
    'spatialOrder',spatialOrder,'reconstructionMode',spatialOrder);
sol=[]; geometry=[];
if abs(freestream.alpha)>1e-10 || abs(freestream.beta)>1e-10
    run.Enable='off';status.Text='当前有攻角或侧滑角。轴对称求解不适用，请在飞行条件页设为零后重新计算。';
end

    function calculate(~,~)
        run.Enable='off';status.Text='正在构建贴体网格并迭代……';drawnow;
        cleanup=onCleanup(@restoreButton); %#ok<NASGU>
        try
            [contour,geometry]=extractAxisymmetricEulerContour(mesh);
            options=struct('axialCells',nx.Value,'radialCells',nr.Value, ...
                'wallRefinement',refine.Value, ...
                'outerRadius',far.Value*(contour(end,1)-contour(1,1)),'progressCallback',@progress, ...
                'fluxScheme',scheme.Value,'spatialOrder',min(spatialOrder.Value,2),'maxIterations',maxIter.Value);
            if spatialOrder.Value>=3,options.boundaryReconstruction='slipGhost';end
            if spatialOrder.Value==4,options.reconstructionVariables='enthalpy';end
            if strcmp(gridType.Value,'noseBlend')
                options.normalBlendLength=min(.15*(contour(end,1)-contour(1,1)),2*max(contour(:,2)));
                grid=buildNormalOffsetEulerGrid(contour,options);
            elseif strcmp(gridType.Value,'normalOffset')
                grid=buildNormalOffsetEulerGrid(contour,options);
            else
                grid=buildAxisymmetricEulerGrid(contour,options);
            end
            % 只沿用同一几何网格的状态；改了网格就重新使用自由来流。
            % 切换通量后也重新检验残差，不继承上次的收敛标记。
            if reuse.Value && ~isempty(sol) && isequal(sol.grid.X,grid.X) && isequal(sol.grid.R,grid.R)
                options.initialState=sol;
            end
            sol=solveAxisymmetricEuler(grid,freestream,options);
            sol.geometry=geometry;
            % 单独保存诊断候选，不写入 mesh.externalFlow 或现有热流链路。
            sol.edgeCandidate=extractAxisymmetricEulerEdge(sol);
            data=fig.UserData;data.eulerResult=sol;fig.UserData=data;
            plotResult([],[]);
            status.Text=sprintf('%s；总焓最大偏差 %.3g%%；实际 %d×%d 格，%d 步 / %.1f s。网格精度未验证，未接入热流。', ...
                sol.status,100*sol.diagnostics.maximumTotalEnthalpyDeviation, ...
                sol.grid.nx,sol.grid.nr,sol.iterations,sol.elapsedSeconds);
        catch problem
            if isgraphics(status),status.Text=['计算未完成：' problem.message];end
            rethrow(problem);
        end
    end
    function progress(iter,res)
        if ~isgraphics(fig),error('RocketAeroHeat:Euler:Closed','窗口已关闭，停止迭代。');end
        status.Text=sprintf('内部 Euler 迭代 %d，残差 %.3g；当前曲线/云图如有显示，仍为上次结果。',iter,res);
        drawnow limitrate nocallbacks;
    end
    function restoreButton()
        if isgraphics(run),run.Enable='on';end
    end
    function plotResult(~,~)
        if isempty(sol),return;end
        description.Text=[sol.method '。初场：' sol.initialization ...
            '。网格：' sol.gridMethod ...
            '。定比热完全气体；未混入 NASA7。总焓偏差小不代表压力/激波网格精度合格。' ...
            '尾端须超声速流出；无底部尾迹、黏性或高温化学。当前仍未将场结果用于热流。'];
        if contains(metric.Value,'近壁') || strcmp(metric.Value,'壁压估计对比')
            description.Text=[sol.edgeCandidate.method '。来源：' sol.method '；网格：' sol.gridMethod ...
                '。' sol.edgeCandidate.gradientMethod ...
                '。不补造驻点零速度；总压含数值误差，不能仅由总焓判断外缘精度。未接入热流。'];
        end
        colorbar(ax,'off');cla(ax,'reset');g=sol.grid;colormap(ax,turbo(256));
        switch metric.Value
            case {'子午面压力 p / kPa','子午面马赫数'}
                if strcmp(metric.Value,'子午面马赫数'),data=sol.Mach;else,data=sol.pressure/1000;end
                % 每个四边形单独着色，避免把单元中心量误当作节点量。
                vertices=[g.X(:),g.R(:)]; inds=reshape(1:numel(g.X),size(g.X));
                a=inds(1:end-1,1:end-1);b=inds(2:end,1:end-1);
                c=inds(2:end,2:end);d=inds(1:end-1,2:end);
                patch(ax,'Vertices',vertices,'Faces',[a(:),b(:),c(:),d(:)], ...
                    'FaceVertexCData',data(:),'FaceColor','flat','EdgeColor','none');
                hold(ax,'on');plot(ax,g.X(:,1),g.R(:,1),'k-','LineWidth',2);hold(ax,'off');
                view(ax,2);axis(ax,'equal');axis(ax,'tight');colorbar(ax);
                if ~fullField.Value
                    ylim(ax,[0,min(g.outerRadius,max(3*max(g.contour(:,2)),.35*g.length))]);
                end
                xlabel(ax,'X / m（鼻尖 +X）');ylabel(ax,'R / m');
            case '三维壁面压力 p / kPa'
                selected=find(g.bodyWall);xs=sol.wall.x(selected);rs=sol.wall.radius(selected);
                theta=linspace(0,2*pi,65);X=xs*ones(size(theta));
                Y=rs*cos(theta)+geometry.centerYZ(1);Z=rs*sin(theta)+geometry.centerYZ(2);
                surf(ax,X,Y,Z,sol.wall.pressure(selected)/1000*ones(size(theta)), ...
                    'EdgeColor','none','FaceColor','interp');
                axis(ax,'equal');axis(ax,'tight');view(ax,[25 18]);colorbar(ax);
                xlabel(ax,'X / m');ylabel(ax,'Y / m');zlabel(ax,'Z / m');
            case '壁面压力曲线'
                selected=g.bodyWall;
                plot(ax,sol.wall.x(selected),sol.wall.pressure(selected)/1000,'LineWidth',1.5);
                xlabel(ax,'X / m');ylabel(ax,'壁面压力 p / kPa');
            case '收敛残差'
                semilogy(ax,sol.residualHistory(:,1),max(sol.residualHistory(:,2),eps),'LineWidth',1.5);
                xlabel(ax,'局部伪时间迭代步（不是物理时间）');ylabel(ax,'最大归一化守恒残差');
            case {'三维近壁切向速度 / m/s','三维近壁温度 / K'}
                e=sol.edgeCandidate;theta=linspace(0,2*pi,65);
                if contains(metric.Value,'速度'),values=e.signedSpeed;else,values=e.temperature;end
                surf(ax,e.x*ones(size(theta)),e.radius*cos(theta)+geometry.centerYZ(1), ...
                    e.radius*sin(theta)+geometry.centerYZ(2),values*ones(size(theta)), ...
                    'EdgeColor','none','FaceColor','interp');
                axis(ax,'equal');axis(ax,'tight');view(ax,[25 18]);colorbar(ax);
                xlabel(ax,'X / m');ylabel(ax,'Y / m');zlabel(ax,'Z / m');
            case '近壁速度及法向分量'
                e=sol.edgeCandidate;
                normalLabel='原单元法向速度';
                if e.sample.fromReconstructedTrace,normalLabel='重构壁面迹法向速度';end
                plot(ax,e.s,e.signedSpeed,'LineWidth',1.5,'DisplayName','投影切向速度');
                hold(ax,'on');plot(ax,e.s,e.sample.normalVelocity,'--','LineWidth',1.3, ...
                    'DisplayName',normalLabel);hold(ax,'off');legend(ax,'Location','best');
                xlabel(ax,'从鼻端向尾部的母线弧长 s / m');ylabel(ax,'速度 / (m/s)');
            case '近壁速度梯度 / 1/s'
                e=sol.edgeCandidate;plot(ax,e.s,e.velocityGradient,'LineWidth',1.5);
                xlabel(ax,'从鼻端向尾部的母线弧长 s / m');ylabel(ax,'dU_t / ds / (1/s)');
            case '近壁总压恢复系数'
                e=sol.edgeCandidate;plot(ax,e.s,e.diagnostics.totalPressureRatio,'LineWidth',1.5);
                xlabel(ax,'从鼻端向尾部的母线弧长 s / m');ylabel(ax,'局部 p_0 / 来流 p_{0∞}');
            case '壁压估计对比'
                e=sol.edgeCandidate;
                pressureLabel='第一层单元压力';
                if e.sample.fromReconstructedTrace,pressureLabel='重构壁面迹压力';end
                plot(ax,e.s,e.pressure/1000,'LineWidth',1.5,'DisplayName',pressureLabel);
                hold(ax,'on');plot(ax,e.s,e.wallFluxPressure/1000,'--','LineWidth',1.3, ...
                    'DisplayName','反射壁面数值通量压力');hold(ax,'off');legend(ax,'Location','best');
                xlabel(ax,'从鼻端向尾部的母线弧长 s / m');ylabel(ax,'压力 / kPa');
        end
        grid(ax,'on');colormap(ax,turbo(256));
        if contains(metric.Value,'p / kPa') && ~strcmp(metric.Value,'壁面压力曲线') && logPressure.Value
            ax.ColorScale='log';
        end
        if sol.available,tag='已收敛，边界检查通过；网格精度尚未验证';else,tag='未通过收敛或边界检查，仅供诊断';end
        title(ax,{['内部 Euler · ' metric.Value],tag},'Interpreter','none');
    end
    function savePlot(~,~)
        [name,path]=uiputfile('*.png','导出 Euler 验证图');
        if isequal(name,0),return;end
        exportgraphics(ax,fullfile(path,name),'Resolution',180);
    end
end
