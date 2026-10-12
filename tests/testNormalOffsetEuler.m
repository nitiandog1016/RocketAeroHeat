function testNormalOffsetEuler(outputFolder)
%TESTNORMALOFFSETEULER 半球鼻联合几何/网格收敛及新网格统一页显示。
% 驻点压力用法向激波后等熵减速的独立极限；不把有限距离速度当驻点值。
% 同时改变母线分辨率，属于几何/网格联合研究，不能冒充固定几何加密。
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);ref=solveNormalShock(fs);
errors=zeros(2,1);intercepts=errors;firstSpeed=errors;firstS=errors;
for k=1:2
    count=64*k;theta=linspace(pi/2,0,count+1)';
    c=[-1,.1;.1*(cos(theta)-1),.1*sin(theta)];
    if k==1
        % 使用解析截面圆环进入真实 GUI 回调，检验选择器与计算/绘图一致。
        mesh=ringMesh(c);fig=uifigure('Visible','off','Position',[30 30 1580 940]);
        cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
        fig.UserData=struct;tabs=uitabgroup(fig,'Position',[0 0 1580 940]);
        ui=createAxisymmetricEulerTab(tabs,mesh,fs);ui.gridType.Value='normalOffset';
        ui.radialCells.Value=96;ui.maxIterations.Value=100000;
        ui.runButton.ButtonPushedFcn(ui.runButton,[]);sol=fig.UserData.eulerResult;
        assert(strcmp(sol.boundary.topology,'axisAtNose'));
        assert(contains(ui.description.Text,'法向偏移'));
        % 编辑下拉框但不重新计算时，算法说明仍应描述实际计算结果。
        ui.gridType.Value='axial';ui.metric.Value=ui.metric.Items{6};
        ui.metric.ValueChangedFcn(ui.metric,[]);
        assert(contains(ui.description.Text,'法向偏移'));
        plotted=findobj(ui.axes,'Type','surface');
        assert(isequal(plotted.CData(:,1),sol.edgeCandidate.signedSpeed));
        ui.gridType.Value='normalOffset'; % 预览恢复为实际采用的选项，避免读图歧义。
        if nargin>0
            if ~isfolder(outputFolder),mkdir(outputFolder);end
            exportapp(fig,fullfile(outputFolder,'normal_grid_velocity.png'));
            ui.metric.Value=ui.metric.Items{1};ui.metric.ValueChangedFcn(ui.metric,[]);
            exportapp(fig,fullfile(outputFolder,'normal_grid_pressure.png'));
        end
        clear cleanup;
    else
        grid=buildNormalOffsetEulerGrid(c,struct('axialCells',120,'radialCells',96));
        sol=solveAxisymmetricEuler(grid,fs,struct('maxIterations',100000));
    end
    edge=extractAxisymmetricEulerEdge(sol);fit=polyfit(edge.s(1:6),edge.signedSpeed(1:6),2);
    errors(k)=abs(sol.wall.pressure(end)/ref.stagnationPressure-1);
    intercepts(k)=abs(fit(end));firstSpeed(k)=edge.signedSpeed(1);firstS(k)=edge.s(1);
    assert(sol.available && sol.conservation.identityError<1e-10);
    assert(sol.conservation.wallMassLeak<1e-12 && sol.conservation.wallEnergyLeak<1e-12);
    assert(sol.diagnostics.maximumTotalEnthalpyDeviation<1e-5);
    assert(isnan(sol.boundary.upstreamStateDisturbance) && ~sol.edgeStateAvailable);
    fprintf('Normal-grid sphere %d: s1=%g U1=%g nose intercept=%g pressure error=%g\n', ...
        count,firstS(k),firstSpeed(k),intercepts(k),errors(k));
end
assert(firstS(2)<.51*firstS(1));
assert(firstSpeed(2)<.6*firstSpeed(1) && intercepts(2)<.6*intercepts(1));
assert(errors(2)<errors(1) && errors(2)<.04);
fprintf('Normal-offset Euler passed: nose velocity convergence, independent pressure and GUI; heat coupling remains disabled.\n');
end
function mesh=ringMesh(c)
% 解析回转体截面夹具；使用正式几何提取接口，不向 GUI 注入已算好的答案。
theta=(0:31)'*2*pi/32;mesh.V=zeros(0,3);ids=cell(size(c,1),1);
for k=1:size(c,1)
    start=size(mesh.V,1);
    if c(k,2)==0,points=[c(k,1),0,0];
    else,points=[repmat(c(k,1),32,1),c(k,2)*cos(theta),c(k,2)*sin(theta)];end
    mesh.V=[mesh.V;points]; %#ok<AGROW>
    ids{k}=(start+1:size(mesh.V,1))';
end
mesh.param.sectionX=c(:,1);mesh.param.sectionNodeIDs=ids;
end
