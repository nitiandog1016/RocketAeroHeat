function testPressureProfiles(outputFolder)
% 检查柱段曲线确实反映攻角变化，并验证局部色标只改变显示范围。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
evalc('raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);');
[raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
evalc('mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
results=zeros(3,4);
for k=1:3
    angle=[0 20 30]; flight=cfg.flight; flight.alpha=angle(k);
    m=modifiedNewtonianPressure(mesh,calcFreestream(flight),cfg);
    profile=computePressureProfiles(m);
    assert(profile.available && nnz(profile.cylinderMask)>1000);
    assert(nnz(isfinite(profile.angularCp))>20);
    assert(nnz(isfinite(profile.windwardCp))>50);
    assert(all(profile.angularCount(isfinite(profile.angularCp))>0));
    assert(all(profile.windwardCount(isfinite(profile.windwardCp))>0));
    results(k,:)=[angle(k),max(profile.windwardCp),max(profile.angularCp),max(m.pressure.Cp)];
    if angle(k)==20
        cfg.pressureGUI.visible='off';
        fig=pressureGUI(m,cfg); cleaner=onCleanup(@()delete(fig)); %#ok<NASGU>
        ui=fig.UserData;
        assert(ui.colorScale.Value==2 && ui.profiles.available);
        assert(abs(ui.axes.CLim(2)-profile.cylinderCpRange(2))<1e-12);
        cylinderLimits=ui.axes.CLim;
        ui.colorScale.Value=1;
        ui.colorScale.ValueChangedFcn(ui.colorScale,[]);
        assert(ui.axes.CLim(2)>cylinderLimits(2));
        ui.colorScale.Value=2;
        ui.colorScale.ValueChangedFcn(ui.colorScale,[]);
        assert(abs(ui.axes.CLim(2)-cylinderLimits(2))<1e-12);
        ui.tabs.SelectedTab=ui.curveTab;
        assert(numel(findall(ui.angularAxes,'Type','line'))==1);
        assert(isgraphics(ui.sectionSurface));
        middle=ui.selectedProfile;
        % 滑块从中部移到鼻部，周向曲线应来自不同的实际环带。
        ui.sectionSlider.Value=max(m.V(:,1));
        ui.sectionSlider.ValueChangedFcn(ui.sectionSlider,[]);
        nose=fig.UserData.selectedProfile;
        assert(nose.stationX>middle.stationX && nose.faceCount>0);
        assert(~isequaln(nose.Cp,middle.Cp));
        if nargin>0
            drawnow; exportapp(fig,fullfile(outputFolder,'nose_section.png'));
        end
        % 模拟三维面元点击：交点 X 决定取样环带，而非点击的周向角。
        ui.sectionSurface.ButtonDownFcn(ui.sectionSurface, ...
            struct('IntersectionPoint',[middle.stationX 0 0]));
        clicked=fig.UserData.selectedProfile;
        assert(abs(clicked.stationX-middle.stationX)<1e-10);
        % 同一站位可切换为绝对压力，数据按该站位面元统计。
        ui.sectionMetric.Value=2;
        ui.sectionMetric.ValueChangedFcn(ui.sectionMetric,[]);
        assert(isequaln(ui.sectionLine.YData(:),clicked.staticPressurePa(:)/1000));
        ui.sectionMetric.Value=1;
        ui.sectionMetric.ValueChangedFcn(ui.sectionMetric,[]);
        assert(isequaln(ui.sectionLine.YData(:),clicked.Cp(:)));
        if nargin>0
            drawnow;
            exportapp(fig,fullfile(outputFolder,'pressure_curves.png'));
            ui.tabs.SelectedTab=ui.mapTab; drawnow;
            exportapp(fig,fullfile(outputFolder,'cylinder_color_scale.png'));
        end
        delete(fig); clear cleaner;
    end
end
assert(results(1,2)<1e-3);
assert(results(2,2)>0.1 && results(3,2)>results(2,2));
assert(results(2,4)>results(2,3)); % 鼻部峰值高于柱段曲线范围。
fprintf('Pressure profiles passed: alpha 0/20/30 windward Cp %.4g / %.4g / %.4g.\n', ...
    results(:,2));
end
