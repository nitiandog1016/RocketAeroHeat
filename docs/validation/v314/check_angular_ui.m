function check_angular_ui(out)
% 模块一：从正式侧栏读取热敏珠冷态基线的归档及两套端面细网，逐个打开实际对照按钮。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);c=fig.UserData.controls.leadGrid;c.Value=3;c.ValueChangedFcn(c,[]);drawnow;
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
v=fig.UserData.webUI.evidenceView;d=v.getData();g=d.angularSpacingComparison;
press(fig.UserData.controls.outerHeat);u=v.getState();assert(strcmp(u.fieldKey,'outerDeltaHeat'));
press(fig.UserData.controls.tipSpacing);u=v.getState();assert(strcmp(u.fieldKey,'tipDeltaHeat'));
press(fig.UserData.controls.angularSpacing);u=v.getState();
assert(strcmp(u.fieldKey,'angularDeltaHeat') && u.selectedFaceID==g.finePeakFace && contains(fig.UserData.status.Text,'96 档'));
html=fig.UserData.webUI.html;t=tic;
while isempty(html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
rendered=html.UserData.rendered;assert(rendered.faces==10272 && rendered.visiblePixelSamples>0 && rendered.webglError==0);
r=d.comparisonColorRanges;
assert(isequal(r.angular96Heat,r.tipFineHeat,r.heatFlux) && isequal(r.angular96Gauge,r.tipFineGauge,r.pressureGauge) && ...
    isequal(r.angular96First,r.tipFineFirst,r.nearWallCellTemperature));
% 模块二：导出两网的周向截线与全部原面的三个差值，以便独立程序逐行核验。
for name={'angularDeltaHeat','angularDeltaPressure','angularDeltaFirst'}
    key=name{1};v.selectField(key);csv=v.exportTo(fullfile(out,[key '.csv']));assert(csv.rows==10272);
end
x=.000172;
for name={'tipFineHeat','angular96Heat','tipFineGauge','angular96Gauge'}
    key=name{1};v.selectField(key);v.setSectionX(x);t=tic;
    while toc(t)<40
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12);
    v.exportSectionTo(fullfile(out,[key '_section.csv']));
end
page=exportAerothermalWeb(s.result,fullfile(out,'workbench'));
proof=struct('page',page,'faces',rendered.faces,'visibleSamples',rendered.visiblePixelSamples,'webglError',rendered.webglError, ...
    'probeFace',u.selectedFaceID,'sectionXM',x,'baseFrame',g.baseFrame,'fineFrame',g.fineFrame);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);guardFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
