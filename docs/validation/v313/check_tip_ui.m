function check_tip_ui(out)
% 模块一：用正式侧栏读取 20D 侧吹，点选细网热峰并比对已有外域按钮。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);
c=fig.UserData.controls.leadGrid;c.Value=3;c.ValueChangedFcn(c,[]);drawnow;
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
v=fig.UserData.webUI.evidenceView;d=v.getData();g=d.tipSpacingComparison;
press(fig.UserData.controls.outerHeat);u=v.getState();assert(strcmp(u.fieldKey,'outerDeltaHeat'));
press(fig.UserData.controls.tipSpacing);u=v.getState();
assert(strcmp(u.fieldKey,'tipDeltaHeat') && u.selectedFaceID==g.finePeakFace);
assert(g.fineSecondFace==10267 && g.finePeakFace==209 && contains(fig.UserData.status.Text,'第二强面'));
gl=fig.UserData.webUI.html;t=tic;
while isempty(gl.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
draw=gl.UserData.rendered;assert(draw.faces==10272 && draw.webglError==0 && draw.visiblePixelSamples>0);
limits=d.comparisonColorRanges;
assert(isequal(limits.tipFineHeat,limits.heatFlux,limits.outerHeat30) && ...
    isequal(limits.tipFineGauge,limits.pressureGauge,limits.outerGauge30) && ...
    isequal(limits.tipMatched,limits.nearWallCellTemperature,limits.tipFineFirst,limits.outerFirst30));
% 模块二：导出同面热流差与同位置温差，截取细化引线材料段以核对角度与面号。
for name={'tipDeltaHeat','tipDeltaMatched'}
    key=name{1};v.selectField(key);rec=v.exportTo(fullfile(out,[key '.csv']));assert(rec.rows==10272);
end
x=.000864;v.selectField('tipFineHeat');v.setSectionX(x);t=tic;
while toc(t)<40
    drawnow;cut=v.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'tipFineHeat') && abs(cut.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'tipFineHeat') && abs(cut.x-x)<1e-12);
v.exportSectionTo(fullfile(out,'tip_section.csv'));
% 模块三：保存实际离线工作台及可复核的探针、两套帧与渲染回执。
page=exportAerothermalWeb(s.result,fullfile(out,'workbench'));
proof=struct('page',page,'faces',draw.faces,'visibleSamples',draw.visiblePixelSamples,'webglError',draw.webglError, ...
    'probeFace',u.selectedFaceID,'sectionXM',x,'baseFrame',g.baseFrame,'fineFrame',g.fineFrame);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);guardFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
