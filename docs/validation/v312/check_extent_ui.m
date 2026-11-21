function check_extent_ui(out)
% 模块一：从已发布的 20D 侧吹入口打开同工况的 30D 对照，定位最大热流差所在的原面。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,3);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
g=s.result.thermistorOuterDomain;assert(strcmp(g.mode,'sidewindGivenWall') && g.wall20K==299.15 && g.wall30K==299.15);
web=fig.UserData.webUI;v=web.evidenceView;d=v.getData();press(fig.UserData.controls.outerHeat);
[~,id]=max(abs(g.deltaHeatWm2));state=v.getState();
assert(strcmp(state.fieldKey,'outerDeltaHeat') && state.selectedFaceID==id);
t=tic;while isempty(web.html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
gl=web.html.UserData.rendered;assert(gl.faces==10272 && gl.webglError==0 && gl.visiblePixelSamples>0);
limits=d.comparisonColorRanges;assert(isequal(limits.outerHeat30,limits.heatFlux) && ...
    isequal(limits.pressureGauge,limits.outerGauge30) && isequal(limits.outerFirst30,limits.nearWallCellTemperature));
csv=v.exportTo(fullfile(out,'deltaHeat.csv'));assert(csv.rows==10272);
% 模块二：对同一平面分别导出两域的表压和热流，逐段与原面核对。
x=.000172;
for name={'pressureGauge','heatFlux','outerGauge30','outerHeat30'}
    key=name{1};v.selectField(key);v.setSectionX(x);t=tic;
    while toc(t)<40
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12);
    path=fullfile(out,[key '_section.csv']);v.exportSectionTo(path);rows=readtable(path,'Encoding','UTF-8');
    f=d.fields(strcmp({d.fields.key},key));assert(height(rows)==96 && max(abs(rows.value-f.values(rows.face_id)))<1e-9);
end
page=exportAerothermalWeb(s.result,fullfile(out,'workbench'));
proof=struct('page',page,'faces',gl.faces,'visibleSamples',gl.visiblePixelSamples,'webglError',gl.webglError, ...
    'probeFace',id,'sectionXM',x,'source20',g.baseFrame,'source30',g.expandedFrame, ...
    'outward20W',g.outward20W,'outward30W',g.outward30W,'peak20Face',g.peak20Face,'peak30Face',g.peak30Face);
% 模块三：相同试件切回轴向第二档，确认它仍然使用自己的 30D 参照包。
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.leadGrid,1);s=fig.UserData.workflowState;
assert(isempty(s.mesh) && isempty(s.result) && isempty(s.cfg.thermistorBoundary));
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);
press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;assert(~isempty(s.result),fig.UserData.status.Text);
axial=s.result.thermistorOuterDomain;
assert(strcmp(axial.mode,'givenWall') && axial.wall30K==305.15 && ~strcmp(axial.expandedFrame,proof.source30));
proof.axialReferenceFrame=axial.expandedFrame;
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
