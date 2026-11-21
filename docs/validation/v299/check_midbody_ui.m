function check_midbody_ui(out)
% 模块一：通过生产侧栏生成体网、确认工况并回读 6,624 面的归档结果。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.midbodyBeadModel);press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;assert(s.mesh.thermistorPrepared.cells==22752,fig.UserData.status.Text);
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
r=s.result;assert(strcmp(r.thermistorPaired.kind,'midbody') && size(r.F,1)==6624);
web=fig.UserData.webUI;t=tic;
while isempty(web.html.UserData.rendered) && toc(t)<35,drawnow;pause(.1);end
gl=web.html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==6624);
v=web.evidenceView;d=v.getData();p=r.thermistorPaired;
assert(isequal(d.comparisonColorRanges.heatFlux,d.comparisonColorRanges.pairedSmoothHeat));
fields={'heatFlux',r.hybridWall.heatFluxWm2/1000;'pressure',r.hybridWall.pressurePa/1000; ...
    'nearWallCellTemperature',r.hybridWall.nearWallCellTemperatureK;'leadReverse',r.thermistorLeadFlow.reverse; ...
    'pairedDeltaHeat',p.deltaHeatWm2;'pairedDeltaPressure',p.deltaPressurePa;'pairedReverseChange',p.reverseChange};
for j=1:size(fields,1)
    field=d.fields(strcmp({d.fields.key},fields{j,1}));
    assert(all(field.status(:)==1) && max(abs(field.values(:)-fields{j,2}))<1e-9,fields{j,1});
end
% 模块二：验证按钮拾取、固定色标、截面与全表面导出均来自同一份差值。
press(fig.UserData.controls.pairedHeat);[~,peak]=max(abs(p.deltaHeatWm2));
state=v.getState();assert(strcmp(state.fieldKey,'pairedDeltaHeat') && state.selectedFaceID==peak);
bound=max(abs(p.deltaHeatWm2));assert(isequal(d.comparisonColorRanges.pairedDeltaHeat,[-bound bound]));
v.exportTo(fullfile(out,'faces.csv'));x=.000172;v.setSectionX(x);t=tic;
while toc(t)<35
    drawnow;sec=v.getSection();
    if isstruct(sec) && isfield(sec,'csv') && strcmp(sec.fieldKey,'pairedDeltaHeat') && abs(sec.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(sec) && isfield(sec,'csv') && strcmp(sec.fieldKey,'pairedDeltaHeat') && abs(sec.x-x)<1e-12);
v.exportSectionTo(fullfile(out,'section.csv'));cut=readtable(fullfile(out,'section.csv'),'Encoding','UTF-8');
assert(height(cut)>0 && max(abs(cut.value-p.deltaHeatWm2(cut.face_id)))<1e-9);
page=exportAerothermalWeb(r,fullfile(out,'workbench'));id=r.evidenceIdentity.sourceFrame;
proof=struct('faces',size(r.F,1),'webglError',gl.webglError,'visibleSamples',gl.visiblePixelSamples,'sourceFrame',id.sha256, ...
    'peakDeltaFaceID',peak,'sectionX',x,'sectionRows',height(cut),'page',page,'meshSHA256',r.evidenceIdentity.meshSHA256, ...
    'viewHeatRangeKWm2',d.comparisonColorRanges.heatFlux,'deltaRangeWm2',[-bound bound]);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
