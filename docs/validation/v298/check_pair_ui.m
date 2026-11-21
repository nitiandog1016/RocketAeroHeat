function check_pair_ui(out,kind)
% 模块一：从生产侧栏按实际工况读入主结果，分别检查带线几何与原光滑几何。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
wired=strcmp(kind,'wired');
if wired
    press(fig.UserData.controls.wiredBeadModel);
else
    press(fig.UserData.controls.beadModel);choose(fig.UserData.controls.beadCells,4);
end
press(fig.UserData.controls.generate);assert(~isempty(fig.UserData.workflowState.mesh),fig.UserData.status.Text);
press(fig.UserData.controls.nextBoundary);
if ~wired,choose(fig.UserData.controls.beadCondition,'published');end
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
pair=s.result.thermistorPaired;assert(strcmp(pair.kind,kind));web=fig.UserData.webUI;t=tic;
while isempty(web.html.UserData.rendered) && toc(t)<35,drawnow;pause(.1);end
gl=web.html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==size(s.result.F,1));
v=web.evidenceView;d=v.getData();mask=isfinite(pair.heatFluxWm2);
if wired,assert(nnz(mask)==6144 && nnz(~mask)==4128);else,assert(all(mask) && numel(mask)==5216);end
fields={'pairedSmoothHeat',pair.heatFluxWm2/1000;'pairedDeltaHeat',pair.deltaHeatWm2; ...
    'pairedDeltaPressure',pair.deltaPressurePa;'pairedSmoothReverse',pair.reverse;'pairedReverseChange',pair.reverseChange};
for k=1:size(fields,1)
    field=d.fields(strcmp({d.fields.key},fields{k,1}));
    assert(isequaln(field.values(:),fields{k,2}) && isequal(field.status(:)>0,mask));
end
assert(isequal(d.comparisonColorRanges.heatFlux,d.comparisonColorRanges.pairedSmoothHeat) && ...
    isequal(d.comparisonColorRanges.pairedReverseChange,[-1 1]));
% 模块二：检查未知面的标记，真实截面按原编号导出，同一字段也进入离线页。
press(fig.UserData.controls.pairedHeat);
good=find(mask);[~,peak]=max(abs(pair.deltaHeatWm2(good)));
assert(web.evidenceTab.Parent.SelectedTab==web.evidenceTab && strcmp(v.getState().fieldKey,'pairedDeltaHeat') && ...
    v.getState().selectedFaceID==good(peak));
if wired
    missing=find(~mask,1);v.inspectFace(missing,[]);assert(v.getState().selectedFaceID==missing);
end
face=find(mask,1);v.inspectFace(face,[]);v.exportTo(fullfile(out,'faces.csv'));
x=.000172;v.setSectionX(x);t=tic;
while toc(t)<35
    drawnow;section=v.getSection();
    if isstruct(section) && isfield(section,'csv') && strcmp(section.fieldKey,'pairedDeltaHeat') && abs(section.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(section) && isfield(section,'csv') && strcmp(section.fieldKey,'pairedDeltaHeat') && abs(section.x-x)<1e-12);
path=fullfile(out,'section.csv');v.exportSectionTo(path);cut=readtable(path,'Encoding','UTF-8');
assert(height(cut)>0 && all(mask(cut.face_id)) && max(abs(cut.value-pair.deltaHeatWm2(cut.face_id)))<1e-9);
page=exportAerothermalWeb(s.result,fullfile(out,'workbench'));
proof=struct('kind',kind,'faces',numel(mask),'paired',nnz(mask),'unmatched',nnz(~mask),'sectionX',x, ...
    'sectionRows',height(cut),'page',page,'sourceFrame',s.result.evidenceIdentity.sourceFrame.sha256, ...
    'buttonPeakFaceID',good(peak),'webglError',gl.webglError);
if ~wired
    % 改变壁温后，明确撤销 305.15 K 的冻结对照。
    press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.beadWall,299.15);
    press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
    s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isfield(s.result,'thermistorPaired'));
    proof.changedTemperatureExcluded=true;
end
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
function choose(control,value)
control.Value=value;control.ValueChangedFcn(control,[]);drawnow;
end
