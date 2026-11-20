function check_view(out,reference)
% 模块一：从 MATLAB 主侧栏进入同刻对照，读取 WebGL 实际显示的差场和图像。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.thermistorMode.Value='timeSteps6';
c.thermistorMode.ValueChangedFcn(c.thermistorMode,[]);c.thermistor.ButtonPushedFcn(c.thermistor,[]);drawnow;
view=fig.UserData.thermistorReview;assert(~isempty(view),fig.UserData.status.Text);
t=tic;while isempty(view.html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
gl=view.html.UserData.rendered;assert(~isempty(gl) && gl.webglError==0 && gl.visiblePixelSamples>0);
if isfield(gl,'observerError'),error('RAH:Test:Observed','%s',gl.observerError);end
assert(gl.faces==5216 && gl.surfaceTimeRows==5216 && contains(gl.csvName,'timeHeatDifference-t') && ...
    contains(gl.pngTitle,'入壁热流差') && contains(gl.pngTitle,'7.8093609 µs') && ...
    contains(gl.snapshotCaption,'主值为一个整步，对照为两个半步') && contains(gl.probeFace,'#'));
png=matlab.net.base64decode(extractAfter(string(gl.snapshotPNG),'base64,'));
assert(isequal(reshape(png(1:4),1,[]),uint8([137 80 78 71])));
fid=fopen(fullfile(out,'heat_difference.png'),'wb');assert(fid>=0);fwrite(fid,png,'uint8');fclose(fid);

% 模块二：按原面号独立重建差值，核对四项对照、共享色标、探针及求解覆盖分母。
v=view.evidenceView;d=v.getData();f=readtable(fullfile(reference,'wall.csv'));
[coarse,cleanup]=openThermistorArchive('startup6'); %#ok<NASGU>
b=readtable(fullfile(coarse,'wall.csv'));assert(isequal(f.faceID,b.faceID));
expected=f.heatFluxWm2-b.heatFluxWm2;
match('timeHalfPressure',f.pressurePa/1000);match('timeHalfGauge',f.pressureGaugePa);
match('timeHalfHeat',f.heatFluxWm2/1000);match('timeHalfNearWall',f.nearWallCellTemperatureK);
match('timePressureDifference',f.pressureGaugePa-b.pressureGaugePa);match('timeHeatDifference',expected);
match('heatFlux',b.heatFluxWm2/1000);match('pressure',b.pressurePa/1000);
assert(isequal(d.comparisonColorRanges.timeHalfHeat,d.comparisonColorRanges.heatFlux) && ...
    isequal(d.comparisonColorRanges.timeHalfGauge,d.comparisonColorRanges.pressureGauge));
assert(contains(gl.pngName,'timeHeatDifference-t') && d.timeStepComparison.originalFaces==gl.faces);
view.evidenceTab.Parent.SelectedTab=view.evidenceTab;
joint=d.evidenceLedger.joint;s=d.timeStepComparison;
assert(abs(joint.pressureHeatFluxAvailableArea-s.areaM2)<1e-15 && joint.pressureHeatFluxQualifiedArea==0);
v.selectField('timeHeatDifference');v.inspectFace(s.heatFluxWm2.maxDifferenceFaceID+1,[]);
assert(v.getState().selectedFaceID==s.heatFluxWm2.maxDifferenceFaceID+1);
x=.000172;v.setSectionX(x);t=tic;
while toc(t)<35
    drawnow;sec=v.getSection();
    if isstruct(sec) && isfield(sec,'csv') && strcmp(sec.fieldKey,'timeHeatDifference') && abs(sec.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(sec) && strcmp(sec.fieldKey,'timeHeatDifference') && abs(sec.x-x)<1e-12);
cutfile=fullfile(out,'difference_section.csv');v.exportSectionTo(cutfile);
cut=readtable(cutfile,'Delimiter',',','NumHeaderLines',0,'ReadVariableNames',true,'Encoding','UTF-8');
assert(height(cut)>0 && max(abs(cut.value-expected(cut.face_id)))<1e-9 && all(abs(cut.elapsed_physical_time_s-s.elapsedPhysicalTimeS)<1e-19));

% 模块三：离线包与主界面同源。核对数据脚本真实内容及原整步模式。
[~,report,page]=runThermistorBeadCase(fullfile(out,'export'),'timeSteps6');
assert(isfile(page) && report.timeStepComparison.elapsedPhysicalTimeS==d.resultCondition.elapsedPhysicalTimeS);
raw=strtrim(fileread(fullfile(fileparts(page),'case-data.js')));start='window.ROCKET_CASE=';
assert(startsWith(raw,start) && endsWith(raw,';'));
payload=jsondecode(raw(numel(start)+1:end-1));field=payload.fields(strcmp({payload.fields.key},'timeHeatDifference'));
assert(max(abs(field.values-expected))<1e-9 && strcmp(payload.defaultField,'timeHeatDifference'));
[old,baseline]=runThermistorBeadCase([],'startup6');
assert(~isfield(old,'beadStepComparison') && ~isfield(baseline,'timeStepComparison') && baseline.fullSurfaceComputed);
receipt=struct('faces',gl.faces,'visiblePixelSamples',gl.visiblePixelSamples,'webglError',gl.webglError, ...
    'physicalTimeS',s.elapsedPhysicalTimeS,'mainAndComparisonAreaM2',s.areaM2,'qualifiedM2',joint.pressureHeatFluxQualifiedArea, ...
    'caption',gl.snapshotCaption,'pngTitle',gl.pngTitle,'pngName',gl.pngName,'csvName',gl.csvName,'surfaceRows',gl.surfaceTimeRows, ...
    'actualPickedFace',gl.probeFace,'maxHeatDifferenceUIFace',s.heatFluxWm2.maxDifferenceFaceID+1,'sectionRows',height(cut), ...
    'sectionXM',x,'page',page,'baselineModeRechecked',true);
fid=fopen(fullfile(out,'comparison_ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(receipt,'PrettyPrint',true));fclose(fid);
    function match(key,values)
        field=d.fields(strcmp({d.fields.key},key));assert(isscalar(field) && all(field.status==1) && max(abs(field.values-values))<1e-9);
    end
end
