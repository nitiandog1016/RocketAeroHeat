function check_near_wall_view(out,tableDir)
% 模块一：从主侧栏读取近壁对照，核对实际 WebGL 绘制、鼠标拾取、图题和下载。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.thermistorMode.Value='nearWall6';
c.thermistorMode.ValueChangedFcn(c.thermistorMode,[]);c.thermistor.ButtonPushedFcn(c.thermistor,[]);drawnow;
view=fig.UserData.thermistorReview;assert(~isempty(view),fig.UserData.status.Text);
t=tic;while isempty(view.html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
gl=view.html.UserData.rendered;assert(~isempty(gl) && gl.webglError==0 && gl.visiblePixelSamples>0);
if isfield(gl,'observerError'),error('RAH:Test:Observed','%s',gl.observerError);end
assert(gl.faces==5216 && gl.surfaceTimeRows==5216 && contains(gl.csvName,'gridHeatDifference-t') && ...
    contains(gl.pngName,'gridHeatDifference-t') && contains(gl.pngTitle,'入壁热流差') && ...
    contains(gl.pngTitle,'7.8093609 µs') && contains(gl.snapshotCaption,'主值为 56 层，对照为 48 层') && contains(gl.probeFace,'#'));
png=matlab.net.base64decode(extractAfter(string(gl.snapshotPNG),'base64,'));
assert(isequal(reshape(png(1:4),1,[]),uint8([137 80 78 71])));
fid=fopen(fullfile(out,'near_wall_heat_difference.png'),'wb');assert(fid>=0);fwrite(fid,png,'uint8');fclose(fid);

% 模块二：从两份原表复核主值、同位置温差、面积、共享色标与截面导出。
v=view.evidenceView;d=v.getData();s=d.nearWallComparison;
f=readtable(fullfile(tableDir,'wall.csv'));m=readtable(fullfile(tableDir,'matched_temperature.csv'));
[folder,cleanup]=openThermistorArchive('halfsteps'); %#ok<NASGU>
b=readtable(fullfile(folder,'wall.csv'));assert(isequal(f.faceID,b.faceID,m.faceID));
dq=f.heatFluxWm2-b.heatFluxWm2;matched=m.temperatureAtOriginalFirstLayerK;
match('pressure',f.pressurePa/1000);match('pressureGauge',f.pressureGaugePa);match('heatFlux',f.heatFluxWm2/1000);
match('nearWallCellTemperature',f.nearWallCellTemperatureK);match('gridMatchedTemperature',matched);
match('gridRefPressure',b.pressurePa/1000);match('gridRefGauge',b.pressureGaugePa);
match('gridRefTemperature',b.nearWallCellTemperatureK);match('gridRefHeat',b.heatFluxWm2/1000);
match('gridHeatDifference',dq);match('gridPressureDifference',f.pressureGaugePa-b.pressureGaugePa);
match('gridTemperatureDifference',matched-b.nearWallCellTemperatureK);
given=d.fields(strcmp({d.fields.key},'wallInput'));
solid=d.fields(strcmp({d.fields.key},'wallSolved'));
assert(all(given.status==2 & given.values==303.15) && all(solid.status==0 & isnan(solid.values)));
ranges=d.comparisonColorRanges;
assert(isequal(ranges.gridRefHeat,ranges.heatFlux) && isequal(ranges.gridRefGauge,ranges.pressureGauge) && ...
    isequal(ranges.gridRefTemperature,ranges.gridMatchedTemperature) && ...
    isequal(ranges.gridHeatDifference,[-max(abs(dq)),max(abs(dq))]));
j=d.evidenceLedger.joint;
assert(abs(j.pressureHeatFluxAvailableArea-s.areaM2)<1e-15 && j.pressureHeatFluxQualifiedArea==0);
assert(strcmp(d.evidenceLedger.sourceFrame.sha256,s.fineSurface.energy.acceptedFrameSHA256));
view.evidenceTab.Parent.SelectedTab=view.evidenceTab;
face=s.heatFluxWm2.maxDifferenceFaceID+1;v.selectField('gridHeatDifference');v.inspectFace(face,[]);
assert(v.getState().selectedFaceID==face);
x=.000172;v.setSectionX(x);t=tic;
while toc(t)<35
    drawnow;sec=v.getSection();
    if isstruct(sec) && isfield(sec,'csv') && strcmp(sec.fieldKey,'gridHeatDifference') && abs(sec.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(sec) && strcmp(sec.fieldKey,'gridHeatDifference') && abs(sec.x-x)<1e-12);
cutfile=fullfile(out,'difference_section.csv');v.exportSectionTo(cutfile);
cut=readtable(cutfile,'Delimiter',',','NumHeaderLines',0,'ReadVariableNames',true,'Encoding','UTF-8');
assert(height(cut)>0 && max(abs(cut.value-dq(cut.face_id)))<1e-9 && ...
    all(abs(cut.elapsed_physical_time_s-s.elapsedPhysicalTimeS)<1e-19));

% 模块三：核对离线页的主场及差场，回查原时间步对照仍读取原来的两份快照。
[~,report,page]=runThermistorBeadCase(fullfile(out,'export'),'nearWall6');
assert(isfile(page) && abs(report.wholeSurfaceHeatPowerW+s.fineSurface.energy.wallReleaseW)<1e-12);
raw=strtrim(fileread(fullfile(fileparts(page),'case-data.js')));start='window.ROCKET_CASE=';
assert(startsWith(raw,start) && endsWith(raw,';'));
p=jsondecode(raw(numel(start)+1:end-1));field=p.fields(strcmp({p.fields.key},'gridHeatDifference'));
assert(max(abs(field.values-dq))<1e-9 && strcmp(p.defaultField,'gridHeatDifference'));
[old,oldReport]=runThermistorBeadCase([],'timeSteps6');
assert(~isfield(old,'beadNearWallComparison') && oldReport.fullSurfaceComputed && ...
    strcmp(oldReport.timeStepComparison.fineFrameSHA256,s.referenceFrameSHA256));
receipt=struct('faces',gl.faces,'visiblePixelSamples',gl.visiblePixelSamples,'webglError',gl.webglError, ...
    'physicalTimeS',s.elapsedPhysicalTimeS,'mainAndComparisonAreaM2',s.areaM2,'qualifiedM2',j.pressureHeatFluxQualifiedArea, ...
    'caption',gl.snapshotCaption,'pngTitle',gl.pngTitle,'pngName',gl.pngName,'csvName',gl.csvName,'surfaceRows',gl.surfaceTimeRows, ...
    'actualPickedFace',gl.probeFace,'maxHeatDifferenceUIFace',face,'sectionRows',height(cut), ...
    'sectionXM',x,'mainHeatPowerIntoWallW',report.wholeSurfaceHeatPowerW,'page',page,'timeStepModeRechecked',true);
fid=fopen(fullfile(out,'near_wall_ui.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(receipt,'PrettyPrint',true));fclose(fid);
    function match(key,values)
        field=d.fields(strcmp({d.fields.key},key));
        assert(isscalar(field) && all(field.status==1) && max(abs(field.values-values))<1e-9,key);
    end
end
