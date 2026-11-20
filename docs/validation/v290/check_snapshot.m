function check_snapshot(out,reference)
% 模块一：从主侧栏读取 6 m/s 启动快照，核对实际 WebGL 场和全部编号。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.thermistorMode.Value='startup6';
c.thermistorMode.ValueChangedFcn(c.thermistorMode,[]);c.thermistor.ButtonPushedFcn(c.thermistor,[]);drawnow;
view=fig.UserData.thermistorReview;assert(~isempty(view),fig.UserData.status.Text);
start=tic;while isempty(view.html.UserData.rendered) && toc(start)<35,drawnow;pause(.1);end
gl=view.html.UserData.rendered;
assert(~isempty(gl) && gl.faces==5216 && gl.webglError==0 && gl.visiblePixelSamples>0);
assert(contains(gl.snapshotCaption,'启动瞬态') && contains(gl.snapshotCaption,'7.8093609 µs'));
checkTime(gl);
assert(all(isfinite(str2double(string(gl.legend)))) && diff(str2double(string(gl.legend)))>0);
png=matlab.net.base64decode(extractAfter(string(gl.snapshotPNG),'base64,'));
assert(isequal(reshape(png(1:4),1,[]),uint8([137 80 78 71])));
fid=fopen(fullfile(out,'startup_heat.png'),'wb');assert(fid>=0);fwrite(fid,png,'uint8');fclose(fid);
v=view.evidenceView;d=v.getData();rows=readtable(fullfile(reference,'wall.csv'));
s=d.thermistorStartup;e=s.energy;assert(s.fullSurfaceComputed && d.resultCondition.elapsedPhysicalTimeS==e.dtSeconds);
check('pressure',rows.pressurePa/1000);check('pressureGauge',rows.pressureGaugePa);
q=check('heatFlux',rows.heatFluxWm2/1000);check('nearWallCellTemperature',rows.nearWallCellTemperatureK);
assert(contains(q.stateNote,'空气储热率') && all(contains(string(q.methodNames),'BDF1')));
tw=d.fields(strcmp({d.fields.key},'wallInput'));ts=d.fields(strcmp({d.fields.key},'wallSolved'));
assert(all(tw.values==303.15 & tw.status==2) && all(isnan(ts.values) & ts.status==0));

% 模块二：选取热峰，逐面回读成对 CSV，并用导出的面积和有符号热流重算整珠功率。
v.selectField('heatFlux');v.inspectFace(s.peakOutwardFaceID+1,[]);
assert(v.getState().selectedFaceID==s.peakOutwardFaceID+1);
pairPath=fullfile(out,'wall_pair.csv');pair=v.exportWallPairTo(pairPath);
tab=readtable(pairPath,'Delimiter',',','NumHeaderLines',0,'ReadVariableNames',true,'Encoding','UTF-8');
assert(isequal(tab.face_id,rows.faceID+1) && max(abs(tab.pressure*1000-rows.pressurePa))<1e-8 && ...
    max(abs(tab.heat_flux*1000-rows.heatFluxWm2))<1e-8 && abs(sum(tab.area_m2.*tab.heat_flux)*1000+e.wallReleaseW)<1e-13);
assert(pair.joint.pressureHeatFluxQualifiedArea==0 && pair.identity.sourceTraceable && ...
    all(contains(string(tab.condition),'启动瞬态')) && all(contains(string(tab.heat_flux_method),'BDF1')));
% 模块三：真实截面按其单个原面编号回查，验证固定时刻的离线出口和旧模式。
view.evidenceTab.Parent.SelectedTab=view.evidenceTab;x=.000172;v.setSectionX(x);start=tic;
while toc(start)<35
    drawnow;sec=v.getSection();
    if isstruct(sec) && isfield(sec,'csv') && strcmp(sec.fieldKey,'heatFlux') && abs(sec.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(sec) && isfield(sec,'csv') && abs(sec.x-x)<1e-12);
cutPath=fullfile(out,'section.csv');v.exportSectionTo(cutPath);cut=readtable(cutPath);
assert(height(cut)>0 && max(abs(cut.value-q.values(cut.face_id)))<1e-10);
assert(all(string(cut.availability)=="研究计算"));
assert(all(abs(cut.elapsed_physical_time_s-e.dtSeconds)<1e-19) && all(contains(string(cut.condition),gl.fieldTime)));
[~,report,page]=runThermistorBeadCase(fullfile(out,'export'),'startup6');
assert(isfile(page) && report.fullSurfaceComputed && report.elapsedPhysicalTimeS==e.dtSeconds);
% 模块四：离线页仅从发布的数据脚本载入，启动时不从 MATLAB 接收任何模型。
folder=fullfile(out,'offline_check');copyfile(fileparts(page),folder);
offlinePage=fullfile(folder,'index.html');doc=fileread(offlinePage);
copyfile(fullfile(fileparts(mfilename('fullpath')),'observe_display.js'),folder);
fid=fopen(offlinePage,'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s',strrep(doc,'</body>','<script src="observe_display.js"></script></body>'));fclose(fid);
offlineFig=uifigure('Position',[50 50 1060 720]);finish=onCleanup(@()delete(offlineFig)); %#ok<NASGU>
layout=uigridlayout(offlineFig,[1 1]);browser=uihtml(layout,'HTMLSource',offlinePage);
start=tic;
while toc(start)<35
    drawnow;standalone=browser.Data;
    if isstruct(standalone) && isfield(standalone,'kind') && strcmp(standalone.kind,'offline-proof'),break;end
    pause(.1);
end
checkTime(standalone);
old=runThermistorBeadCase([],'archive7p5');assert(~isfield(old,'thermistorStartup') && old.hybridCondition.wallTemperatureK==305.15);
try
    runThermistorBeadCase([],'startup6',[],[],299.15);error('RAH:Test:Accepted','固定时间快照接收了另一壁温');
catch ex
    assert(strcmp(ex.identifier,'RAH:Bead:Startup'));
end
proof=struct('faces',gl.faces,'visiblePixelSamples',gl.visiblePixelSamples,'webglError',gl.webglError, ...
    'elapsedPhysicalTimeS',e.dtSeconds,'caption',gl.snapshotCaption,'legend',{gl.legend},'selectedUIFace',s.peakOutwardFaceID+1, ...
    'sectionXM',x,'sectionRows',height(cut),'wallPowerIntoW',report.wholeSurfaceHeatPowerW, ...
    'computedAreaM2',report.jointComputedAreaM2,'jointExperimentVerifiedM2',pair.joint.pressureHeatFluxQualifiedArea, ...
    'wallCsvSHA256',cm011FileSHA256(pairPath),'sectionCsvSHA256',cm011FileSHA256(cutPath),'page',page, ...
    'fieldTime',gl.fieldTime,'pngTitle',gl.pngTitle,'pngName',gl.pngName,'csvName',gl.csvName, ...
    'pickedProbe',gl.probeFace,'probeTimes',{gl.probeTimes},'surfaceTimeRows',gl.surfaceTimeRows, ...
    'offlineSurfaceTimeRows',standalone.surfaceTimeRows,'offlineProbe',standalone.probeFace, ...
    'offlineSectionTimeRows',standalone.sectionTimeRows);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
    function field=check(key,expected)
        field=d.fields(strcmp({d.fields.key},key));
        assert(isscalar(field) && all(field.status==1) && max(abs(field.values-expected))<1e-9);
    end
    function checkTime(r)
        assert(isstruct(r) && isfield(r,'fieldTime'),'没有收到页面时刻回执。');
        if isfield(r,'observerError'),error('RAH:Test:Web','%s',r.observerError);end
        assert(strcmp(r.fieldTime,'距构造初场 7.8093609 µs') && contains(r.pngTitle,r.fieldTime));
        assert(contains(r.probeFace,'#') && ~isempty(r.probeTimes) && all(contains(string(r.probeTimes),r.fieldTime)));
        name='RocketAeroHeat-heatFlux-t0.00000780936086544745s';
        assert(strcmp(r.csvName,[name '.csv']) && strcmp(r.pngName,[name '.png']));
        assert(r.surfaceTimeRows==5216 && r.sectionRows>0 && r.sectionTimeRows==r.sectionRows);
    end
end
