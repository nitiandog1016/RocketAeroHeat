function check_angle_view(out)
% 模块一：从主窗口原菜单打开角向细网，在同一工作台读取云图和原面台账。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
controls=fig.UserData.controls;controls.thermistorMode.Value='solve7p5';
controls.thermistorMode.ValueChangedFcn(controls.thermistorMode,[]);
assert(contains(controls.thermistor.Text,'计算') && isempty(fig.UserData.thermistorReview));
controls.thermistorMode.Value='archive7p5';controls.thermistorMode.ValueChangedFcn(controls.thermistorMode,[]);
controls.thermistor.ButtonPushedFcn(controls.thermistor,[]);drawnow;
assert(~isempty(fig.UserData.thermistorReview),fig.UserData.status.Text);
assert(strcmp(fig.UserData.controls.thermistorMode.Value,'archive7p5') && contains(fig.UserData.controls.thermistor.Text,'打开'));
view=fig.UserData.thermistorReview;start=tic;
while isempty(view.html.UserData.rendered) && toc(start)<35,drawnow;pause(.1);end
render=view.html.UserData.rendered;
assert(~isempty(render) && render.faces==5216 && render.webglError==0 && render.visiblePixelSamples>0);
d=view.evidenceView.getData();assert(d.thermistorBead.cellAngleDeg==7.5);
old=readtable(fullfile(root,'outputs','v282_validation','result_outer01','wall.csv'));
new=readtable(fullfile(root,'outputs','v283_validation','result01','wall.csv'));
deltaQ=new.heatFluxWm2-old.heatFluxWm2;deltaP=new.pressureGaugePa-old.pressureGaugePa;
checkField('angularHeatDifference',deltaQ,1e-9);
checkField('angularPressureDifference',deltaP,1e-11);
checkField('heatFlux',new.heatFluxWm2/1000,1e-11);
checkField('pressure',new.pressurePa/1000,1e-10);

% 模块二：选取差异最大的原面，复查选择、差值表和主压力／热流表的来源。
[~,id]=max(abs(deltaQ));v=view.evidenceView;
v.selectField('angularHeatDifference');v.inspectFace(id,[]);
state=v.getState();assert(state.selectedFaceID==id);
receipt=v.exportTo(fullfile(out,'delta_heat.csv'));
export=readtable(receipt.csvPath);
assert(isequal(export.face_id,new.faceID+1) && max(abs(export.value-deltaQ))<1e-9);
assert(all(export.availability_code==2) && all(export.qualification_code==0));
pair=v.exportWallPairTo(fullfile(out,'wall_pair.csv'));tab=readtable(pair.csvPath);
assert(pair.rows==5216 && pair.identity.sourceTraceable && ...
    pair.joint.pressureHeatFluxQualifiedArea==0 && pair.joint.pressureHeatFluxEngineeringArea==0);
assert(max(abs(tab.heat_flux*1000-new.heatFluxWm2))<1e-9 && max(abs(tab.pressure*1000-new.pressurePa))<1e-8);

% 模块三：验证公开版页面实际画出的横截面，并按其原面号反查每个曲线点。
view.evidenceTab.Parent.SelectedTab=view.evidenceTab;
sectionX=0.000172;v.setSectionX(sectionX);start=tic;
while toc(start)<30
    drawnow;section=v.getSection();
    if isstruct(section) && isfield(section,'fieldKey') && strcmp(section.fieldKey,'angularHeatDifference') ...
            && isfield(section,'x') && abs(section.x-sectionX)<1e-11 && isfield(section,'csv'),break;end
    pause(.1);
end
assert(isstruct(section) && isfield(section,'csv') && abs(section.x-sectionX)<1e-11);
sectionPath=fullfile(out,'section.csv');v.exportSectionTo(sectionPath);
sectionRows=readtable(sectionPath);assert(height(sectionRows)>0);
assert(max(abs(sectionRows.value-deltaQ(sectionRows.face_id)))<1e-9);
[~,report,page]=runThermistorBeadCase(fullfile(out,'export'),'archive7p5');
assert(isfile(page) && report.fullSurfaceComputed);
status=struct('faces',height(export),'render',render,'selectedNativeID',id-1,...
    'pressureHeatQualifiedAreaM2',pair.joint.pressureHeatFluxQualifiedArea,'comparedEveryFace',true,...
    'heatDifferenceAtSelectedFaceWm2',deltaQ(id),'sectionX',sectionX,'sectionRows',height(sectionRows),...
    'sourceFrame',d.evidenceIdentity.sourceFrame,'powerW',report.wholeSurfaceHeatPowerW,'page',page);
fid=fopen(fullfile(out,'view.json'),'w','n','UTF-8');assert(fid>=0);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(status,'PrettyPrint',true));
    function checkField(key,expected,tolerance)
        field=d.fields(strcmp({d.fields.key},key));
        assert(isscalar(field) && all(field.status==1) && max(abs(field.values-expected))<tolerance);
    end
end
