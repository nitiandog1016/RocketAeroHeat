function check_bead_ui(out)
% 模块一：通过主窗口现有控件打开热敏珠，核对三维渲染回执和原面探针。
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
s=fig.UserData.workflowState;button=fig.UserData.controls.thermistor;
button.ButtonPushedFcn(button,[]);drawnow;
assert(isfield(fig.UserData,'thermistorReview') && ~isempty(fig.UserData.thermistorReview),fig.UserData.status.Text);
view=fig.UserData.thermistorReview;
start=tic;
while isempty(view.html.UserData.rendered) && toc(start)<40,drawnow;pause(.1);end
render=view.html.UserData.rendered;
assert(~isempty(render) && render.faces==5216 && render.webglError==0 && render.visiblePixelSamples>0);
assert(fig.UserData.workflowState.revision==s.revision && isequal(fig.UserData.workflowState.result,s.result));
v=view.evidenceView;v.selectField('heatFlux');v.inspectFace(2610,[]);
vs=v.getState();assert(vs.selectedFaceID==2610);
pair=v.exportWallPairTo(fullfile(out,'wall_pair.csv'));
assert(pair.rows==5216 && pair.identity.sourceTraceable && strcmp(pair.identity.origin,'research') && ...
    pair.joint.pressureHeatFluxQualifiedArea==0 && pair.joint.pressureHeatFluxEngineeringArea==0);

% 模块二：逐面复查导出单位、原编号和温度身份，再单独生成无需服务的云图目录。
[cur,report,page]=runThermistorBeadCase(fullfile(out,'export'),'archive');
d=view.html.UserData.payload;w=cur.hybridWall;t=readtable(pair.csvPath);
assert(isequal(t.face_id,w.faceID+1) && all(abs(t.pressure*1000-w.pressurePa)<1e-8) && ...
    all(abs(t.heat_flux*1000-w.heatFluxWm2)<1e-9));
assert(max(abs(d.fields(strcmp({d.fields.key},'nearWallCellTemperature')).values-w.nearWallCellTemperatureK))<1e-10);
assert(all(d.fields(strcmp({d.fields.key},'wallInput')).status==2));
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    a=d.fields(strcmp({d.fields.key},key{1}));assert(all(a.status==0 & isnan(a.values)));
end
assert(isfile(page) && contains(d.fields(strcmp({d.fields.key},'heatFlux')).stateNote,'图 7'));
% 新计算按钮使用的解释器须能从 MATLAB 启动并载入全链依赖；这里仅查询入口。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
runCM011Python({fullfile(root,'scripts','runThermistorCase.py'),'--help'},fullfile(out,'solver_entry.log'));
assert(contains(fileread(fullfile(out,'solver_entry.log')),'new_output'));
b=fig.UserData.workflowButtons(1);b.ButtonPushedFcn(b,[]);
assert(isempty(fig.UserData.thermistorReview) && fig.UserData.workflowState.revision==s.revision);
proof=struct('render',render,'probeDisplayID',2610,'probeNativeID',w.faceID(2610), ...
    'areaM2',sum(t.area_m2),'powerW',report.wholeSurfaceHeatPowerW,'verifiedJointAreaM2',0, ...
    'sourceFrame',cur.evidenceIdentity.sourceFrame,'pairedRows',height(t),'sourceTraceable',pair.identity.sourceTraceable, ...
    'comparedAllRows',true,'returnedToOriginalInput',true,'solverEntryLoadedFromMatlab',true,'page',page);
fid=fopen(fullfile(out,'view.json'),'w','n','UTF-8');assert(fid>=0);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
