function check_power_ui(out,mode,job)
% 模块一：从实际 STL 完成同一侧栏的网格、热边界和计算结果读取。
fig=rocketAeroHeatWorkbench();closeFig=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(out,'input.stl');c.scale.Value=1;press(c.load);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'shortSupport65');press(fig.UserData.controls.generate);
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.shortHeatBoundary,'power');press(fig.UserData.controls.confirm);
c=fig.UserData.controls;ready=fig.UserData.workflowState;
frozen=aeroHeatWorkflowState('beginSolve',ready,struct('shortSupportAction','solve'));
assert(strcmp(frozen.frozenInputs.cfg.shortSupportBoundary.branch,'power') && ~isempty(frozen.frozenInputs.cfg.shortSupportJob));
if strcmp(mode,'archive'),press(c.readShort);else,c.shortDirectory.Value=job;press(c.reopenShort);end
s=fig.UserData.workflowState;assert(strcmp(s.status,'solved'),'%s',fig.UserData.status.Text);
m=s.result;r=fig.UserData.caseReport;ui=fig.UserData.webUI;
assert(r.jointComputedFaceCount==32480 && r.jointExperimentalAreaM2==0);
timer=tic;while isempty(ui.html.UserData.rendered) && toc(timer)<45,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.faces==32480 && gl.visiblePixelSamples>0 && gl.webglError==0);
ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;ev=ui.evidenceView;d=ev.getData();sphere=m.shortSupportImport.regionID==1;
get=@(key)d.fields(strcmp({d.fields.key},key));a=get('wallSolved');b=get('wallInput');near=get('nearWallCellTemperature');
assert(isequal(isfinite(a.values(:)),sphere) && isequal(isfinite(b.values(:)),~sphere) && all(isfinite(near.values)));
for name={'edgeTemperature','recoveryTemperature'},f=get(name{1});assert(all(isnan(f.values)));end
assert(all(a.values(sphere)==m.shortSupportPower.temperatureK) && max(abs(b.values(~sphere)-m.shortSupportPower.temperatureK))<1e-8);
assert(all(a.status(sphere)==1) && all(b.status(~sphere)==2));
ledger=d.evidenceLedger.fields;sa=ledger(strcmp({ledger.key},'wallSolved'));sb=ledger(strcmp({ledger.key},'wallInput'));
assert(all(ismember(sa.availability(sphere),[1 2])) && all(sa.availability(~sphere)==0) && ...
    all(sb.availability(~sphere)==3) && all(sb.availability(sphere)==0));
% 模块二：点选球面、支管和端盖，读取真实截面及完整原面 CSV，离线页复用当前快照。
for id=1:5
    face=find(m.shortSupportImport.regionID==id,1);
    if id==1,ev.selectField('wallSolved');else,ev.selectField('wallInput');end
    ev.inspectFace(face,[]);probe=ev.getState();assert(probe.selectedFaceID==face);
end
ev.selectField('heatFlux');x=m.shortSupportImport.translationM(1)+.00137;ev.setSectionX(x);timer=tic;
while toc(timer)<45
    drawnow;cut=ev.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-13,break;end
    pause(.12);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-13);
ev.exportSectionTo(fullfile(out,'section.csv'));rec=ev.exportWallPairTo(fullfile(out,'wall.csv'));assert(rec.rows==32480);
exportAerothermalWeb(m,fullfile(out,'workbench'));
proof=struct('mode',mode,'frame',m.evidenceIdentity.sourceFrame.sha256,'targetW',s.freestream.targetSphereW, ...
    'temperatureK',m.shortSupportPower.temperatureK,'outW',-r.wholeSurfaceHeatPowerW,'regionsW',[r.regions.outW], ...
    'faces',gl.faces,'pixels',gl.visiblePixelSamples,'peakFile',r.peakOutFace,'peakReference',r.peakReferenceFaceID);
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.shortHeatBoundary,'temperature');press(fig.UserData.controls.confirm);
assert(isempty(fig.UserData.workflowState.result));ready=fig.UserData.workflowState;
rejected=false;try,aeroHeatWorkflowState('beginSolve',ready,struct('shortSupportAction','reopen','runDirectory',job));
catch ex,rejected=strcmp(ex.identifier,'RAH:Short:Incomplete') || strcmp(ex.identifier,'RAH:Short:Job');end
assert(rejected);press(fig.UserData.controls.readShort);assert(strcmp(fig.UserData.workflowState.status,'solved'));
assert(~isfield(fig.UserData.workflowState.result,'shortSupportPower'));proof.oldBoundaryRestored=true;
f=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(f>=0);close=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,val)
c.Value=val;c.ValueChangedFcn(c,[]);drawnow;
end
