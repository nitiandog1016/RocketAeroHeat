function check_power_ui(out,nativeDir)
% 模块一：从实际侧栏选择带线中段加密网，确认整件功率并取得全表面解。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,1);
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(s.mesh.thermistorPrepared.cells==23224,fig.UserData.status.Text);
press(fig.UserData.controls.nextBoundary);c=fig.UserData.controls;
assert(strcmp(c.beadCondition.Value,'published') && strcmp(c.confirm.Enable,'on'));
choose(c.beadCondition,'power');assert(strcmp(c.beadWall.Enable,'off') && strcmp(c.beadPower.Enable,'on'));
assert(c.beadPower.Value==.004);press(c.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
r=s.result;t=r.thermistorLead.thermal;n=size(r.F,1);assert(n==10272);
[folder,g]=openThermistorArchive('leadsMidbody3'); %#ok<NASGU>
w=readtable(fullfile(folder,'wall.csv'));ref=jsondecode(fileread(fullfile(folder,'summary.json')));
checkResult(r,.004,w,ref);gl=fig.UserData.webUI.html;t0=tic;
while isempty(gl.UserData.rendered) && toc(t0)<35,drawnow;pause(.1);end
render=gl.UserData.rendered;assert(render.webglError==0 && render.faces==n && render.visiblePixelSamples>0);
% 模块二：三维导出和面元探针应使用同一个功率输入与计算温度。
v=fig.UserData.webUI.evidenceView;d=v.getData();
tw=d.fields(strcmp({d.fields.key},'wallSolved'));given=d.fields(strcmp({d.fields.key},'wallInput'));
assert(all(tw.status==1 & tw.values==t.wallTemperatureK) && all(given.status==0 & isnan(given.values)));
assert(strcmp(d.resultCondition.wallInputType,'netHeatingPower') && d.resultCondition.wallKnownFaces==0 && ...
    contains(d.caseNote,'整件净对流功率 4 mW') && all(contains(string(tw.methodNames),'热平衡')));
v.selectField('wallSolved');v.inspectFace(207,[]);assert(v.getState().selectedFaceID==207);
v.exportTo(fullfile(out,'temperature.csv'));receipt=v.exportWallPairTo(fullfile(out,'wall.csv'));assert(receipt.joint.pressureHeatFluxQualifiedArea==0);
x=.000172;keys={'wallSolved','heatFlux'};
for j=1:numel(keys)
    key=keys{j};v.selectField(key);v.setSectionX(x);t0=tic;
    while toc(t0)<35
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,key) && abs(cut.x-x)<1e-12);
    file=fullfile(out,[key '_section.csv']);v.exportSectionTo(file);table=readtable(file,'Encoding','UTF-8');
    f=d.fields(strcmp({d.fields.key},key));assert(height(table)>0 && max(abs(table.value-f.values(table.face_id)))<1e-9);
end
page=exportAerothermalWeb(r,fullfile(out,'workbench'));
% 模块三：按同一求解入口检验端点、超域、独立半温差功率反演及每个面元。
cfg=s.cfg;cfg.thermistorResultSource='archive';cfg.thermistorRunDirectory='';
limit=thermistorNetPowerLimit(4,cfg);assert(limit==ref.totalOutwardW);
half=jsondecode(fileread(fullfile(nativeDir,'half.json')));direct=readtable(fullfile(nativeDir,'direct_wall.csv'));
for watts=[0 half.wallReleaseW limit]
    cfg.thermistorBoundary.netHeatingPowerW=watts;
    [trial,~,report]=runAeroHeatSingleCase(s.mesh,cfg);
    assert(report.fullSurfaceComputed);checkResult(trial,watts,w,ref);
    if watts==half.wallReleaseW
        inverseErrorK=abs(trial.hybridSolvedWallTemperatureK(1)-299.15);
        err=abs(trial.hybridWall.heatFluxWm2./direct.heatFluxWm2-1);
        assert(inverseErrorK<.001 && max(err)<.0005 && max(abs(trial.hybridWall.pressurePa-direct.pressurePa))<.001);
        assert(max(abs(trial.hybridWall.nearWallCellTemperatureK-direct.firstLayerK))<.001);
    elseif watts==0
        zero=buildAerothermalViewerData(trial);
        assert(zero.uniformFieldValues.heatFlux.value==0 && zero.uniformFieldValues.wallSolved.value==293.15);
    end
end
for value={-eps,limit+eps(limit),NaN,Inf,[0 .001],1i}
    cfg.thermistorBoundary.netHeatingPowerW=value{1};
    try
        runAeroHeatSingleCase(s.mesh,cfg);error('RAH:Test:Accepted','错误输入通过');
    catch ex,assert(strcmp(ex.identifier,'RAH:Bead:Boundary'),ex.message);end
end
% 模块四：改为零功率并重新确认；再切回给定温度，核对场量身份和原来的网距对照。
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;choose(c.beadPower,0);
s=fig.UserData.workflowState;assert(isempty(s.result) && ~s.boundaryConfirmed && isempty(fig.UserData.webUI));
press(c.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;checkResult(s.result,0,w,ref);
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;choose(c.beadCondition,'published');
press(c.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && isfield(s.result,'thermistorLeadGrid') && ~isfield(s.result,'hybridSolvedWallTemperatureK') && all(s.result.hybridWallTemperatureK==305.15));
proof=struct('faces',n,'areaM2',sum(w.areaM2),'visibleSamples',render.visiblePixelSamples,'webglError',render.webglError, ...
    'thermal',t,'regionW',[r.thermistorLead.regions.outwardW],'regionNames',{ {r.thermistorLead.regions.name} }, ...
    'probeFaceID',207,'sectionX',x,'inverseNativeErrorK',inverseErrorK,'conditionSHA256',receipt.identity.conditionSHA256,'page',page);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function checkResult(r,power,w,ref)
% 校验所有面与参照场的对应关系，并对每个暴露区域独立积分。
t=r.thermistorLead.thermal;assert(t.netHeatingPowerW==power && isequaln(r.thermistorLead.referenceNative,ref));
assert(~isfield(r,'hybridWallTemperatureK') && ~isfield(r.hybridCondition,'wallTemperatureK') && ...
    ~isfield(r,'thermistorLeadGrid') && ~isfield(r,'thermistorPaired') && isequal(r.evidenceIdentity.sourceFrame.thermalConversion,t));
assert(r.hybridCondition.netHeatingPowerW==power && isequal(r.hybridCondition,r.evidenceIdentity.physicalCondition));
scale=power/ref.totalOutwardW;
assert(all(r.hybridSolvedWallTemperatureK==293.15+12*scale) && isequal(w.pressurePa,r.hybridWall.pressurePa) && ...
    isequal(r.hybridWall.faceID,w.faceID) && isequal(r.hybridWall.heatFluxWm2,scale*w.heatFluxWm2));
assert(max(abs(r.hybridWall.nearWallCellTemperatureK-(293.15+scale*(w.nearWallCellTemperatureK-293.15))))<1e-12);
for j=1:7
    use=w.regionID==j;
    assert(abs(r.thermistorLead.regions(j).outwardW+sum(w.areaM2(use).*r.hybridWall.heatFluxWm2(use)))<=1e-11*power);
end
assert(abs(r.thermistorLead.totalOutwardW-power)<=1e-11*power);
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,val)
c.Value=val;c.ValueChangedFcn(c,[]);drawnow;
end
