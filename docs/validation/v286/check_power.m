function check_power(out,cells)
% 模块一：用绑定原面表独立核对端点、输入准入和同功率网格对照。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
[folder,cleanup]=openThermistorArchive(cells); %#ok<NASGU>
rows=readtable(fullfile(folder,'wall.csv'));
stats=jsondecode(fileread(fullfile(folder,'summary.json')));
limit=-stats.powerW;mode='archive';if cells==4,mode='archive7p5';end
assert(thermistorNetPowerLimit(cells)==limit);
proof=struct('cellAngleDeg',30/cells,'referenceFrameSHA256',stats.sourceFrameSHA256, ...
    'nativeTableSHA256',stats.wallTableSHA256,'limitW',limit,'cases',[]);
nativePath=fullfile(root,'docs','validation','v285','evidence','native.json');
assert(strcmp(cm011FileSHA256(nativePath),'41cec8d1f80bd116687b0c7911ab09f1921b867da2994d8949dce546cb547d09'));
direct=jsondecode(fileread(nativePath));powers=[0 .004 limit];
if cells==2,powers(end+1)=-direct.half.wallPowerIntoW;
else,powers(end+1)=(thermistorNetPowerLimit(2)+limit)/2;end
for power=powers
    [r,report]=runThermistorBeadCase([],mode,[],[],struct('netHeatingPowerW',power));
    t=r.thermistorBead.thermal;w=r.hybridWall;scale=power/limit;
    release=-sum(rows.areaM2.*w.heatFluxWm2);
    assert(report.fullSurfaceComputed && height(w)==5216 && isequal(w.faceID,rows.faceID));
    assert(isequal(w.pressurePa,rows.pressurePa) && isequal(w.heatFluxWm2,scale*rows.heatFluxWm2));
    assert(max(abs(w.nearWallCellTemperatureK-(293.15+scale*(rows.nearWallCellTemperatureK-293.15))))<1e-12);
    assert(abs(release-power)<=1e-11*power && abs(report.wholeSurfaceHeatPowerW+power)<=1e-11*power);
    assert(~isfield(r,'hybridWallTemperatureK') && ~isfield(r.hybridCondition,'wallTemperatureK') && ...
        all(r.hybridSolvedWallTemperatureK==293.15+12*scale) && ...
        r.evidenceIdentity.physicalCondition.netHeatingPowerW==power && ~isfield(r.thermistorBead,'experiment'));
    entry=struct('powerW',power,'solvedWallK',t.wallTemperatureK,'releaseW',release, ...
        'balanceErrorW',release-power,'areaM2',sum(rows.areaM2),'faces',height(w), ...
        'peakOutwardFluxWm2',max(-w.heatFluxWm2),'conductanceWK',t.conductanceWK);
    if cells==2 && power==-direct.half.wallPowerIntoW
        proof.nativeInverseErrorK=t.wallTemperatureK-direct.half.wallTemperatureK;
        assert(abs(proof.nativeInverseErrorK)<.001);
    end
    if cells==4
        if power<=thermistorNetPowerLimit(2)
            coarse=runThermistorBeadCase([],'archive',[],[],struct('netHeatingPowerW',power));
            assert(r.thermistorBead.baselineSolvedWallTemperatureK==coarse.hybridSolvedWallTemperatureK(1));
            assert(isequal(r.beadAngularDelta.heatFluxWm2,w.heatFluxWm2-coarse.hybridWall.heatFluxWm2));
            if power==.004,proof.samePowerWallDifferenceK=t.wallTemperatureK-coarse.hybridSolvedWallTemperatureK(1);end
        else
            assert(~isfield(r,'beadAngularDelta') && isfield(r.thermistorBead,'comparisonUnavailable'));
            data=buildAerothermalViewerData(r);
            assert(contains(data.fields(strcmp({data.fields.key},'heatFlux')).stateNote,'15°'));
        end
    end
    proof.cases=[proof.cases;entry];
end
b=thermistorWorkbenchBoundary();b.netHeatingPowerW=0;reject(b);
b=rmfield(b,{'netHeatingPowerW','wallTemperatureK'});reject(b);
for invalid={NaN,Inf,-eps,limit+eps(limit),[0 1],complex(1,1)}
    reject(struct('netHeatingPowerW',invalid{1}));
end
old=runThermistorBeadCase([],mode,[],[],299.15);
assert(all(old.hybridWallTemperatureK==299.15) && ~isfield(old,'hybridSolvedWallTemperatureK'));

% 模块二：走真实侧栏、WebGL、面元探针和 CSV，核对冻结身份及上游失效。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.beadModel);
if cells==4,choose(fig.UserData.controls.beadCells,4);end
press(fig.UserData.controls.generate);assert(~isempty(fig.UserData.workflowState.mesh),fig.UserData.status.Text);
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;choose(c.beadCondition,'power');
assert(strcmp(c.beadWall.Enable,'off') && strcmp(c.beadPower.Enable,'on') && c.beadPower.Value==.004);
press(c.confirm);assert(fig.UserData.workflowState.currentStep==4,fig.UserData.status.Text);
press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
view=fig.UserData.webUI;start=tic;
while isempty(view.html.UserData.rendered) && toc(start)<35,drawnow;pause(.1);end
gl=view.html.UserData.rendered;
assert(~isempty(gl) && gl.webglError==0 && gl.faces==5216 && gl.visiblePixelSamples>0);
v=view.evidenceView;data=v.getData();solved=s.result.hybridSolvedWallTemperatureK(1);
tw=data.fields(strcmp({data.fields.key},'wallSolved'));input=data.fields(strcmp({data.fields.key},'wallInput'));
assert(all(tw.values==solved & tw.status==1) && all(isnan(input.values) & input.status==0));
assert(all(contains(string(tw.methodNames),'热平衡')) && ...
    strcmp(data.resultCondition.wallInputType,'netHeatingPower') && data.resultCondition.wallKnownFaces==0 && ...
    data.resultCondition.solvedWallTemperatureK==solved && contains(data.resultCondition.caption,'净加热功率 4 mW'));
v.selectField('wallSolved');v.inspectFace(5204,[]);assert(v.getState().selectedFaceID==5204);
v.exportTo(fullfile(out,'temperature.csv'));csv=readtable(fullfile(out,'temperature.csv'));
assert(isequal(csv.face_id,rows.faceID+1) && max(abs(csv.value-solved))<1e-12 && all(string(csv.availability)=="研究计算"));
pairPath=fullfile(out,'wall.csv');receipt=v.exportWallPairTo(pairPath);
csv=readtable(pairPath,'Delimiter',',','NumHeaderLines',0,'ReadVariableNames',true,'Encoding','UTF-8');
assert(isequal(csv.face_id,rows.faceID+1) && max(abs(csv.pressure*1000-rows.pressurePa))<1e-9 && ...
    abs(sum(csv.area_m2.*csv.heat_flux)*1000+.004)<1e-13 && receipt.joint.pressureHeatFluxQualifiedArea==0);
exportAerothermalWeb(s.result,fullfile(out,'workbench'));
proof.workbench=struct('faces',gl.faces,'visiblePixelSamples',gl.visiblePixelSamples,'webglError',gl.webglError, ...
    'probeFaceOneBased',5204,'temperatureK',solved,'conditionSHA256',receipt.identity.conditionSHA256, ...
    'pressureHeatVerifiedM2',receipt.joint.pressureHeatFluxQualifiedArea, ...
    'wallCsvSHA256',cm011FileSHA256(pairPath),'temperatureCsvSHA256',cm011FileSHA256(fullfile(out,'temperature.csv')));
press(fig.UserData.workflowButtons(3));c=fig.UserData.controls;choose(c.beadPower,0);
s=fig.UserData.workflowState;assert(isempty(s.result) && isempty(fig.UserData.webUI) && ~s.boundaryConfirmed);
choose(c.beadCondition,'published');press(c.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && all(s.result.hybridWallTemperatureK==305.15));
v=fig.UserData.webUI.evidenceView;data=v.getData();
assert(all(data.fields(strcmp({data.fields.key},'wallSolved')).status==0) && all(data.fields(strcmp({data.fields.key},'wallInput')).status==2));
proof.boundaryInvalidationChecked=true;
fid=fopen(fullfile(out,'result.json'),'w','n','UTF-8');assert(fid>=0);file=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
    function reject(value)
        % 无效输入必须在分配输出目录和启动冷求解器前拒绝。
        dest=fullfile(out,'invalid');
        try,runThermistorBeadCase(dest,strrep(mode,'archive','solve'),[],[],value);error('RAH:Test:Accepted','无效边界通过。');
        catch ex,assert(strcmp(ex.identifier,'RAH:Bead:Boundary') && ~isfolder(dest),ex.message);end
    end
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
