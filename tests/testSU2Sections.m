function testSU2Sections(out,capture)
%TESTSU2SECTIONS 从主工作台重开真实已有场，核对任意X/点击联动/未知断开和导出。
% 模块一：复用v107小型真实快照，不运行数值引擎，不把旧状态当新三维解。
% 模块二：截图、网页绘制回执、原面CSV和状态版本绑定同一当前快照。
if nargin<2,capture=false;end
assert(~isfolder(out));mkdir(out);source='docs/validation/v107/replay_seed/mapped_research.mat';
input=load('docs/validation/v107/replay_seed/input_snapshot.mat','frozen');
fig=rocketAeroHeatWorkbench(input.frozen.cfg);guard=onCleanup(@()close(fig));
fig.UserData.actions.loadNative();fig.UserData.actions.generateMesh();fig.UserData.actions.selectStep(3);
fig.UserData.actions.confirmBoundary();fig.UserData.actions.readResearch(source);
before=fig.UserData.workflowState;ui=fig.UserData.researchUI;assert(isempty(fig.UserData.engineJob));
records=cell(0,1);keys={'pressureRaw','heatFluxRaw','fluidNodeTemperatureMeanRaw','wallInput','wallSolved'};
locations=[.21 .7 1.25];
for k=1:numel(keys)
    ui.select.Value=keys{k};ui.select.ValueChangedFcn(ui.select,[]);
    for j=1:numel(locations)
        x=locations(j);ui.setSectionX(x);section=waitSection(ui,x);
        assert(~isempty(section.samples) && numel(section.components)==1 && section.components.closed);
        if strcmp(keys{k},'wallSolved'),assert(all(arrayfun(@(s)isempty(s.value),section.samples)));
        elseif strcmp(keys{k},'wallInput'),assert(all([section.samples.value]==300) && all([section.samples.status]==2));
        else
            values=before.researchResult.(keys{k});ids=[section.samples.face]+1;
            assert(isequal([section.samples.value]',values(ids)) && all([section.samples.status]==1));
        end
        file=fullfile(out,sprintf('%s_x%02d.csv',keys{k},j));ui.exportSectionTo(file);
        records{end+1}=section; %#ok<AGROW>
    end
end
% 点击使用真实空间位置的X；程序化面探针使用该原三角面的中心，不限柱段。
ui.select.Value='pressureRaw';ui.select.ValueChangedFcn(ui.select,[]);
ui.inspectFace(1235);x=mean(before.mesh.V(before.mesh.F(1235,:),1));picked=waitSection(ui,x);
assert(contains(ui.probe.Value{1},'1235') && picked.x==x);
ui.low.Value=0;ui.high.Value=100000;ui.apply.ButtonPushedFcn(ui.apply,[]);waitSection(ui,x);
assert(isequal(ui.geometry.html.UserData.payload.faceField.range,[0 100000]));
write(fullfile(out,'section_records.json'),records);
write(fullfile(out,'case_for_node.json'),struct('V',before.mesh.V,'F',before.mesh.F-1,'pressureRaw',before.researchResult.pressureRaw, ...
    'heatFluxRaw',before.researchResult.heatFluxRaw,'fluidNodeTemperatureMeanRaw',before.researchResult.fluidNodeTemperatureMeanRaw));
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('figureName',fig.Name));fprintf('V109_SECTION_READY\n');
    start=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(start)<60,drawnow;pause(.1);end
end
after=fig.UserData.workflowState;assert(after.revision==before.revision && after.solveId==before.solveId && isempty(fig.UserData.engineJob));
write(fullfile(out,'report.json'),struct('mainWorkbench',true,'solverRun',false,'sectionsChecked',numel(records), ...
    'unknownPreserved',true,'givenTemperatureSeparate',true,'probeSelectsOriginalFaceX',true,'inputVersionUnchanged',true, ...
    'newQualifiedArea',0,'rendered',ui.geometry.html.UserData.rendered));
fprintf('Main workbench actual WebGL/15 sections/current CSV/unknown and given roles passed; no numerical solver.\n');
end
function s=waitSection(ui,x)
start=tic;
while toc(start)<40
    drawnow;s=ui.getSection();web=ui.geometry.html.UserData;
    if ~isempty(s) && s.revision==web.revision && abs(s.x-x)<1e-12 && ~isempty(web.rendered)
        assert(web.rendered.webglError==0 && web.rendered.visiblePixelSamples>0 && web.rendered.sectionSegments==numel(s.samples));return;
    end
    pause(.1);
end
error('当前截面未完成真实网页绘制。');
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
