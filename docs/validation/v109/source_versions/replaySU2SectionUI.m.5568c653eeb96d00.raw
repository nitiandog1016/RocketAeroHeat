function replaySU2SectionUI(out)
%REPLAYSU2SECTIONUI 在最终主窗口核对说明与截面交互；复用v107快照，不求解流场。
assert(~isfolder(out));mkdir(out);
input=load('docs/validation/v107/replay_seed/input_snapshot.mat','frozen');
fig=rocketAeroHeatWorkbench(input.frozen.cfg);guard=onCleanup(@()close(fig));
fig.UserData.actions.loadNative();fig.UserData.actions.generateMesh();fig.UserData.actions.selectStep(3);
fig.UserData.actions.confirmBoundary();fig.UserData.actions.readResearch('docs/validation/v107/replay_seed/mapped_research.mat');
ui=fig.UserData.researchUI;ui.setSectionX(.7);start=tic;
while isempty(ui.getSection()) && toc(start)<40,drawnow;pause(.1);end
assert(~isempty(ui.getSection()));
labels=findall(fig.UserData.editorPanel,'Type','uilabel');assert(any(arrayfun(@(x)contains(x.Text,'点击三维表面或拖动X滑块'),labels)));
write(fullfile(out,'ready_for_capture.json'),struct('figureName',fig.Name));fprintf('V109_FINAL_SECTION_READY\n');
start=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(start)<60,drawnow;pause(.1);end
ui.exportSectionTo(fullfile(out,'final_current_section.csv'));
write(fullfile(out,'report.json'),struct('finalSidebarTextChecked',true,'noSolverRun',true, ...
    'rendered',ui.geometry.html.UserData.rendered,'section',ui.getSection(),'newQualifiedArea',0));
fprintf('Final main-window section replay passed.\n');
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
