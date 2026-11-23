function check_cold_ui(out)
% 模块一：从实际侧栏创建独立冷态任务，验证计算、端点检查和三维页自动接续。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
folder=fullfile(root,'outputs');old=dir(fullfile(folder,'short65_supported_*'));
fig=rocketAeroHeatWorkbench();closeFig=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.shortSphereMode.Value='supportSolve';c.shortSphere.ButtonPushedFcn(c.shortSphere,[]);drawnow;
v=fig.UserData.shortSphereReview;assert(isstruct(v)&&isfield(v,'web'),'%s',fig.UserData.status.Text);
jobs=dir(fullfile(folder,'short65_supported_*'));jobs=jobs([jobs.isdir] & ~ismember({jobs.name},{old.name}));
assert(numel(jobs)==1);job=fullfile(folder,jobs.name);
ack=jsondecode(fileread(fullfile(job,'completed.json')));s=v.report;
assert(strcmp(ack.frameSHA256,s.frameSHA256) && s.faces==32480 && s.thermal.passed && s.coverage.fullSurfaceComputed);
html=v.web.html;t=tic;while isempty(html.UserData.rendered)&&toc(t)<55,drawnow;pause(.15);end
gl=html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==s.faces);
proof=struct('job',job,'frame',s.frameSHA256,'faces',s.faces,'pixels',gl.visiblePixelSamples, ...
    'gapFraction',s.thermal.early.physicalGapFraction,'archive',ack.archiveSHA256);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>0);finish=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
