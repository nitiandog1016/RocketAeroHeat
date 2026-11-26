function check_point_pick(out)
% 模块用途：在真实 uihtml 页面驱动生产鼠标处理函数，核对同像素的前后点和背景。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
web=fullfile(root,'web');pageDir=fullfile(out,'pick-page');mkdir(pageDir);
copyfile(fullfile(web,'geometry.js'),pageDir);copyfile(fullfile(web,'geometry.css'),pageDir);
body=fileread(fullfile(web,'geometry.html'));
shim=fileread(fullfile(fileparts(mfilename('fullpath')),'pointer_fixture.js'));
f=fopen(fullfile(pageDir,'pick.html'),'w','n','UTF-8');assert(f>=0);
fprintf(f,'%s',strrep(body,'</body>',['<script>' shim '</script></body>']));fclose(f);
fig=uifigure('Position',[100 100 780 580]);guard=onCleanup(@()delete(fig)); %#ok<NASGU>
grid=uigridlayout(fig,[1 1]);h=uihtml(grid);seen=[];
h.DataChangedFcn=@receive;h.HTMLSource=fullfile(pageDir,'pick.html');
wait('ready');
v=[0 -1 0;0 1 0;0 1 0;-2 0 0;2 0 0];
field=struct('context','published-samples','productionQualified',false,'selectedIndex',0, ...
    'name','测试点','unit','1','range',[-1 4],'values',(-1:3)');
h.Data=struct('kind','aerothermal-geometry','schemaVersion',1,'revision',1,'vertices',v,'faces',[], ...
    'pointField',field,'resetCamera',true,'initialCameraAzimuthRad',pi/2,'initialCameraElevationRad',0);
wait('geometry-rendered');assert(seen.sampleCount==5 && seen.faces==0 && seen.webglError==0);
h.Data=struct('kind','fixture-click','x',.5,'y',.5);wait('geometry-point-probe');
assert(seen.revision==1 && seen.sampleIndex==1);picked=seen.sampleIndex+1;
h.Data=struct('kind','fixture-click','x',.14,'y',.14);wait('fixture-background');
f=fopen(fullfile(out,'picking.json'),'w');assert(f>=0);
fprintf(f,'%s\n',jsonencode(struct('frontmostSampleID',picked,'equalDepthKeptFirst',true,'backgroundIgnored',true)));fclose(f);
    function receive(src,~)
        if isstruct(src.Data) && isfield(src.Data,'kind'),seen=src.Data;end
    end
    function wait(kind)
        t=tic;
        while toc(t)<35
            drawnow;
            if isstruct(seen) && strcmp(seen.kind,kind),return;end
            if isstruct(seen) && strcmp(seen.kind,'geometry-error'),error('%s',seen.message);end
            pause(.05);
        end
        error('等待事件 %s 超时，最近收到 %s。',kind,jsonencode(seen));
    end
end
