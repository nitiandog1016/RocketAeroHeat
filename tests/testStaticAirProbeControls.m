function testStaticAirProbeControls(source,out)
%TESTSTATICAIRPROBECONTROLS 在独立隐藏组件内读取表格、触发真实下载并比对字节。
assert(~isfolder(out));mkdir(out);
actual=fileread(fullfile(source,'index.html'));
js=fileread(fullfile(fileparts(mfilename('fullpath')),'staticAirProbeControls.js'));
page=strrep(actual,'</body>',['<script>' js '</script></body>']);
fig=uifigure('Visible','off');guard=onCleanup(@()delete(fig));
html=uihtml(fig,'HTMLSource',page);
t=tic;
while ~isstruct(html.Data) && toc(t)<20,drawnow;pause(.05);end
x=html.Data;assert(isstruct(x),'离线页面未回传点击结果。');
if ~x.ok,error('StaticAir:DOM','%s',x.error);end
assert(strcmp(x.selected,'62') && contains(x.name,'62') && contains(x.method,'STABLE_PIECEWISE'));
assert(contains(x.position,'内部节点') && contains(x.sources,'stable/terminal.csv'));
assert(contains(x.sources,'stable/case.cfg') && contains(x.sources,'stable/initialization.json'));
expected={'流体节点温度','300.123','K','初始化回读';'流体节点压力','100000','Pa','初始化回读'; ...
    '给定壁温','—','K','不适用';'壁面热流','—','W/m²','不适用';'求解固体壁温','—','K','不适用'};
cells=string(jsondecode(x.cells));
assert(isequal(cells(:),reshape(string(expected.'),[],1)),'实际显示的物理量、数值、单位或身份不一致。');
assert(strcmp(x.filename,'static-air-stable-terminal-node-62.csv'));
fid=fopen(fullfile(source,x.filename),'rb');assert(fid>=0);reader=onCleanup(@()fclose(fid));
bytes=fread(fid,inf,'*uint8');assert(isequal(bytes,uint8(x.bytes(:))));
x=rmfield(x,'bytes');x.matchedDownloadBytes=numel(bytes);x.externalBrowserStarted=false;
fidOut=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fidOut>=0);writer=onCleanup(@()fclose(fidOut));
fprintf(fidOut,'%s\n',jsonencode(x));
end
