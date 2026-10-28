function testRocketTestGeometryView()
%TESTROCKETTESTGEOMETRYVIEW 核对离线几何导出的真实面元与明确物理身份。
% 模块一：使用两个不同解析几何，防止模型选择器重复显示同一个负载。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
a=buildCanonicalBody('coneCylinder',struct('resolution',16));
b=buildCanonicalBody('hemisphereCylinder',struct('resolution',20));
folder=tempname;mkdir(folder);cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
path=exportRocketTestGeometryView({a,b},folder,{'光顺对照','第二模型'});
assert(isfile(path));
script=fileread(fullfile(folder,'rocket-geometry-data.js'));
prefix='window.ROCKET_GEOMETRY_CASES=';
assert(startsWith(script,prefix)&&endsWith(script,';'));
payload=jsondecode(script(numel(prefix)+1:end-1));
assert(numel(payload)==2&&isequal(payload(1).vertices,a.V));
assert(isequal(payload(2).faces,b.F-1)&&payload(1).revision~=payload(2).revision);
assert(~isfield(payload,'pressure')&&~isfield(payload,'heatFlux'));
assert(contains(payload(1).coordinateNote,'尚无物理场'));
html=fileread(path);assert(contains(html,'fixtureChoice'));
for name={'geometry.js','geometry.css'}
    assert(isequal(fileread(fullfile(folder,name{1})),fileread(fullfile(root,'web',name{1}))));
end
% 模块二：非法编号与退化三角形必须明确拒绝，不能生成误导预览。
bad=a;bad.F(1,1)=size(a.V,1)+1;refused=false;
try,exportRocketTestGeometryView({bad},folder);catch,refused=true;end
assert(refused);
bad=a;bad.F(1,2)=bad.F(1,1);refused=false;
try,exportRocketTestGeometryView({bad},folder);catch err
    refused=strcmp(err.identifier,'RocketAeroHeat:RocketPreview:Degenerate');
end
assert(refused);
fprintf('Rocket offline geometry export identity and refusal checks passed.\n');
end
