function path=exportRocketTestGeometryView(meshes,folder,captions,surfaceFlows)
%EXPORTROCKETTESTGEOMETRYVIEW 导出可切换对照模型的离线三维几何页。
% 模块用途：复用正式工作台的几何WebGL管线，检查自绘箭体及突起。
% 只输出真实V/F，不生成压力、热流或温度；可选surfaceFlows叠加已保存的条件流线。
% 流线附带来源和非生产身份，不把输入p/T变成几何云图中的求解场。
% meshes 为非空cell，每项为含V/F的网格；captions为等长中文标题cell。
% 示例：exportRocketTestGeometryView({m0,m1},'outputs/rocket_view', ...
%       {'光顺基线','带整流罩箭体'});

%% 模块一：核验完整网格和对照标题，拒绝非法数值及三角形编号
if nargin<3,captions=arrayfun(@(k)sprintf('箭体几何 %d',k), ...
        1:numel(meshes),'UniformOutput',false);end
if nargin<4,surfaceFlows={};end
assert(iscell(surfaceFlows)&&(isempty(surfaceFlows)||numel(surfaceFlows)==numel(meshes)), ...
    'RocketAeroHeat:RocketPreview:Flows','可选流线负载须与模型数一致。');
if ~iscell(meshes)||isempty(meshes)||~iscell(captions)||numel(meshes)~=numel(captions)
    error('RocketAeroHeat:RocketPreview:Cases','模型与标题必须为等长非空cell。');
end
validateattributes(folder,{'char','string'},{'nonempty'});
if isstring(folder)&&~isscalar(folder)
    error('RocketAeroHeat:RocketPreview:Folder','输出目录必须为单一文本。');
end
payloads=cell(numel(meshes),1);
for k=1:numel(meshes)
    mesh=meshes{k};label=captions{k};
    if ~isstruct(mesh)||~isscalar(mesh)||~isfield(mesh,'V')||~isfield(mesh,'F')
        error('RocketAeroHeat:RocketPreview:Mesh','每个模型必须含V/F。');
    end
    validateattributes(mesh.V,{'numeric'},{'real','finite','nonempty','ncols',3});
    validateattributes(mesh.F,{'numeric'},{'real','finite','nonempty','ncols',3, ...
        'integer','>=',1,'<=',size(mesh.V,1)});
    if ~((ischar(label)&&isrow(label))||(isstring(label)&&isscalar(label)))
        error('RocketAeroHeat:RocketPreview:Caption','标题必须为单行文本。');
    end
    if max(max(mesh.V,[],1)-min(mesh.V,[],1))<=0
        error('RocketAeroHeat:RocketPreview:Extent','模型尺寸必须为正。');
    end
    a=cross(mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,1),:), ...
        mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,1),:),2);
    if any(sum(a.^2,2)<=0)
        error('RocketAeroHeat:RocketPreview:Degenerate','不能输出退化三角形。');
    end
    payloads{k}=struct('kind','aerothermal-geometry','schemaVersion',1, ...
        'vertices',double(mesh.V),'faces',double(mesh.F)-1, ...
        'caption',char(label),'revision',k,'showWire',false,'resetCamera',true, ...
        'coordinateNote','几何测试，尚无物理场 · SI / m · 主轴X，鼻尖+X · 突起不得视为已解区域');
    if ~isempty(surfaceFlows)
        flow=surfaceFlows{k};
        assert(isstruct(flow)&&isscalar(flow)&&isfield(flow,'kind')&& ...
            strcmp(flow.kind,'pressure-surface-flow-research')&&isfield(flow,'productionQualified')&& ...
            isequal(flow.productionQualified,false)&&isfield(flow,'curves'), ...
            'RocketAeroHeat:RocketPreview:Flows','流线须具有明确研究身份和曲线。');
        curves=flow.curves;if isstruct(curves),curves=num2cell(curves);end
        assert(iscell(curves));
        for j=1:numel(curves)
            c=curves{j};validateattributes(c.points,{'double'},{'real','finite','ncols',3});
            validateattributes(c.normals,{'double'},{'real','finite','size',size(c.points)});
            assert(all(vecnorm(c.normals,2,2)>0));
        end
        flow.curves=curves;payloads{k}.surfaceFlowLines=flow;
        % 模块用途：默认科研视角属于模型，隐藏流线查看网格时也保持同一方向。
        for cameraKey={'initialCameraAzimuthRad','initialCameraElevationRad'}
            key=cameraKey{1};
            if isfield(flow,key)
                validateattributes(flow.(key),{'double'},{'real','finite','scalar'});
                payloads{k}.(key)=flow.(key);
            end
        end
        payloads{k}.coordinateNote='SI / m · X鼻尖正向 · 条件流线数值研究；压力/热状态为输入，热流未求解';
        if isfield(flow,'thermodynamicStateIsConditionallyRestored')&& ...
                isequal(flow.thermodynamicStateIsConditionallyRestored,true)
            payloads{k}.coordinateNote='SI / m · X鼻尖正向 · 压力为流线输入，Te/rho/U为同物性条件反算；Tr/Tw/热流未求';
        end
    end
end

%% 模块二：复制既有WebGL查看器，不另建不一致的几何显示算法
if ~exist(folder,'dir'),mkdir(folder);end
root=fileparts(mfilename('fullpath'));
for name={'geometry.html','geometry.css','geometry.js'}
    copyfile(fullfile(root,'web',name{1}),fullfile(folder,name{1}));
end
source=fileread(fullfile(folder,'geometry.html'));
source=strrep(source,'<div class="controls">', ...
    ['<div class="controls"><label>对照模型 <select id="fixtureChoice" ' ...
     'aria-label="选择光顺或带突起箭体"></select></label>']);
source=strrep(source,'</body>', ...
    '<script src="rocket-geometry-data.js"></script><script src="rocket-geometry-bootstrap.js"></script></body>');
writeUtf8(fullfile(folder,'index.html'),source);
writeUtf8(fullfile(folder,'rocket-geometry-data.js'), ...
    ['window.ROCKET_GEOMETRY_CASES=' jsonencode(payloads) ';']);

%% 模块三：最小本地桥接，记录真实WebGL回执并保持几何与物理结果分离
% 使用与MATLAB uihtml相同的setup入口。确认消息只记录，不再次触发加载。
bootstrap=strjoin({ ...
    '''use strict'';', ...
    '(function(){', ...
    ' const cases=window.ROCKET_GEOMETRY_CASES,select=document.getElementById("fixtureChoice");', ...
    ' const listeners=[];let value=null;', ...
    ' const bridge={addEventListener:(name,fn)=>{if(name==="DataChanged")listeners.push(fn);}};', ...
    ' Object.defineProperty(bridge,"Data",{get:()=>value,set:v=>{value=v;window.ROCKET_GEOMETRY_RECEIPT=v;}});', ...
    ' cases.forEach((c,i)=>{const o=document.createElement("option");o.value=i;o.textContent=c.caption;select.append(o);});', ...
    ' function load(){value=cases[Number(select.value)];listeners.forEach(fn=>fn());}', ...
    ' window.setup(bridge);select.addEventListener("change",load);load();', ...
    '})();'},newline);
writeUtf8(fullfile(folder,'rocket-geometry-bootstrap.js'),bootstrap);
path=fullfile(folder,'index.html');
end

function writeUtf8(path,value)
% 模块四：完整UTF-8文件写入，保证中文说明可以离线打开
fid=fopen(path,'w','n','UTF-8');
assert(fid>=0,'RocketAeroHeat:RocketPreview:Write','不能写入文件：%s',path);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',value);
end
