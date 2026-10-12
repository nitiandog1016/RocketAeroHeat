function grid=buildLayeredWallGrid(material)
%BUILDLAYEREDWALLGRID 建立从外表面到内表面的多层一维控制体网格。
% material.layers 为结构数组或元胞数组；每层有 thickness、density、
% specificHeat、thermalConductivity，可选 propertyTable 和 thermalNodes。
% contactResistance(i) 为第 i 与 i+1 层之间的面接触热阻 [m² K/W]。
% 零热阻共用界面节点，两侧半控制体分别储热；正热阻保留同深度的
% 两个温度节点，以无质量的热阻连接。不得用平均材料替代实际分层。
if ~isfield(material,'layers') || isempty(material.layers)
    error('RocketAeroHeat:LayeredWall:Layers','需要非空 material.layers。');
end
layers=material.layers;
if isfield(material,'propertyTable') && ~isempty(material.propertyTable)
    error('RocketAeroHeat:LayeredWall:Ambiguous','多层模式不能同时提供顶层 propertyTable。');
end
if isstruct(layers), layers=num2cell(layers); end
if ~iscell(layers) || ~isvector(layers) || numel(layers)>10
    error('RocketAeroHeat:LayeredWall:Layers','layers 须为不超过 10 层的结构/元胞向量。');
end
layers=layers(:)'; nl=numel(layers);
resistance=zeros(1,nl-1);
if isfield(material,'contactResistance')
    validateattributes(material.contactResistance,{'numeric'},{'real','finite','nonnegative'});
    resistance=double(material.contactResistance(:)');
end
if ~isreal(resistance) || numel(resistance)~=nl-1 || any(~isfinite(resistance)|resistance<0)
    error('RocketAeroHeat:LayeredWall:Contact','接触热阻须为层数减一的非负有限向量。');
end
x=[]; ids=cell(1,nl); mass=ids; links=ids; widths=zeros(1,nl);
contactLinks=zeros(1,nl-1); interfaceNodes=zeros(nl-1,2);
offset=0;
for j=1:nl
    layer=layers{j};
    if ~isstruct(layer) || ~isscalar(layer)
        error('RocketAeroHeat:LayeredWall:Layer','每层须为标量材料结构。');
    end
    for key={'thickness','density'}
        validateattributes(layer.(key{1}),{'numeric'},{'real','finite','scalar','positive'});
    end
    n=8; if isfield(layer,'thermalNodes'), n=layer.thermalNodes; end
    validateattributes(n,{'numeric'},{'real','finite','scalar','integer','>=',3,'<=',100});
    layers{j}.thermalNodes=n;
    wallMaterialProperties(layer,NaN); % 同时检查常数/温度表格式，不猜测物性。
    dx=layer.thickness/(n-1); widths(j)=dx;
    begin=numel(x)+1;
    if j>1 && resistance(j-1)==0, begin=begin-1; end
    ids{j}=begin:begin+n-1; links{j}=ids{j}(1:end-1);
    x(ids{j})=offset+linspace(0,layer.thickness,n); %#ok<AGROW>
    mass{j}=layer.density*dx*[.5 ones(1,n-2) .5];
    if j>1
        interfaceNodes(j-1,:)=[ids{j-1}(end) ids{j}(1)];
        if resistance(j-1)>0, contactLinks(j-1)=ids{j-1}(end); end
    end
    offset=offset+layer.thickness;
end
if numel(x)>100
    error('RocketAeroHeat:LayeredWall:Nodes','多层总温度节点数超过 100，请减少每层节点数。');
end
lo=zeros(1,numel(x)); hi=inf(size(lo));
for j=1:nl
    p=wallMaterialProperties(layers{j},NaN);
    lo(ids{j})=max(lo(ids{j}),p.temperatureRange(1));
    hi(ids{j})=min(hi(ids{j}),p.temperatureRange(2));
end
if any(lo>=hi)
    error('RocketAeroHeat:LayeredWall:Range','完美接触界面的两种材料温度适用范围没有重叠。');
end
grid=struct('layers',{layers},'layerNodes',{ids},'layerLinks',{links}, ...
    'layerMass',{mass},'spacing',widths,'depth',x,'nodeCount',numel(x), ...
    'contactResistance',resistance,'contactLinks',contactLinks, ...
    'interfaceNodes',interfaceNodes,'temperatureMinimum',lo,'temperatureMaximum',hi, ...
    'method','逐层节点有限体积；零热阻共温度节点，有限热阻双温度节点');
end
