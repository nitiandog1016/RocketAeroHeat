function result=readSU2SurfaceResult(filename)
%READSU2SURFACERESULT 读取SU2实际输出的ASCII legacy VTK壁面量。
% 模块一：严格核对坐标、单元连接和节点字段，拒绝二进制及不支持格式。
% 模块二：保留原始压力/温度/热流；单位和热流方向须由算例合同核验。
% 模块三：Temperature为流体节点状态，与给定Tw可能有离散偏差，不是壁温解。
fid=fopen(filename,'r');assert(fid>=0,'RocketAeroHeat:SU2:Read','无法读取壁面文件。');
guard=onCleanup(@()fclose(fid));
assert(startsWith(fgetl(fid),'# vtk DataFile Version'),'RocketAeroHeat:SU2:VTK','仅支持legacy VTK。');
fgetl(fid);assert(strcmp(strtrim(fgetl(fid)),'ASCII'),'RocketAeroHeat:SU2:VTK','仅支持ASCII输出。');
assert(strcmp(strtrim(fgetl(fid)),'DATASET UNSTRUCTURED_GRID'),'RocketAeroHeat:SU2:VTK','需要非结构网格。');
expect(fid,'POINTS');n=number(fid);fscanf(fid,'%s',1);V=reshape(numbers(fid,3*n),3,[])';
expect(fid,'CELLS');nc=number(fid);count=number(fid);packed=numbers(fid,count);
connectivity=zeros(count-nc,1);offsets=zeros(nc,1);cursor=1;next=1;
for j=1:nc
    width=packed(cursor);assert(width>=2 && width<=4 && width==fix(width));
    assert(cursor+width<=count);connectivity(next:next+width-1)=packed(cursor+(1:width));
    cursor=cursor+width+1;next=next+width;offsets(j)=next-1;
end
assert(cursor==count+1 && next==numel(connectivity)+1);
expect(fid,'CELL_TYPES');assert(number(fid)==nc);types=numbers(fid,nc);
widths=diff([0;offsets]);assert(all((types==3 & widths==2)|(types==5 & widths==3)|(types==9 & widths==4)));
expect(fid,'POINT_DATA');assert(number(fid)==n);values=struct;
while ~feof(fid)
    token=fscanf(fid,'%s',1);if isempty(token),break;end
    name=fscanf(fid,'%s',1);fscanf(fid,'%s',1);
    if strcmp(token,'SCALARS')
        components=number(fid);expect(fid,'LOOKUP_TABLE');expect(fid,'default');
    elseif strcmp(token,'VECTORS'),components=3;
    else,error('RocketAeroHeat:SU2:VTK','不支持的字段类型：%s。',token);end
    assert(components>=1 && components<=3);
    key=matlab.lang.makeValidName(name);assert(~isfield(values,key),'重复字段。');
    values.(key)=reshape(numbers(fid,n*components),components,[])';
end
assert(all(isfield(values,{'Pressure','Temperature','Heat_Flux'})), ...
    'RocketAeroHeat:SU2:Fields','缺少p/T/q，不能以保守量代替热流。');
assert(all([numel(values.Pressure),numel(values.Temperature),numel(values.Heat_Flux)]==n));
assert(all(isfinite(V(:))) && all(isfinite(values.Pressure) & values.Pressure>0 & ...
    isfinite(values.Temperature) & values.Temperature>0 & isfinite(values.Heat_Flux)));
assert(all(connectivity>=0 & connectivity<n & connectivity==fix(connectivity)));
result=struct('vertices',V,'connectivity',connectivity,'offsets',offsets,'cellTypes',types, ...
    'pressureRaw',values.Pressure,'temperatureRaw',values.Temperature,'heatFluxRaw',values.Heat_Flux, ...
    'arrays',values,'sourceFile',filename,'method','SU2 8.6.0 compressible Navier-Stokes wall output', ...
    'temperatureIdentity','raw-fluid-node-state; distinct-from-prescribed-wall-temperature', ...
    'unitConversionApplied',false,'heatFluxSignConfirmed',false,'productionQualified',false);
end

function expect(fid,wanted)
actual=fscanf(fid,'%s',1);assert(strcmp(actual,wanted),'RocketAeroHeat:SU2:VTK','预期%s，得到%s。',wanted,actual);
end
function value=number(fid)
value=numbers(fid,1);assert(isfinite(value) && value>=1 && value==fix(value));
end
function value=numbers(fid,count)
value=fscanf(fid,'%f',count);assert(numel(value)==count,'RocketAeroHeat:SU2:VTK','数据被截断或格式错误。');
end
