function sub=readThermistorSubfaceFields(input,result)
% 模块一：将本次边界、源帧及对照绑定到实际的壁反力和首层节点。
b=input.thermistorImport;p=input.thermistorPrepared;id=result.evidenceIdentity;
key=id.sourceFrame;folder=id.runDirectory;stage='';baseStage='';
if ~isempty(folder)
    task=folder;
    if isfile(fullfile(task,'run','completed.json')),task=fullfile(task,'run');end
    done=jsondecode(fileread(fullfile(task,'completed.json')));stage=fullfile(task,done.stage);
    profile=fullfile(task,'prepared','profile.json');
    if isfile(profile),pr=jsondecode(fileread(profile));baseStage=pr.base;end
end
sources.main=struct('frame',key.sha256,'reaction',key.reactionSHA256,'directory',stage);
factor=1;controlFactor=1;
if isfield(result,'thermistorWire')
    m=result.thermistorWire;sources.control=struct('frame',m.baseFrame,'reaction','','directory',baseStage);
    if isfield(m,'thermometry')
        t=m.thermometry;factor=t.fieldFactor;controlFactor=t.scale;
        sources.original=struct('frame',t.referenceFrame,'reaction','','directory','');
    end
else
    factor=key.ratio;
end
ids=b.childIDs;req=struct('prepared',p.directory,'mesh',p.meshSHA256,'nodeMap',p.nodeMapSHA256, ...
    'body',p.bodySHA256,'parents',b.referenceFaceID(ids),'childID',ids,'bary',b.childBarycentric, ...
    'area',b.areaM2(ids),'sources',sources,'stlSHA256',b.source.sha256);
root=fileparts(mfilename('fullpath'));folder=tempname(fullfile(root,'outputs'));mkdir(folder);
source=fullfile(folder,'input.json');out=fullfile(folder,'subfaces.json');
fid=fopen(source,'w','n','UTF-8');assert(fid>=0);
closeFile=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(req));clear closeFile;
runCM011Python({fullfile(root,'scripts','thermistorSubfaces.py'),source,out},fullfile(folder,'read.log'),[],false);
raw=jsondecode(fileread(out));assert(isequal(raw.childID,ids) && strcmp(raw.stlSHA256,b.source.sha256));
% 模块二：沿用本次计算的温差比例和符号。壁温输入和首层温度各自保留取值位置。
main=raw.fields.main;air=result.hybridCondition.ambientTemperatureK;
sub=struct('id',ids,'gaugePa',main.gauge,'fluxWm2',main.flux*factor,'wallK',air+main.wall*factor, ...
    'firstK',air+main.first*factor,'firstXYZ',raw.firstXYZ+b.translationM,'sources',raw.sources, ...
    'shapeError',raw.shapeError,'overlaps',raw.overlaps,'receipt',out);
if isfield(raw.fields,'control'),sub.wireDelta=sub.fluxWm2-controlFactor*raw.fields.control.flux;end
if isfield(raw.fields,'original'),sub.thermometryDelta=sub.fluxWm2-raw.fields.original.flux;end
end
