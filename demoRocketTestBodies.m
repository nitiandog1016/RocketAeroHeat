function summary = demoRocketTestBodies(outputDir)
%DEMOROCKETTESTBODIES 导出完整箭体及突起配对几何，供工作台和保真审核使用。
% 文件含 MAT、ASCII STL、JSON 与审核摘要；均为几何夹具，不含虚构气动热结果。
if nargin<1,outputDir=fullfile(fileparts(mfilename('fullpath')),'outputs','rocket_fixture_2026_10_17');end
if ~exist(outputDir,'dir'),mkdir(outputDir);end
names={'smooth','fairing','ridge','two_protrusions'};records=cell(0,1);
for k=1:numel(names)
    cfg=struct('resolution',48);
    if k==1,cfg.protrusionType='none';elseif k<4,cfg.protrusionType=names{k};
    else
        cfg.protrusions(1)=struct('type','fairing','height',0.0192,'length',0.24, ...
            'centerFraction',0.35,'azimuthDeg',35,'halfWidthDeg',22);
        cfg.protrusions(2)=struct('type','ridge','height',0.0144,'length',0.30, ...
            'centerFraction',0.68,'azimuthDeg',205,'halfWidthDeg',10);
    end
    [mesh,meta]=buildRocketTestBody(cfg);stem=fullfile(outputDir,['rocket_',names{k}]);
    save([stem,'.mat'],'mesh','meta','cfg','-v7');writeStl([stem,'.stl'],mesh,names{k});
    data=struct('name',names{k},'V',mesh.V,'F',mesh.F,'faceRegion',{meta.faceRegion}, ...
        'faceBaseRegion',{meta.faceBaseRegion},'faceFeatureMask',meta.faceFeatureMask, ...
        'features',meta.features,'audit',meta.geometryAudit,'geometryOnly',true,'physicalValidation',false);
    fid=fopen([stem,'.json'],'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'%s',jsonencode(data));clear cleanup;
    records{end+1}=struct('name',names{k},'vertices',size(mesh.V,1),'faces',size(mesh.F,1), ...
        'audit',meta.geometryAudit,'exactVolume',meta.exactVolume, ...
        'matFile',[stem,'.mat'],'stlFile',[stem,'.stl'],'jsonFile',[stem,'.json']); %#ok<AGROW>
end
summary=struct('status','generated','geometryOnly',true,'physicalValidation',false,'cases',[records{:}]);
fid=fopen(fullfile(outputDir,'geometry_summary.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(summary));
fprintf('Generated %d full rocket geometry fixtures in %s\n',numel(names),outputDir);
end

function writeStl(path,m,name)
% 每一面使用真实单位法向，保持单一闭合曲面的三角面拓扑。
fid=fopen(path,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'solid %s\n',name);
for j=1:size(m.F,1)
    fprintf(fid,' facet normal %.16g %.16g %.16g\n  outer loop\n',m.normal(j,:));
    for k=1:3,fprintf(fid,'   vertex %.16g %.16g %.16g\n',m.V(m.F(j,k),:));end
    fprintf(fid,'  endloop\n endfacet\n');
end
fprintf(fid,'endsolid %s\n',name);
end
