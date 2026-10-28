function report=testNasaBoattailBenchmark(folder)
%TESTNASABOATTAILBENCHMARK 核对原图数据身份、坐标、闭合几何及原字节篡改拒绝。
% 本测试证明来源读取和几何实现，不把来源校核或几何收敛当成气动热实验精度。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'tests'));
assert(nargin==1&&~exist(folder,'dir'),'NasaBoattailTest:Overwrite','须指定新的来源测试目录。');mkdir(folder);
b=loadNasaBoattailBenchmark();assert(numel(b.points)==15&&~b.solverValidationComplete);
assert(max(abs([b.points.projectXFromFigureM]-(1.18*.0254-3.375*.0254*[b.points.sourceFigureXprimeOverD])))<2e-15);
source=fullfile(root,'tests','data','nasa_tnd2761');

%% 模块一：原始PDF、清单及读图结果篡改，分别保留失败副本和异常身份
labels={'raw','manifest','digitization'};rejections=cell(3,1);
for j=1:3
    trial=fullfile(folder,['bad_' labels{j}]);copyfile(source,trial);
    if j==1
        p=fullfile(trial,'raw','NASA_TN_D_2761_19650013004.pdf');fid=fopen(p,'r+');assert(fid>=0);fwrite(fid,uint8(0),'uint8');fclose(fid);
        expected='NasaBoattail:Source';
    elseif j==2
        p=fullfile(trial,'source_manifest.json');fid=fopen(p,'a');fprintf(fid,' ');fclose(fid);expected='NasaBoattail:Manifest';
    else
        p=fullfile(trial,'pressure_digitization.json');d=jsondecode(fileread(p));d.points(1).pressureRatio=0;
        fid=fopen(p,'w','n','UTF-8');fprintf(fid,'%s',jsonencode(d));fclose(fid);expected='NasaBoattail:Source';
    end
    actual='';try,loadNasaBoattailBenchmark(trial);catch ex,actual=ex.identifier;end
    assert(strcmp(actual,expected),'NasaBoattailTest:Rejection','篡改未按来源身份拒绝。');rejections{j}=actual;
end

%% 模块二：三组闭合整模型，独立解析体积/面积与周向细化，而非只检查有无三角形
rows=struct([]);
for angle=[6 12 18]
    errors=zeros(2,1);
    for k=1:2
        n=48*2^(k-1);[m,t]=buildNasaBoattailBody(angle,n);
        edges=[m.F(:,[1 2]);m.F(:,[2 3]);m.F(:,[3 1])];[~,~,id]=unique(sort(edges,2),'rows');
        assert(all(accumarray(id,1)==2)&&all(accumarray(id,sign(edges(:,2)-edges(:,1)))==0));
        assert(max(conncomp(graph(edges(:,1),edges(:,2))))==1);
        assert(all(abs(vecnorm(m.normal,2,2)-1)<2e-13)&&all(m.area>0));
        cylinder=strcmp(t.faceRegion,'cylinder');base=strcmp(t.faceRegion,'base');tail=strcmp(t.faceRegion,'boattail');
        assert(all(m.normal(base,1)<-1+2e-13)&&all(abs(m.normal(cylinder,1))<2e-13));
        assert(all(m.normal(tail,1)<0)&&all(m.center(tail,1)>0)&&all(m.center(tail,1)<t.boattailLengthM));
        R=3.375*.0254/2;Lbt=1.18*.0254;Ln=R/tand(15);Lc=14*.0254-Lbt;Rb=R-Lbt*tand(angle);
        volume=pi*(Lbt*(R*R+R*Rb+Rb*Rb)/3+R*R*Lc+R*R*Ln/3);
        assert(abs(volume-t.exactVolumeM3)<1e-16&&abs(t.facetVolumeM3/volume-sin(2*pi/n)/(2*pi/n))<2e-13);
        assert(abs(max(m.V(:,1))-(14*.0254+Ln))<1e-14&&min(m.V(:,1))==0);
        errors(k)=abs(t.facetAreaM2/t.exactAreaM2-1);
        row=struct('angleDeg',angle,'resolution',n,'faces',size(m.F,1),'closedOrientedConnected',true, ...
            'facetVolumeM3',t.facetVolumeM3,'analyticVolumeM3',volume,'relativeAreaError',errors(k), ...
            'stingIncluded',t.stingIncluded,'gritBandIncluded',t.gritBandGeometryIncluded);
        if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    end
    assert(errors(2)<errors(1)&&errors(2)<.002);
end
report=struct('status','passed-source-calibration-three-tamper-rejections-and-three-closed-nominal-geometries', ...
    'sourceManifestSHA256',b.sourceManifestSHA256,'experimentalSymbols',15,'geometry',rows, ...
    'rejections',{rejections},'figureCoordinatesNotForcedToOrificeIDs',true, ...
    'geometryOrSourceCheckIsSolverValidation',false,'experimentalFlowSolved',false, ...
    'engineProcessID',feature('getpid'),'newProductionSolvedAreaM2',0);
save(fullfile(folder,'source_geometry_evidence.mat'),'b','report');
fid=fopen(fullfile(folder,'source_geometry_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
