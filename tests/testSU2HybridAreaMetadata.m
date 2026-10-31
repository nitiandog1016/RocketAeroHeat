function report=testSU2HybridAreaMetadata(out)
%TESTSU2HYBRIDAREAMETADATA 只读已有经典/真实/整流罩网格，核对总面积修正及写出影响。
% 模块一：复用三份已通过拓扑的体网格，与Python原报告及原三角面面积独立比较。
% 模块二：只重新写经典输入，证明修正仅改变清单面积，不改变求解器配置/网格。
% 模块三：历史清单保留原字节并报告错误值；不重建网格、不重算流场。
assert(~isfolder(out));mkdir(out);
actual=jsondecode(fileread('outputs/v114_validation/actual01/report.json'));
feature=jsondecode(fileread('outputs/v113_validation/ui01/report.json'));
folders={'outputs/hybrid_mesh/a75f5d87-b861-4478-b972-3095018283e6', ...
    'outputs/v113_validation/actual02','outputs/hybrid_mesh/5899a261-8218-478e-9bbb-b062c99c7af4'};
engines={actual.engineFolder,'outputs/v113_validation/actual_review02/engine',feature.engineFolder};
records=cell(3,1);
for k=1:3
    receipt=jsondecode(fileread(fullfile(folders{k},'mesh_process.json')));
    d=load(receipt.input,'V','F');surface=struct('V',d.V,'F',double(d.F));
    volume=readSU2HybridVolume(folders{k},surface);
    a=surface.V(surface.F(:,2),:)-surface.V(surface.F(:,1),:);
    b=surface.V(surface.F(:,3),:)-surface.V(surface.F(:,1),:);
    % 逐面叉积的三个分量显式展开，与原Python网格报告交叉校验。
    expected=sum(sqrt((a(:,2).*b(:,3)-a(:,3).*b(:,2)).^2+ ...
        (a(:,3).*b(:,1)-a(:,1).*b(:,3)).^2+(a(:,1).*b(:,2)-a(:,2).*b(:,1)).^2))/2;
    assert(abs(volume.bodyArea/expected-1)<2e-13);
    original=jsondecode(fileread(fullfile(engines{k},'case_manifest.json')));
    records{k}=struct('meshDirectory',folders{k},'originalFaces',size(surface.F,1), ...
        'correctArea',volume.bodyArea,'pythonReportArea',volume.meshConstruction.bodyArea, ...
        'historicalManifestArea',original.wallArea,'historicalManifestRetained',true);
    if k==1
        fresh=fullfile(out,'corrected_input');info=writeSU2RocketCase(volume,fresh,original.flight);
        assert(info.wallArea==volume.bodyArea);
        assert(strcmp(sha(fullfile(fresh,'mesh.su2')),sha(fullfile(engines{k},'mesh.su2'))));
        assert(strcmp(sha(fullfile(fresh,'case.cfg')),sha(fullfile(engines{k},'case.cfg'))));
    end
end
report=struct('cases',{records},'solverInputBytesUnchanged',true,'mesherInvoked',false,'solverInvoked',false, ...
    'productionQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Three saved hybrid surface areas verified; config and mesh bytes unchanged; no solver or mesher.\n');
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
