function mesh=attachThermistorLeadNearWall(mesh,runFolder,outDir)
% 模块一：已算归档由保存的同帧原面表恢复，新算从自己的速度和壁剪切输出提取。
root=fileparts(mfilename('fullpath'));
if isempty(runFolder) && isfield(mesh,'thermistorMidbody')
    % 光滑加密归档携带同一终帧导出的首层速度、壁剪切和全表面流向表。
    stats=jsondecode(fileread(fullfile(outDir,'flow.json')));
    csv=fullfile(outDir,'flow.csv');digest=stats.tableSHA256;
elseif isempty(runFolder)
    proof=fullfile(root,'docs','validation','v297','evidence');
    stats=jsondecode(fileread(fullfile(proof,'nearwall.json')));stats=stats.fine;
    target=tempname;mkdir(target);guard=onCleanup(@()clearTable(target)); %#ok<NASGU>
    files=gunzip(fullfile(proof,'nearwall_faces.csv.gz'),target);csv=files{1};
    digest='2d92881abecb315b0ac2f10a0cf92f231ec86a07c3716db8324649e6bcc43c4e';
else
    if nargin<3,outDir=fullfile(runFolder,'surface');end
    runCM011Python({fullfile(root,'scripts','thermistorLeadNearWall.py'),runFolder,'--out',outDir}, ...
        fullfile(outDir,'nearwall.log'));
    stats=jsondecode(fileread(fullfile(outDir,'flow.json')));csv=fullfile(outDir,'flow.csv');digest=stats.tableSHA256;
    assert(strcmp(cm011FileSHA256(fullfile(runFolder,'surface.vtk')),stats.surfaceVTKSHA256));
end
% 模块二：逐面核对标号、流向判据、面积和来源，以实际首层距离和原生摩擦系数估算 y+。
w=readtable(csv);n=size(mesh.F,1);
assert(strcmp(cm011FileSHA256(csv),digest) && isequal(w.face_id,(1:n)') && ...
    strcmp(stats.sourceFrameSHA256,mesh.evidenceIdentity.sourceFrame.sha256) && ...
    all(isfinite(w.downstream_Cf) & isfinite(w.downstream_first_U_over_Uinf) & isfinite(w.estimated_first_yplus)) && ...
    all(w.estimated_first_yplus>=0) && isequal(logical(w.reverse_both),w.downstream_Cf<0 & w.downstream_first_U_over_Uinf<0) && ...
    nnz(w.reverse_both)==stats.reverseFaces,'RAH:Lead:Flow','流向或近壁数据与求解帧不符。');
mesh.thermistorLeadFlow=struct('summary',stats,'downstreamCf',w.downstream_Cf, ...
    'downstreamU',w.downstream_first_U_over_Uinf,'estimatedYPlus',w.estimated_first_yplus,'reverse',w.reverse_both);
end

function clearTable(folder)
file=fullfile(folder,'nearwall_faces.csv');if isfile(file),delete(file);end
if isfolder(folder),rmdir(folder);end
end
