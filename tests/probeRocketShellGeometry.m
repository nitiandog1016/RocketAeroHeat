function report=probeRocketShellGeometry(folder,options)
%PROBEROCKETSHELLGEOMETRY 记录真实箭体层网格的失败候选，不改变候选几何或接受判据。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
assert(~exist(folder,'dir'),'RocketShellProbe:Overwrite','不得覆盖失败或成功的独立试验。');mkdir(folder);
p=fullfile(root,'outputs','v56_validation','state_validity','state_validity_actual_rocket.mat');
if ~exist(p,'file'),p=fullfile(root,'validation','inherited_v56','v56','state_validity','state_validity_actual_rocket.mat');end
d=load(p,'result');options.observer=@observe;captured=false;
try
    domain=buildRocketTriangularShellDomain(d.result.V,d.result.F,options);
    report=struct('status','geometry-built-after-original-checks','quality',domain.quality,'captured',captured);
    save(fullfile(folder,'probe_accepted_domain.mat'),'domain','-v7');
catch ex
    report=struct('status','rejected-original-geometry-candidate','identifier',ex.identifier,'message',ex.message,'captured',captured);
end
report.engineProcessID=feature('getpid');fid=fopen(fullfile(folder,'probe_summary.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function observe(probe)
        captured=true;save(fullfile(folder,'volume_candidate.mat'),'probe','-v7');
        nf=probe.sourceFaces;bad=reshape(probe.volume<=0,nf,[]);fprintf('SHELL_GEOMETRY_NEGATIVE_COUNTS_BY_LAYER=%s\n',mat2str(sum(bad,1)));
    end
end
