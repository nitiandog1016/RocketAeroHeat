function delta=readShortSupportTimeDifference(job,cases,branch,p,s)
% 模块一：差分文件随本次两条时间链封存，核对同一分支的早、晚帧与原面顺序。
if isempty(job),folder=fullfile(cases,'later1');else,folder=fullfile(job,'export');end
file=fullfile(folder,'timeDifference.zip');c=p.timeComparison;ref=c.(branch);
assert(strcmp(cm011FileSHA256(file),c.differenceSHA256) && strcmp(ref.laterFrame,s.frameSHA256) && ...
    s.elapsedSeconds==c.laterTimeS && c.laterTimeS==4*c.earlyTimeS,'RAH:Short:TimeDifference','同面差分包或物理时刻不符。');
scratch=tempname;mkdir(scratch);clean=onCleanup(@()discard(scratch)); %#ok<NASGU>
unzip(file,scratch);t=readtable(fullfile(scratch,[branch '.csv']));
assert(isequal(t.referenceFaceID,(1:s.faces)') && all(isfinite(t{:,2:end}),'all'));
delta=struct('earlyTimeS',c.earlyTimeS,'laterTimeS',c.laterTimeS,'dtS',s.dtSeconds,'baseFrame',ref.earlyFrame, ...
    'terminalFrame',ref.laterFrame,'powerChangeW',ref.wallChangeW,'fluxChangeWm2',t.deltaFluxInWm2, ...
    'gaugeChangePa',t.deltaGaugePa,'firstChangeK',t.deltaFirstK);
% 模块二：两个步长共用以零为中心的色标，直接比较差值的大小和符号。
a=c.full.fields;b=c.half.fields;
delta.sharedMax=struct('timeFluxChange',max(a.fluxInWm2.maxAbs,b.fluxInWm2.maxAbs), ...
    'timeGaugeChange',max(a.gaugePa.maxAbs,b.gaugePa.maxAbs),'timeFirstChange',max(a.firstK.maxAbs,b.firstK.maxAbs));
end
function discard(folder)
for name={'full.csv','half.csv'}
    file=fullfile(folder,name{1});if isfile(file),delete(file);end
end
rmdir(folder);
end
