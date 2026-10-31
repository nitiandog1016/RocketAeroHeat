function report=studySU2DNSPlate(out)
%STUDYSU2DNSPLATE 运行公开M6冷壁条件下的内嵌SST发展平板。
% 模块一：匹配来流速度、压力、温度、壁温和R/Pr，保留-X流向。
% 模块二：显式记录均匀全湍流入口假设，只运行一个有界任务。
% 模块三：读取真实终态并独立对比相同Re_theta；不声称复现DNS入口或箭体资格。
assert(~isfolder(out));mkdir(out);folder=fullfile(out,'plate');
root=fileparts(fileparts(mfilename('fullpath')));python='D:/anacondapython/python.exe';
o=struct('nx',80,'ny',80,'length',.35,'upstreamLength',.02,'height',.04, ...
    'firstCellHeight',.5e-6,'Mach',869.1/sqrt(1.4*287*55.2), ...
    'temperature',55.2,'pressure',.044*287*55.2,'wallTemperature',97.5, ...
    'gasConstant',287,'Prandtl',.71,'iterations',10000,'outputEvery',5000);
info=buildSU2FlatPlateCase(folder,o);assert(info.R==287 && info.Prandtl==.71);
[code,text]=system(sprintf('"%s" -X utf8 "%s" --folder "%s"', ...
    python,fullfile(root,'scripts','prepareSU2DNSPlate.py'),folder));assert(code==0,text);
runEmbeddedSU2Case(folder,struct('pythonExecutable',python,'threads',2,'maximumSeconds',720));
[code,text]=system(sprintf('"%s" -X utf8 "%s" --folder "%s" --out "%s"', ...
    python,fullfile(root,'scripts','compareSU2DNSPlate.py'),folder,fullfile(out,'comparison')));
assert(code==0,text);report=reviewSU2DNSPlate(out,fullfile(out,'comparison','comparison.json'));
end
