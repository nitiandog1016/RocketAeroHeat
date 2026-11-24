function [b,pair]=shortSupportImportBoundary(branch,radiusD)
% 模块用途：从已验收的时间层和试验原表生成唯一的银球边界与物理时刻。
assert(any(strcmp(branch,{'full','half'})),'RAH:Short:Branch','请选择一个整步或两个半步。');
folder=fullfile(fileparts(mfilename('fullpath')),'cases','short_sphere');
if nargin<2,radiusD=14;end
pair=readShortSupportTimePair(folder,radiusD);file=fullfile(folder,'test65.json');
assert(strcmp(cm011FileSHA256(file),pair.caseSHA256));
s=jsondecode(fileread(file));c=pair.condition;e=pair.([branch 'Energy']);
b=struct('branch',branch,'velocityMps',[-c.speedMps 0 0],'airK',c.airK,'airPa',c.airPa, ...
    'reynolds',c.re,'angleDeg',s.angleDeg,'wallProfileK',(s.surfaceF-32)*5/9+273.15, ...
    'supportWallK',(s.surfaceF(s.angleDeg==90)-32)*5/9+273.15,'elapsedS',pair.elapsedS, ...
    'dtS',e.dt*c.diameterM/c.speedMps,'initialFrame',pair.seed,'caseSHA256',pair.caseSHA256);
end
