function [p,archive,seal]=readShortSupportPowerCase(job,propertyModel)
% 模块一：区分版本归档和本次安装生成的量热任务，核对两端功率、目标、候选及实际复算。
root=fullfile(fileparts(mfilename('fullpath')),'cases','short_sphere');
if nargin<2,propertyModel='fixedFilm';end
film=strcmp(propertyModel,'candidateFilm');assert(film || strcmp(propertyModel,'fixedFilm'));
stem='calorimetry65';expected='94097a483b00142615e272ec6da6baa7f061d106cb9f530f03346c7ff3dcf642';schema='short65-sphere-power-v1';
if film,stem='calorimetryFilm65';expected='cca04cdb314df2ae3f1860b40162209f5c3a460f34263497f933832679b06fb1';schema='short65-sphere-power-film-v1';end
if isempty(job)
    file=fullfile(root,[stem '.json']);archive=fullfile(root,[stem '.zip']);
else
    file=fullfile(job,'completed.json');archive=fullfile(job,'power65.zip');
end
assert(isfile(file) && isfile(archive),'RAH:Short:Job','所选目录缺少银球量热工况。');
seal=cm011FileSHA256(file);
if isempty(job),assert(strcmp(seal,expected),'RAH:Short:Job','量热归档与发布版本不符。');end
p=jsondecode(fileread(file));[b,reference]=shortSupportImportBoundary('power',14);
assert(strcmp(p.schema,schema) && strcmp(p.source.caseSHA256,b.caseSHA256) && ...
    p.faces==32480 && abs(p.areaM2-.000700423853595)<1e-15 && p.experimentalJointAreaM2==0 && ...
    p.targetSphereW==b.targetSphereW && strcmp(cm011FileSHA256(archive),p.archiveSHA256),'RAH:Short:Job','量热任务身份不符。');
goal=p.targetSphereW;
if film
    validateFilm(p,b,reference.condition);dt=p.balance.excessK;
else
    lo=p.low.earlyRegionW(1);hi=p.high.earlyRegionW(1);
    dt=p.low.excessK+(goal-lo)*(p.high.excessK-p.low.excessK)/(hi-lo);
    assert(p.low.excessK==20 && p.high.excessK==40 && lo<goal && goal<hi && abs(hi/2/lo-1)<.002 && abs(dt-p.balance.excessK)<1e-10);
end
assert(abs(b.airK+dt-p.temperatureK)<1e-10 && ...
    p.temperatureK==p.balance.temperatureK && abs(p.balance.earlyRegionW(1)/goal-1)<.0002 && ...
    abs(p.balance.lateRegionW(1)/goal-1)<.0002,'RAH:Short:Job','球面功率或独立壁温终点不符。');
assert(strcmp(cm011FileSHA256(file),seal),'RAH:Short:Job','读取期间量热任务回执已改变。');
end

function validateFilm(p,b,c)
% 模块二：独立重建每次括区、父场和候选物性，热功率取第三次及后续的实际端点。
keys=cellstr(string(p.order(:)));assert(numel(keys)>=3 && numel(keys)<=5 && isequal(keys(1:2),{'low';'high'}),'RAH:Short:Job','候选顺序错误。');
lo=p.trials.low;hi=p.trials.high;previous=p.source.frameSHA256;
assert(lo.excessK==20 && hi.excessK==40 && isequal(p.balance,p.trials.(p.balanceName)) && strcmp(keys{end},p.balanceName));
for j=1:numel(keys)
    row=p.trials.(keys{j});f=row.filmTransport;
    temp=(b.airK+row.temperatureK)/2;k=.0241*(temp/273.15)^1.5*(273.15+194)/(temp+194);
    assert(strcmp(row.parentFrame,previous) && strcmp(f.model,'air-conductivity-at-candidate-film') && ...
        abs(row.temperatureK-b.airK-row.excessK)<1e-10 && abs(f.filmK-temp)<1e-10 && f.wallK==row.temperatureK && f.airK==b.airK && ...
        abs(k/f.kWmK-1)<1e-12 && abs(c.rho*c.cp*c.speedMps*c.diameterM/k/f.pe-1)<1e-12 && ...
        abs(f.diffusionND*f.pe-1)<1e-12,'RAH:Short:Job','源帧或膜温物性不符。');
    previous=row.earlyFrame;
    if j<=2,continue;end
    low=lo.earlyRegionW(1);high=hi.earlyRegionW(1);
    x=lo.excessK+(p.targetSphereW-low)*(hi.excessK-lo.excessK)/(high-low);
    assert(low<p.targetSphereW && p.targetSphereW<high && abs(x-row.excessK)<1e-10 && strcmp(keys{j},sprintf('root%02d',j-2)));
    if row.earlyRegionW(1)<p.targetSphereW,lo=row;else,hi=row;end
end
end
