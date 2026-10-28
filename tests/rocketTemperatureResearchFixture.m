function file=rocketTemperatureResearchFixture(name,legacyRelativePath)
%ROCKETTEMPERATURERESEARCHFIXTURE 优先使用Git仓库内小型验证输入，旧本地结果路径保留兼容。
% 模块一：默认研究入口无需解压大型工程包或访问被忽略的整箭输出目录。
% 模块二：随库文件与原已测试输入按SHA逐字节一致，不改变制造场或数值验收要求。
root=fileparts(fileparts(mfilename('fullpath')));
file=fullfile(root,'tests','data','temperature_research',name);
if ~isfile(file),file=fullfile(root,legacyRelativePath);end
assert(isfile(file),'TemperatureResearch:Fixture','缺少小型验证输入%s；请检查仓库测试数据。',name);
end
