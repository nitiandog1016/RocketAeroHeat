function testSU2FlatPlateProperties(out)
%TESTSU2FLATPLATEPROPERTIES 验证可选物性与旧默认兼容，不启动流场求解。
% 模块一：默认R/Pr和旧配置字面值保持，防止严格重启动身份无意变化。
% 模块二：显式DNS参数能传入；非法值在建目录/启动前被拒绝。
assert(~isfolder(out));mkdir(out);
a=buildSU2FlatPlateCase(fullfile(out,'default'),struct('nx',24,'ny',24));
assert(a.R==287.05287 && a.Prandtl==.72);
text=fileread(fullfile(out,'default','case.cfg'));assert(contains(text,'PRANDTL_LAM= 0.72'));
b=buildSU2FlatPlateCase(fullfile(out,'matched'),struct('nx',24,'ny',24,'gasConstant',287,'Prandtl',.71));
assert(b.R==287 && b.Prandtl==.71);
for key={'gasConstant','Prandtl'}
    o=struct;o.(key{1})=-1;rejected=false;
    try,buildSU2FlatPlateCase(fullfile(out,['invalid_' key{1}]),o);catch,rejected=true;end
    assert(rejected && ~isfolder(fullfile(out,['invalid_' key{1}])));
end
fid=fopen(fullfile(out,'checks.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));
fprintf(fid,'{"defaultPreserved":true,"explicitProperties":true,"negativeGasConstantRejected":true,"negativePrandtlRejected":true,"nativeRuns":0}\n');
end
