function report=verifyPublishedShockTable()
%VERIFYPUBLISHEDSHOCKTABLE 对比独立转录的 NACA 1135 表 II。
% 模块用途：检查正激波求解数值；这是理论表格校核，不是实验验证。
root=fileparts(mfilename('fullpath'));
source=jsondecode(fileread(fullfile(root,'docs','validation','naca1135_normal_shock.json')));
evalc('cfg=createCaseConfig();');
reference=source.rows(:,2:end);computed=zeros(size(reference));
for j=1:size(reference,1)
    cfg.flight.Mach=source.rows(j,1);fs=calcFreestream(cfg.flight);
    assert(abs(fs.gamma-source.gamma)<1e-12,'参考表与求解器 gamma 不一致');
    s=solveNormalShock(fs);
    computed(j,:)=[s.postShockMach,s.pressureRatio,s.densityRatio, ...
        s.temperatureRatio,s.totalPressureRatio];
end
% 印刷精度：允许半个末位单位，加浮点舍入裕度，不以宽百分比掩盖误差。
tolerance=0.5*source.lastDigitUnits+1e-12;
err=abs(computed-reference);
report=struct('source',source,'computed',computed,'reference',reference, ...
    'absoluteError',err,'tolerance',tolerance,'passed',all(err<=tolerance,'all'), ...
    'maximumRelativeError',max(err./abs(reference),[],'all'), ...
    'experimentalValidation',false);
end
