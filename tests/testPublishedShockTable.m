function testPublishedShockTable()
% 模块用途：以手工转录的公开表格检查求解器，避免用同一算法产生期望值。
addpath(fileparts(fileparts(mfilename('fullpath'))));
r=verifyPublishedShockTable();
assert(r.passed,'NACA 1135 表格超出印刷精度：最大相对偏差 %.6g',r.maximumRelativeError);
assert(~r.experimentalValidation);
fprintf('NACA 1135: 20 values passed at half printed last-digit tolerance; max relative error %.6g.\n',r.maximumRelativeError);
end
