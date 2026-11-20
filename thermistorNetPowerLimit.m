function limit=thermistorNetPowerLimit(angularCells)
%THERMISTORNETPOWERLIMIT 从所选网格的已核验参照场读取 12 K 温差下整珠放热功率。
[folder,cleanup]=openThermistorArchive(angularCells); %#ok<NASGU>
stats=jsondecode(fileread(fullfile(folder,'summary.json')));
limit=-stats.powerW;
assert(isfinite(limit) && limit>0,'RAH:Bead:Power','参照场的放热功率须为正。');
end
