function limit=thermistorNetPowerLimit(angularCells,cfg)
%THERMISTORNETPOWERLIMIT 从所选网格的已核验参照场读取 12 K 温差下整珠放热功率。
wired=nargin>1 && isfield(cfg,'thermistorVariant') && strcmp(cfg.thermistorVariant,'leads');
archive=angularCells;
if wired
    assert(angularCells==4 && isfield(cfg,'thermistorRefineLeads') && cfg.thermistorRefineLeads, ...
        'RAH:Lead:Boundary','整件净对流功率使用带线中段加密网。');
    archive='leadsMidbody3';
end
[folder,cleanup]=openThermistorArchive(archive); %#ok<NASGU>
stats=jsondecode(fileread(fullfile(folder,'summary.json')));
if wired,limit=stats.totalOutwardW;else,limit=-stats.powerW;end
assert(isfinite(limit) && limit>0,'RAH:Bead:Power','参照场的放热功率须为正。');
end
