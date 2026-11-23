function pair=readShortSupportAngularPair(folder)
% 模块用途：绑定 128 扇区的独立终场、周向交叠及 64 扇区对照源。
file=fullfile(folder,'supportedAngular65.json');
assert(strcmp(cm011FileSHA256(file),'08396677f33935fb5e40ee34658a26d788cb386b995d6a69cb85375b4cbcaeb8'));
pair=jsondecode(fileread(file));original=readShortSupportGridPair(folder);
assert(strcmp(pair.schema,'short-65-angular-64-128-v1') && pair.sectors==128 && pair.baseSectors==64 && ...
    strcmp(pair.coarseFrame,original.fineFrame) && strcmp(pair.coarseMesh,original.fineMesh) && ...
    strcmp(pair.coarseArchiveSHA256,original.fineArchiveSHA256) && strcmp(pair.bodySHA256,original.bodySHA256) && ...
    pair.faces==original.faces && strcmp(pair.geometry.publicTransferSHA256,pair.transferSHA256));
end
