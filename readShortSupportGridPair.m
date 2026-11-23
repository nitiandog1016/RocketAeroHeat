function pair=readShortSupportGridPair(folder)
% 绑定本次粗细网、各自独立终帧和同面比较小表。
file=fullfile(folder,'supportedNearWall65.json');
assert(strcmp(cm011FileSHA256(file),'d5c97ea022d968533a7e0b98c0b7cd599b14d6aa106860c575fedd8952df82d3'));
pair=jsondecode(fileread(file));
end
