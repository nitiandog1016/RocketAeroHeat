function [folder,cleanup]=openThermistorArchive(angularCells)
%OPENTHERMISTORARCHIVE 校验并展开随工程发布的珠体数据；调用方持有临时文件的清理句柄。
caseDir=fullfile(fileparts(mfilename('fullpath')),'cases','thermistor_b');
switch angularCells
    case 2
        archive=fullfile(caseDir,'surface.zip');
        expected='11ba9c504f8df41a1be7c317d1047e77ab76800d3eae7b569d0b34333ffbb563';
    case 4
        archive=fullfile(caseDir,'surface_7p5.zip');
        expected='d02c8ac7fb41a7072e416e251efd2ccf8d19101ec388c98a46b3ed0a6e72449b';
    otherwise
        error('RAH:Bead:Angle','周向单元数须为 2 或 4。');
end
assert(strcmp(cm011FileSHA256(archive),expected),'RAH:Bead:Archive','B 型珠体的归档指纹不符。');
folder=tempname;mkdir(folder);cleanup=onCleanup(@()removeExtracted(folder));
unzip(archive,folder);
end

function removeExtracted(folder)
% 模块四：释放此调用独有的四个临时文件。新计算的现场和用户导出保存在项目中。
for name={'case.json','body.json','wall.csv','summary.json'}
    path=fullfile(folder,name{1});if isfile(path),delete(path);end
end
rmdir(folder);
end
