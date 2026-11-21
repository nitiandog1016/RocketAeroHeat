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
    case 'paired3'
        archive=fullfile(caseDir,'paired3.zip');
        manifest=jsondecode(fileread(fullfile(caseDir,'paired3.json')));expected=manifest.archiveSHA256;
    case {'leads3','midbody3','leadsMidbody3','leadsThird3','outer30','sidewind3','sidewindOuter30','sidewindTipSpacing','sidewindAngular96','windPoint12'}
        contract=jsondecode(fileread(fullfile(caseDir,[angularCells '.json'])));
        archive=fullfile(caseDir,[angularCells '.zip']);expected=contract.archiveSHA256;
    case 'startup6'
        archive=fullfile(caseDir,'startup_6ms.zip');
        expected='25bdbfa8f11b6a12a6488b81b41002e7c3bc9808f8490f4b63987ee782c54b61';
    case 'halfsteps'
        archive=fullfile(caseDir,'startup_6ms_halfsteps.zip');
        expected='9a1ce5800c51a19b845877b83c735cd56d082a873c04024194f1f666eb1f84f3';
    case 'nearWall6'
        archive=fullfile(caseDir,'startup_6ms_nearwall56.zip');
        expected='9bd48452a2516bd42b4aab6e644d101a883146d0e6e225322e109f66d15dab9e';
    otherwise
        error('RAH:Bead:Angle','请选择已登记的 3 m/s 扇区或 6 m/s 全周启动快照。');
end
assert(strcmp(cm011FileSHA256(archive),expected),'RAH:Bead:Archive','B 型珠体的归档指纹不符。');
folder=tempname;mkdir(folder);cleanup=onCleanup(@()removeExtracted(folder));
unzip(archive,folder);
end

function removeExtracted(folder)
% 模块四：释放本次从已校验数据包展开的临时文件。
for name={'case.json','body.json','wall.csv','summary.json','matched_temperature.csv','comparison.json','historical.csv','shared.csv','flow.csv','flow.json','six_case.json','six_wall.csv','six_summary.json'}
    path=fullfile(folder,name{1});if isfile(path),delete(path);end
end
rmdir(folder);
end
