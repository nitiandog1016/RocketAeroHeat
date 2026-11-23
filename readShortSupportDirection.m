function [wall,summary]=readShortSupportDirection(folder,pair,geometry,angular)
% 模块一：核对两套终帧、几何与传递表，再读取同面剪切、首层速度及拒判标志。
if nargin<4,angular=false;end
if angular
    summary=pair;archive=fullfile(folder,'supportedAngularDirection65.zip');digest=summary.directionArchiveSHA256;
    grids={'angular'};
else
    file=fullfile(folder,'supportedDirection65.json');
    assert(strcmp(cm011FileSHA256(file),'ca1dcaf627465ca2bb9e5467dc31c799ea0f71629ef597f773c39f03ded13af9'));
    summary=jsondecode(fileread(file));
    assert(strcmp(summary.coarse.frame,pair.coarseFrame) && strcmp(summary.fine.frame,pair.fineFrame) && ...
        strcmp(summary.caseSHA256,pair.caseSHA256) && strcmp(summary.bodySHA256,pair.bodySHA256) && ...
        strcmp(summary.transferSHA256,pair.transferSHA256) && summary.faces==pair.faces);
    archive=fullfile(folder,'supportedDirection65.zip');
    digest=summary.archiveSHA256;grids={'coarse','fine'};
end
assert(strcmp(cm011FileSHA256(archive),digest));
temp=tempname;mkdir(temp);done=onCleanup(@()discard(temp)); %#ok<NASGU>
unzip(archive,temp);wall=readtable(fullfile(temp,'direction.csv'));
assert(isequal(wall.faceID,(1:summary.faces)') && isequal(wall.regionID,geometry.regionID));
if angular
    suffix={'ShearPa','ShearMagnitudePa','SpeedMps','SpeedMagnitudeMps','Status','Flags'};
    wall.Properties.VariableNames(3:8)=strcat('angular',suffix);
end
wall=shortSupportDirectionReasons(wall,summary.flagNames,grids);
end

function discard(folder)
% 仅清除本次校验过的压缩表的临时展开文件。
file=fullfile(folder,'direction.csv');
if isfile(file),delete(file);end
rmdir(folder);
end
