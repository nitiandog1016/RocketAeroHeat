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
% 模块二：各原因允许并存；云图用缺值显示无法判别的面，数值分量保留自身的定义域。
names={'axis','noDirection','targetEdge','sourceEdge','smallShear','smallSpeed','disagreement','mixedShear','mixedSpeed'};
labels={'球面迎背流轴点邻面','来流切向投影退化','原面接触分区棱边','参与插值的壁片接触棱边', ...
    '壁剪切分量在舍入带内','首层速度分量在舍入带内','剪切与速度异号','不同壁片的剪切分量异号','不同壁片的速度分量异号'};
assert(isequal(summary.flagNames(:),names(:)));
for grid=grids
    key=grid{1}; flags=wall.([key 'Flags']);status=wall.([key 'Status']);
    assert(all(ismember(status,[-1 0 1])) && isequal(flags==0,status~=0));
    [codes,~,index]=unique(flags);reasons=cell(size(codes));
    for i=1:numel(codes)
        mask=logical(bitget(uint16(codes(i)),1:numel(labels)));
        reasons{i}=strjoin(labels(mask),'、');
    end
    wall.([key 'Reason'])=reasons(index);
end
end

function discard(folder)
% 仅清除本次校验过的压缩表的临时展开文件。
file=fullfile(folder,'direction.csv');
if isfile(file),delete(file);end
rmdir(folder);
end
