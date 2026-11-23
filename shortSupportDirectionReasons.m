function wall=shortSupportDirectionReasons(wall,flagNames,grids)
% 模块二：各原因允许并存；云图用缺值显示无法判别的面，数值分量保留自身的定义域。
names={'axis','noDirection','targetEdge','sourceEdge','smallShear','smallSpeed','disagreement','mixedShear','mixedSpeed'};
labels={'球面迎背流轴点邻面','来流切向投影退化','原面接触分区棱边','参与插值的壁片接触棱边', ...
    '壁剪切分量在舍入带内','首层速度分量在舍入带内','剪切与速度异号','不同壁片的剪切分量异号','不同壁片的速度分量异号'};
assert(isequal(flagNames(:),names(:)));
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
