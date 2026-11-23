function [xyz,facets]=readSTLTextCorners(path)
% 按文本 STL 原始记录顺序读取十进制角点，用于精确坐标恢复及来源核对。
contents=fileread(path);
if nargout>1,facets=numel(regexp(contents,'(?mi)^\s*facet\s+normal\s','start'));end
tokens=regexp(contents,'(?mi)^\s*vertex\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)[^\S\n]*$','tokens');
if isempty(tokens),xyz=zeros(0,3);return;end
xyz=str2double(vertcat(tokens{:}));
end
