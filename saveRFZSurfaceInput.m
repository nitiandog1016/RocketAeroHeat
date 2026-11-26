function saveRFZSurfaceInput(file,input)
% 模块用途：在目标目录完成写入并回读后，再替换保存的工况文件。
folder=fileparts(file);if isempty(folder),folder=pwd;end
tmp=[tempname(folder) '.json'];
clean=onCleanup(@()removeTemporary(tmp)); %#ok<NASGU>
[f,msg]=fopen(tmp,'w','n','UTF-8');assert(f>=0,'RAH:RFZ:Write','%s',msg);
closer=onCleanup(@()fclose(f));
fprintf(f,'%s\n',jsonencode(input,'PrettyPrint',true));clear closer;
assert(isequaln(input,jsondecode(fileread(tmp))),'RAH:RFZ:Write','写入的输入与当前确认内容不一致。');
[ok,msg]=movefile(tmp,file,'f');assert(ok,'RAH:RFZ:Write','%s',msg);
end
function removeTemporary(tmp)
if isfile(tmp),delete(tmp);end
end
