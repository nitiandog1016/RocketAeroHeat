function p=shortAzimuthPayload(az)
% 模块用途：按同一实际方位传递轻量曲线，未解读数在离线 JSON 中保留 null。
t=az.points;states=["sphere","tube","cap","seam","mesh-disagreement","missing","multiple"];
[ok,id]=ismember(string(t.status),states);assert(all(ok));
% 把连续的缺项区间按极角合并，几何边界点本身的身份仍保存在逐点表。
masked=t(t.kind=="interval" & t.status~="sphere",:);gaps=zeros(0,3);
for j=1:height(masked)
    if ~isempty(gaps) && gaps(end,1)==masked.angleDeg(j) && abs(gaps(end,3)-masked.loDeg(j))<1e-8
        gaps(end,3)=masked.hiDeg(j);
    else
        gaps(end+1,:)=[masked.angleDeg(j),masked.loDeg(j),masked.hiDeg(j)]; %#ok<AGROW>
    end
end
p=struct('kind','short65-azimuth','frame',az.summary.frameSHA256,'timeMs',az.summary.elapsedS*1e3,'stepMs',az.summary.dtS*1e3, ...
    'angles',[az.summary.stations.angleDeg],'experiment',[az.summary.stations.experimentNu], ...
    'gaps',gaps,'rear',az.summary.rear, ...
    'states',{cellstr(states)},'rows',[t.angleDeg,t.phiDeg,t.nu,t.qInWm2,t.wallK,t.fileFaceID,t.referenceFaceID,id,t.faceCount]);
end
