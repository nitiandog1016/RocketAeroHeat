function report=testShockNormalBodyIntersections(folder,inletFile)
%TESTSHOCKNORMALBODYINTERSECTIONS 在真实箭体原面元上校核激波内法线几何接口。
% 模块一：只读取已求候选激波和原V/F，不重算跃迁或修改已有物理解。
% 模块二：保存原面号/重心权，核对刚体变换及内部/不相交射线拒绝。
% 模块三：此结果是待匹配厚度的几何目标，不是已求激波层厚度或热流。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupEngineeringShockLayer;
assert(~exist(folder,'dir'),'ShockLayer:Overwrite','几何接口验证必须写新目录。');
source=load(inletFile,'geometry','cases');V=source.geometry.V;F=source.geometry.F;
inlet=source.cases{2};positions=cell2mat(cellfun(@(s)s.position',inlet.states,'UniformOutput',false));
directions=cell2mat(cellfun(@(s)s.normal',inlet.states,'UniformOutput',false));
mapping=intersectShockNormalsWithBody(V,F,positions,directions);
assert(all(mapping.hit)&&all(mapping.distanceM>0)&&all(mapping.originalFaceID>=1));
assert(mapping.inputFaceCount==19752&&~mapping.flowValuesAssignedToBody&&~mapping.bodyMatchedByFlowSolution);
for j=1:size(positions,1)
    id=mapping.originalFaceID(j);b=mapping.barycentric(j,:);
    assert(abs(sum(b)-1)<1e-14&&min(b)>=-mapping.barycentricTolerance);
    assert(norm(b*V(F(id,:),:)-mapping.bodyPosition(j,:))<1e-13);
end
angle=.51;rotation=[1 0 0;0 cos(angle) -sin(angle);0 sin(angle) cos(angle)];shift=[2.3 -.4 1.1];
changed=intersectShockNormalsWithBody(V*rotation'+shift,F,positions*rotation'+shift,directions*rotation');
assert(all(changed.hit)&&max(abs(changed.distanceM-mapping.distanceM))<1e-11);
assert(max(vecnorm(changed.bodyPosition-(mapping.bodyPosition*rotation'+shift),2,2))<1e-11);
for j=1:size(positions,1)
    if numel(mapping.coincidentOriginalFaces{j})==1
        assert(changed.originalFaceID(j)==mapping.originalFaceID(j));
    end
end
miss=intersectShockNormalsWithBody(V,F, ...
    [[max(V(:,1))+.1 2 2];[mean(V(:,1)) 0 0]],repmat([1 0 0],2,1));
assert(~any(miss.hit)&&all(isnan(miss.originalFaceID))&&all(isnan(miss.distanceM)));
assert(all(isnan(miss.bodyPosition),'all')&&all(strlength(string(miss.missingReason))>0));
refused(@()intersectShockNormalsWithBody(V,F(2:end,:),positions,directions),'ShockLayer:BodyMesh');
refused(@()intersectShockNormalsWithBody(V,F(:,[1 3 2]),positions,directions),'ShockLayer:BodyMesh');
report=struct('status','tested-original-rocket-shock-normal-geometric-interface-not-flow-matching', ...
    'originalFaces',19752,'normalRays',size(positions,1),'hitCount',nnz(mapping.hit), ...
    'rayMissAndInsideRemainNaN',true,'closedTopologyAndOutwardOrientationChecked',true, ...
    'originalGeometryChanged',false,'shockLayerThicknessSolved',false, ...
    'flowValuesAssignedToBody',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'shock_body_geometry_evidence.mat'), ...
    'mapping','changed','miss','rotation','shift','report','inletFile','-v7');
fid=fopen(fullfile(folder,'shock_body_geometry_report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(report);
end
function refused(action,expected)
name='';try,action();catch ex,name=ex.identifier;end
assert(strcmp(name,expected),'ShockLayer:TestRejection','实际拒绝%s，期望%s。',name,expected);
end
