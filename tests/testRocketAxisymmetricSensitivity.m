function report=testRocketAxisymmetricSensitivity(folder)
%TESTROCKETAXISYMMETRICSENSITIVITY 单参数、几何合同和五量初猜的独立重建检查。
% 模块一：四档真实完整箭体的墙、源面、边界与体积度量不被静默更换。
% 模块二：逐格独立复算实际距离、供体权重和五量初猜，域外检查自由流。
% 模块三：未知类型/选项/未绑定根拒绝；这些检查不等于敏感性求解已完成。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v76_validation','factory_test');end
assert(~exist(folder,'dir'),'SensitivityTest:Overwrite','不能覆盖测试证据。');parent=createRocketAxisymmetricRootCase;rows=struct([]);
for kinds={'baseline','radial48','nearwall_half','farfield1p5'}
    kind=kinds{1};c=createRocketAxisymmetricSensitivityCase(struct('kind',kind));
    assert(isequaln(c,createRocketAxisymmetricSensitivityCase(c.options))&& ...
        isequaln(buildAxisymmetricFVGeometry(c.mesh.X,c.mesh.R),c.mesh));
    assert(isequaln(parent.surfaceMesh,c.surfaceMesh)&&isequal(parent.domain.wallPolylineXR,c.domain.wallPolylineXR)&& ...
        isequal(parent.boundary,c.boundary)&&isequal(parent.gas,c.gas));
    oldDomain=parent.domain.options;newDomain=c.domain.options;
    oldCollar=parent.domain.collar.options;newCollar=c.domain.collar.options;
    switch kind
        case 'radial48',oldDomain.radialLayers=48;
        case 'nearwall_half',oldCollar.firstNormalSpacing=1e-5;
        case 'farfield1p5',oldDomain.outerRadius=1.5*oldDomain.outerRadius;
    end
    assert(isequal(oldDomain,newDomain)&&isequal(oldCollar,newCollar));
    a=c.initialization.mapping;inside=a.insideParentDomain;ids=a.donorCellIDs;w=a.donorWeights;
    assert(all(w(inside,:)>=0,'all')&&all(sum(w(inside,:),2)==1)&&all(ids(inside,:)>=1,'all')&& ...
        ~a.stateClipped&&~a.conservativeRemappingClaimed&&~a.initialGuessIsBoundaryCondition);
    primitive=c.inputContract.freestreamPrimitiveRhoUxUrT;
    free=axisymmetricFVConservative(primitive,c.gas);free=[free,primitive(1)*c.inputContract.farfieldNuTildeInputM2s];
    expected=repmat(free,c.mesh.cellCount,1);
    expected(inside,:)=w(inside,1).*parent.Q(ids(inside,1),:)+w(inside,2).*parent.Q(ids(inside,2),:);
    assert(isequal(expected,c.Q));
    % 独立距离验证：内插格须在两个真实供体距离之间，流向条带编号一致。
    ni=c.mesh.shape(1);wall=parent.mesh.faces.centroid(parent.domain.wallBoundaryFaceID,:);
    for j=find(inside&~a.endpointConstantExtension)'
        i=mod(j-1,ni)+1;d=vecnorm(parent.mesh.volumeCentroid(ids(j,:),:)-wall(i,:),2,2);
        x=norm(c.mesh.volumeCentroid(j,:)-wall(i,:));t=(x-d(1))/(d(2)-d(1));
        assert(all(mod(ids(j,:)-1,ni)+1==i)&&w(j,2)==t&&t>=0&&t<=1);
    end
    if strcmp(kind,'baseline'),assert(isequal(c.Q,parent.Q));end
    if strcmp(kind,'farfield1p5'),assert(any(~inside));else,assert(all(inside));end
    assert(all(c.mesh.volumePerRadian>0)&&all(c.domain.wallNormalDistanceM>0)&& ...
        isequal(sum(c.Q.*c.mesh.volumePerRadian,1),c.initialization.initialConservedIntegralPerRadian));
    row=struct('kind',kind,'cellCount',c.mesh.cellCount,'interpolatedCells',nnz(inside), ...
        'outsideParentCells',nnz(~inside),'actualWallDistanceM',[min(c.domain.wallNormalDistanceM),max(c.domain.wallNormalDistanceM)]);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
reject(@()createRocketAxisymmetricSensitivityCase(struct('kind','unknown')),'RocketSensitivity:Kind');
reject(@()createRocketAxisymmetricSensitivityCase(struct('kind','radial48','alpha',20)),'RocketSensitivity:Options');
reject(@()createRocketAxisymmetricSensitivityCase(struct('rootEvidenceSHA256','bad')),'RocketSensitivity:Evidence');
report=struct('status','passed-real-geometry-single-parameter-and-independent-five-quantity-initial-guess-checks', ...
    'cases',rows,'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'factory_test_evidence.mat'),'report');
fid=fopen(fullfile(folder,'factory_test_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('%s\n',jsonencode(report));
end
function reject(action,id)
try,action();error('SensitivityTest:ExpectedRefusal','非法输入未被拒绝。');catch ex,assert(strcmp(ex.identifier,id));end
end
