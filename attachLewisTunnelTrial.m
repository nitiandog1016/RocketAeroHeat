function mesh=attachLewisTunnelTrial(mesh,scheme)
%ATTACHLEWISTUNNELTRIAL 给匹配的风洞工况附上已存三维诊断，主壁压与热流仍取原工程结果。
% 模块一：选择同一物理步的计算归档。沿用单参数调用时读取固定 CFL 对照。
if nargin<2,scheme='fixed';end
switch scheme
    case 'fixed'
        version='v279';name='固定 CFL';
    case 'adaptive'
        version='v280';name='自适应 CFL';
    otherwise
        error('RAH:Lewis:TrialScheme','无法识别圆柱试算方案。');
end
folder=fullfile(fileparts(mfilename('fullpath')),'docs','validation',version,'evidence');
receipt=jsondecode(fileread(fullfile(folder,'trial.json')));
assert(isfield(mesh,'lewisCondition') && ~receipt.passesNineResiduals && ~receipt.productionQualified, ...
    'RAH:Lewis:Trial','试算身份或工况不符合归档。');
names=fieldnames(receipt.condition);
for j=1:numel(names)
    ref=receipt.condition.(names{j});
    assert(isfield(mesh.lewisCondition,names{j}) && abs(mesh.lewisCondition.(names{j})-ref)<=4e-14*max(1,abs(ref)), ...
        'RAH:Lewis:TrialCondition','已存三维试算的边界与当前工况不同。');
end
% 模块二：核对边界、数据校验值与每个原面的位置、面积，避免载入其他工况或几何。
file=fullfile(folder,'trial_faces.csv');
assert(strcmpi(cm011FileSHA256(file),receipt.faceSHA256),'RAH:Lewis:Trial','试算逐面数据与来源不符。');
rows=readtable(file);f=mesh.F;v=mesh.V;
center=(v(f(:,1),:)+v(f(:,2),:)+v(f(:,3),:))/3;
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(isequal(rows.faceID,(1:size(f,1))')&&max(abs(rows{:,{'xM','yM','zM'}}-center),[],'all')<1e-13 && ...
    max(abs(area-rows.areaM2))<1e-16,'RAH:Lewis:TrialGeometry','原面编号、位置或面积不匹配。');
mesh.lewisTrial=struct('pressurePa',rows.trialPressurePa,'heatFluxWm2',rows.trialHeatFluxWm2, ...
    'sourceFrameSHA256',receipt.sourceFrameSHA256,'physicalDtS',receipt.physicalDtS, ...
    'innerIterations',receipt.innerIterations,'schemeName',name);
end
