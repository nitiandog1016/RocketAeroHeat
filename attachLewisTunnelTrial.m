function mesh=attachLewisTunnelTrial(mesh)
%ATTACHLEWISTUNNELTRIAL 给匹配的风洞工况附上已存三维诊断，主壁压与热流仍取原工程结果。
folder=fullfile(fileparts(mfilename('fullpath')),'docs','validation','v279','evidence');
receipt=jsondecode(fileread(fullfile(folder,'trial.json')));
assert(isfield(mesh,'lewisCondition') && ~receipt.passesNineResiduals && ~receipt.productionQualified, ...
    'RAH:Lewis:Trial','试算身份或工况不符合归档。');
names=fieldnames(receipt.condition);
for j=1:numel(names)
    ref=receipt.condition.(names{j});
    assert(isfield(mesh.lewisCondition,names{j}) && abs(mesh.lewisCondition.(names{j})-ref)<=4e-14*max(1,abs(ref)), ...
        'RAH:Lewis:TrialCondition','已存三维试算的边界与当前工况不同。');
end
file=fullfile(folder,'trial_faces.csv');
assert(strcmpi(cm011FileSHA256(file),receipt.faceSHA256),'RAH:Lewis:Trial','试算逐面数据与来源不符。');
rows=readtable(file);f=mesh.F;v=mesh.V;
center=(v(f(:,1),:)+v(f(:,2),:)+v(f(:,3),:))/3;
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
assert(isequal(rows.faceID,(1:size(f,1))')&&max(abs(rows{:,{'xM','yM','zM'}}-center),[],'all')<1e-13 && ...
    max(abs(area-rows.areaM2))<1e-16,'RAH:Lewis:TrialGeometry','原面编号、位置或面积不匹配。');
mesh.lewisTrial=struct('pressurePa',rows.trialPressurePa,'heatFluxWm2',rows.trialHeatFluxWm2, ...
    'sourceFrameSHA256',receipt.sourceFrameSHA256,'physicalDtS',receipt.physicalDtS);
end
