function report=writeSU2HybridInitialState(volume,folder,thickness)
%WRITESU2HYBRIDINITIALSTATE 在原混合体网格上生成正状态、满足无滑移/等温壁的守恒初猜。
% 模块一：仅接受完整法向棱柱区和给定等温壁，外层必须超过显式初始化厚度。
% 模块二：速度及总焓用C1函数连接壁面与自由流；静温由总焓扣除动能求得。
% 模块三：记录初猜身份，不当成边界层解、激波解或正式温度/热流结果。
validateattributes(thickness,{'numeric'},{'scalar','finite','real','positive'});
assert(strcmp(volume.meshConstruction.method,'original-triangle-normal-prisms-and-gmsh-exterior-tetrahedra'));
file=fullfile(folder,'initial_flow.csv');assert(~isfile(file) && ~isfile(fullfile(folder,'engine_process.json')));
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));f=manifest.flight;
assert(strcmp(manifest.thermalBoundary.mode,'isothermal'),'当前初猜仅支持显式等温壁。');
Tw=manifest.thermalBoundary.specification.temperature;U=manifest.expectedVelocity(:)';
nv=size(volume.sourceVertices,1);nf=size(volume.sourceFaces,1);
layers=volume.meshConstruction.prisms/nf;assert(layers==fix(layers));
np=size(volume.vertices,1);assert((layers+1)*nv<=np);
assert(isequal(volume.vertices(1:nv,:),volume.sourceVertices));
b=ones(np,1);distance=zeros((layers+1)*nv,1);
for k=0:layers
    ids=k*nv+(1:nv);distance(ids)=vecnorm(volume.vertices(ids,:)-volume.sourceVertices,2,2);
end
assert(all(distance(end-nv+1:end)>=thickness),'棱柱区不足以完整容纳初猜厚度，禁止在四面体中突跳。');
eta=min(distance/thickness,1);b(1:numel(distance))=eta.^2.*(3-2*eta);
assert(all(b(1:nv)==0) && all(b(layers*nv+1:end)==1));
far=unique([volume.farfield{:}]);assert(all(b(far)==1),'远场不满足精确自由流。');
R=287.05287;gamma=1.4;cp=gamma*R/(gamma-1);K=sum(U.^2)/2;
velocity=b.*U;
% 此恒等式明确保证T>=min(Tw,Tinf)，无温压裁剪；总焓沿b线性变化。
T=Tw+(f.temperature-Tw)*b+(K/cp)*b.*(1-b);
rho=f.pressure./(R*T);momentum=rho.*velocity;
energy=rho.*(R*T/(gamma-1)+sum(velocity.^2,2)/2);
assert(all(isfinite([rho;momentum(:);energy])) && all(rho>0) && all(T>0));
values=[(0:np-1)',volume.vertices,rho,momentum,energy];
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));
fprintf(fid,'"PointID","x","y","z","Density","Momentum_x","Momentum_y","Momentum_z","Energy"\n');
fprintf(fid,'%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n',values');clear g;
cfg=fileread(fullfile(folder,'case.cfg'));assert(count(cfg,'RESTART_SOL= NO')==1);
cfg=strrep(cfg,'RESTART_SOL= NO','RESTART_SOL= YES');
cfg=[cfg newline 'SOLUTION_FILENAME= initial_flow' newline 'READ_BINARY_RESTART= NO' newline];
writeText(fullfile(folder,'case.cfg'),cfg);
report=struct('method','C1 velocity/total-enthalpy wall-to-freestream initial guess', ...
    'thicknessM',thickness,'normalLayers',layers,'originalWallNodes',nv,'nodes',np, ...
    'temperatureRangeK',[min(T),max(T)],'minimumDensity',min(rho),'wallNoSlip',true, ...
    'exactFarfield',true,'initialPressurePa',f.pressure,'givenWallTemperatureK',Tw, ...
    'constantCp',cp,'initialGuessOnly',true,'boundaryLayerSolved',false,'physicsQualified',false);
manifest.initialization=report;writeText(fullfile(folder,'initialization_manifest.json'),jsonencode(report));
writeText(fullfile(folder,'case_manifest.json'),jsonencode(manifest));
end
function writeText(file,text)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',text);
end
