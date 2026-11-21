function mesh=attachThermistorThreeGrids(mesh)
% 模块一：从两次已核验的冷态解取回同一试件的原面，并绑定本次第三档结果。
c=mesh.hybridCondition;w=mesh.hybridWall;s=mesh.thermistorLead;
if isfield(s,'referenceNative'),s=s.referenceNative;end
v=mesh.V;f=mesh.F;n=size(f,1);
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
heat=nan(n,3);heat(:,3)=w.heatFluxWm2;
powerMode=isfield(c,'netHeatingPowerW');
info=struct('mode','givenWall','boundaryValue',305.15,'boundaryUnit','K','available',true);
if powerMode
    info.mode='givenPower';info.boundaryValue=c.netHeatingPowerW;info.boundaryUnit='W';
end
regionNames={s.regions.name};
for level=1:3
    if level==3
        src=s;wallK=305.15;
        if powerMode,wallK=mesh.hybridSolvedWallTemperatureK(1);end
    else
        [archive,~,~]=thermistorLeadGridChoice(level-1);
        [folder,guard]=openThermistorArchive(archive);
        src=jsondecode(fileread(fullfile(folder,'summary.json')));
        body=jsondecode(fileread(fullfile(folder,'body.json')));
        spec=jsondecode(fileread(fullfile(folder,'case.json')));
        path=fullfile(folder,'wall.csv');r=readtable(path);
        assert(strcmp(src.bodySHA256,s.bodySHA256) && strcmp(src.caseSHA256,s.caseSHA256) && ...
            isequal(body.V,v) && isequal(body.F,f) && isequal(r.faceID,w.faceID) && isequal(r.regionID,body.regionID) && ...
            all(abs(r.areaM2-area)<1e-18) && strcmp(cm011FileSHA256(path),src.wallSHA256) && ...
            abs(area'*r.heatFluxWm2+src.totalOutwardW)<1e-12 && isequal({src.regions.name},regionNames));
        clear guard;
        % 模块二：各网用本域热导反算同功率壁温，检查输入是否在该参照场范围内。
        wallK=NaN;
        if ~powerMode || c.netHeatingPowerW<=src.totalOutwardW
            [r,~,~,wallK]=scaleThermistorWallField(r,c,spec.airK,spec.wallK,src.totalOutwardW);
            heat(:,level)=r.heatFluxWm2;
        else
            info.available=false;
        end
    end
    power=-area'*heat(:,level);
    if powerMode && isfinite(power),assert(abs(power-c.netHeatingPowerW)<=1e-11*c.netHeatingPowerW);end
    info.grids(level)=struct('level',level,'meshSHA256',src.meshSHA256,'frameSHA256',src.sourceFrameSHA256, ...
        'reactionSHA256',src.sourceReactionSHA256,'wallSHA256',src.wallSHA256,'conductanceWK',src.totalOutwardW/12, ...
        'wallK',wallK,'totalOutwardW',power);
end
% 模块三：按暴露原面的名称和面积累计七区放热；正热流入壁，零功率保留为物理零。
for j=1:numel(regionNames)
    use=body.regionID==j;
    info.regions(j)=struct('name',regionNames{j},'faces',nnz(use),'areaM2',sum(area(use)),'outwardW',-area(use)'*heat(use,:));
end
assert(strcmp(info.grids(3).frameSHA256,mesh.evidenceIdentity.sourceFrame.sha256));
mesh.evidenceIdentity.sourceFrame.midbodyGrids=info.grids;
mesh.thermistorThreeGrids=struct('summary',info,'heat',heat);
end
