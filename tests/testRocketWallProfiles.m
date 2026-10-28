function report=testRocketWallProfiles(folder,sourceFolder)
%TESTROCKETWALLPROFILES 从实际四档Q读取剖面，与独立Python逐格/逐内面审计对照。
% 温度/速度/距离为实际控制体代理；给定墙值及未知Te/Tr/求解Tw明确分离。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v77_validation','profile_test');end
if nargin<2,sourceFolder=fullfile(root,'outputs','v77_validation');end
assert(~exist(folder,'dir'),'RocketWallProfileTest:Overwrite','不能覆盖剖面核对。');
a=jsondecode(fileread(fullfile(sourceFolder,'wall_resolution_independent.json')));rows=struct([]);
versions=[75 76 77 77];names={'spatial_root_double','radial48','collar72','nested144'};
for k=1:4
    p=fullfile(root,'outputs',sprintf('v%d_validation',versions(k)),names{k},'spatial_root_evidence.mat');
    if ~exist(p,'file'),p=fullfile(root,'validation',sprintf('v%d',versions(k)),names{k},'spatial_root_evidence.mat');end
    d=load(p,'c','solution');profiles=buildRocketWallNormalProfiles(d.c,d.solution);oracle=a.cases(k).wallNormalProfiles;
    assert(numel(profiles)==115&&numel(oracle)==115);maxError=0;
    numeric={'cellIDs','cellCentroidXR','distanceM','temperatureK','densityKgM3','pressurePa','speedMps', ...
        'tangentialVelocityMps','muTOverMu','radialFaceID','radialFaceNormalDistanceM', ...
        'EulerEnergyFluxMinusArithmeticPhysicalWm2','viscousTotalEnergyFluxWm2'};
    for i=1:115
        b=profiles(i);o=oracle(i);assert(strcmp(b.region,o.region)&&strcmp(b.wallBoundaryValueKind,'prescribed')&& ...
            ~b.edgeTemperatureSolved&&~b.recoveryTemperatureSolved&&~b.wallTemperatureSolved&&b.wallTemperatureInputK==300);
        for fields=numeric
            key=fields{1};x=b.(key);y=o.(key);assert(isequal(size(x),size(y)));
            scale=max(1,max([abs(x(:));abs(y(:))]));
            if contains(key,'Flux'),scale=max(scale,d.c.residualScale(4)*d.c.freestream.speed);end
            error=max(abs(x-y),[],'all')/scale;maxError=max(maxError,error);
            assert(error<=512*eps,'RocketWallProfileTest:Independent','剖面真实字段%s与独立审计不符。',key);
        end
    end
    row=struct('kind',a.cases(k).kind,'wallProfiles',115,'radialLayers',d.c.mesh.shape(2), ...
        'maximumScaledDifference',maxError,'independentEOSCellAndFaceComparisonPassed',true);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
end
report=struct('status','passed-four-real-states-all-wall-normal-cell-and-energy-flux-profiles', ...
    'cases',rows,'engineProcessID',feature('getpid'),'physicalSolution',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'profile_test_evidence.mat'),'report');
fid=fopen(fullfile(folder,'profile_test_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
