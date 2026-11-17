function testCoupledEntropyV254(folder)
%TESTCOUPLEDENTROPYV254 实际46入口的联合积分/固定熵对照/预定细化/标签协变。
% 保存各完整分支或失败前缀，不重算v252三条实际剖面。
source=jsondecode(fileread('outputs/v254_validation/prepare01/source.json'));
base=solvePhaseAwareEntropyIntegral(source);write('base.json',base);assert(base.available,base.failure);
fine=solvePhaseAwareEntropyIntegral(source,struct('stepDeg',.25));write('fine.json',fine);assert(fine.available,fine.failure);
reference=solvePhaseAwareEntropyIntegral(source,struct('mode','normal-shock-reference'));write('reference.json',reference);assert(reference.available,reference.failure);
scaled=solvePhaseAwareEntropyIntegral(source,struct('metricScale',3));write('scaled.json',scaled);assert(scaled.available,scaled.failure);
ids=base.targetIndices;fineIDs=fine.targetIndices;assert(numel(ids)==47&&numel(fineIDs)==47);
keys={'heatFluxWm2','massThicknessM','edgeTemperatureK','edgeSpeedMps'};refinement=0;label=0;
for name=keys
    x=base.(name{1});y=fine.(name{1});z=scaled.(name{1});
    refinement=max(refinement,max(abs(x(ids)-y(fineIDs))./max(1e-30,abs(y(fineIDs)))));
    label=max(label,max(abs(x-z)./max(1e-30,abs(x))));
end
assert(refinement<=.02,'RAH:EntropyIntegral:Refinement','预定细化未过，不自动再加密。');
assert(label<=1e-8,'RAH:EntropyIntegral:Label','标签尺度改变物理输出。');
assert(max(abs(base.integralI*9-scaled.integralI)./max(1e-30,abs(scaled.integralI)))<=1e-8);
assert(max(abs(base.massCoordinate*3-scaled.massCoordinate)./max(1e-30,abs(scaled.massCoordinate)))<=1e-8);
rejections=0;
bad=source;bad.sourceRole='accepted-edge';reject(@()solvePhaseAwareEntropyIntegral(bad));
bad=source;bad.targetAnglesRad(end)=pi/2;reject(@()solvePhaseAwareEntropyIntegral(bad));
bad=source;bad.sourceMassCoordinate(1)=1;reject(@()solvePhaseAwareEntropyIntegral(bad));
reject(@()solvePhaseAwareEntropyIntegral(source,struct('wallTemperatureK',400)));
bad=source;bad.totalEnthalpyJkg=bad.totalEnthalpyJkg+100;reject(@()solvePhaseAwareEntropyIntegral(bad));
report=struct('passed',true,'actualRingPositions',46,'allRingIntegralResultsAvailable',true, ...
    'maximumStepRefinementRelativeDifference',refinement,'maximumLabelPhysicalRelativeDifference',label, ...
    'maximumImplicitResidual',max([base.maximumImplicitResidual fine.maximumImplicitResidual scaled.maximumImplicitResidual]), ...
    'jointLeadingGradient',base.jointLeadingVelocityGradient,'referenceLeadingGradient',base.referenceVelocityGradient, ...
    'oldMissingPositionsIntegralResultsAvailable',5,'actualProfilesFilled',0,'profileReadyForCFD',false, ...
    'rejections',rejections,'physicalQualification',false,'newQualifiedAreaM2',0,'MatlabVersion',version);
write('report.json',report);disp(report);
    function write(name,value)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');assert(fid>=0);clean=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value));
    end
    function reject(fn)
        passed=false;try,fn();catch err,passed=startsWith(err.identifier,'RAH:EntropyIntegral:');end
        assert(passed);rejections=rejections+1;
    end
end
