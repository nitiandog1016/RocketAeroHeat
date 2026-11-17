function testRegularizedEntropyV255(folder)
%TESTREGULARIZEDENTROPYV255 全角域、原面、半种子和标签控制；失败不丢完整输出。
source=jsondecode(fileread('outputs/v255_validation/prepare01/source.json'));
base=run('base',struct);fine=run('fine',struct('relativeTolerance',1e-10,'absoluteTolerance',1e-12));
seed=run('seed',struct('startAngleRad',5e-5));scaled=run('scaled',struct('metricScale',3));
keys={'heatFluxWm2','massThicknessM','edgeTemperatureK','edgeSpeedMps'};metrics=struct;
for pair={'fine','seed','scaled'}
    name=pair{1};other=eval(name);ids=base.targetIndices;target=other.targetIndices;difference=0;
    for key=keys
        x=base.(key{1});y=other.(key{1});relative=abs(x(ids)-y(target))./max(abs(y(target)),1e-30);
        difference=max(difference,max(relative));
    end
    metrics.(name)=difference;
end
write('comparison.json',metrics);
assert(metrics.fine<=.002&&metrics.seed<=.002,'RAH:RegularIntegral:Refinement','全角域容差/有限种子差超0.2%，不自动续试。');
assert(metrics.scaled<=1e-8,'RAH:RegularIntegral:Label','标签尺度改变物理结果。');
assert(max(abs(base.massCoordinate*3-scaled.massCoordinate)./max(abs(scaled.massCoordinate),1e-30))<=1e-8);
assert(max(abs(base.integralI*9-scaled.integralI)./max(abs(scaled.integralI),1e-30))<=1e-8);
for item={base,fine,seed,scaled}
    r=item{1};assert(r.maximumEnergyRelativeError<=1e-8&&r.maximumEntropyRelativeError<=1e-8);
    assert(abs(r.heatFluxWm2(2)/r.heatFluxWm2(1)-1)<.002);
end
% 模块二：新差量表达式与既有Cp绝对积分同源核对；小量按绝对舍入界。
T0=source.stagnationTemperatureK;drops=[0 1e-10 1e-8 1e-6 1 25 100 500 T0-60];c0=phaseAwareAirCaloric(T0);
[dh,ds]=phaseAwareAirCaloricDrop(T0,drops);ct=phaseAwareAirCaloric(T0-drops);
assert(max(abs(dh-(c0.h-ct.h)))<=32*eps*max(abs([c0.h ct.h])));
assert(max(abs(ds-(c0.s0-ct.s0)))<=32*eps*max(abs([c0.s0 ct.s0])));
rejections=0;bad=source;bad.sourceRole='accepted-edge';reject(@()solveRegularizedEntropyIntegral(bad));
bad=source;bad.targetAnglesRad(end)=pi/2;reject(@()solveRegularizedEntropyIntegral(bad));
bad=source;bad.targetAnglesRad(end)=1e-5;reject(@()solveRegularizedEntropyIntegral(bad));
bad=source;bad.totalEnthalpyJkg=bad.totalEnthalpyJkg+100;reject(@()solveRegularizedEntropyIntegral(bad));
report=struct('passed',true,'fullTargetSamples',numel(source.targetAnglesRad),'actualRingPositions',46, ...
    'fullAngleToleranceDifference',metrics.fine,'finiteSeedDifference',metrics.seed,'labelDifference',metrics.scaled, ...
    'nativeSteps',[base.acceptedSteps fine.acceptedSteps seed.acceptedSteps scaled.acceptedSteps], ...
    'maximumEnergyRelativeError',max([base.maximumEnergyRelativeError fine.maximumEnergyRelativeError]), ...
    'maximumEntropyRelativeError',max([base.maximumEntropyRelativeError fine.maximumEntropyRelativeError]), ...
    'maximumDenseODERelativeDefect',base.maximumDenseODERelativeDefect,'rejections',rejections, ...
    'physicalQualification',false,'profileReadyForCFD',false,'newQualifiedAreaM2',0,'MatlabVersion',version);
write('report.json',report);disp(report);
    function r=run(name,options)
        [r,solution]=solveRegularizedEntropyIntegral(source,options);
        save(fullfile(folder,[name '_continuous.mat']),'solution','-v7');write([name '.json'],r);
        assert(r.available,r.failure);
    end
    function write(name,obj)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');assert(fid>=0);clean=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(obj));
    end
    function reject(fn)
        refused=false;try,fn();catch err,refused=startsWith(err.identifier,'RAH:RegularIntegral:');end
        assert(refused);rejections=rejections+1;
    end
end
