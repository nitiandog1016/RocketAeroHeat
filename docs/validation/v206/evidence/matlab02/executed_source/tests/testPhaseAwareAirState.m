function report=testPhaseAwareAirState(outputName)
%TESTPHASEAWAREAIRSTATE 跨语言数据、独立积分及导数核对新工程热接口。
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1, outputName='matlab02'; end
fixture=jsondecode(fileread(fullfile(root,'outputs','v206_validation','python_tests01','matlab_fixture.json')));
p=fixture.pressurePa(:); T=fixture.temperatureK(:); actual=phaseAwareAirState('PT',p,T);
keys=fieldnames(fixture.states); maximum=0;
for k=1:numel(keys)
    key=keys{k}; reference=fixture.states.(key)(:); observed=actual.(key)(:);
    maximum=max(maximum,max(abs(observed-reference)./max(1,abs(reference))));
end
assert(maximum<2e-9,'RAH:PhaseTest:State','MATLAB/Python状态不一致。');
fromH=phaseAwareAirState('Ph',p,actual.enthalpyJkg); fromE=phaseAwareAirState('Pe',p,actual.internalEnergyJkg);
inverseError=max([abs(fromH.temperatureK-T);abs(fromE.temperatureK-T)]);
assert(inverseError<2e-8);
row=phaseAwareAirState('PT',p.',T.'); assert(isequal(size(row.temperatureK),size(T.')));
assert(max(abs(row.enthalpyJkg(:)-actual.enthalpyJkg))<1e-7);
inputs=fixture.referenceInputs;
ref=phaseAwareAirReferenceState(inputs(:,1),inputs(:,2),inputs(:,3),inputs(:,4),inputs(:,5));
keys=fieldnames(fixture.referenceOutputs); refMaximum=0;
for k=1:numel(keys)
    key=keys{k};
    % Python附带平均自由程；MATLAB参考接口只返回所需输运，不伪造重复字段。
    if ~isfield(ref,key), assert(strcmp(key,'meanFreePathViscosityDefinitionM')); continue; end
    expected=fixture.referenceOutputs.(key)(:); observed=ref.(key)(:);
    refMaximum=max(refMaximum,max(abs(observed-expected)./max(1,abs(expected))));
end
assert(refMaximum<2e-9);
step=.01;
plus=phaseAwareAirReferenceState(inputs(:,1),inputs(:,2),inputs(:,3),inputs(:,4)+step,inputs(:,5));
minus=phaseAwareAirReferenceState(inputs(:,1),inputs(:,2),inputs(:,3),inputs(:,4)-step,inputs(:,5));
derivative=max(abs((plus.referenceTemperatureK-minus.referenceTemperatureK)/(2*step)-ref.referenceTemperatureDerivativeByWallTemperature));
assert(derivative<1e-7);
assert(ref.thermalDrivingEnthalpyJkg(end)<0 && isempty(ref.wallHeatFluxWm2) && ~ref.engineeringHeatReady);
rejects={@()phaseAwareAirState('PT',1,59.9),@()phaseAwareAirState('PT',1,2000.1), ...
    @()phaseAwareAirState('PT',10000,60),@()phaseAwareAirState('PT',0,300), ...
    @()phaseAwareAirState('PT',1,NaN),@()phaseAwareAirState('PT',[1;2],[300 400]), ...
    @()phaseAwareAirState('Ph',1,-1e10),@()phaseAwareAirState('Pe',1,1e10), ...
    @()phaseAwareAirState('Ph',100000,actual.enthalpyJkg(1)), ...
    @()phaseAwareAirReferenceState(1000,300,5000,300,.9), ...
    @()phaseAwareAirReferenceState(1000,300,1000,300,1.1), ...
    @()phaseAwareAirReferenceState([1;2],[300 400],1000,300,.9)};
for k=1:numel(rejects)
    failed=false; try, rejects{k}(); catch, failed=true; end
    assert(failed,'RAH:PhaseTest:Reject','非法输入未拒绝。');
end
report=struct('passed',true,'stateMaxScaledDifference',maximum,'referenceMaxScaledDifference',refMaximum, ...
    'inverseMaxErrorK',inverseError,'wallDerivativeMaxAbs',derivative,'rejections',numel(rejects), ...
    'dataSHA256',actual.sourceSHA256,'physicalQualification',false,'newQualifiedHeatAreaM2',0);
out=fullfile(root,'outputs','v206_validation',outputName);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8'); closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report)); disp(report);
end
