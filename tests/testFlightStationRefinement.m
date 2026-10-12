function testFlightStationRefinement()
% 时间站线性插值：保留输入站、限制实际步长并核对每段物理输入。
addpath(fileparts(fileparts(mfilename('fullpath'))));
flight=struct('time',[0 2 5], ...
    'altitude',[30000 20000 40000], ...
    'Mach',[5 4 6], ...
    'alpha',[0 10 -5], ...
    'beta',[0 0 5]);
[dense,info]=refineFlightStations(flight,1);
assert(isequal(dense.time,0:5));
assert(isequal(info.inputStationIndex,[1 3 6]));
assert(isequal(info.segmentSubsteps,[2 3]));
assert(info.inputStationCount==3 && info.evaluationStationCount==6);
assert(max(diff(dense.time))<=1);
assert(isequal(dense.altitude(info.inputStationIndex),flight.altitude));
assert(isequal(dense.Mach(info.inputStationIndex),flight.Mach));
assert(isequal(dense.alpha(info.inputStationIndex),flight.alpha));
assert(dense.altitude(2)==25000 && dense.Mach(2)==4.5);
assert(abs(dense.alpha(4)-5)<1e-12);
assert(abs(dense.beta(4)-5/3)<1e-12);
assert(contains(info.method,'线性插值'));
assert(isequal(flight.time,[0 2 5]));
failed=false;
try, refineFlightStations(flight,0); catch, failed=true; end
assert(failed);
bad=flight; bad.time=[0 2 2];
failed=false;
try, refineFlightStations(bad,1); catch, failed=true; end
assert(failed);
fprintf('Flight station refinement passed: %d input -> %d evaluation stations.\n', ...
    info.inputStationCount,info.evaluationStationCount);
end
