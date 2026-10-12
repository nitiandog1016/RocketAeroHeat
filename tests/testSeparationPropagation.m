function testSeparationPropagation()
% 分叉准流线：层流分离后整条下游路径留空，旁支保持独立。
addpath(fileparts(fileparts(mfilename('fullpath'))));
path=struct;
path.predecessorFace=[0;1;2;2;3;4;0;7];
path.pathCost=[0;1;2;2;3;3;0;1];
path.connectedFace=true(8,1);
onset=false(8,1); onset(2)=true;
result=propagateSurfaceSeparation(path,onset);
assert(isequal(result.unresolvedFace, ...
    logical([0;1;1;1;1;1;0;0])));
assert(result.summary.onsetFaces==1);
assert(result.summary.downstreamFaces==4);
assert(contains(result.scope,'再附'));
bad=path; bad.pathCost(3)=1;
failed=false;
try, propagateSurfaceSeparation(bad,onset); catch, failed=true; end
assert(failed);
fprintf('Separation propagation passed: 1 onset, 4 downstream unresolved faces.\n');
end
