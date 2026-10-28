function testNSOpenRegionAdvance(folder)
%TESTNSOPENREGIONADVANCE 开放特征边界与真实RK2推进的联合数值校核。
% 模块用途：检查声学/熵扰动离域、回流储库身份和五量边界积分。
% 这是小网格平直开放域，不是箭体尾迹实验或真实分离/底热验证。
addpath(fileparts(mfilename('fullpath')),'-begin');
gas=struct('gamma',1.4,'R',287,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]);
shape=[10 2 2];grid=struct('spacing',[.01 .01 .01]);solid=false(shape);
periodic=struct('type','periodic');cv=gas.R/(gas.gamma-1);
record=struct('experimentalValidation',false,'actualRocketSolved',false,'newSolvedSurfaceAreaM2',0);
for direction=[1 -1]
    given=[1.2,100*direction,0,0,300];
    far=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',given);
    outlet=far;outlet.mode='pressureOutlet';outlet.pressure=given(1)*gas.R*given(5);
    if direction<0,outlet.pressure=1e8;end % 回流时该压力必须明确不被使用。
    boundary=struct('xmin',far,'xmax',outlet,'ymin',periodic,'ymax',periodic,'zmin',periodic,'zmax',periodic);
    primitive=repmat(reshape(given,1,1,1,5),[shape,1]);
    Q=conserved(primitive);uniform=advanceCompressibleNavierStokes3D(Q,grid,solid,gas,boundary, ...
        struct('targetTime',5e-5,'maxSteps',200,'conservationTolerance',1e-11));
    assert(uniform.completedRequestedTime && ~uniform.converged && ~uniform.physicalSolution);
    assert(max(abs(uniform.Q(:)-Q(:)))/max(abs(Q(:)))<1e-12);
    assert(uniform.globalBalanceRelative<1e-11);
    % 从偏离均匀流的温度/压力扰动开始；必须经过真实时间步和开放通量输运。
    x=((1:shape(1))-.5)*grid.spacing(1);perturb=sin(2*pi*x/(shape(1)*grid.spacing(1)));
    primitive(:,:,:,5)=primitive(:,:,:,5)+reshape(3*perturb,[],1,1);
    Q=conserved(primitive);result=advanceCompressibleNavierStokes3D(Q,grid,solid,gas,boundary, ...
        struct('targetTime',.002,'maxSteps',3000,'conservationTolerance',1e-10));
    assert(result.completedRequestedTime && result.acceptedSteps>100 && ~result.converged);
    assert(result.globalBalanceRelative<1e-10 && result.positivityClippingCount==0);
    initialT=primitive(:,:,:,5);final=result.finalSpatialInfo.primitive;
    finalT=final(:,:,:,5);
    initialError=mean((initialT(:)-300).^2);finalError=mean((finalT(:)-300).^2);
    assert(finalError<.25*initialError,'扰动未按明确开放边界离域，不能称正确开放输运。');
    faces=result.finalSpatialInfo.openFaces;right=faces([faces.axis]==1);
    right=right(arrayfun(@(f)f.normal(1)>0,right));
    assert(numel(right)==4,'Expected four right boundary faces; got %d of %d, axes=%s, normals=%s', ...
        numel(right),numel(faces),mat2str([faces.axis]),mat2str(reshape([faces.normal],3,[])'));
    if direction<0
        assert(all(arrayfun(@(f)f.boundaryState.isInflow,right)));
        assert(~any(arrayfun(@(f)f.boundaryState.pressureInputApplied,right)));
    else
        assert(all(arrayfun(@(f)f.boundaryState.pressureInputApplied,right)));
    end
    key='outflow';if direction<0,key='backflow';end
    record.(key)=struct('steps',result.acceptedSteps,'time',result.time, ...
        'initialTemperatureMeanSquareK2',initialError,'finalTemperatureMeanSquareK2',finalError, ...
        'globalBalanceRelative',result.globalBalanceRelative,'rejectedAttempts',result.rejectedAttempts);
    fprintf('Open NS %s passed: %d real RK2 steps, T disturbance %.6g -> %.6g K^2; balance %.3g.\n', ...
        key,result.acceptedSteps,initialError,finalError,result.globalBalanceRelative);
end
% 算子程序错误不属于正性失败；必须保留原异常并停止，不能反复减dt。
badCalls=0;lastBoundary=boundary;
failed=false;
try
    advanceCompressibleNavierStokes3D(Q,grid,solid,gas,lastBoundary, ...
        struct('targetTime',.001,'residualFunction',@brokenOperator));
catch exception
    failed=strcmp(exception.identifier,'RocketAeroHeat:Test:OperatorContract');
end
assert(failed && badCalls==2,'算子合同错误被错误重试或丢失原始identifier。');
record.unrecoverableOperatorCalls=badCalls;
if nargin>0
    save(fullfile(folder,'open_advance_evidence.mat'),'record','result','uniform');
    fid=fopen(fullfile(folder,'open_advance_summary.json'),'w','n','UTF-8');assert(fid>=0);
    close=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(record),'char'); %#ok<NASGU>
end
    function q=conserved(p)
        q=p;rho=p(:,:,:,1);velocity=p(:,:,:,2:4);
        q(:,:,:,2:4)=rho.*velocity;
        q(:,:,:,5)=rho.*(cv*p(:,:,:,5)+.5*sum(velocity.^2,4));
    end
    function [rate,meta]=brokenOperator(q,g,s,air,bc)
        badCalls=badCalls+1;
        if badCalls>1,error('RocketAeroHeat:Test:OperatorContract','明确的非物理算子程序错误校核。');end
        [rate,meta]=compressibleNavierStokesResidual3D(q,g,s,air,bc);
    end
end
