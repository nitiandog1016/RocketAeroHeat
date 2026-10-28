function report=verifyRocketShellRecomputed(newRoot,frozenRoot)
%VERIFYROCKETSHELLRECOMPUTED 独立工程中重算后逐字段核对真实三维域、全算子和RK2。
% 模块一：完整数值对象严格isequaln，不用摘要通过替代全状态身份。
% 模块二：运行耗时/PID不要求相同，三组原面结果JSON须逐字节一致。
files={'domain02/shell_domain_evidence.mat','initial01/initial_residual_evidence.mat', ...
    'thermal02/wall_thermal_evidence.mat','advance01/advance_evidence.mat'};
variables={{'domain'},{'Q','rate','info','mapped','gas','boundary'}, ...
    {'rateNew','infoNew','infoPerturbed','perturbation','heatGeometryMatrix','boundary','gas','audit'}, ...
    {'result','mapped','gas','boundary'}};
checked=0;
for j=1:numel(files)
    a=load(fullfile(newRoot,files{j}));b=load(fullfile(frozenRoot,files{j}));
    for name=variables{j}
        key=name{1};assert(isequaln(a.(key),b.(key)),'RocketShellReplay:Field','重算完整字段不同：%s / %s',files{j},key);checked=checked+1;
    end
    clear a b;
end
for name={'legacy-initial','dirichlet-initial','dirichlet-startup'}
    file=fullfile('cases',name{1},'results.json');
    assert(strcmp(fileread(fullfile(newRoot,'workbench02',file)),fileread(fullfile(frozenRoot,'workbench02',file))), ...
        'RocketShellReplay:View','实际原面结果JSON未逐字节复现。');
end
a=jsondecode(fileread(fullfile(newRoot,'workbench02','shell-study-data.json')));
b=jsondecode(fileread(fullfile(frozenRoot,'workbench02','shell-study-data.json')));
for name={'kind','cases','history','sourceFaces','cellCount','layerCount','wallNormalDistanceRangeM','gas','commonRanges', ...
          'geometryURL','domainURL','readOnlyFrozenStudy','productionQualified','newProductionSolvedAreaM2','browserVisualOrWebGLRuntimeVerified'}
    key=name{1};assert(isequaln(a.(key),b.(key)),'RocketShellReplay:Study','研究数值/身份未复现：%s',key);checked=checked+1;
end
report=struct('status','passed-isolated-recomputed-full-domain-spatial-thermal-RK2-and-three-JSON-identity', ...
    'completeNumericalFieldsChecked',checked,'resultJSONExactlyReproduced',3, ...
    'fullIndependentNavierStokesSolver',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
fid=fopen(fullfile(newRoot,'recomputed_identity.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('%s\n',jsonencode(report));
end
