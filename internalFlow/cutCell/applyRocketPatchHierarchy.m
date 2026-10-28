function [value,report]=applyRocketPatchHierarchy(hierarchy,rhs)
%APPLYROCKETPATCHHIERARCHY 一次固定线性V循环：局部五量逆，再递归校正跨组剩余。
% 模块一：层数和作用次数固定，不按右端/剩余早停或裁掉物理分量。
% 模块二：每层重新使用该层完整C；最粗解仍只是预条件的一部分。
assert(strcmp(hierarchy.kind,'fixed-all-five-field-unsmoothed-strong-connection-V-cycle')&& ...
    isequal(size(rhs),[hierarchy.cellCount 5])&&isreal(rhs)&&all(isfinite(rhs),'all'),'PatchHierarchy:Input','固定层次或原五量右端非法。');
nl=numel(hierarchy.levels);actions=zeros(nl,1);localMaps=actions;terminalMaps=actions;value=cycle(1,rhs);
report=struct('CApplicationsByLevel',actions,'localMapApplicationsByLevel',localMaps,'terminalApplicationsByLevel',terminalMaps, ...
    'fixedLinearMap',true,'originalNSAcceptanceClaimed',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function x=cycle(level,b)
        entry=hierarchy.levels{level};n=entry.packed.cellCount;
        if entry.terminal
            x=reshape(entry.inverse*b(:),n,5);terminalMaps(level)=terminalMaps(level)+1;
        else
            x=solveRocketCoupledPatches(entry.factor,b);localMaps(level)=localMaps(level)+1;
            remaining=b-applyRocketPackedGraph(entry.packed,x);actions(level)=actions(level)+1;
            coarse=restrictRocketPatchResidual(entry.factor.partition,remaining);correction=cycle(level+1,coarse);
            x=x+prolongRocketPatchCorrection(entry.factor.partition,correction);
        end
        assert(all(isfinite(x),'all'),'PatchHierarchy:Nonfinite','层次作用非有限，拒绝使用。');
    end
end
