function [value,report]=rocketShellResidualDirectionalDerivative(Q,directionU,baseRate,cut,gas,boundary,scale,options)
%ROCKETSHELLRESIDUALDIRECTIONALDERIVATIVE 直接调用原完整三维NS残差形成矩阵自由方向导数。
% 模块一：方向是U=Q/scale；按无穷范数设置明确扰动，不以预条件通量替代原方程。
% 模块二：前向/中心差分及实际调用数单列；原状态/温域/HLLC/几何护栏保持。
% 模块三：差分导数是数值近似，不宣称完整解析Jacobian或物理时间推进。
if nargin<8,options=struct;end
defaults=struct('differenceStep',1e-6,'scheme','forward','saveEvidence',false,'spatialConservationTolerance',1e-9);
assert(isstruct(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'ShellJvp:Options','未知方向导数选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.differenceStep,{'numeric'},{'real','finite','positive','scalar'});
assert(ismember(options.scheme,{'forward','central'}),'ShellJvp:Options','差分形式必须明确。');
scale=reshape(scale,1,5);assert(isequal(size(directionU),size(Q))&&isequal(size(baseRate),size(Q))&&all(isfinite(directionU),'all'));
normDirection=max(abs(directionU),[],'all');assert(normDirection>0,'ShellJvp:ZeroDirection','零方向应由Krylov调用方直接处理。');
h=options.differenceStep/normDirection;plus=Q+h*directionU.*scale;
[ratePlus,infoPlus]=cutCellNavierStokesResidual(plus,cut,gas,boundary);calls=1;check(infoPlus);
if strcmp(options.scheme,'central')
    minus=Q-h*directionU.*scale;[rateMinus,infoMinus]=cutCellNavierStokesResidual(minus,cut,gas,boundary);calls=2;check(infoMinus);
    value=(ratePlus-rateMinus)./scale/(2*h);
else
    rateMinus=baseRate;value=(ratePlus-baseRate)./scale/h;
end
assert(all(isfinite(value),'all'),'ShellJvp:Nonfinite','原残差方向导数非有限。');
report=struct('scheme',options.scheme,'differenceStepInScaledState',options.differenceStep,'actualScalarStep',h, ...
    'directionInfinityNorm',normDirection,'actualOriginalResidualEvaluations',calls, ...
    'fullThreeDimensionalResidualCalled',true,'physicalTimeAdvanced',0,'originalEquationReplaced',false);
if options.saveEvidence,report.directionU=directionU;report.ratePlus=ratePlus;report.rateMinusOrBase=rateMinus;end
    function check(info)
        f=cut.openFaces;w=cut.walls;mask=f.cellPlus==0;
        denominator=max(1,sum(abs(info.faceFlux(mask,:)).*f.area(mask),1)+sum(abs(info.wallFlux).*w.area,1));
        assert(max(abs(info.balanceResidual)./denominator)<=options.spatialConservationTolerance, ...
            'ShellJvp:Conservation','差分实际空间算子未通过原五量账。');
    end
end
