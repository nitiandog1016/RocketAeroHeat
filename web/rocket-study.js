/* 模块用途：只切换冻结研究快照的五个实际查看入口，始终复用同一iframe。 */
(function(){'use strict';
  const d=window.ROCKET_STUDY,$=id=>document.getElementById(id),buttons=[...document.querySelectorAll('[data-step]')];
  if(!d||d.kind!=='complete-rocket-frozen-research-workbench-v1')throw Error('完整箭体研究快照无效');
  const num=(x,n=6)=>Number(x).toPrecision(n),f=d.flight,b=d.boundary;
  /* 模块用途：按真实工况的两个独立策略说明来源；旧v68快照缺字段时保持原方法身份。 */
  const strategies=d.numericalStrategies||{},wallPolicy=strategies.wallVelocityGradientPolicy||'legacy',workPolicy=strategies.viscousWorkVelocityPolicy||'legacy-wls';
  if(!['legacy','dirichlet-tangent'].includes(wallPolicy)||!['legacy-wls','arithmetic'].includes(workPolicy))throw Error('未知离散策略，不能猜测算法来源');
  const wallLabel=wallPolicy==='dirichlet-tangent'?'静止壁Dirichlet切向约束：两速度分量沿壁导数为0':'原WLS与位移修正墙梯度';
  const workLabel=workPolicy==='arithmetic'?'内部面两侧单元速度代理的算术平均；边界速度迹保持':'原内部面WLS速度迹；边界速度迹保持';
  const pressurePolicy=strategies.pressureSourceQuadraturePolicy||'linear-planar';
  if(!['linear-planar','cell-constant'].includes(pressurePolicy))throw Error('未知压力源积分策略');
  const pressureLabel=pressurePolicy==='cell-constant'?'单元压力分片常数 × 平面面积':'WLS压力在平面质心积分';
  const pressureNote='该选项仅作用于压力几何源；通量由独立算法决定；不是裁零或给定边界';
  const interiorPolicy=strategies.interiorConvectiveFluxPolicy||'hllc';
  if(!['hllc','hlle-davis','hlle-einfeldt'].includes(interiorPolicy))throw Error('未知内部Euler通量策略');
  const interiorLabel=interiorPolicy==='hllc'?'HLLC':interiorPolicy==='hlle-davis'?'HLLE（两侧声学波速）':'HLLE（Einfeldt / Roe波速）';
  const interiorNote='仅内部Euler通量；开放面仍用原HLLE；SA采用实际质量通量。接触/剪切耗散需单独验证';
  /* 模块用途：数值Jacobian与内部初猜独立呈现，不能把暖初猜标成均匀新链或给定入口。 */
  const jacobianPolicy=strategies.jacobianDifferencePolicy||'legacy-normalized-unit-floor';
  if(!['legacy-normalized-unit-floor','local-conservative-acoustic-SA-roundoff-floor-v1'].includes(jacobianPolicy))throw Error('未知Jacobian扰动策略');
  const localJacobian=jacobianPolicy==='local-conservative-acoustic-SA-roundoff-floor-v1';
  const jacobianLabel=localJacobian?'局部密度/声学动量/能量/SA尺度，SA舍入底限；相对扰动 '+num(strategies.jacobianRelativeStep):'原归一化量与单位底限的差分扰动';
  const seed=d.initialization?.seed;
  if(localJacobian&&(!seed||!['uniform','v73-uniform23','v75-spatial-root','v76-sensitivity-spatial-root','v77-wall-resolution-root','v78-HLLC-spatial-root','v78-HLLC-refined-spatial-root'].includes(seed.kind)||seed.initialGuessIsBoundaryCondition!==false))throw Error('缺少真实内部初猜身份');
  const seedLabel=seed?.kind==='v78-HLLC-spatial-root'||seed?.kind==='v78-HLLC-refined-spatial-root'?'原HLLC内部定常根，'+seed.priorNewtonIterations+' 次Newton准备；父根初猜历史单列；本链从0计BE步；源SHA256 '+seed.sourceEvidenceSHA256:seed?.kind==='v77-wall-resolution-root'?'近壁分辨率网格内部定常根，'+seed.priorNewtonIterations+' 次Newton；本链从0计BE步；源SHA256 '+seed.sourceEvidenceSHA256:seed?.kind==='v76-sensitivity-spatial-root'?'同一箭体新网格的内部定常根，'+seed.priorNewtonIterations+' 次Newton；父v75根仅用于数值初猜；新根SHA256 '+seed.sourceEvidenceSHA256+'；本链从0计BE步':seed?.kind==='v75-spatial-root'?'内部定常Newton解，'+seed.priorNewtonIterations+' 次Newton；原空间残差 '+num(seed.spatialRootResidualNorm)+' s⁻¹；本链从0计BE步，不是给定入口。来源SHA256 '+seed.sourceEvidenceSHA256:seed?.kind==='v73-uniform23'?'内部冻结Q23，来源SHA256 '+seed.sourceEvidenceSHA256+'；此前23步单列，本链从0计步；不是给定入口':seed?.kind==='uniform'?'完整外域均匀来流初猜，本链从0计步':'历史快照未单列初猜元数据';
  const arithmetic=strategies.interiorFluxArithmetic||'double';
  if(!['double','two-double-euler'].includes(arithmetic))throw Error('未知内部通量算术');
  const arithmeticLabel=arithmetic==='double'?'普通double；黏性、边界与源同为double':'内部Euler形成使用双double余项，最终通量为double；黏性、边界与源仍为double';
  const wallScope=wallPolicy==='dirichlet-tangent'?'仅改选中墙面的Newton应力梯度；内场WLS、SA与Fourier保持':'沿用原墙面Newton应力梯度';
  const workScope=workPolicy==='arithmetic'?'仅改内部面黏性功；牵引与动量通量保持。一阶算术平均，未证明一般偏斜网格正性或二阶精度':'沿用原黏性功速度；不启用内部面算术平均';
  const gridPolicy=strategies.gridPolicy||'star-ray-reference',collar=strategies.effectiveCollarOptions||{};
  if(!['normal-collar','star-ray-reference'].includes(gridPolicy))throw Error('未知网格构造策略');
  if(gridPolicy==='normal-collar'&&!(Number.isInteger(collar.normalLayers)&&collar.normalLayers>0&&Number.isFinite(collar.firstNormalSpacing)&&collar.firstNormalSpacing>0&&Number.isFinite(collar.normalGrowth)&&collar.normalGrowth>=1))throw Error('缺少真实近壁法向层配置');
  const gridLabel=gridPolicy==='normal-collar'?'分段等法向偏移近壁层与外域连接层':'完整星形射线网格';
  const gridNote=collar.radialRefinementFactor===2?collar.normalLayers+' 个近壁层，由父网格逐层二分；首节点偏移 '+num(collar.firstNormalSpacing*1e6)+' μm；相邻层厚比交替为1和父增长率 '+num(collar.normalGrowth)+'，不是单一几何级数；节点偏移不等于体质心墙距。':gridPolicy==='normal-collar'?collar.normalLayers+' 个近壁层；首节点法向偏移 '+num(collar.firstNormalSpacing*1e6)+' μm，增长率 '+num(collar.normalGrowth)+'。节点偏移不等于体质心墙距':'未启用法向偏移层；射线长度不等于墙法向距离';
  const temperatureNote=(d.wallAdjacentGasProxyIncluded?'邻壁气体温度另列为单元代理；':'')+'Tw为给定；外缘温度、恢复温度、求解Tw仍未知';
  const panels=[
    {title:'输入模型 · 已冻结',description:'原生闭合零攻角箭体，含球鼻、锥段、柱段、收缩、后体与真实底面。',rows:[['总长',num(d.totalLengthM)+' m'],['原模型',d.sourceVertices+' 顶点 / '+d.sourceFaces+' 三角面'],['支持工况','轴对称、无周向突起；当前不代表20°横流解']]},
    {title:'生成网格 · 已构建',description:'完整外域由程序内部建立；几何页可切换原模型和细分墙网格。',rows:[['流场控制体',d.cellCount+' 格 / '+d.radialLayers+' 层'],['网格构造',gridLabel],['近壁节点设置',gridNote],['真实墙段',d.wallSegments+' 段，包含真实底面'],['墙首层法向距离',d.nearWallNormalDistanceM.map(x=>num(x)).join(' ～ ')+' m'],['曲面 / 原弦面面积',num(d.curvedAreaM2)+' / '+num(d.nativeFacetAreaM2)+' m²'],['质量状态','近壁与外域充分性尚未验收']]},
    {title:'添加边界条件 · 已给定',description:'只给定远场、静止等温壁和SA远场量；前体与后体不导入外部CFD。',rows:[['飞行条件','M '+f.Mach+' / '+f.altitude+' m / α '+f.alpha+'°'],['壁温输入',b.wallTemperatureInputK+' K，非求解壁温'],['远场原始量',b.freestreamPrimitiveRhoUxUrT.map(x=>num(x)).join(' / ')+' （ρ, ux, ur, T）'],['湍流输入','ν̃/ν∞ = '+b.nuTildeRatioInput],['边界方法','完整远场Riemann储库；无滑移等温墙；镜像奇偶轴对称'],['壁速度梯度',wallLabel],['内面黏性功速度',workLabel]]},
    {title:'计算 · '+d.stage,description:'显示每个真实接受端点；完整五量BE方程和守恒账已经逐步重放。',rows:[['接受步数',String(d.acceptedSteps)],['初猜来源',seedLabel],['Jacobian扰动',jacobianLabel],['内部通量算术',arithmeticLabel],['原五量残差',num(d.residual)+' s⁻¹'],['稳态门槛',num(d.steadyTolerance)+' s⁻¹，连续3个接受端点'],['累计五量账',num(d.fiveQuantityBalanceRelative)],['求解方法','五量SA-noft2 RANS / '+interiorLabel+'、WLS、Sutherland、Fourier；全局隐式BE'],['内部Euler通量',interiorLabel+'。'+interiorNote],['压力源积分',pressureLabel+'。'+pressureNote],['壁应力梯度策略',wallLabel+'。'+wallScope],['黏性功速度策略',workLabel+'。'+workScope],['结果身份','数值研究；稳态、物理精度和生产资格分别验收；离散策略不代表物理精度已验证']]},
    {title:'查看结果 · 研究云图',description:'在同一三维工作台切换压力、热流及可用的温度场，拾取面元查看算法；点击模型或滑块选择任意截面周向曲线。',rows:[['结果状态',d.stage],['显示三角面',String(d.displayFaces)],['原面保守投影',d.sourceFaces+' 原三角面，真曲面积热积分保持'],['温度说明',temperatureNote],['当前范围','零攻角研究端点，新增合格生产面积0']]}
  ];
  function select(index){
    const panel=panels[index];$('stepTitle').textContent=panel.title;$('description').textContent=panel.description;
    $('details').replaceChildren();panel.rows.forEach(([name,value])=>{const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=name;dd.textContent=value;$('details').append(dt,dd);});
    buttons.forEach((button,i)=>button.setAttribute('aria-selected',String(i===index)));
    $('iterationPanel').hidden=index!==3;
    const url=index===4?d.resultsURL:d.geometryURL;if($('viewer').getAttribute('src')!==url)$('viewer').setAttribute('src',url);
  }
  $('stage').textContent=d.stage;$('footer').textContent='单工况冻结研究 · 物理解资格未验收 · 新增生产面积0 · 伪时间不是飞行时间';
  $('rootIterationGroup').hidden=!d.rootIterationHistory;
  if(d.rootIterationHistory)d.rootIterationHistory.forEach(s=>{const tr=document.createElement('tr');[s.iteration,num(s.residual),num(s.boundaryBalance)].forEach(value=>{const td=document.createElement('td');td.textContent=String(value);tr.append(td);});$('rootIterationRows').append(tr);});
  d.steps.forEach(s=>{const tr=document.createElement('tr');[s.number,num(s.residual),num(s.dt),num(s.balance)].forEach(value=>{const td=document.createElement('td');td.textContent=String(value);tr.append(td);});$('iterationRows').append(tr);});
  buttons.forEach((button,i)=>button.addEventListener('click',()=>select(i)));select(4);
})();
