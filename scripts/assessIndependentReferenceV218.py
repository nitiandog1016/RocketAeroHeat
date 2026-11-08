"""模块职责：盘点完整独立参考的实际缺项，绑定准确来源并关闭未满足数值门槛的生产分支。"""
from pathlib import Path
from collections import Counter
import csv,hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v218_validation'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(rel):return json.loads((ROOT/rel).read_text(encoding='utf8'))
def save(p,data):p.write_bytes((json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode('utf8'))

def main():
    start=time.monotonic();out=BASE/'assessment01';out.mkdir(exist_ok=False)
    base='outputs/v214_validation/control01/'
    config=read(base+'lmrsim/config');nk=read(base+'lmrsim/nk-config')
    reconstruct=read('outputs/v218_validation/reconstruct01/report.json')
    assert reconstruct['passed'] and reconstruct['checkedComponents']==96
    # 模块一：确认实际分派选项；只读原记录，不重新生成或复跑算子。
    expected=dict(dimensions=3,viscous=True,viscous_least_squares_type='weighted_qr',
        spatial_deriv_locn='cells',include_boundary_faces_in_spatial_deriv_correction=False,
        flux_calculator='ausmdv',apply_entropy_fix=True,thermo_interpolator='rhou',
        turbulence_model='spalart_allmaras_edwards',apply_shear_stress_relative_limit=False,
        mass_diffusion_model='none',adjust_invalid_cell_data=False)
    assert all(config[k]==v for k,v in expected.items())
    phase=nk['NewtonKrylovPhase_0']
    assert phase['jacobian_interpolation_order']==1 and phase['residual_interpolation_order']==2
    assert not phase['frozen_limiter_for_jacobian'] and not phase['frozen_limiter_for_residual']
    stencil=list(csv.DictReader((ROOT/'outputs/v217_validation/control01/stencil.csv').open()))
    u0=sorted(set((int(r['regime']),int(r['cell'])) for r in stencil if r['kind']=='U0' and r['sample']=='0'))
    assert len(u0)==16
    regimes=[];raw=list(csv.DictReader((ROOT/base/'nk_samples.csv').open()))
    for regime in range(3):
        kinds={}
        for r in raw:
            if int(r['regime'])!=regime:continue
            a=float(r['frozenDerivative']);b=(float(r['frozenPlus'])-float(r['frozenMinus']))/(2e-4)
            err=abs(a-b)/max(1,abs(a),abs(b))
            s=kinds.setdefault(r['kind'],dict(count=0,failed=0,maxScaledError=0.))
            s['count']+=1;s['failed']+=int(err>1e-5);s['maxScaledError']=max(s['maxScaledError'],err)
        assert sum(x['count'] for x in kinds.values())==15918
        assert kinds['residual']['maxScaledError']>1e-5
        regimes.append(dict(regime=regime,oldFiniteDifference=kinds))
    # 模块二：来源只绑定已运行构建；摘录是公式依据，不充当独立数值结果。
    source=read('outputs/v218_validation/source01/report.json');assert source['passed']
    ranges={
        'fluxcalc.d':[(504,515),(553,562),(600,614),(657,691)],
        'fvinterface.d':[(450,464),(500,518),(539,546),(638,648),(784,800),(836,848)],
        'flowstate.d':[],
        'turbulence.d':[(636,660),(679,699),(1021,1031),(1069,1100)],
        'flowgradients.d':[(509,519),(550,574),(620,656),(679,699)],
        'ufluidblock.d':[(621,636),(794,805),(828,846),(1039,1052)],
        'gas_model.d':[],
        'boundary_interface_effect.d':[],
    }
    excerpts=[]
    for row in source['sources']:
        p=ROOT/row['path'];assert sha(p)==row['sha256']
        lines=p.read_text(encoding='utf8').splitlines()
        excerpts.append(dict(path=row['path'],origin=row['origin'],sha256=row['sha256'],
            excerpts=[dict(first=a,last=b,text='\n'.join('%d: %s'%(i+1,lines[i]) for i in range(a-1,b)))
                      for a,b in ranges[p.name]]))
    save(out/'source_contract.json',dict(sources=excerpts,config=expected,
        jacobianInterpolationOrder=1,residualInterpolationOrder=2,
        fixedSupersonicBoundaryUnchangedAcrossRegimes=True,
        evidenceScope='Fixed-source inspection, not a high-precision full residual evaluation'))
    # 模块三：区分“有原始输入、原则上可重建”与“已有完整独立参考”。
    gaps=[
        dict(id='face-state-chain',available='原配置、源代码、供体几何和制造态定义',
             missing='独立的完整一阶左右态、实际边界动作、面平均和rhou物性更新实现及验证'),
        dict(id='convective-branches',available='准确AUSMDV实现与全部原生积分通量/切向输出',
             missing='独立局部坐标通量、平滑声速max、上风/压力切换和熵修正的全量方向及分支界'),
        dict(id='viscous-chain',available='原网格/原生中心法向壁距、局部16工况单元的实际权重',
             missing='全部单元点云/权重的独立重建验证、面梯度修正、湍黏性/热导和黏性通量完整参考'),
        dict(id='full-error-bound',available='独立局部解码、原生梯度/通量驱动的SA源项及散度审核',
             missing='不依赖原生通量/梯度切向的完整残差方向与同h截断/上游误差界'),
    ]
    inputs=[base+'directions.csv',base+'nk_samples.csv',base+'donors.json',base+'lmrsim/config',base+'lmrsim/nk-config',
        base+'lmrsim/grid/grid-0000','outputs/v213_validation/control03/geometry.csv',
        'outputs/v217_validation/control01/stencil.csv','outputs/v217_validation/audit01/report.json',
        'outputs/v218_validation/reconstruct01/report.json','outputs/v218_validation/source01/report.json',
        'data/thermo/phase_aware_dry_air_v1.json','scripts/auditNativeNKAssemblyV214.py',
        'docs/SST_ENERGY_CONSTRAINT_V167.md','docs/HYBRID_OPENFOAM_ENTRY_V198.md',
        'docs/HYBRID_REGIONAL_LAYER_COMPATIBILITY_V208.md',Path(__file__).relative_to(ROOT).as_posix()]
    report=dict(passed=True,meaningOfPassed='Evidence inventory and closure decision completed; numerical gate remains false',
        branch='Eilmer5 complex NK SA-Edwards frozen mapped sequential seven-region production',
        decision='closed-to-production; return to H1 complete mature-combination review',
        knownDefectInSolverProven=False,rawInputsMathematicallyInsufficientProven=False,
        fullIndependentReferenceImplemented=False,fullUpstreamAndTruncationBoundAvailable=False,
        implementationDecision='Do not duplicate the full residual graph as an open-ended project-specific kernel',
        fullStatesReconstructed=1248,observedU0States=16,observedU0Components=96,
        originalScalarOutputs=len(raw),actualLocalRecordKinds=dict(Counter(r['kind'] for r in stencil)),
        unchangedPerturbation=1e-4,unchangedGate=1e-5,regimes=regimes,missingIndependentEvidence=gaps,
        sourceFilesBound=len(excerpts),sourceBytes=sum(r['bytes'] for r in source['sources']),
        nextReviewTarget='Existing SU2 8.6 real NS/RANS path: source-only SA/thermo/physical-boundary/resource combination; not admitted',
        historicalFullDomainSU2PeakBytes=2239016960,historicalPeakIsNotNewConfigurationGuarantee=True,
        fullNumericalAcceptance=False,physicalQualification=False,productionAdmitted=False,
        newNativeQueries=0,newBuilds=0,newMeshes=0,CFDAdvanceRuns=0,newQualifiedAreaM2=0,goalComplete=False,
        inputSHA256={p:sha(ROOT/p) for p in inputs},seconds=time.monotonic()-start)
    save(out/'report.json',report)
    print(json.dumps({k:report[k] for k in ('passed','decision','observedU0Components','originalScalarOutputs','sourceFilesBound','fullNumericalAcceptance')}))

if __name__ == '__main__':
    raise ValueError('v218评估已完成、冻结生产分支关闭；请按NEXT_DECISION登记新任务。')
