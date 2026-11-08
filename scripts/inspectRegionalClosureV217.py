"""模块职责：绑定最小观察的准确源，并核算外包络几何/拓扑设计；不生成网格或原生结果。"""
from pathlib import Path
import hashlib,json,math,time
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v217_validation'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d): p.write_bytes((json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n').encode('utf8'))

def main():
    start=time.monotonic(); out=BASE/'source01';out.mkdir(parents=True,exist_ok=False)
    build=json.loads((ROOT/'outputs/v214_validation/build01/report.json').read_text())
    definitions=[
      ('outputs/v214_validation/prepared01/nk_probe.d','/mnt/e/RocketAeroHeat/outputs/v214_validation/prepared01/nk_probe.d',[(113,158),(179,212)]),
      ('outputs/v213_validation/source02/source/src/lmr/newtonkrylovsolver.d','/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/newtonkrylovsolver.d',[(3514,3539),(3557,3635)]),
      ('outputs/v211_validation/source01/source/lmr/flowgradients.d','/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/flowgradients.d',[(34,48),(509,578),(626,699)]),
      ('outputs/v212_validation/source01/source/src/lmr/ufluidblock.d','/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/ufluidblock.d',[(621,637),(726,737)]),
      ('outputs/v212_validation/source01/source/src/lmr/bc/ghost_cell_effect/extrapolate_copy.d','/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/bc/ghost_cell_effect/extrapolate_copy.d',[(43,86)]),
      ('outputs/v212_validation/source01/source/src/lmr/bc/ghost_cell_effect/flow_state_copy.d','/opt/rocket-aeroheat/eilmer-v5.0.0/source/src/lmr/bc/ghost_cell_effect/flow_state_copy.d',[(60,94)]),
    ]
    rows=[]
    for rel,origin,ranges in definitions:
        p=ROOT/rel;assert sha(p)==build['inputSHA256'][origin],rel
        lines=p.read_text(encoding='utf8').splitlines()
        rows.append(dict(path=rel,buildSource=origin,sha256=sha(p),excerpts=[dict(first=a,last=b,text='\n'.join('%d: %s'%(j+1,lines[j]) for j in range(a-1,b))) for a,b in ranges]))
    save(out/'source_excerpts.json',dict(passed=True,buildSHA256=sha(ROOT/'outputs/v214_validation/build01/report.json'),sources=rows))
    inv=json.loads((ROOT/'outputs/v216_validation/inventory.json').read_text())
    decision=json.loads((ROOT/'outputs/v216_validation/final01/decision.json').read_text())['preferredBoundaryDesign']
    a,b=inv['domainBounds'];xmin,xmax=a[0],b[0];half=b[1];length=inv['sourceBounds'][1][0]-inv['sourceBounds'][0][0]
    h=xmax-xmin;cube=h*(2*half)**2
    cf,cb=decision['frontRadiusM'],decision['backRadiusM'];ch=decision['frontXM']-decision['backXM']
    cone=math.pi*ch*(cf*cf+cf*cb+cb*cb)/3
    # 同40度法向条件的方截锥：六平面可与旧方盒标签逐片对应，前后均留正壳厚。
    front=xmax+.1*length;back=xmin-.1*length;front_half=1.05*half
    back_half=front_half+(front-back)*math.tan(math.radians(40))
    square=4*(front-back)*(front_half**2+front_half*back_half+back_half**2)/3
    rows_mn=[-5*math.cos(math.radians(20)),5*math.cos(math.radians(20)),
             -5*math.sin(math.radians(40))*math.cos(math.radians(20)),
             -5*math.sin(math.radians(60)),-5*math.sin(math.radians(20))]
    # 壳拓扑仅组合设计：每个旧外三角对应一棱柱，闭合三角面E=3F/2、V=F/2+2。
    F=6572;E=3*F//2;V=F//2+2
    counts=dict(cells=F,faces=2*F+E,vertices=2*V,ghosts=2*F,internalFaces=E,
                cellVertexLinks=6*F,cellFaceLinks=5*F,faceVertexLinks=6*F+4*E)
    from eilmerNKSAResourceBudget import lower_bound,EXPECTED_SIZES
    payload=lower_bound(counts,EXPECTED_SIZES)['payloadLowerBoundBytes']
    report=dict(passed=True,sourceFilesBound=len(rows),savedStencilAndWeightsAvailable=False,
        formalJvClearsPrimitiveImaginaryBeforeReturn=True,
        proposedObservation='isolated helper: same native settings/U1/decode/evalResidual; capture before cleanup; compare selected outputs to saved v214; no resident rerun',
        originalCubeVolumeM3=cube,circularConeVolumeM3=cone,circularAddedVolumeM3=cone-cube,
        planarCandidate=dict(status='continuous-design-only',frontXM=front,backXM=back,frontHalfWidthM=front_half,backHalfWidthM=back_half,
                             halfAngleDeg=40,volumeM3=square,addedVolumeM3=square-cube,initialNormalMachValues=rows_mn,
                             proposedSingleShellCounts=counts,singleShellPayloadLowerBoundBytes=payload,
                             singleShellPayloadLowerBoundGiB=payload/1024**3,
                             shellCellsGenerated=False,positiveJacobiansProven=False,meshQualityProven=False,
                             proposedLayerCountIsNotVerifiedResolution=True,peakMemoryProven=False),
        physicalBoundaryAccepted=False,fullNumericalAcceptance=False,CFDAdmitted=False,seconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(report))

if __name__ == '__main__':
    raise ValueError('v217观察与审核已结束；当前无构建/查询/网格/CFD准入。')
