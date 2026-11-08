"""核对v209源码身份、引用与资产保护；审核通过不等于区域求解准入。"""
from pathlib import Path
import hashlib,json,time,shutil

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v209_validation'
COMMIT='f53f4609a0331d48efee69a4e4f3c3598378cc03'

def sha(path):
    """模块一：分块读取大资产，仅核对摘要，不加载或修改网格。"""
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def save(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def main():
    start=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    # 模块二：逐字绑定所有成功提取；同一路径必须始终指向同一固定源码。
    sources={};elapsed=0.
    for name in ['source01','source02','source04']:
        r=json.loads((BASE/name/'report.json').read_text(encoding='utf-8'))
        assert r['passed'] and r['sourceCommit']==COMMIT
        elapsed+=r['elapsedSeconds']
        for row in r['sources']:
            rel=row['relativePath'];p=BASE/name/'source'/rel
            assert sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
            if rel in sources:assert sources[rel]['sha256']==row['sha256']
            sources[rel]=dict(path=p,sha256=row['sha256'],bytes=row['bytes'])
    overlay=BASE/'source04/geometric_length_fluidfvcell.d'
    r=json.loads((BASE/'source04/report.json').read_text())['overlay']
    assert sha(overlay)==r['sha256']==sha(ROOT/r['path'])
    sources['RAH_overlay/fluidfvcell.d']=dict(path=overlay,sha256=r['sha256'],bytes=overlay.stat().st_size)
    # 模块三：保存有准确行号的最小引用，避免再次公开完整第三方源码和大样例。
    selections=[
        ('src/geom/elements/properties.d',499,580,'wedge_properties'),
        ('src/lmr/ufluidblock.d',609,635,'not implemented for unstructured grid'),
        ('src/lmr/ufluidblock.d',726,736,'set_up_workspace_leastsq_via_qr_factorization'),
        ('src/lmr/flowgradients.d',377,408,'dx*dx+dy*dy+dz*dz'),
        ('src/lmr/flowgradients.d',504,658,'r33 = sqrt'),
        ('RAH_overlay/fluidfvcell.d',469,499,'wedge_properties'),
        ('RAH_overlay/fluidfvcell.d',997,1064,'first_face'),
        ('src/lmr/newtonkrylovsolver.d',1707,1725,'foreach (f; cell.iface)'),
        ('src/lmr/newtonkrylovsolver.d',3787,3834,'catch (FlowSolverException e)'),
        ('RAH_overlay/fluidfvcell.d',774,837,'catch (GasModelException err)'),
        ('src/gas/thermo/therm_perf_gas_mix.d',424,470,'MAX_ITERATIONS = 10'),
        ('src/gas/thermo/cea_thermo_curves.d',80,193,'h_low - Cp_low'),
        ('src/gas/uniform_lut.d',194,292,'EXTRAP_MARGIN = 0.2'),
        ('src/gas/uniform_lut_plus_ideal.d',115,128,'Q.u <= 0.0'),
        ('src/gas/gas_model.d',332,344,'Q.T = T_MIN'),
        ('src/lmr/lua-modules/bc.lua',821,871,'UpdateThermoTransCoeffs'),
        ('src/lmr/lua-modules/bc.lua',1570,1641,'OutFlowBC_Simple = OutFlowBC_SimpleFlux'),
        ('src/lmr/bc/boundary_flux_effect.d',575,644,'face.F[cqi.mass] = 0.0'),
        ('src/lmr/bc/ghost_cell_effect/fixed_p.d',40,79,'ghost0.fs.gas.p = p_outside'),
        ('src/lmr/bc/ghost_cell_effect/from_upwind.d',47,72,'dot(inside.fs.vel, f.n)'),
        ('src/lmr/turbulence.d',1021,1102,'class saeTurbulenceModel'),
        ('src/nm/brent.d',38,82,'Root must be bracketed'),
    ]
    refs=[]
    for index,(rel,lo,hi,needle) in enumerate(selections,1):
        row=sources[rel];lines=row['path'].read_text(encoding='utf-8').splitlines()
        excerpt='\n'.join(lines[lo-1:hi])+'\n'
        assert needle in excerpt,(rel,lo,hi,needle)
        name='%02d_%s.txt'%(index,Path(rel).name)
        (out/name).write_text(excerpt,encoding='utf-8')
        refs.append(dict(id=index,path=rel,startLine=lo,endLine=hi,sourceSHA256=row['sha256'],
            excerptFile=name,excerptSHA256=sha(out/name),
            url=None if rel.startswith('RAH_overlay/') else 'https://github.com/gdtk-uq/gdtk/blob/'+COMMIT+'/'+rel+'#L'+str(lo)))
    save(out/'source_references.json',refs)
    # 模块四：保留原物性、真实层和并发结果，不以引用核对授予CFD/热资格。
    contract=json.loads((ROOT/'docs/validation/v208/route_contract.json').read_text())
    assert sha(ROOT/contract['regionLayerMesh'])==contract['regionLayerMeshSHA256']
    assert contract['originalFullFaces']==19752 and contract['layers']==32
    thermo=ROOT/'data/thermo/phase_aware_dry_air_v1.json'
    assert sha(thermo)=='7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29
    assert all(sha(ROOT/p)==digest for p,digest in protected.items())
    assert json.loads((BASE/'source03/failure.json').read_text())['passed'] is False
    assert (BASE/'source03/executed_entry.py').is_file()
    report=dict(passed=True,scope='source identity, citations and preservation only',sourceCommit=COMMIT,
        uniqueUpstreamFiles=len(sources)-1,overlayFiles=1,referenceCount=len(refs),
        successfulExtractionSeconds=elapsed,failedExtractionEntryWallSeconds=1.1719234,
        reviewTotalWallSeconds=None,elapsedSeconds=time.monotonic()-start,
        protectedFilesUnchanged=29,originalLayerMeshSHA256=contract['regionLayerMeshSHA256'],
        caloricTableSHA256=sha(thermo),existingCombinationAccepted=False,
        builds=0,CFDRuns=0,newQualifiedHeatAreaM2=0,physicalQualification=False,
        sourceManifest=[dict(path=k,sha256=v['sha256'],bytes=v['bytes']) for k,v in sorted(sources.items())])
    assert report['elapsedSeconds']<60
    save(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='sourceManifest'}))

if __name__=='__main__':main()
