"""一次绑定固定成熟内核的各向异性/气体/边界源码，不启动求解器。"""
from pathlib import Path
import json,re,shutil,subprocess,time,sys
from prepareOpenFoamControl import ROOT,save,sha
from runEilmerProcess import PREFIX,COMMIT

def main():
    start=time.monotonic();out=ROOT/'outputs/v209_validation/source01';out.mkdir(parents=True,exist_ok=False);src=PREFIX/'source'
    assert subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()==COMMIT
    paths=[];index=[]
    for folder in ['src/gas','src/lmr/bc','src/lmr/turbulence']:
        for p in (src/folder).rglob('*'):
            if not p.is_file() or p.suffix not in ['.d','.lua']:continue
            text=p.read_text(errors='replace')
            hits=[line.strip()[:180] for line in text.splitlines() if re.search(r'class .*Gas|class .*Thermo|tabulat|look.up|UserDefinedGas|Spalart|Edwards|class .*OutFlow|class .*InFlow|Extrapolate',line,re.I)]
            if hits:index.append(dict(path=p.relative_to(src).as_posix(),hits=hits[:12]))
    for p in [src/'COPYING',src/'src/lmr/ufluidblock.d',src/'src/lmr/fluidfvcell.d',src/'src/lmr/viscousflux.d',src/'src/lmr/flowstate.d']:
        if p.is_file():paths.append(p)
    for row in index:
        p=src/row['path']
        if any(k in p.name.lower() for k in ['gas_model','therm','lookup','lut','user_defined','composite','edwards','outflow','extrapolate']):paths.append(p)
    for folder in ['src/lmr','src/nm']:
        for p in (src/folder).glob('*.d'):
            if any(k in p.name.lower() for k in ['wlsq','least','interp','deriv','turb','boundary']):paths.append(p)
    for rel in ['examples/lmr/2D/turbulent-flat-plate/README.adoc','examples/lmr/2D/turbulent-flat-plate/unstructured/fp.lua']:
        paths.append(src/rel)
    rows=[]
    for p in sorted(set(paths)):
        assert time.monotonic()-start<60
        rel=p.relative_to(src);raw=p.read_bytes();fixed=subprocess.check_output(['git','-C',str(src),'show',COMMIT+':'+rel.as_posix()])
        assert raw==fixed,rel
        q=out/'source'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(raw)
        rows.append(dict(path=str(p),relativePath=rel.as_posix(),sha256=sha(p),bytes=len(raw),matchesFixedCommit=True))
    save(out/'interface_index.json',index);shutil.copyfile(__file__,out/'executed_entry.py')
    assert sum(r['bytes'] for r in rows)<16*1024**2
    save(out/'report.json',dict(passed=True,sourceCommit=COMMIT,sources=rows,elapsedSeconds=time.monotonic()-start,builds=0,CFDRuns=0,physicalQualification=False))
    print(json.dumps(dict(files=len(rows),bytes=sum(r['bytes'] for r in rows),seconds=time.monotonic()-start)))

def supplement(closure=False):
    """补齐同一次来源审核的真实黏性梯度、边界和气体接受实现。"""
    start=time.monotonic();out=ROOT/('outputs/v209_validation/source04' if closure else 'outputs/v209_validation/source02');out.mkdir(parents=True,exist_ok=False)
    src=PREFIX/'source'
    assert subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()==COMMIT
    paths=['COPYING','src/lmr/flowgradients.d','src/lmr/fvinterface.d',
           'src/lmr/lua-modules/bc.lua','src/lmr/bc/boundary_interface_effect.d',
           'src/lmr/bc/boundary_condition.d','src/lmr/bc/ghost_cell_effect/fixed_p.d',
           'src/lmr/bc/ghost_cell_effect/fixed_pt.d','src/lmr/bc/ghost_cell_effect/from_upwind.d',
           'src/nm/brent.d','src/gas/thermo/therm_perf_gas_mix.d']
    if closure:
        paths=['src/lmr/newtonkrylovsolver.d','src/lmr/bc/boundary_flux_effect.d',
               'src/geom/elements/properties.d','src/gas/gas_state.d']
    rows=[]
    for rel in paths:
        p=src/rel;raw=p.read_bytes()
        fixed=subprocess.check_output(['git','-C',str(src),'show',COMMIT+':'+rel],timeout=10)
        assert raw==fixed,rel
        q=out/'source'/rel;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(raw)
        rows.append(dict(relativePath=rel,sha256=sha(p),bytes=len(raw),matchesFixedCommit=True))
        assert time.monotonic()-start<60
    # 已交付几何覆盖层另列身份，不能误称未经修改的上游源码。
    p=ROOT/'tools/runtime/eilmer-geometric-length-v1/overlay/src/lmr/fluidfvcell.d'
    shutil.copyfile(p,out/'geometric_length_fluidfvcell.d')
    overlay=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,matchesFixedCommit=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    elapsed=time.monotonic()-start
    prior=json.loads((out.parent/'source01/report.json').read_text())['elapsedSeconds']
    if closure:prior+=json.loads((out.parent/'source02/report.json').read_text())['elapsedSeconds']
    assert elapsed+prior<120 and sum(r['bytes'] for r in rows)+overlay['bytes']<16*1024**2
    save(out/'report.json',dict(passed=True,scope='source identity only',sourceCommit=COMMIT,sources=rows,
        overlay=overlay,elapsedSeconds=elapsed,combinedExtractionSeconds=prior+elapsed,
        builds=0,CFDRuns=0,physicalQualification=False))
    print(json.dumps(dict(files=len(rows),seconds=elapsed)))

if __name__=='__main__':
    raise ValueError('v209定向来源审核已完成；现成组合关闭，不重提取/构建/网格/CFD。下一项见NEXT_DECISION。')
