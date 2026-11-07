"""核对真实四步小域与独立几何/热流真值；不重复启动CFD。"""
from pathlib import Path
import argparse
import copy
import csv
import json
import shutil
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import audit_prepared,read_field,export_loads,sha


def verify(case,run_receipt,output):
    started=time.monotonic();case=Path(case).resolve();out=Path(output).resolve()
    out.mkdir(parents=True,exist_ok=False)
    report,grid,mapping=audit_prepared(case)
    f0=read_field(case,0,4);f1=read_field(case,1,4)
    # 模块一：解析体积真值，分别核对四类单元，不能只比较总和。
    expected=np.array([1,1/3,1/6,5/8])
    assert np.max(np.abs(f0['vol']-expected))<2e-15
    assert np.all(f1['vol']>0) and np.isfinite(f1['dt_local']).all()
    assert np.all(f1['T']>0) and np.all(f1['rho']>0) and np.all(f1['p']>0)
    old=read_field(case.parent/'native02',0,4)
    assert old['vol'][3]==-.625 and np.isnan(old['dt_local']).all()
    # 本小例早于自动输出绑定实现；首次只读封存明示时间点，原回执不改写、不补跑。
    receipt=json.loads(Path(run_receipt).read_text())
    receipt['outputSHA256']={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()}
    receipt['outputBinding']='first-read-after-run; original run receipt preserved; no earlier hash claim'
    receipt['originalRunReceiptSHA256']=sha(run_receipt)
    seal=out/'observed_run.json';seal.write_text(json.dumps(receipt,indent=2)+'\n')
    loads=case/'lmrsim/loads/0000/blk-0000-bndry-0-wall.dat'
    actual=export_loads(case,loads,out/'actual-wall.csv',seal)
    # 模块二：不调用生产符号公式，按Fourier向外流体域法向直接形成真值。
    a=np.genfromtxt(loads,names=True,deletechars='');rows=list(csv.DictReader((out/'actual-wall.csv').open()))
    independent_power=0.;signs=set()
    for face,row in enumerate(rows):
        j=int(row['nativeBoundaryRow']);sign=int(a['outsign'][j]);signs.add(sign)
        expected_heat=a['q_total'][j] if sign==-1 else -a['q_total'][j]
        assert abs(float(row['heatFluxIntoWall_W_m2'])-expected_heat)<2e-12
        assert int(row['originalFaceId'])==face+1
        assert row['physicsQualified']=='False' and row['status']=='research'
        independent_power+=float(row['area_m2'])*expected_heat
    assert signs=={-1,1} and abs(independent_power-actual['heatPowerIntoWallW'])<2e-12
    # 模块三：在独立小副本内制造线性温度梯度，覆盖双向法向、正/负/零热流。
    synth=out/'synthetic';shutil.copytree(case,synth)
    sl=synth/'lmrsim/loads/0000/blk-0000-bndry-0-wall.dat'
    data=np.loadtxt(sl,skiprows=1);names=sl.read_text().splitlines()[0].split();col={k:i for i,k in enumerate(names)}
    grad=np.array([2.,-3.,5.]);conductivity=.1
    normals=data[:,[col[k] for k in ('n.x','n.y','n.z')]]
    raw=conductivity*(normals@grad);raw[0]=0
    data[:,col['q_total']]=raw;data[:,col['q_cond']]=raw;data[:,col['q_diff']]=0
    np.savetxt(sl,data,header=' '.join(names),comments='',fmt='%.17e')
    sr=copy.deepcopy(receipt);sr['case']=str(synth)
    sr['outputSHA256']={str(p):sha(p) for p in (synth/'lmrsim').rglob('*') if p.is_file()}
    sr['outputBinding']='synthetic test only';sp=out/'synthetic_run.json';sp.write_text(json.dumps(sr))
    signed=export_loads(synth,sl,out/'synthetic-wall.csv',sp)
    sn=-normals*data[:,col['outsign'],None] # 固体外向单位法向，与流体域外向相反。
    fourier_heat=conductivity*(sn@grad);fourier_heat[0]=0
    synthrows=list(csv.DictReader((out/'synthetic-wall.csv').open()))
    assert max(abs(float(r['heatFluxIntoWall_W_m2'])-fourier_heat[int(r['nativeBoundaryRow'])]) for r in synthrows)<2e-15
    assert min(signed[k] for k in ('positiveHeatFaces','negativeHeatFaces','zeroHeatFaces'))>0
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError,FileNotFoundError):rejected.append(name)
        else:raise AssertionError('未拒绝: '+name)
    reject('missing_run',lambda:export_loads(case,loads,out/'no.csv',None))
    for key,value,label in [('passed',False,'failed_run'),('case','/wrong/case','wrong_case'),
                            ('caseManifestSHA256','0'*64,'wrong_manifest'),('outputSHA256',{},'unbound_loads')]:
        bad=copy.deepcopy(receipt);bad[key]=value;p=out/(label+'.json');p.write_text(json.dumps(bad))
        reject(label,lambda p=p:export_loads(case,loads,out/'no.csv',p))
    reject('old_wedge_case',lambda:audit_prepared(case.parent/'native02'))
    reject('overwrite_csv',lambda:export_loads(case,loads,out/'actual-wall.csv',seal))
    # 合法来源也不能绕过壁温、法向、分项和未知值检查。
    original=sl.read_bytes()
    for key,value in [('T',301.),('outsign',0.),('area',-1.),('q_diff',100.),('p',float('nan'))]:
        changed=data.copy();changed[0,col[key]]=value
        np.savetxt(sl,changed,header=' '.join(names),comments='',fmt='%.17e')
        sr['outputSHA256'][str(sl)]=sha(sl);sp.write_text(json.dumps(sr))
        reject('invalid_'+key,lambda:export_loads(synth,sl,out/'no.csv',sp))
    sl.write_bytes(original);sr['outputSHA256'][str(sl)]=sha(sl);sp.write_text(json.dumps(sr))
    sl.write_bytes(original+b'\n')
    reject('changed_native_output',lambda:export_loads(synth,sl,out/'no.csv',sp))
    sl.write_bytes(original)
    r=dict(passed=True,actual=actual,synthetic=signed,rejected=rejected,
           oldNegativeWedgeVolumeM3=float(old['vol'][3]),correctedCellVolumesM3=f0['vol'].tolist(),
           finalTemperatureRangeK=[float(f1['T'].min()),float(f1['T'].max())],
           nativeRuns=0,physicsQualified=False,newQualifiedAreaM2=0,
           elapsedSeconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');return r


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--run-receipt',required=True)
    p.add_argument('--output',required=True);a=p.parse_args()
    print(json.dumps(verify(a.case,a.run_receipt,a.output)))
