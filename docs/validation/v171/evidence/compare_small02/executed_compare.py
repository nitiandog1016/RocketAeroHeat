"""模块用途：比较实际四种边界组合，逐列核对旧数值路径及非均匀壁温实际施加。"""
import argparse
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2ThermalRuntime import read,sha
from su2NACALaminarInterface import load_state
from validateIdealAir import save


def compare(base,out):
    out.mkdir(parents=True,exist_ok=False);(out/'executed_compare.py').write_bytes(Path(__file__).read_bytes())
    names=('observed_scalar01','spatial_scalar01','spatial_constant01','spatial_profile01')
    rows=[];locks={str(Path(__file__).resolve()):sha(Path(__file__).resolve())}
    for i,name in enumerate(names):
        folder=base/name;r=read(folder/'engine_process.json');a=read(base/('audit_small%02d'%(i+1))/'report.json')
        assert r['exitCode']==0 and all(r['inputsStable'].values()) and a['interfacePassed']
        assert [p['localIteration'] for p in r['checkpoints']['complete']]==[0,2,3]
        locks[str(folder/'engine_process.json')]=sha(folder/'engine_process.json')
        locks[str(base/('audit_small%02d'%(i+1))/'report.json')]=sha(base/('audit_small%02d'%(i+1))/'report.json')
        rows.append(dict(variant=name,nativeSeconds=r['elapsedSeconds'],
            wallHeatMaximumScaledError=max(f['nativeHeatMaximumScaledError'] for f in a['frames']),
            energyMaximumRelativeError=max(f['energyMaximumRelativeError'] for f in a['frames']),
            eosMaximumRelativeError=max(f['eosMaximumRelativeError'] for f in a['frames']),
            powerMaximumScaledError=max(f['powerScaledDifference'] for f in a['frames']),
            appliedTemperatureCheck=a['appliedTemperatureCheck']))
    comparisons=[]
    for target in names[1:3]:
        for step in (0,2,3):
            filename='restart_flow_%06d.csv'%step
            pa=base/names[0]/filename;pb=base/target/filename
            a,_=load_state(pa);b,_=load_state(pb)
            err=float(np.max(abs(a-b)/np.maximum(1,abs(a))))
            assert err<2e-13
            comparisons.append(dict(target=target,localIteration=step,maximumScaledDifference=err,byteEqual=pa.read_bytes()==pb.read_bytes()))
            locks[str(pa)]=sha(pa);locks[str(pb)]=sha(pb)
    a,_=load_state(base/names[0]/'restart_flow_000003.csv');b,_=load_state(base/names[3]/'restart_flow_000003.csv')
    difference=float(np.max(abs(a-b)));assert difference>0
    assert rows[2]['appliedTemperatureCheck']['maximumGivenTwDifferenceK']<1e-10
    assert rows[3]['appliedTemperatureCheck']['maximumGivenTwDifferenceK']<1e-10
    assert all(sha(Path(p))==h for p,h in locks.items())
    report=dict(passed=True,variants=rows,comparisons=comparisons,spatialProfileActualFieldDifference=difference,
        inputSHA256=locks,inputsStable=True,physicalQualified=False,readyForRegisteredValidationCandidate=True)
    save(out/'report.json',report);print(dict(passed=True,comparisons=len(comparisons),allComparedBytesEqual=all(r['byteEqual'] for r in comparisons)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--base',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();compare(Path(args.base).resolve(),Path(args.output).resolve())
