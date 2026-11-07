# -*- coding: utf-8 -*-
"""独立原二进制/CSV复核与错误来源拒绝；不写原始场或准备输入。"""
from pathlib import Path
import argparse
import copy
import json
import sys
import shutil
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from readHybridEulerNative import read_wall_frame
from readEilmerNative import sha


def verify(base, final, out):
    start=time.monotonic();base,final,out=map(Path,(base,final,out));out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_verify.py')
    case=base/'control01';prep=base/'preparation01/report.json';run=base/'run01/report.json'
    source=Path(__file__).resolve().parents[1]/'scripts/readHybridEulerNative.py'
    shutil.copyfile(source,out/source.name)
    result=dict(passed=False);initial=sha(source)
    try:
        # 模块一：直接按已审核布局读原始double；独立网格索引，不复用映射器求owner。
        raw=np.fromfile(case/'lmrsim/snapshots/0002/fluid-0000',dtype='<f8').reshape(15,2400)
        owners=np.repeat(np.arange(120)*20,2)
        rows=np.genfromtxt(final/'wall_comparison.csv',delimiter=',',names=True,dtype=None,encoding='ascii')
        checks={}
        for csvkey,index in [('p_Pa',8),('T_K',13),('Ux_mps',5),('Uy_mps',6),('Uz_mps',7)]:
            diff=float(np.max(abs(rows[csvkey]-raw[index,owners])));checks[csvkey]=diff;assert diff==0
        assert np.array_equal(rows['faceID'],np.arange(1,241))
        assert np.max(abs(rows['Uz_mps']))>150 and np.max(abs(rows['Uy_mps']))==0
        surf=json.loads((case/'surface.json').read_text());xyz=np.asarray(surf['vertices'])[np.asarray(surf['facesOneBased'])-1]
        a=np.linalg.norm(xyz[:,1]-xyz[:,0],axis=1);b=np.linalg.norm(xyz[:,2]-xyz[:,1],axis=1);c=np.linalg.norm(xyz[:,0]-xyz[:,2],axis=1)
        semi=(a+b+c)/2;area=np.sqrt(semi*(semi-a)*(semi-b)*(semi-c))
        assert np.max(abs(area-rows['area_m2']))<1e-16
        h0=raw[12,owners]+raw[8,owners]/raw[4,owners]+np.sum(raw[5:8,owners]**2,axis=0)/2
        assert np.max(abs(h0-rows['h0_Jkg'])/np.maximum(1,abs(h0)))<2e-10
        # 模块二：只改小回执夹具，拒绝错误/过期/失败来源，不修改任何原始场。
        rejected=[]
        def reject(name,fn):
            try:fn()
            except (ValueError,FileNotFoundError):rejected.append(name)
            else:raise AssertionError('未拒绝: '+name)
        reject('missing-run-receipt',lambda:read_wall_frame(case,2,prep))
        reject('negative-snapshot',lambda:read_wall_frame(case,-1,prep,run))
        reject('boolean-snapshot',lambda:read_wall_frame(case,True,prep,run))
        reject('missing-snapshot',lambda:read_wall_frame(case,9,prep,run))
        for name,change in [('wrong-preparation-case',lambda d:d.update(case=d['case']+'-wrong')),
                            ('failed-preparation',lambda d:d.update(passed=False)),
                            ('changed-preparation-input',lambda d:d['nativeInputSHA256'].update({next(iter(d['nativeInputSHA256'])):'0'*64}))]:
            d=json.loads(prep.read_text());change(d);p=out/(name+'.json');p.write_text(json.dumps(d))
            reject(name,lambda p=p:read_wall_frame(case,2,p,run))
        for name,target,change in [
            ('changed-output','report',lambda d:d['outputSHA256'].update({str(case.resolve()).replace('\\','/').replace('E:/','/mnt/e/',1)+'/lmrsim/snapshots/0002/fluid-0000':'0'*64})),
            ('failed-run','report',lambda d:d.update(processPassed=False)),
            ('wrong-binary','native',lambda d:d['argv'].__setitem__(0,'wrong-runner')),
            ('unstable-run-input','native',lambda d:d.update(inputsStable=False)),
            ('live-process','native',lambda d:d.update(liveGroupAfterCleanup={'1':{}}))]:
            folder=out/name;(folder/'native').mkdir(parents=True)
            d=json.loads(run.read_text());n=json.loads((run.parent/'native/report.json').read_text())
            change(d if target=='report' else n)
            p=folder/'report.json';p.write_text(json.dumps(d));(folder/'native/report.json').write_text(json.dumps(n))
            reject(name,lambda p=p:read_wall_frame(case,2,prep,p))
        valid,_=read_wall_frame(case,2,prep,run)
        assert valid['algebraicInterfacePassed'] and not valid['engineeringHeatReady']
        result.update(passed=True,rejections=rejected,csvDirectBinaryMaxAbsoluteDifferences=checks,
            independentHeronMaxAreaDifferenceM2=float(np.max(abs(area-rows['area_m2']))),
            independentTotalEnthalpyMaxRelativeError=float(np.max(abs(h0-rows['h0_Jkg'])/abs(h0))),
            physicalQualification=False)
    except Exception as error:result['error']=repr(error)
    result.update(elapsedSeconds=time.monotonic()-start,readerStable=sha(source)==initial,readerSHA256=initial)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
    return result['passed']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--final',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();raise SystemExit(0 if verify(a.base,a.final,a.out) else 1)
