"""区分前缘全局热峰和扩张段热峰；独立记录固定50步与最终短窗稳定性。"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from compareEilmerRun14Control import read_wall
from readEilmerRun14Control import need
from runEilmerProcess import sha


def assess(case,receipt,report):
    start=time.monotonic();case=Path(case).resolve();dest=Path(report)
    need(not dest.exists(),'拒绝覆盖')
    run=json.loads(Path(receipt).read_text())
    need(run['case']==str(case) and run['stage']=='run','来源错配')
    for path,digest in run['outputSHA256'].items():need(sha(path)==digest,'源场改变')
    meta=np.atleast_2d(np.loadtxt(case/'lmrsim/loads/loads-metadata',skiprows=1))
    regular=meta[meta[:,2]%50==0]
    pairs=[('last-saved',meta[-2:]),('last-full-50-step',regular[-2:])]
    rows=[]
    for label,pair in pairs:
        need(len(pair)==2,'保存窗不足')
        old,new=[read_wall(case,int(row[0])) for row in pair]
        need(np.array_equal(old[:,:4],new[:,:4]),'几何身份变化')
        regions={}
        for name,mask in [('all-wall',np.ones(len(new),dtype=bool)),('cylinder',new[:,2]==0),('flare',new[:,2]==1)]:
            a=new[mask];b=old[mask]
            regions[name]=dict(qPeakWm2=float(a[:,5].max()),qPeakXM=float(a[np.argmax(a[:,5]),0]),
                qPeakChange=float(abs(a[:,5].max()-b[:,5].max())/max(1e-300,abs(a[:,5].max()))),
                qAreaL1Change=float(np.dot(abs(a[:,5]-b[:,5]),a[:,3])/np.dot(abs(a[:,5]),a[:,3])),
                pAreaL1Change=float(np.dot(abs(a[:,4]-b[:,4]),a[:,3])/np.dot(abs(a[:,4]),a[:,3])),
                heatPowerW=float(np.dot(a[:,5],a[:,3])))
        allr=regions['all-wall'];flare=regions['flare']
        rows.append(dict(window=label,steps=pair[:,2].astype(int).tolist(),regions=regions,
            originalHeatStable=bool(allr['qAreaL1Change']<=.01),originalPressureStable=bool(allr['pAreaL1Change']<=.005),
            flarePeakWithin2Percent=bool(flare['qPeakChange']<=.02)))
    result=dict(windows=rows,physicsQualified=False,newQualifiedAreaM2=0,
        warning='Flare maximum is a regional numerical maximum, not independently validated reattachment heating. Leading-edge maximum must not substitute for it.',
        runReceiptSHA256=sha(receipt),executedSourceSHA256=sha(Path(__file__)),elapsedSeconds=time.monotonic()-start)
    dest.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--receipt',required=True);p.add_argument('--report',required=True)
    a=p.parse_args();print(json.dumps(assess(a.case,a.receipt,a.report)))
