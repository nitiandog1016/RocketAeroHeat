# -*- coding: utf-8 -*-
"""独立几何体积/边界及交换拒绝检查；合成场不冒充原生Euler结果。"""
from pathlib import Path
import argparse,copy,hashlib,json,math,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepareHybridEulerCase import prepare
from hybridExternalState import exchange,fingerprint


def test(out):
    started=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    manifest=prepare(out/'control')
    # 顶边来流法向Mach为0，不能因总Mach3就给它超声速全状态入口。
    assert 'outer=InOutFlowBC_Ambient:new{flowState=inflow}' in (out/'control/job.lua').read_text(encoding='utf-8')
    # 模块一：独立读取生成网格，六四面体分解检查正体积及完整配对。
    lines=iter((out/'control/mesh.su2').read_text().splitlines())
    def head(name):
        key,val=next(lines).split('=');assert key.strip()==name;return int(val)
    assert head('NDIME')==3
    n=head('NELEM');cells=np.array([[int(s) for s in next(lines).split()][1:] for _ in range(n)])
    count=head('NPOIN');vertices=np.array([[float(s) for s in next(lines).split()[:3]] for _ in range(count)])
    sides={};marker_count=head('NMARK')
    for _ in range(marker_count):
        key,name=next(lines).split('=');assert key=='MARKER_TAG';name=name.strip()
        count=head('MARKER_ELEMS');sides[name]=[[int(s) for s in next(lines).split()[1:]] for _ in range(count)]
    assert next(lines,None) is None
    volume=0.;minimum=float('inf')
    for tet in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
        a,b,c,d=[vertices[cells[:,i]] for i in tet]
        volumes=np.einsum('ij,ij->i',b-a,np.cross(c-a,d-a))/6
        minimum=min(minimum,float(volumes.min()));volume+=volumes.sum()
    assert minimum>0
    exact=.04*(.6-(.5*.4*.4+.4*.4)*math.tan(math.radians(10)))
    assert abs(volume-exact)<1e-14
    assert abs(manifest['geometry']['analyticalVolumeM3']-exact)<1e-14
    owners={}
    face_orders=((0,1,2,3),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
    for c in cells:
        for face in face_orders:
            key=tuple(sorted(c[list(face)]));owners[key]=owners.get(key,0)+1
    assert all(n in (1,2) for n in owners.values())
    boundary=[tuple(sorted(face)) for group in sides.values() for face in group]
    assert len(boundary)==len(set(boundary))
    assert set(boundary)=={k for k,n in owners.items() if n==1}
    assert not manifest['executionAuthorized'] and not manifest['nativePrepared']
    # 模块二：非等面积原面，非零两个横向速度，同一理想气体独立真值。
    surface=dict(vertices=[[0,0,0],[1,0,0],[0,1,0],[2,0,0],[0,3,0]],
                 facesOneBased=[[1,2,3],[1,4,3],[1,2,5]],faceIDs=[1,2,3])
    condition=dict(case='synthetic-only',alphaDeg=20)
    model=dict(name='synthetic-caloric-control',R=287.,cp=1004.5,temperatureRangeK=[200,2000],
               energyReference='h=cp*T',entropyReference='cp*log(T/300)-R*log(p/100000)')
    p=np.array([2000,5000,9000.]);T=np.array([300,500,800.]);u=np.array([[-900,20,100],[-600,-40,200],[-450,80,250.]])
    def thermo(pressure,temperature):
        return dict(densityKgm3=pressure/(287*temperature),internalEnergyJkg=(1004.5-287)*temperature,
                    enthalpyJkg=1004.5*temperature,entropyJkgK=1004.5*np.log(temperature/300)-287*np.log(pressure/1e5))
    fields={k:v.tolist() for k,v in thermo(p,T).items()}
    fields.update(pressurePa=p.tolist(),temperatureK=T.tolist(),velocityMps=u.tolist(),
                  totalEnthalpyJkg=(1004.5*T+(u*u).sum(axis=1)/2).tolist())
    source=out/'synthetic_source.json';source.write_text(json.dumps(fields)+'\n')
    sources={str(source.resolve()):hashlib.sha256(source.read_bytes()).hexdigest()}
    xyz=np.array(surface['vertices'])[np.array(surface['facesOneBased'])-1]
    frame=dict(schema='RAH.HybridEulerFrame.v1',solver='Euler',turbulenceModel='none',snapshot=2,
               surfaceSHA256=fingerprint(surface),conditionSHA256=fingerprint(condition),modelSHA256=fingerprint(model),
               sourceFilesSHA256=sources,faceIDs=[1,2,3],faceCentersM=xyz.mean(axis=1).tolist(),
               originalAreaM2=3.,fields=fields)
    result=exchange(frame,surface,condition,model,sources,thermo)
    assert result['algebraicInterfacePassed'] and result['wallHeatFluxWm2']==[None]*3
    assert not result['engineeringHeatReady'] and not result['physicalQualification']
    rejections=[]
    def reject(name,modify):
        altered=copy.deepcopy(frame);modify(altered)
        try: exchange(altered,surface,condition,model,sources,thermo)
        except ValueError:rejections.append(name);return
        raise AssertionError('not rejected: '+name)
    for key in ('surfaceSHA256','conditionSHA256','modelSHA256'):
        reject(key,lambda f,k=key:f.update({k:'0'*64}))
    reject('reorderedFaces',lambda f:f.update(faceIDs=[2,1,3]))
    reject('duplicateFaces',lambda f:f.update(faceIDs=[1,1,3]))
    reject('missingFace',lambda f:f.update(faceIDs=[1,2]))
    reject('changedArea',lambda f:f.update(originalAreaM2=2.9))
    reject('wrongSolver',lambda f:f.update(solver='RANS'))
    reject('wrongTurbulence',lambda f:f.update(turbulenceModel='SST'))
    reject('wrongSnapshot',lambda f:f.update(snapshot=-1))
    reject('wrongDensity',lambda f:f['fields']['densityKgm3'].__setitem__(0,1.))
    reject('missingLateralKineticEnergy',lambda f:f['fields']['totalEnthalpyJkg'].__setitem__(0,1004.5*T[0]+u[0,0]**2/2))
    reject('wrongEnthalpyReference',lambda f:f['fields']['enthalpyJkg'].__setitem__(0,1.))
    reject('temperatureOutsideModel',lambda f:f['fields']['temperatureK'].__setitem__(0,150.))
    reject('nonfinite',lambda f:f['fields']['pressurePa'].__setitem__(0,float('nan')))
    reject('wrongSource',lambda f:f.update(sourceFilesSHA256={}))
    source.write_text('changed fixture\n')
    try: exchange(frame,surface,condition,model,sources,thermo)
    except ValueError: rejections.append('changedSourceBytes')
    else: raise AssertionError('changed source accepted')
    # 保留修改前后的夹具字节；它们不是真实计算数据。
    (out/'synthetic_source_original.json').write_text(json.dumps(fields)+'\n')
    (out/'exchanged_synthetic.json').write_text(json.dumps(result,indent=2)+'\n')
    report=dict(passed=True,cells=len(cells),vertices=len(vertices),uniqueFaces=len(owners),
                boundaryFaces=len(boundary),minimumTetVolumeM3=minimum,volumeM3=float(volume),
                exactVolumeM3=exact,volumeErrorM3=abs(float(volume)-exact),
                rejectionCases=rejections,actualEulerRun=False,nativePreparation=False,
                newPhysicalQualification=False,elapsedSeconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    test(parser.parse_args().out)
