"""模块用途。独立状态构造、真实CVariable/CNSVariable夹具及原物性参照，有限小检查。"""
import csv
import json
from pathlib import Path
import sys
import numpy as np
from idealAirProperties import IdealAir
from prepareSU2SSTEnergyPatch import ROOT,SOURCE,BUILD,overlay
from su2ThermalRuntime import read,sha
from validateIdealAir import BIN,run,save,write_rows,read_rows

OUT=ROOT/'outputs/v167_validation/fixture01'
CONFIG=ROOT/'outputs/v164_validation/small01/case.cfg'


def main():
    OUT.mkdir(parents=True,exist_ok=False)
    assert read(ROOT/'outputs/v167_validation/path01/report.json')['passed']
    core=IdealAir();rows=[]
    def add(rho,t,k,dk,dw=1.,alpha=.75,kind=0,dim=3,velocity=(500.,-20.,75.)):
        u=np.array(velocity);u[2]=u[2] if dim==3 else 0
        e=float(core.state('PT',rho*287.0491267875601*t,t)['e'][0])
        state=[rho,*list(rho*u),rho*(e+np.sum(u*u)/2+k)]
        rows.append([len(rows),0,dim,*state,k,100.,rho*dk,rho*dw,1e-10,1e-4,1e10,1e15,alpha,kind])
    for rho in (.0061265027371194313,.02,1.):
        add(rho,300,100,20)
        add(rho,60.001,100,100,7,kind=1)
        add(rho,1999.999,100,-20,-3,kind=1)
        add(rho,300,100,-1e6,-1e9) # 原k/omega下限裁限不改变。
        add(rho,60.001,100,100,7,kind=1,dim=2)
    add(1,60,1024,1000,kind=2,alpha=1,velocity=(0,0,0))
    add(1,2000,1024,-100,kind=2,alpha=1,velocity=(0,0,0))
    events=list(csv.DictReader((ROOT/'outputs/v166_validation/real02/admissibility.csv').open()))
    a,b=events[335],events[337];flow=[float(a['old_'+n]) for n in ('rho','rhou','rhov','rhow','rhoE')]
    rows.append([len(rows),0,3,*flow,float(a['k']),100.,flow[0]*(float(b['k'])-float(a['k'])),flow[0]*9,1e-10,1e-4,1e10,1e15,1.,1])
    bad=list(rows[-1]);bad[0]=len(rows);bad[1]=1;bad[8]=float(b['k']);rows.append(bad)
    for index,value in ((3,0.),(3,-1.),(10,float('nan')),(10,float('inf')),(16,-1.),(2,4.)):
        row=list(rows[0]);row[0]=len(rows);row[1]=1;row[index]=value;rows.append(row)
    with (OUT/'inputs.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['case','reject','dim','rho','rhou','rhov','rhow','rhoE','oldk','oldw','drhok','drhow','lowk','loww','highk','highw','alpha','kind']);writer.writerows(rows)
    command=[str(BIN/'g++.exe'),'-std=c++17','-O3','-fopenmp','-ffast-math','-fno-finite-math-only',
        '-I'+str(SOURCE/'SU2_CFD/include'),'-I'+str(BUILD/'Common/src'),'-I'+str(SOURCE/'externals/eigen'),'-I'+str(SOURCE/'externals/mel'),
        str(ROOT/'tests/testSU2SSTEnergy.cpp'),str(SOURCE/'SU2_CFD/src/variables/CNSVariable.cpp'),
        '-Wl,--start-group',str(BUILD/'SU2_CFD/src/libSU2core.a'),str(BUILD/'Common/src/libSU2Common.a'),'-Wl,--end-group','-o',str(OUT/'fixture.exe')]
    save(OUT/'commands.json',dict(compile=command))
    before={str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}
    with overlay(OUT):
        run(OUT,'compile',command,30)
        run(OUT,'test',[str(OUT/'fixture.exe'),str(CONFIG),str(OUT/'inputs.csv'),str(OUT/'results.csv')],30)
    assert all(sha(BUILD/p)==digest for p,digest in before.items())
    data=read_rows(OUT/'results.csv');assert len(data)==len(rows)
    maximum=0.;oracle_rows=[];energies=[]
    for row,result in zip(rows,data):
        if row[1]:assert result['accepted']==0;continue
        assert result['accepted']==1
        alpha=float(result['alpha']);rho=row[3]
        values=[np.clip((row[8+j]*rho+alpha*row[10+j])/rho,row[12+j],row[14+j]) for j in range(2)]
        err=max(abs(values[j]-result[n])/max(1,abs(values[j])) for j,n in enumerate(('k','omega')))
        maximum=max(maximum,err);assert err<2e-13
        u=np.array(row[4:7])/rho;e=row[7]/rho-np.dot(u,u)/2-float(result['k'])
        oracle_rows.append(('rhoT',rho,float(result['T'])));energies.append(e)
    helper=ROOT/'outputs/v153_validation/native03/idealAirProbe.exe';assert sha(helper)==read(helper.parent/'runtime.json')['sha256']
    write_rows(OUT/'oracle_input.csv',oracle_rows)
    run(OUT,'oracle',[str(helper),'--oracle',str(OUT/'oracle_input.csv'),str(OUT/'oracle.csv')],30)
    oracle=read_rows(OUT/'oracle.csv');assert np.all(oracle['accepted']==1)
    error=float(max(abs(oracle['e']-energies)/np.maximum(1,abs(np.array(energies)))))
    assert error<2e-11
    sources=['scripts/verifySU2SSTEnergy.py','scripts/prepareSU2SSTEnergyPatch.py','tests/testSU2SSTEnergy.cpp','third_party/ideal_air/SSTEnergyUpdate.hpp','docs/validation/v167/PLAN.md']
    for name in sources:
        target=OUT/'project_sources'/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes((ROOT/name).read_bytes())
    save(OUT/'report.json',dict(passed=True,accepted=int(sum(data['accepted'])),rejected=int(sum(data['accepted']==0)),limited=int(sum(data['limited'])),
        zeroSteps=int(sum((data['accepted']==1)&(data['alpha']==0))),parallelCases=32,
        nativeUpdateRelativeError=maximum,independentEnergyError=error,engineObjectsUnchanged=len(before),sourceRestored=True,
        sourceSHA256={str(ROOT/n):sha(ROOT/n) for n in sources},physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(read(OUT/'report.json')))


if __name__=='__main__':main()
