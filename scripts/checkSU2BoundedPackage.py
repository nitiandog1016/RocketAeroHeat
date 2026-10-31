"""v131原生小型兼容/有界接口；固定同一均匀初态，两步各最多30秒。"""
from pathlib import Path
import json,subprocess,sys
import numpy as np
from prepareSU2SSTBoundedStudy import ROOT,BASE,prepare,sha,write
from auditSU2SSTState import audit,vtk_tokens
from su2SSTInterface import state
def main():
    target=BASE/'small_comparison.json';assert not target.exists()
    for name in ['official','compat','bounded']:
        prepare(name);folder=BASE/('small_'+name)
        subprocess.run([sys.executable,str(ROOT/'scripts/runEmbeddedSU2.py'),'--request',str(folder/'engine_request.json')],check=True,timeout=90)
        write(folder/'sst_state_audit.json',audit(folder))
    errors={};arrays={}
    for name in ['official','compat','bounded']:
        folder=BASE/('small_'+name);a,cols,nd,_,_,_=state(folder/'restart_flow.csv',True,True)
        arrays[name]=a
        r=json.loads((folder/'engine_process.json').read_text())
        assert r['exitCode']==0 and not r['timedOut'] and all(r['inputsStable'].values())
    assert arrays['official'].shape==arrays['compat'].shape
    errors['allRestartColumnsRelativeMax']=float(np.max(abs(arrays['official']-arrays['compat'])/np.maximum(np.maximum(abs(arrays['official']),abs(arrays['compat'])),1)))
    for key in ['Pressure','Heat_Flux','Temperature','Y_Plus']:
        _,a=vtk_tokens(BASE/'small_official/surface_flow.vtk');_,b=vtk_tokens(BASE/'small_compat/surface_flow.vtk')
        errors[key]=float(np.max(abs(a[key]-b[key])/np.maximum(np.maximum(abs(a[key]),abs(b[key])),1)))
    passed=max(errors.values())<=2e-6
    report=dict(passed=passed,compatibilityErrors=errors,compatibilityTolerance=2e-6,
        sameUniformInitialState=True,updatesPerRun=2,boundedInterfacePassed=True,physicsQualified=False,
        newQualifiedArea=0,sourceSHA256={str(p):sha(p) for name in ['official','compat','bounded']
            for p in (BASE/('small_'+name)).iterdir() if p.is_file()},
        recipeSHA256=sha(Path(__file__).resolve()))
    write(target,report);print(json.dumps(dict(passed=passed,errors=errors)))
    assert passed,'未改路径兼容失败；停止三维对照'
if __name__=='__main__':main()
