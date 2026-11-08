"""模块职责：审核实际原生导出，准备语义损坏用例，并独立核对原生导入逐位结果。"""
from pathlib import Path
import json,sys,copy,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import nativeFrozenSnapshot as api
def save(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf8')
def mutate(d,path,value):
    for key in path[:-1]:d=d[key]
    d[path[-1]]=value
def main(case):
    start=time.monotonic();case=Path(case);meta=api.load(case/'binding.json');doc=api.load(case/'native-snapshot.json')
    donors=api.load(case/'donors.json')['donors'];mother={tuple(d['motherFaceNodes']):d['motherCellID'] for d in donors}
    assert len(mother)==len(donors)==len(doc['records'])==288
    mapping={}
    for r in doc['records']:
        key=tuple(r['motherFaceNodes']);assert mother[key]==r['motherCellID']
        mapping[r['receiverFaceID']]=dict(motherFaceNodes=list(key),motherCellID=mother[key],geometry=r['values'][:10])
    assert len(mapping)==288
    originalSHA=api.sha(case/'native-snapshot.json');a=api.read_snapshot(case/'native-snapshot.json',originalSHA,meta,mapping)
    assert np.all(a[:,4,0]==0) and np.isnan(a).any() and np.any(a[:,21,0]!=a[:,10,0])
    report=dict(schema='RAH.NativeFileBridgeAudit.v1',nativeFileSHA256=originalSHA,records=len(a),scalarsPerRecord=140,
        areaxyZeroPreserved=True,unusedNaNPreserved=True,nativeRhoSNotRepaired=True,CFDAdvanceRuns=0,
        newQualifiedAreaM2=0,fullNumericalAcceptance=False)
    if (case/'import-results.jsonl').exists():
        requests=api.load(case/'requests.json');results=[json.loads(s) for s in (case/'import-results.jsonl').read_text().splitlines()]
        assert len(results)==len(requests['positive'])+len(requests['negative'])+4 and all(r['passed'] for r in results)
        for p in requests['positive']:
            wanted=api.load(case/p['path']);got=api.load(case/p['output'])
            assert got['metadata']==wanted['metadata']
            assert {r['receiverFaceID']:r for r in got['records']}=={r['receiverFaceID']:r for r in wanted['records']}
            api.read_snapshot(case/p['output'],api.sha(case/p['output']),meta,mapping)
        report.update(passed=True,bitwiseNativeRoundtrip=True,nativePositiveCount=len(requests['positive']),
            nativeRejectedCount=len(requests['negative'])+4,nativeResidentDonorReferences=0)
        save(case/'final-audit.json',report);print(json.dumps(report));return
    reverse=copy.deepcopy(doc);reverse['records'].reverse();save(case/'reordered-snapshot.json',reverse)
    assert np.array_equal(api.validate(reverse,meta,mapping).view(np.uint64),a.view(np.uint64))
    requests=dict(positive=[dict(name='original',path='native-snapshot.json',sha256=originalSHA,output='roundtrip.json'),
        dict(name='reordered',path='reordered-snapshot.json',sha256=api.sha(case/'reordered-snapshot.json'),output='reordered-roundtrip.json')],negative=[])
    def add(name,path,value):requests['negative'].append(dict(name=name,path=path,replacement=value))
    for key in meta:
        v=meta[key];v=v+1 if type(v) in (int,float) else v+'-wrong'
        add('identity-'+key,['metadata',key],v)
    add('duplicate-face',['records',1,'receiverFaceID'],doc['records'][0]['receiverFaceID'])
    add('wrong-face',['records',0,'receiverFaceID'],10000000)
    add('wrong-mother-cell',['records',0,'motherCellID'],10000000)
    add('wrong-mother-face',['records',0,'motherFaceNodes',0],10000000)
    add('short-state',['records',0,'values'],doc['records'][0]['values'][:-1])
    add('short-records',['records'],doc['records'][:-1])
    add('malformed-hex',['records',0,'values',0],'x'*32)
    for name,index,value in [('zero-volume',3,0),('wrong-position',0,1),('negative-length',5,-1),
        ('low-temperature',14,59),('bad-density',10,1),('bad-sound-speed',13,1),('bad-energy',15,1),
        ('bad-viscosity',16,.001),('bad-conductivity',17,1),('negative-SA',25,-1),('negative-eddy-mu',27,-1),
        ('negative-eddy-k',28,-1),('bad-mass-fraction',20,.5),('active-nan',22,np.nan),
        ('active-inf',23,np.inf),('unused-inf',31,np.inf),('active-gradient-nan',35,np.nan),
        ('inactive-active-rhou',35+3*6+3,-1),('reversed-extrema',35+4,1e10),
        ('active-viscous-nan',119,np.nan),('mass-gradient',131,1)]:
        add(name,['records',0,'values',index],api.scalar(value))
    add('active-complex',['records',0,'values',22],api.scalar(1,1e-30))
    add('unused-complex',['records',0,'values',31],api.scalar(0,1e-30))
    requests['negative'] += [dict(name='wrong-sha',wrongSHA=True),dict(name='truncated-json',truncate=True)]
    rejected=[]
    for q in requests['negative']:
        d=copy.deepcopy(doc)
        if 'path' in q:mutate(d,q['path'],q['replacement'])
        try:
            if 'wrongSHA' in q:api.read_snapshot(case/'native-snapshot.json','0'*64,meta,mapping)
            elif 'truncate' in q:json.loads(json.dumps(d)[:100])
            else:api.validate(d,meta,mapping)
        except (ValueError,TypeError,KeyError):rejected.append(q['name'])
        else:raise AssertionError('not rejected: '+q['name'])
    save(case/'requests.json',requests)
    report.update(passed=True,pythonRejectedCount=len(rejected),rejected=rejected,seconds=time.monotonic()-start)
    save(case/'export-audit.json',report);print(json.dumps({k:report[k] for k in ('passed','records','pythonRejectedCount','seconds')}))
if __name__=='__main__':main(sys.argv[1])
