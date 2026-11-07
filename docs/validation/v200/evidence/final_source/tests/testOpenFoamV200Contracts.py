"""重启来源、显式/隐式SU2编号与已完成任务关闭的定向拒绝测试。"""
from pathlib import Path
import copy,json,shutil,sys,tempfile,time,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runOpenFoamRestart import validate_source
from prepareOpenFoamControl import read_tetra


def main():
    start=time.monotonic();out=ROOT/'outputs/v200_validation/contract_tests01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_test.py')
    source=json.loads((ROOT/'outputs/v200_validation/restart_review02/report.json').read_text());assert validate_source(source)
    rows=[]
    def reject(name,action):
        try:action()
        except (AssertionError,ValueError) as error:rows.append(dict(name=name,rejected=True,error=str(error)));return
        raise AssertionError('未拒绝 '+name)
    edits=[('failed_review',lambda r:r.update(passed=False)),('bitwise_identity',lambda r:r.update(continuousBitwiseEquivalent=True)),
        ('wrong_source',lambda r:r.update(sourceFrame='outputs/other/0.002')),('wrong_index',lambda r:r['time'].update(index=0)),
        ('wrong_time',lambda r:r['time'].update(value=.003)),('missing_energy',lambda r:r['sourceFieldsSHA256'].pop('e')),
        ('changed_field',lambda r:r['sourceFieldsSHA256'].update(p='0'*64)),('changed_mesh',lambda r:r.update(sourceGeometrySHA256='0'*64)),
        ('changed_clock',lambda r:r.update(sourceTimeSHA256='0'*64))]
    for name,edit in edits:
        r=copy.deepcopy(source);edit(r);reject(name,lambda r=r:validate_source(r))
    text='NDIME= 3\nNELEM= 1\n10 0 1 2 3{ID}\nNPOIN= 4\n0 0 0 0\n1 0 0 1\n0 1 0 2\n0 0 1 3\nNMARK= 1\nMARKER_TAG= wall\nMARKER_ELEMS= 4\n5 1 2 3\n5 0 3 2\n5 0 1 3\n5 0 2 1\n'
    with tempfile.TemporaryDirectory(dir=out) as directory:
        path=Path(directory)/'mesh.su2';path.write_text(text.replace('{ID}',' 0'))
        a=read_tetra(path);reject('explicit_as_implicit',lambda:read_tetra(path,True))
        path.write_text(text.replace('{ID}',''));b=read_tetra(path,True)
        assert np.array_equal(a[0],b[0]) and np.array_equal(a[1],b[1]) and a[2]==b[2]
        reject('implicit_without_named_mode',lambda:read_tetra(path))
        path.write_text(text.replace('{ID}',' 1'));reject('wrong_explicit_element_id',lambda:read_tetra(path))
        path.write_text(text.replace('{ID}','').replace('0 0 1 3','0 0 1 2'));reject('wrong_point_id',lambda:read_tetra(path,True))
    p=subprocess.run([sys.executable,str(ROOT/'scripts/runOpenFoamRestart.py')],capture_output=True,timeout=10)
    assert p.returncode!=0 and b'ValueError' in p.stderr and b'ModuleNotFoundError' not in p.stderr
    rows.append(dict(name='completed_restart_cli_closed',rejected=True))
    result=dict(passed=True,rejectionCount=len(rows),rejections=rows,positiveChecks=2,CFDRuns=0,engineBuilds=0,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['passed','rejectionCount','positiveChecks','elapsedSeconds']}))


if __name__=='__main__':main()
