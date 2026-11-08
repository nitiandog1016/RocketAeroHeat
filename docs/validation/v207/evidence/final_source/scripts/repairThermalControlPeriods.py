"""仅重定严格一维制造控制的周期长度；不改物理壁距或重算BVP。"""
from pathlib import Path
import json,shutil,time
from prepareOpenFoamControl import ROOT,save,sha

def main():
    start=time.monotonic();base=ROOT/'outputs/v207_validation';out=base/'period_repair01';out.mkdir(exist_ok=False)
    old=base/'thermal01';new=base/'thermal02'
    before={str(p.relative_to(ROOT)):sha(p) for p in old.rglob('*') if p.is_file()}
    report=json.loads((base/'mesh01/report.json').read_text())
    if report['passed'] or not report['block']['passed']:raise ValueError('不是已记录的单一质量失败')
    text=(base/'mesh01/check/stdout.txt').read_text()
    if 'Failed 1 mesh checks.' not in text or 'number of cells: 32' not in text:raise ValueError('原质量失败不符')
    shutil.copytree(old,new,ignore=shutil.ignore_patterns('polyMesh'))
    path=new/'system/blockMeshDict';text=path.read_text()
    if text.count('.02')!=8:raise ValueError('周期顶点格式改变')
    path.write_text(text.replace('.02','.00003125'))
    for n in ['0/p','0/T','0/U','constant/physicalProperties','constant/momentumTransport','constant/thermophysicalTransport','system/fvSchemes','system/fvSolution','system/controlDict']:
        if sha(old/n)!=sha(new/n):raise ValueError('周期长度以外输入改变')
    oracle=json.loads((base/'prepare01/oracle.json').read_text());oracle['wallAreaM2']=(.001/32)**2
    oracle['periodChangeIsOneDimensionalInvariance']=True
    save(out/'oracle.json',oracle);shutil.copyfile(__file__,out/'executed_entry.py')
    if any(sha(ROOT/n)!=h for n,h in before.items()):raise ValueError('原失败数据改变')
    save(out/'report.json',dict(passed=True,originalMeshPassed=False,oldDataUnchanged=True,BVPReruns=0,CFDRuns=0,
        periodicLengthM=.001/32,wallAreaM2=oracle['wallAreaM2'],oldInputSHA256=before,
        newInputSHA256={str(p.relative_to(ROOT)):sha(p) for p in new.rglob('*') if p.is_file()},
        elapsedSeconds=time.monotonic()-start,physicalQualification=False))
    print('Period-only repair prepared; BVP and original failed mesh retained.')

if __name__=='__main__':main()
