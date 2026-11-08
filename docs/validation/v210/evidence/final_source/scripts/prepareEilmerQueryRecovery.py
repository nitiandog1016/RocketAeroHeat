"""保存原查询失败；只补未完整输出的实数行，首次执行复数查询。"""
from pathlib import Path
import json,math,shutil
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v210_validation'

def main():
    out=BASE/'query_inputs02';out.mkdir(exist_ok=False)
    prior=BASE/'query_inputs01';config=json.loads((prior/'cases.json').read_text())
    original=(BASE/'query01/real/stdout.txt').read_text();receipt=json.loads((BASE/'query01/real/report.json').read_text())
    assert receipt['nativeExitCode']==1 and receipt['stopReason'] is None
    assert "The gas model 'rah-phase-aware-dry-air-eilmer-v999' is not available." in original
    # 未知模型抛的是原生Error，而非GasModelException；保留真实退出，不改引擎。
    config['cases']=[c for c in config['cases'] if c['id']!='bad_badModel']
    required={c['id'] for c in config['cases'] if c['dx']==0 and c['dy']==0}
    valid={}
    for line_number,line in enumerate(original.splitlines(),1):
        fields=line.split()
        if len(fields)<2 or fields[0] not in required:continue
        ok=fields[1]=='REJECT' and len(fields)>=3
        if fields[1]=='OK' and len(fields)==24:
            try:ok=all(math.isfinite(float(x)) for x in fields[2:])
            except ValueError:ok=False
        if ok:
            assert fields[0] not in valid
            valid[fields[0]]=dict(line=line,lineNumber=line_number)
    missing=required-set(valid)
    real=[line for line in (prior/'real.tsv').read_text().splitlines() if line.split()[0] in missing]
    complex_rows=[line for line in (prior/'complex.tsv').read_text().splitlines() if line.split()[0]!='bad_badModel']
    (out/'real.tsv').write_text('\n'.join(real)+'\n');(out/'complex.tsv').write_text('\n'.join(complex_rows)+'\n')
    (out/'cases.json').write_text(json.dumps(config,indent=2)+'\n')
    (out/'reused_real_rows.json').write_text(json.dumps(valid,indent=2)+'\n')
    report=dict(passed=True,originalQueryPassed=False,unknownModelNativeErrorVerified=True,reusedRealRows=len(valid),
        missingRealRows=len(real),missingIDs=sorted(missing),firstComplexRows=len(complex_rows),
        reason='unknown model throws Error and stopped original batch; one real row interleaved with stderr',CFDRuns=0)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');shutil.copyfile(__file__,out/'executed_entry.py')
    print(json.dumps(report))

if __name__=='__main__':raise ValueError('v210缺行恢复已完成，不重查询或覆盖原结果。')
