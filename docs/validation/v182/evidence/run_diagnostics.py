"""一次四组合原生算子检查：两数字类型读取相同初态/末态，不发展CFD。"""
from pathlib import Path
import json,shutil,sys,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import run,sha
base=ROOT/'outputs/v182_validation';build=base/'build01';source=ROOT/'outputs/v181_validation/small01'
assert json.loads((build/'recovery04_report.json').read_text())['passed']
start=time.monotonic();results=[]
for number,stem in [('real','lmr-run'),('complex','lmrZ-run')]:
    for snapshot in [0,2]:
        name=number+'-state'+str(snapshot);case=base/name;shutil.copytree(source,case)
        # 这些是复制的历史场，只把指定CSV认作本轮诊断输出；没有原生时间推进。
        exe=build/(stem+'-diagnostic');target=case/'operator.csv'
        inputs=[p for p in source.rglob('*') if p.is_file()]
        r=run([str(exe),str(snapshot),str(target)],case,base/(name+'-process'),min(15,60-(time.monotonic()-start)),128,
              inputs+[exe,ROOT/'scripts/native/eilmerOperatorDiagnostic.d',Path(__file__)])
        text=(base/(name+'-process')/'stdout.txt').read_text(errors='replace')
        good=r['passed'] and text.count('RAH_NATIVE_OPERATOR_V182_ENTRY')==1 and text.count('RAH_OPERATOR_COMPLETED')==1 and ('RAH_NUMBER_TYPE '+number) in text and target.is_file()
        results.append(dict(name=name,passed=good,process=r,outputSHA256=sha(target) if target.exists() else None))
        if not good:break
    if not results[-1]['passed']:break
report=dict(passed=len(results)==4 and all(x['passed'] for x in results),results=results,
            elapsedSeconds=time.monotonic()-start,newCFDSteps=0,physicalQualification=False)
(base/'diagnostics_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['passed','elapsedSeconds','newCFDSteps']}))
