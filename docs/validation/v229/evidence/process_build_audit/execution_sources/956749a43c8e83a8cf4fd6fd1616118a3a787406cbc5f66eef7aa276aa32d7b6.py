"""模块职责：保留构建输入守卫失败，准确审核唯一调度入口变动与真实编译依赖。"""
import ast,json
from pathlib import Path
from runPhaseEulerCompatibilityV229 import ROOT,BASE,OLD_BUILD,read,save,sha

def main():
    out=BASE/'build_audit01';out.mkdir(exist_ok=False)
    process=BASE/'process_recover_build';report=read(process/'report.json')
    assert report['exitCode']==0 and report['stopReason'] is None and not report['liveGroupAfterCleanup']
    changed=[n for n,h in report['inputSHA256'].items() if sha(n)!=h]
    target=ROOT/'scripts/runPhaseEulerCompatibilityV229.py';assert changed==[str(target)],changed
    before=process/'execution_sources'/(report['inputSHA256'][str(target)]+'.py')
    old=ast.parse(before.read_text(encoding='utf8'));new=ast.parse(target.read_text(encoding='utf8'))
    # 构建只导入常量/size/locks等；改动只在未调用的另一个主入口main。
    def exclude_main(tree):
        tree.body=[n for n in tree.body if not isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) or n.name!='main']
        return ast.dump(tree)
    assert exclude_main(old)==exclude_main(new)
    assert report['inputsStable'] is False and report['passed'] is False
    actual={}
    for name in ('CEulerSolver.obj.d','compatibility.d'):
        text=(BASE/'build02'/name).read_text(encoding='utf8').replace('\\\n',' ')
        names=text.split(': ',1)[1].split()
        for value in names:
            p=Path(value);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
            assert report['inputSHA256'][str(p)]==h,str(p);actual[str(p)]=h
    assert read(BASE/'process_build/report.json')['passed']
    members=read(BASE/'build02/members.json')
    for rows in members.values():
        for row in rows:assert sha(row['path'])==row['sha256']
    save(out/'actual_dependencies.json',actual)
    save(out/'report.json',dict(passed=True,originalProcessPassed=False,originalInputsStable=False,
      changedFile=str(target),changeScope='unused main dispatcher only; imported constants/functions AST identical',
      actualBuildInputsStable=True,actualDependencyFiles=len(actual),noBuildRepeated=True,
      scientificQualification=False,nativeControlAdmission=True,originalFailurePreserved=True))
    print(json.dumps(dict(passed=True,actualDependencyFiles=len(actual),originalProcessPassed=False)))

if __name__=='__main__':main()
