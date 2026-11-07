"""只读恢复原生错误流格式核验；原失败/输出/数值不变。"""
from pathlib import Path
import json,re,shutil,time
from prepareOpenFoamControl import ROOT,sha,save


def main():
    start=time.monotonic();base=ROOT/'outputs/v203_validation';out=base/'scalar_audit02';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');old=json.loads((base/'scalar01/report.json').read_text())
    assert old['passed'] is False;process=old['process'];assert process['nativeExitCode']==1 and not process['liveGroupAfterCleanup'] and all(process['inputStable'].values())
    rows=[]
    for p in Path('/opt/openfoam13/src/OpenFOAM/db/IOstreams').rglob('*'):
        if p.suffix not in ['.C','.H'] or 'lnInclude' in p.parts:continue
        lines=p.read_text().splitlines();found=[dict(line=i,text=s) for i,s in enumerate(lines,1) if ('precision' in s.lower() and ('6' in s or 'default' in s.lower()))]
        if found:
            q=out/'source'/p.name;shutil.copyfile(p,q) if q.parent.exists() else (q.parent.mkdir(),shutil.copyfile(p,q))
            rows.append(dict(path=str(p),sha256=sha(p),matches=found))
    # 读取实际原生六有效数字字符串；std::cout的17位设置不影响Foam错误流。
    text=(base/'scalar01/native/stdout.txt').read_text();printed=re.search(r'([\deE.+-]+) out of range 100 to 2500',text).group(1)
    expected=format(old['predictedTrialTemperatureK'],'.6g');assert printed==expected
    assert 'ROW 1197.0279007383137 100.01000000000001 ' in text and 'of NonUniformTable Cp' in text
    assert any('6' in m['text'] for row in rows for m in row['matches'])
    result=dict(passed=True,originalAuditPassed=False,postRunReadOnlyRecovery=True,CFDRuns=0,nativeQueries=0,
        originalReportSHA256=sha(base/'scalar01/report.json'),nativeOutputSHA256=sha(base/'scalar01/native/stdout.txt'),
        printedTemperature=printed,expectedAtActualPrecision=expected,knownTrueTemperatureK=100.01,seedK=180,
        falseRejectionOfInRangeEnergy=True,exactRocketRootCauseEstablished=False,sourceFiles=rows,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','printedTemperature','falseRejectionOfInRangeEnergy','elapsedSeconds']}))


if __name__=='__main__':main()
