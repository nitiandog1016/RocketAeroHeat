from pathlib import Path
import json,subprocess,sys,time
root=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
name=sys.argv[1];assert name in ('analysis01','analysis02')
command=[sys.executable,str(root/'scripts/analyzeHybridRegionInputs.py'),'--out',str(out/name)]
start=time.monotonic()
with (out/(name+'.stdout.txt')).open('wb') as stream:
    result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=60)
record=dict(argv=command,exitCode=result.returncode,elapsedSeconds=time.monotonic()-start,parentLimitSeconds=60)
(out/(name+'_parent.json')).write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record));raise SystemExit(result.returncode)
