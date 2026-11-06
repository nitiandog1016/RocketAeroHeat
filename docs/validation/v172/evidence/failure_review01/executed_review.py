"""模块用途：有限只读异常符号映射与实际反馈摘要；不读大场、不启动CFD。"""
import argparse
import bisect
import csv
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
PACKAGE=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
BIN=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf8')


def main(folder,out):
    """准确PE基址、实际记录偏移及nm地址相互绑定；残缺调用栈不补造物理调用者。"""
    started=time.monotonic();folder=Path(folder).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    paths=[folder/name for name in ('engine_process.json','engine_stdout.log','history.csv')]
    runtime_path=PACKAGE/'runtime.json';exe=PACKAGE/'bin/SU2_CFD.exe';nm=BIN/'nm.exe';addr=BIN/'addr2line.exe'
    paths += [runtime_path,exe,nm,addr,Path(__file__).resolve()]
    locks={str(path):sha(path) for path in paths}
    receipt=json.loads(paths[0].read_text());runtime=json.loads(runtime_path.read_text())
    assert receipt['status']=='terminal' and receipt['exitCode']==134 and not receipt['timedOut']
    assert all(value is True for value in receipt['inputsStable'].values())
    assert locks[str(exe)]==receipt['inputSHA256'][str(exe)]==next(x['sha256'] for x in runtime['files'] if x['path']=='bin/SU2_CFD.exe')
    assert locks[str(runtime_path)]==receipt['inputSHA256'][str(runtime_path)]
    raw=exe.read_bytes();pe=struct.unpack_from('<I',raw,0x3c)[0];assert raw[pe:pe+4]==b'PE\0\0'
    optional=pe+24;assert struct.unpack_from('<H',raw,optional)[0]==0x20b
    image_base=struct.unpack_from('<Q',raw,optional+24)[0]
    stdout=paths[1].read_text(encoding='utf8',errors='replace');frames=[]
    for number,address,module,base,offset in re.findall(
            r'frame=(\d+) address=([0-9a-fA-F]+) module=(.*?) base=([0-9a-fA-F]+) offset=(0x[0-9a-fA-F]+)',stdout):
        assert int(address,16)-int(base,16)==int(offset,16)
        frames.append(dict(frame=int(number),actualAddress='0x'+address,module=module,
            loadedBase='0x'+base,offset=offset,isExactExecutable=Path(module).resolve()==exe,
            preferredAddress=hex(image_base+int(offset,16)) if Path(module).resolve()==exe else None))
    own=[frame for frame in frames if frame['isExactExecutable']];assert own
    commands=[]
    def tool(name,command):
        before=time.monotonic();result=subprocess.run(command,capture_output=True,timeout=8,creationflags=subprocess.CREATE_NO_WINDOW)
        commands.append(dict(name=name,command=list(map(str,command)),exitCode=result.returncode,
            elapsedSeconds=time.monotonic()-before,stdoutSHA256=hashlib.sha256(result.stdout).hexdigest(),
            stderr=result.stderr.decode('utf8',errors='replace')))
        assert result.returncode==0
        return result.stdout.decode('utf8',errors='replace')
    symbols=tool('nm',[str(nm),'--defined-only','--numeric-sort','--demangle',str(exe)])
    text_symbols=[]
    for address,kind,name in re.findall(r'^([0-9a-fA-F]+) ([Tt]) (.+)$',symbols,re.M):
        text_symbols.append((int(address,16),name))
    addresses=[row[0] for row in text_symbols];assert addresses==sorted(addresses)
    selected=[]
    for frame in own:
        address=int(frame['preferredAddress'],16);index=bisect.bisect_right(addresses,address)-1
        assert index>=0
        start,name=text_symbols[index];stop=text_symbols[index+1][0] if index+1<len(text_symbols) else None
        frame.update(nearestTextSymbol=name,symbolStart=hex(start),offsetWithinSymbol=address-start,
            nextTextSymbolAddress=hex(stop) if stop is not None else None,
            symbolMappingScope='nearest preceding text symbol; exact source line only when addr2line provides one')
        selected.extend(['%016x %s'%(a,n) for a,n in text_symbols[max(0,index-1):index+2]])
    (out/'selected_text_symbols.txt').write_text('\n'.join(selected)+'\n',encoding='utf8')
    resolved=tool('addr2line',[str(addr),'-f','-C','-e',str(exe)]+[frame['preferredAddress'] for frame in own])
    (out/'addr2line.txt').write_text(resolved,encoding='utf8');save(out/'symbol_commands.json',commands)
    with paths[2].open(encoding='utf8') as stream:rows=[{k.strip():float(v) for k,v in row.items()} for row in csv.DictReader(stream,skipinitialspace=True)]
    assert [int(row['Inner_Iter']) for row in rows]==list(range(46))
    selected_keys=('Inner_Iter','Min CFL','Max CFL','Avg CFL','LinSolRes','Linear_Solver_Iterations','Min DT','Max DT','rms[Rho]','rms[RhoE]','HF')
    records=[{key:row[key] for key in selected_keys} for row in rows]
    save(out/'actual_history_summary_rows.json',records)
    missed=[int(row['Inner_Iter']) for row in rows if row['LinSolRes']>-5]
    decrease=[int(row['Inner_Iter']) for row in rows if 10**row['LinSolRes']>.0012]
    reset=[int(row['Inner_Iter']) for row in rows if 10**row['LinSolRes']>.99]
    can_grow=[int(row['Inner_Iter']) for row in rows if 10**row['LinSolRes']<.001]
    stable={path:sha(Path(path))==value for path,value in locks.items()};assert all(stable.values())
    result=dict(reviewCompleted=True,exceptionType=re.findall(r'^exception_type=(.+)$',stdout,re.M),
        exceptionWhat=re.findall(r'^what=(.+)$',stdout,re.M),peImageBase=hex(image_base),actualRecordedFrames=frames,
        stackEndMarkerPresent='RAH_NATIVE_TERMINATION_END' in stdout,lastStdoutCharacters=stdout[-80:],
        physicalCallSiteIdentified=False,failingNodeKnown=False,failedTemperatureDirectionKnown=False,
        failedFlowStateAvailable=False,linearTargetMissedLocalIterations=missed,
        linearFeedbackReductionPredicateLocalIterations=decrease,linearFeedbackResetPredicateLocalIterations=reset,
        linearFeedbackCanIncreasePredicateLocalIterations=can_grow,
        maximumReportedCFL=max(row['Max CFL'] for row in rows),first=records[0],last=records[-1],
        nonlinearHundredUpdateWindowReached=False,underRelaxationInferred=False,
        actualTerminationCFLKnown=False,strictUpperCFL100HitInFiniteHistory=any(row['Max CFL']>=100 for row in rows),
        inputSHA256=locks,inputsStable=stable,symbolCommands=commands,elapsedSeconds=time.monotonic()-started,
        CFDLaunches=0,builds=0,largeFieldsRead=0,physicalQualification=False)
    save(out/'report.json',result)
    print(json.dumps(dict(elapsedSeconds=result['elapsedSeconds'],frames=len(frames),
        symbols=[f['nearestTextSymbol'] for f in own],maximumReportedCFL=result['maximumReportedCFL'],
        linearTargetMissedLocalIterations=missed,linearReductionPredicate=decrease,
        last=result['last']),ensure_ascii=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();main(args.folder,args.output)
