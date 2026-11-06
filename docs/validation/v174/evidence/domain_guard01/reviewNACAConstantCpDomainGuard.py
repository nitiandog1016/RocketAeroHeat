"""有预算的小型保存场守卫审核；实际只读回放与合成拒绝分列。"""
import argparse
import copy
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy.spatial import cKDTree
import nacaConstantCpDomainGuard as guard
from nacaConstantCpModel import CV,R

ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/'docs/validation/v174/DOMAIN_GUARD_PLAN.md'
BUILD=ROOT/'outputs/v174_validation/build01'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
REAL=ROOT/'outputs/v174_validation/small02'
BASE=ROOT/'outputs/v174_validation/small03'


def source_review():
    tree_path=BUILD/'source_tree.json';tree=json.loads(tree_path.read_text(encoding='utf8'))
    entries=[
        ('SU2_CFD/src/output/COutput.cpp',True,[434,637,807,814,880,907]),
        ('SU2_CFD/src/output/filewriter/CSU2FileWriter.cpp',True,[40,72]),
        ('SU2_CFD/src/output/filewriter/CParaviewFileWriter.cpp',False,[59,168]),
        ('SU2_CFD/src/output/filewriter/CParallelFileWriter.cpp',False,[264,305]),
        ('SU2_CFD/include/output/filewriter/CFileWriter.hpp',False,[192]),
        ('SU2_CFD/src/iteration/CFluidIteration.cpp',False,[437,447,453]),
        ('SU2_CFD/src/iteration/CIteration.cpp',False,[210,213,214]),
        ('SU2_CFD/src/drivers/CSinglezoneDriver.cpp',False,[94,204,205,206,207]),
        ('Common/src/CConfig.cpp',True,[9042,9043,9046])]
    rows=[];locks={str(tree_path):guard.digest(tree_path)}
    for rel,executed,lines in entries:
        path=(BUILD/'executed_sources' if executed else SOURCE)/rel
        expected=tree[str(SOURCE/rel)];actual=guard.digest(path);assert actual==expected
        text=path.read_text(encoding='utf8');split=text.splitlines()
        rows.append(dict(path=str(path),sourceTreePath=str(SOURCE/rel),SHA256=actual,
            evidence=[dict(line=line,text=split[line-1]) for line in lines]))
        locks[str(path)]=actual
    assert 'append ? "ab" : "wb"' in (SOURCE/entries[3][0]).read_text()
    assert 'append = false' in (SOURCE/entries[4][0]).read_text()
    assert 'GetMultizone_Problem() ? curOuterIter : curInnerIter' in (BUILD/'executed_sources'/entries[-1][0]).read_text()
    return dict(passed=True,files=rows,sourceTreeSHA256=locks[str(tree_path)],
        synchronousSingleProcessWriters=True,atomicRenamePublication=False,
        mostRecentNamedPairMayBeOverwrittenAtTerminal=True,
        runtimeEligibility='continuous finite history maximum > named local iteration',
        terminalEligibility='native process confirmed terminated; continuous finite history includes frame',
        unsupportedContexts=['MPI','multizone','unsteady','restart','changed writer source']),locks


def request(folder,manifest,step,terminal=True):
    return dict(folder=str(folder),nodes=manifest['nodes'],contract=manifest['idealGas'],step=step,
        nativeTerminalAtSubmission=terminal,token='independent-small-check',
        signatures=[guard.signature(folder/('restart_flow_%06d.csv'%step)),
                    guard.signature(folder/('surface_flow_%06d.vtk'%step))])


def vtk_bytes(xyz,faces,fields):
    lines=['# vtk DataFile Version 3.0','synthetic guard fixture','ASCII','DATASET UNSTRUCTURED_GRID',
           'POINTS %s double'%len(xyz),' '.join('%.17g'%v for v in xyz.ravel()),
           'CELLS %s %s'%(len(faces),sum(len(f)+1 for f in faces)),
           ' '.join(str(v) for f in faces for v in [len(f)]+f),'CELL_TYPES %s'%len(faces),
           ' '.join('5' if len(f)==3 else '9' for f in faces),'POINT_DATA %s'%len(xyz)]
    for name,value in fields.items():
        width=value.shape[1]
        lines+=['SCALARS %s double 1'%name,'LOOKUP_TABLE default'] if width==1 else ['VECTORS %s double'%name]
        lines+=[' '.join('%.17g'%v for v in value.ravel())]
    return ('\n'.join(lines)+'\n').encode('ascii')


def review(folder):
    started=time.monotonic();folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    assert not folder.exists();folder.mkdir(parents=True)
    for path in [Path(__file__),Path(guard.__file__),PLAN]:shutil.copy2(str(path),str(folder/path.name))
    report=dict(passed=False,newCFDRuns=0,newBuilds=0,physicalQualified=False,tests=[])
    sources=[Path(__file__),Path(guard.__file__),PLAN,BASE/'case_manifest.json',BASE/'mesh.su2',BASE/'history.csv',
        BASE/'restart_flow_000002.csv',BASE/'surface_flow_000002.vtk',REAL/'case_manifest.json',REAL/'mesh.su2',REAL/'history.csv']
    sources += [REAL/('%s_%06d.%s'%(prefix,step,extension)) for step in (2,3)
                for prefix,extension in [('restart_flow','csv'),('surface_flow','vtk')]]
    locks={str(path):guard.digest(path) for path in sources}
    handles=[]
    try:
        source,source_locks=source_review();locks.update(source_locks);report['writerSourceReview']=source
        guard.write_json(folder/'writer_source_review.json',source)
        def check(name,condition,**values):
            assert condition,name
            report['tests'].append(dict(name=name,passed=True,**values))
            assert time.monotonic()-started<60,'独立测试累计60秒预算耗尽'
        def rejects(name,fn):
            try:fn()
            except (AssertionError,ValueError,StopIteration,RuntimeError,OSError) as error:
                check(name,True,reason=type(error).__name__+': '+str(error));return
            raise AssertionError(name+' 未拒绝')
        manifest=json.loads((BASE/'case_manifest.json').read_text(encoding='utf8'))
        manifest['idealGas']['savedFieldDomainGuard']=guard.domain_guard_contract()
        real_manifest=json.loads((REAL/'case_manifest.json').read_text(encoding='utf8'))
        real_manifest['idealGas']['savedFieldDomainGuard']=guard.domain_guard_contract()
        actual=[]
        for step,wanted in [(2,[176.73863274849217,361.938931018377]),(3,[176.7056454871528,367.6974140032273])]:
            result=guard.inspect_pair(request(REAL,real_manifest,step))
            check('actual-small02-read-only-frame-%s'%step,result['status']=='complete' and result['exceeded']
                  and np.max(abs(np.array(result['temperatureRangeK'])-wanted))<1e-11)
            result['evidenceKind']='read-only replay of actual saved frame; not original online stop';actual.append(result)
        report['actualSavedFramesReadOnlyReplay']=actual
        report['originalSmall02MissingFirstFrameFailureUnchanged']=True
        fixture=folder/'synthetic';fixture.mkdir()
        for name in ('mesh.su2','history.csv'):shutil.copy2(str(BASE/name),str(fixture/name))
        a=np.loadtxt(str(BASE/'restart_flow_000002.csv'),delimiter=',',skiprows=1)
        xyz,faces,fields=guard.parse_surface((BASE/'surface_flow_000002.vtk').read_bytes())
        mapped=cKDTree(a[:,1:4]).query(xyz)[1]
        def save_pair(t,step=0):
            data=a.copy();data[:,10]=t;data[:,9]=data[:,4]*R*t
            data[:,8]=data[:,4]*(CV*t+.5*np.sum((data[:,5:8]/data[:,4,None])**2,axis=1))
            stream=io.BytesIO();np.savetxt(stream,data,delimiter=',',fmt='%.17g',comments='',
                header=','.join('"'+v+'"' for v in guard.COLUMNS))
            sf={k:v.copy() for k,v in fields.items()};sf['Temperature'][:,0]=np.asarray(t)[mapped]
            (fixture/('restart_flow_%06d.csv'%step)).write_bytes(stream.getvalue())
            (fixture/('surface_flow_%06d.vtk'%step)).write_bytes(vtk_bytes(xyz,faces,sf))
        t=np.full(len(a),300.);t[0]=150.;t[1]=360.;save_pair(t)
        good=guard.inspect_pair(request(fixture,manifest,0,False))
        check('synthetic-inclusive-endpoints',good['temperatureRangeK']==[150.,360.] and not good['exceeded'])
        for name,index,value in [('below',0,np.nextafter(150.,0)),('above',1,np.nextafter(360.,1000))]:
            bad=t.copy();bad[index]=value;save_pair(bad)
            result=guard.inspect_pair(request(fixture,manifest,0,False))
            check('synthetic-'+name,result['exceeded'] and result['outOfRangeCount']==1)
        save_pair(t);vp=fixture/'restart_flow_000000.csv';sp=fixture/'surface_flow_000000.vtk'
        vr=vp.read_bytes();sr=sp.read_bytes()
        with vp.open('ab'):
            try:guard.inspect_pair(request(fixture,manifest,0,False))
            except guard.WriterBusy:check('existing-writer-blocks-read-handle',True)
            else:raise AssertionError('已有写者未被排斥')
        with guard.no_writer_reader(vp):
            rejects('read-handle-blocks-new-writer',lambda:vp.open('ab'))
        mutations=[('missing-final-newline',vr[:-1]),('wrong-row-count',b'\n'.join(vr.splitlines()[:-1])+b'\n'),
            ('wrong-schema',vr.replace(b'PointID',b'Unknown',1)),
            ('nan',vr.replace(vr.splitlines()[1],vr.splitlines()[1].replace(b'150',b'nan',1),1)),
            ('duplicate-point-id',vr.replace(vr.splitlines()[2],b'0'+vr.splitlines()[2][1:],1))]
        for name,value in mutations:
            vp.write_bytes(value);rejects('synthetic-'+name,lambda:guard.inspect_pair(request(fixture,manifest,0,False)))
        vp.write_bytes(vr);sp.write_bytes(sr[:-40]);rejects('synthetic-truncated-vtk',lambda:guard.inspect_pair(request(fixture,manifest,0,False)))
        sp.write_bytes(sr);sf={k:v.copy() for k,v in fields.items()};sf['Temperature'][:]=999.
        sp.write_bytes(vtk_bytes(xyz,faces,sf));rejects('synthetic-mismatched-pair',lambda:guard.inspect_pair(request(fixture,manifest,0,False)))
        sp.write_bytes(sr)
        history=(fixture/'history.csv').read_bytes();lines=history.splitlines(keepends=True)
        (fixture/'history.csv').write_bytes(b''.join(lines[:2]))
        rejects('live-latest-frame-refused',lambda:guard.inspect_pair(request(fixture,manifest,0,False)))
        late=guard.DomainGuard(fixture,manifest);handles.append(late)
        state=late.poll();check('latest-frame-not-submitted-live',not state['active'] and not state['error'])
        state=late.finish(deadline=time.monotonic()+10)
        check('latest-frame-terminal-only-check',not state['error'] and len(state['checkedFrames'])==1)
        n=len(state['helperAttempts']);late.poll(native_terminal=True)
        check('unchanged-complete-frame-not-reparsed',len(late.snapshot()['helperAttempts'])==n)
        (fixture/'history.csv').write_bytes(b''.join(lines[:3]))
        # 已检查0后最新历史1无保存场，不应凭空要求异常瞬间存在。
        state=late.finish(deadline=time.monotonic()+10)
        check('abnormal-latest-unsaved-stays-unknown',not state['error'] and len(state['checkedFrames'])==1)
        late.close()
        # 新目录防止复用既有helper回执；原生进程并未启动，此处明确为合成调度。
        online=folder/'synthetic_online';online.mkdir()
        for path in [fixture/'mesh.su2',vp,sp]:shutil.copy2(str(path),str(online/path.name))
        (online/'history.csv').write_bytes(history)
        bad=t.copy();bad[1]=361.;save_pair(bad)
        for path in (vp,sp):shutil.copy2(str(path),str(online/path.name))
        live=guard.DomainGuard(online,manifest);handles.append(live)
        before=time.monotonic();state=live.poll();duration=time.monotonic()-before
        check('asynchronous-live-submission',state['active'] and duration<1.,pollSeconds=duration)
        while not state['exceeded'] and not state['error']:
            assert time.monotonic()-before<15;time.sleep(.025);state=live.poll()
        check('live-exceeded-is-not-terminal-detection',state['exceeded'] and not state['error']
              and not state['terminalOnlyDetection'] and state['stopReason']=='model-domain-exceeded')
        check('guard-never-writes-user-cancel',not (online/'engine_cancel.json').exists());live.close()
        terminal=folder/'synthetic_terminal';terminal.mkdir()
        for path in (online/'mesh.su2',online/'restart_flow_000000.csv',online/'surface_flow_000000.vtk'):
            shutil.copy2(str(path),str(terminal/path.name))
        (terminal/'history.csv').write_bytes(b''.join(lines[:2]))
        final=guard.DomainGuard(terminal,manifest);handles.append(final)
        state=final.finish(deadline=time.monotonic()+10)
        check('terminal-only-exceeded-identified',state['exceeded'] and state['terminalOnlyDetection']);final.close()
        for name,elapsed,age in [('single-helper-timeout',0.,16.),('cumulative-helper-timeout',59.,2.)]:
            timed=guard.DomainGuard(online,manifest);handles.append(timed)
            log=(folder/(name+'.log')).open('wb')
            process=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'],stdout=log,stderr=log,
                                     creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            timed.elapsed=elapsed;timed.active=dict(process=process,start=time.monotonic()-age,step=0,log=log)
            state=timed.poll()
            check(name,bool(state['error']) and not state['exceeded'] and process.poll() is not None);timed.close()
        deadline=guard.DomainGuard(online,manifest);handles.append(deadline)
        state=deadline.finish(deadline=time.monotonic()-1)
        check('parent-deadline-refusal',bool(state['error']) and not state['exceeded']);deadline.close()
        # 缺行与部分尾行都不能使latest编号被视为安全的旧帧。
        history_test=folder/'history_gap.csv';history_test.write_bytes(b''.join([lines[0],lines[1],lines[3]]))
        check('history-gap-stops-prefix',guard.history_prefix(history_test,4)==0)
        history_test.write_bytes(b''.join(lines[:2])+lines[2][:-10])
        check('history-partial-tail-ignored',guard.history_prefix(history_test,4)==0)
        report['passed']=True
    except Exception as error:
        report['error']=type(error).__name__+': '+str(error);raise
    finally:
        for instance in handles:instance.close()
        report['elapsedSeconds']=time.monotonic()-started
        report['inputSHA256']=locks;report['inputsStable']={name:guard.digest(name)==value for name,value in locks.items()}
        report['passed']=report['passed'] and all(report['inputsStable'].values())
        guard.write_json(folder/'report.json',report)
        print(json.dumps(dict(passed=report['passed'],tests=len(report['tests']),elapsedSeconds=report['elapsedSeconds'],
                             error=report.get('error')),ensure_ascii=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True)
    review(parser.parse_args().folder)
