"""原生定比热具名保存场温区守卫；异步只读检查，不改变CFD、不发用户取消。"""
import argparse
import csv
import ctypes
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import uuid


def domain_guard_contract():
    """固定开发域和有界观察合同；不保证未保存瞬时状态。"""
    return dict(schema='RAH.NACA.SavedFieldDomainGuard.v1',
        temperatureRangeK=[150.,360.],helperMaximumSeconds=15.,
        totalHelperMaximumSeconds=60.,historyStrictlyAfterForLive=True,
        requireWindowsNoWriterHandles=True,completeVolumeAndSurfacePair=True,
        observesUnsavedInstantaneousState=False,changesNumericalUpdates=False,
        stopReason='model-domain-exceeded',errorStopReason='model-domain-guard-error')


def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()


def write_json(path,value):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.pending')
    temporary.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    os.replace(str(temporary),str(path))


def signature(path):
    value=Path(path).stat()
    return [value.st_size,value.st_mtime_ns,value.st_ctime_ns,value.st_ino]


class WriterBusy(Exception):
    pass


def no_writer_reader(path):
    """Win32双向共享约束：已有写者拒绝，新写者在本句柄关闭前不能打开。"""
    if os.name!='nt':raise RuntimeError('该守卫只审核Windows无写者句柄，不能降级为普通open')
    import msvcrt
    from ctypes import wintypes
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.CreateFileW.argtypes=[wintypes.LPCWSTR,wintypes.DWORD,wintypes.DWORD,
        wintypes.LPVOID,wintypes.DWORD,wintypes.DWORD,wintypes.HANDLE]
    kernel.CreateFileW.restype=wintypes.HANDLE
    kernel.CloseHandle.argtypes=[wintypes.HANDLE];kernel.CloseHandle.restype=wintypes.BOOL
    handle=kernel.CreateFileW(str(Path(path).resolve()),0x80000000,1,None,3,0x80,None)
    if handle==wintypes.HANDLE(-1).value:
        error=ctypes.get_last_error()
        if error in (32,33):raise WriterBusy('具名文件仍有写者或排他锁')
        raise ctypes.WinError(error)
    try:fd=msvcrt.open_osfhandle(handle,os.O_RDONLY|os.O_BINARY)
    except BaseException:
        kernel.CloseHandle(handle);raise
    return os.fdopen(fd,'rb')


def history_prefix(path,maximum):
    """只接受以换行结束、从0连续且每列有限的完整历史前缀。"""
    if not Path(path).is_file():return -1
    raw=Path(path).read_bytes();end=raw.rfind(b'\n')
    if end<0:return -1
    rows=csv.reader(io.StringIO(raw[:end+1].decode('ascii')),skipinitialspace=True)
    keys=[x.strip() for x in next(rows,[])];index=keys.index('Inner_Iter')
    last=-1
    for row in rows:
        if not row:continue
        values=[float(x) for x in row]
        if len(values)!=len(keys) or not all(math.isfinite(x) for x in values):break
        step=values[index]
        if step!=last+1 or step>=maximum:break
        last=int(step)
    return last


COLUMNS=['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy',
    'Pressure','Temperature','Mach','Pressure_Coefficient','Velocity_x','Velocity_y','Velocity_z',
    'Laminar_Viscosity','Skin_Friction_Coefficient_x','Skin_Friction_Coefficient_y',
    'Skin_Friction_Coefficient_z','Heat_Flux','Y_Plus','Residual_Density',
    'Residual_Momentum_x','Residual_Momentum_y','Residual_Momentum_z','Residual_Energy']
WIDTHS={name:(3 if name in ('Momentum','Velocity','Skin_Friction_Coefficient','Residual_Momentum') else 1)
    for name in ('Density','Momentum','Energy','Pressure','Temperature','Mach','Pressure_Coefficient',
                 'Velocity','Laminar_Viscosity','Skin_Friction_Coefficient','Heat_Flux','Y_Plus',
                 'Residual_Density','Residual_Momentum','Residual_Energy')}


def mesh_geometry(path):
    """只读取合同锁定的原网格；不从任意表面子集推断壁身份。"""
    import numpy as np
    with Path(path).open() as stream:
        assert stream.readline().strip()=='NDIME= 3'
        for _ in range(int(stream.readline().split('=')[1])):stream.readline()
        count=int(stream.readline().split('=')[1])
        rows=np.array([list(map(float,stream.readline().split())) for _ in range(count)])
        assert rows.shape==(count,4) and np.array_equal(rows[:,3],np.arange(count))
        walls=[];names=[]
        for _ in range(int(stream.readline().split('=')[1])):
            name=stream.readline().split('=')[1].strip();assert name not in names;names.append(name)
            size=int(stream.readline().split('=')[1])
            for _ in range(size):
                line=list(map(int,stream.readline().split()))
                if name=='wall':walls.append(line)
        assert not stream.read().strip() and walls
    return rows[:,:3],walls


def parse_surface(raw):
    """解析已持无写者句柄的ASCII VTK，拒绝截断、重复字段和不合法单元。"""
    import numpy as np
    stream=io.StringIO(raw.decode('ascii'))
    assert stream.readline().startswith('# vtk');stream.readline()
    assert stream.readline().strip()=='ASCII'
    assert stream.readline().strip()=='DATASET UNSTRUCTURED_GRID'
    words=iter(stream.read().split());assert next(words)=='POINTS'
    count=int(next(words));assert count>0 and next(words) in ('double','float')
    xyz=np.fromiter((float(next(words)) for _ in range(count*3)),float).reshape(count,3)
    assert next(words)=='CELLS';nc=int(next(words));size=int(next(words));cells=[];used=0
    for _ in range(nc):
        width=int(next(words));assert width in (3,4)
        face=[int(next(words)) for _ in range(width)]
        assert len(set(face))==width and all(0<=i<count for i in face)
        cells.append(face);used+=width+1
    assert used==size and next(words)=='CELL_TYPES' and int(next(words))==nc
    kinds=[int(next(words)) for _ in range(nc)]
    assert all(k==(5 if len(f)==3 else 9) for k,f in zip(kinds,cells))
    assert next(words)=='POINT_DATA' and int(next(words))==count
    fields={}
    for tag in words:
        name=next(words);assert next(words) in ('double','float')
        assert name not in fields and name in WIDTHS
        if tag=='SCALARS':
            width=int(next(words));assert next(words)=='LOOKUP_TABLE' and next(words)=='default'
        else:assert tag=='VECTORS';width=3
        assert width==WIDTHS[name]
        fields[name]=np.fromiter((float(next(words)) for _ in range(count*width)),float).reshape(count,width)
    assert set(fields)==set(WIDTHS) and np.isfinite(xyz).all()
    assert all(np.isfinite(value).all() for value in fields.values())
    return xyz,cells,fields


def inspect_pair(request):
    """helper的唯一工作：完整几何/状态配对通过后，比较实际已保存Temperature。"""
    import numpy as np
    from scipy.spatial import cKDTree
    folder=Path(request['folder']);step=request['step'];contract=request['contract']
    assert contract['savedFieldDomainGuard']==domain_guard_contract()
    assert contract['selectedFluidModel']=='IDEAL_GAS' and not contract['turbulent']
    assert contract['energyReference']=='e=cv*T; h=cp*T; no-OTH-offset'
    assert contract['developmentTemperatureRangeK']==[150.,360.]
    maximum=contract['checkpointContract']['maximumUpdates']
    last=history_prefix(folder/'history.csv',maximum)
    assert last>=step and (request['nativeTerminalAtSubmission'] or last>step), '历史尚未越过运行中最新帧'
    assert digest(folder/'mesh.su2')==contract['meshSHA256']
    points,walls=mesh_geometry(folder/'mesh.su2')
    assert len(points)==request['nodes'] and np.isfinite(points).all()
    paths=[folder/('restart_flow_%06d.csv'%step),folder/('surface_flow_%06d.vtk'%step)]
    with no_writer_reader(paths[0]) as volume, no_writer_reader(paths[1]) as surface:
        before=[signature(p) for p in paths]
        if before!=request['signatures']:return dict(status='pending',reason='检查前身份已变')
        vr=volume.read();sr=surface.read()
        assert vr.endswith(b'\n') and sr.endswith(b'\n'), '未完整结束的保存文件'
        keys=next(csv.reader(io.StringIO(vr.split(b'\n',1)[0].decode('ascii'))))
        assert keys==COLUMNS, '保存场列合同不匹配'
        a=np.loadtxt(io.BytesIO(vr),delimiter=',',skiprows=1,ndmin=2)
        assert a.shape==(len(points),len(COLUMNS)), '保存场行数不匹配'
        assert np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(len(points)))
        assert np.allclose(a[:,1:4],points,rtol=2e-14,atol=1e-14), '体节点身份不匹配'
        assert min(a[:,4])>0 and min(a[:,9])>0 and min(a[:,10])>0
        xyz,faces,fields=parse_surface(sr)
        wallids=np.unique([i for face in walls for i in face[1:]])
        assert len(xyz)==len(wallids) and len(faces)==len(walls)
        distance,mapped=cKDTree(points[wallids]).query(xyz)
        assert max(distance)<1e-6 and len(np.unique(mapped))==len(wallids), '表面节点不是完整原壁'
        mapped=wallids[mapped]
        expected=sorted(tuple(sorted(face[1:])) for face in walls)
        actual=sorted(tuple(sorted(mapped[face].tolist())) for face in faces)
        assert actual==expected, '表面单元不是完整原壁'
        # VTK固定6位指数，使用既有2e-6配对门槛，绝不以表面圆整改变体场温区。
        surface_error=float(max(abs(fields['Temperature'][:,0]/a[mapped,10]-1)))
        assert surface_error<2e-6, '体/面温度不属于同一保存态'
        assert [signature(p) for p in paths]==before, '持锁读取期间文件身份变化'
        low,high=contract['savedFieldDomainGuard']['temperatureRangeK']
        outside=np.flatnonzero((a[:,10]<low)|(a[:,10]>high))
        result=dict(status='complete',step=step,temperatureRangeK=[float(a[:,10].min()),float(a[:,10].max())],
            outOfRangeCount=len(outside),firstOutOfRangePointID=int(outside[0]) if len(outside) else None,
            exceeded=bool(len(outside)),volumeSHA256=hashlib.sha256(vr).hexdigest(),
            surfaceSHA256=hashlib.sha256(sr).hexdigest(),surfaceTemperatureRelativeError=surface_error,
            volumeNodes=len(points),wallNodes=len(wallids),wallFaces=len(walls),signatures=before,
            noWriterHandlesAcquired=True,historyMaximumAtInspection=last)
    return result


def worker(request_path,result_path):
    request=json.loads(Path(request_path).read_text(encoding='utf8'));start=time.monotonic()
    try:result=inspect_pair(request)
    except WriterBusy as error:result=dict(status='pending',reason=str(error))
    except Exception as error:result=dict(status='error',reason=type(error).__name__+': '+str(error))
    result.update(requestSHA256=digest(request_path),token=request['token'],elapsedSeconds=time.monotonic()-start)
    write_json(result_path,result)
    return 0


class DomainGuard:
    """主监督器只做小文件轮询；所有网格/大场读取在一个自有异步helper中。"""
    def __init__(self,folder,manifest):
        self.folder=Path(folder).resolve();self.contract=manifest['idealGas'];self.nodes=manifest['nodes']
        assert self.contract['savedFieldDomainGuard']==domain_guard_contract()
        assert manifest['caseKind']=='naca-tn3344-constant-cp-spatial-wall-v1'
        assert self.contract['selectedFluidModel']=='IDEAL_GAS' and not self.contract['turbulent']
        assert self.contract['developmentTemperatureRangeK']==[150.,360.]
        self.maximum=self.contract['checkpointContract']['maximumUpdates']
        self.every=self.contract['checkpointContract']['every']
        self.active=None;self.checked={};self.attempts=[];self.elapsed=0.;self.serial=0
        self.exceeded=False;self.error=None;self.terminal_detection=False;self.closed=False
        self.last_history=-1;self.retry_after=0.;self.pending=[];self.failure_frame=None;self.cleanup_errors=[]

    def snapshot(self):
        elapsed=self.elapsed+(time.monotonic()-self.active['start'] if self.active else 0)
        return dict(exceeded=self.exceeded,error=self.error,
            stopReason='model-domain-exceeded' if self.exceeded else ('model-domain-guard-error' if self.error else None),
            terminalOnlyDetection=self.terminal_detection,active=bool(self.active),
            pending=list(self.pending),helperElapsedSeconds=elapsed,helperAttempts=list(self.attempts),
            checkedFrames=[self.checked[k] for k in sorted(self.checked)],firstExceededFrame=self.failure_frame,
            lastContinuousFiniteHistoryIteration=self.last_history,contract=domain_guard_contract(),
            cleanupErrors=list(self.cleanup_errors),
            observesUnsavedInstantaneousState=False)

    def _stop_helper(self):
        if self.active:
            active=self.active;process=active['process']
            try:
                if process.poll() is None:process.kill()
                process.wait(timeout=5)
            except Exception as error:
                self.cleanup_errors.append(type(error).__name__+': '+str(error))
                self.error='自有守卫helper回收失败：'+self.cleanup_errors[-1]
            finally:
                try:active['log'].close()
                except Exception as error:
                    self.cleanup_errors.append(type(error).__name__+': '+str(error))
                    self.error='守卫日志关闭失败：'+self.cleanup_errors[-1]
                self.elapsed+=time.monotonic()-active['start'];self.active=None

    def poll(self,native_terminal=False):
        if self.closed or self.exceeded or self.error:return self.snapshot()
        try:self._poll(bool(native_terminal))
        except Exception as error:
            self.error=type(error).__name__+': '+str(error);self._stop_helper()
        return self.snapshot()

    def _poll(self,terminal):
        if self.active:
            active=self.active;duration=time.monotonic()-active['start']
            if duration>=15 or self.elapsed+duration>=60:
                self.attempts.append(dict(step=active['step'],status='timeout',elapsedSeconds=duration))
                self.error='保存场helper超过单次15秒或累计60秒';self._stop_helper();return
            if active['process'].poll() is None:return
            exitcode=active['process'].returncode;self._stop_helper()
            assert exitcode==0 and active['result'].is_file(), '守卫helper异常退出或结果缺失'
            result=json.loads(active['result'].read_text(encoding='utf8'))
            assert result['requestSHA256']==active['requestSHA256'] and result['token']==active['token']
            self.attempts.append(dict(step=active['step'],status=result['status'],elapsedSeconds=duration,
                                     request=active['request'].name,result=active['result'].name))
            if result['status']=='error':self.error=result['reason'];return
            if result['status']=='complete':
                self.checked[active['step']]=result
                if result['exceeded']:
                    self.exceeded=True;self.terminal_detection=terminal;self.failure_frame=result;return
            else:
                assert result['status']=='pending';self.retry_after=time.monotonic()+.25;return
        self.last_history=history_prefix(self.folder/'history.csv',self.maximum)
        eligible=[step for step in range(0,self.last_history+1,self.every)
                  if terminal or step<self.last_history]
        if terminal and self.last_history>=0:
            # 异常/取消后最后历史可能从未保存；正常缺末帧另由严格inventory拒绝。
            latest=[self.folder/('restart_flow_%06d.csv'%self.last_history),
                    self.folder/('surface_flow_%06d.vtk'%self.last_history)]
            if any(p.is_file() for p in latest):eligible=sorted(set(eligible+[self.last_history]))
        self.pending=[step for step in eligible if step not in self.checked]
        # 已接收的旧帧若改变，必须拒绝，不能把上次温度结果嫁接新内容。
        for step,result in self.checked.items():
            paths=[self.folder/('restart_flow_%06d.csv'%step),self.folder/('surface_flow_%06d.vtk'%step)]
            assert [signature(p) for p in paths]==result['signatures'], '已检查旧帧身份发生变化'
        if not self.pending or time.monotonic()<self.retry_after:return
        if self.elapsed>=60:self.error='保存场helper累计60秒预算已耗尽';return
        step=self.pending[0];paths=[self.folder/('restart_flow_%06d.csv'%step),self.folder/('surface_flow_%06d.vtk'%step)]
        assert all(p.is_file() for p in paths), '已不可重写的计划帧缺少CSV/VTK配对：%s'%step
        self.serial+=1;stem='domain_guard_%03d_%06d'%(self.serial,step)
        request_path=self.folder/(stem+'_request.json');result_path=self.folder/(stem+'_result.json')
        assert not request_path.exists() and not result_path.exists(), '守卫回执名称已存在'
        token=uuid.uuid4().hex
        request=dict(folder=str(self.folder),step=step,nodes=self.nodes,contract=self.contract,
                     nativeTerminalAtSubmission=terminal,signatures=[signature(p) for p in paths],token=token)
        write_json(request_path,request);log=(self.folder/(stem+'.log')).open('wb');start=time.monotonic()
        try:
            process=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--worker',str(request_path),
                '--result',str(result_path)],stdout=log,stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        except BaseException:log.close();raise
        self.active=dict(process=process,start=start,step=step,log=log,request=request_path,result=result_path,
                         token=token,requestSHA256=digest(request_path))

    def finish(self,deadline=None):
        """原生已结束后的有界排空；deadline为父监督器monotonic绝对截止。"""
        while True:
            if deadline is not None and time.monotonic()>=deadline:
                self.error='父监督器截止前守卫未完成';self._stop_helper();break
            state=self.poll(native_terminal=True)
            if state['exceeded'] or state['error'] or (not state['active'] and not state['pending']):break
            time.sleep(.025)
        return self.snapshot()

    def close(self):
        self._stop_helper();self.closed=True


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker',required=True);parser.add_argument('--result',required=True)
    args=parser.parse_args();sys.exit(worker(args.worker,args.result))
