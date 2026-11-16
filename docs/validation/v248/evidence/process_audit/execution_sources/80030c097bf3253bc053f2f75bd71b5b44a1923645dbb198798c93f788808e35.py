"""模块职责：严格区分完整循环FAS修正与终态空间残差，禁止旧残差字段身份混用。"""
from pathlib import Path
import struct
import numpy as np


def binary(path,magic,points,columns):
    p=Path(path)
    with p.open('rb') as f:
        head=f.read(24)
        if len(head)!=24 or struct.unpack('<8sQQ',head)!=(magic,points,columns):raise ValueError('FRAME_IDENTITY')
        if p.stat().st_size!=24+points*columns*8:raise ValueError('FRAME_LENGTH')
        a=np.fromfile(f,dtype='<f8').reshape(points,columns)
    if not np.array_equal(np.sort(a[:,0]),np.arange(points)):raise ValueError('FRAME_GLOBAL_IDS')
    return a


def state(path,points=335726,initial=False):
    a=binary(path,b'RAH248LE',points,34)
    if not np.isfinite(a[:,:27]).all() or np.any(a[:,33]!=0):raise ValueError('FRAME_STATE')
    if np.any(a[:,4]<=0) or np.any(a[:,5]<=0):raise ValueError('FRAME_DENSITY_VOLUME')
    if initial:
        if not np.isnan(a[:,27:33]).all():raise ValueError('INITIAL_ARRAY_UNKNOWN')
    elif not np.isfinite(a[:,27:33]).all():raise ValueError('FRAME_WORK_ARRAY')
    return a


def residual(path,points=335726):
    a=binary(path,b'RAH248RS',points,7)
    if not np.isfinite(a).all():raise ValueError('RESIDUAL_NONFINITE')
    physical=a[:,1:].copy();physical[:,4]-=300000*physical[:,0]
    return a,physical
