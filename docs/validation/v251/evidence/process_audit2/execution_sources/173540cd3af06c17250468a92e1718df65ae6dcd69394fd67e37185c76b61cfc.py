"""模块职责：严格读取具名单连通状态帧；只解码，不把结构正确等同于物理资格。"""
from pathlib import Path
import struct
import numpy as np

COLUMNS='globalID x y z volume rho rhou rhov rhow rhoE_internal nuTilde T u v w p rho_primitive H_internal a mu mut_flow k Cp mut_SA wallDistance flowRelaxation flowCFL flowRHS0 flowRHS1 flowRHS2 flowRHS3 flowRHS4_internal saRHS reserved0'.split()


def read_frame(path,expected_points,initial=False):
    p=Path(path)
    with p.open('rb') as f:
        head=f.read(24)
        if len(head)!=24:raise ValueError('FRAME_HEADER_TRUNCATED')
        magic,n,k=struct.unpack('<8sQQ',head)
        if magic!=b'RAH235LE' or n!=expected_points or k!=len(COLUMNS):raise ValueError('FRAME_IDENTITY')
        if p.stat().st_size!=24+n*k*8:raise ValueError('FRAME_LENGTH')
        a=np.fromfile(f,dtype='<f8',count=n*k).reshape((n,k))
    if not np.isfinite(a[:,:27]).all() or not np.isfinite(a[:,33]).all():raise ValueError('FRAME_STATE_NONFINITE')
    if initial:
        if not np.isnan(a[:,27:33]).all():raise ValueError('INITIAL_RHS_MUST_BE_UNKNOWN')
    elif not np.isfinite(a[:,27:33]).all():raise ValueError('FRAME_RHS_NONFINITE')
    if not np.array_equal(np.sort(a[:,0]),np.arange(n)):raise ValueError('FRAME_GLOBAL_IDS')
    if np.any(a[:,4]<=0) or np.any(a[:,5]<=0) or np.any(a[:,33]!=0):raise ValueError('FRAME_PHYSICAL_LAYOUT')
    return dict(columns=COLUMNS,values=a,physicalQualification=False,initialRHSUnknown=initial,
      physicalEnergy=a[:,9]-300000*a[:,5],physicalEnergyRHS=None if initial else a[:,31]-300000*a[:,27])
