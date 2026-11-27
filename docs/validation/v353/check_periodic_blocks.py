"""在一个小离散系统上直接组装 360° 矩阵，校核全部角向模态和共享轴点。"""
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0,str(Path(__file__).resolve().parents[3]/"scripts"))
from simpleLauncherAngularHeat import SECTORS,steady_fields


def check():
    # 第 0 个轴点自由，第 5 个轴点给温。第 2 个点的单个扇区升温激励全部角向模态。
    axis=np.array([0,5]);wet=np.array([2,5]);inlet=np.array([4])
    block=np.zeros((6,6))
    for left,right,forward,back in [(0,1,2,2.3),(1,2,1.1,1),(1,3,3,2),
            (3,4,1.5,2),(3,5,.7,.9),(2,3,.6,.6),(1,5,.4,.5)]:
        block[left,left]+=forward;block[left,right]-=forward
        block[right,right]+=back;block[right,left]-=back
    gamma=np.array([0.,.8,.4,.55,.3,0.])
    identifiers=np.empty((SECTORS,6),int)
    identifiers[:,axis]=[0,1]
    identifiers[:,[1,2,3,4]]=np.arange(2,4*SECTORS+2).reshape(SECTORS,4)
    size=identifiers.max()+1
    matrix=np.zeros((size,size))
    for j,idx in enumerate(identifiers):
        matrix[np.ix_(idx,idx)]+=block
        right=identifiers[(j+1)%SECTORS]
        for a,b,g in zip(idx,right,gamma):
            matrix[a,a]+=g;matrix[b,b]+=g;matrix[a,b]-=g;matrix[b,a]-=g
    wall=np.ones((1,SECTORS,2));wall[0,7,0]=1.1
    fixed=np.unique(np.r_[identifiers[:,wet].ravel(),identifiers[:,inlet].ravel()])
    free=np.setdiff1d(np.arange(size),fixed)
    expected=np.zeros(size)
    expected[identifiers[:,wet]]=wall[0]
    expected[free]=np.linalg.solve(matrix[np.ix_(free,free)],-matrix[np.ix_(free,fixed)]@expected[fixed])
    model=dict(axis=axis,wet=wet,inlet=inlet,xy=np.zeros((6,2)),matrix=csr_matrix(block),conductance=gamma)
    computed=steady_fields(model,wall)[0]
    modes=int(np.count_nonzero(abs(np.fft.rfft(wall[0,:,0]))>1e-12))
    error=float(np.abs(computed-expected[identifiers]).max())
    assert modes==25 and error<1e-11 and np.ptp(computed[:,0])==0
    return dict(physicalDOF=int(size),excitedFourierModes=modes,maxTemperatureDifference=error)


if __name__=="__main__":
    print(json.dumps(check()))
