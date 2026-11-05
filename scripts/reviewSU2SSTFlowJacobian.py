"""独立解析核对已存矩阵的缺项与规模；不重跑原生片段或修改求解器。"""
from pathlib import Path
import json, numpy as np
from verifySU2AUSMUp2Flux import ROOT, sha, write


def main():
    folder=ROOT/'outputs/v149_validation/jacobian01'; target=folder/'independent_review.json'
    assert not target.exists()
    cases={r['name']:r for r in json.loads((folder/'cases.json').read_text())}
    req=json.loads((folder/'requests.json').read_text())
    data=np.loadtxt(str(folder/'native_result.csv'),delimiter=',')
    gamma=1.4; uniform=[]; actual=[]
    for r in req:
        c=cases[r['case']];k=c['k'];normal=np.array(c['normal']);U=np.array(c['left']);rho=U[0];u=U[1:4]/rho
        p=(gamma-1)*(U[4]-.5*U[1:4]@U[1:4]/rho-rho*k[0]);H=(U[4]+p)/rho;un=u@normal
        dp=np.r_[(gamma-1)*(.5*u@u-k[0]),-(gamma-1)*u,gamma-1]
        # 独立解析微分物理通量，既不调用复步oracle，也不使用原函数矩阵公式。
        A=np.zeros((5,5));A[0,1:4]=normal
        for i in range(3):
            A[i+1]=u[i]*A[0]+normal[i]*dp
            A[i+1,0]-=u[i]*un;A[i+1,i+1]+=un
        A[4]=un*dp;A[4,0]-=un*H;A[4,1:4]+=H*normal;A[4,4]+=un
        fs=np.array(r['fluxScale']);us=np.array(r['uScale']);J=data[r['base'],6:].reshape(2,5,5)
        if c['kind']=='uniform' and r['mode'] in [0,1]:
            defect=np.zeros((5,5));defect[1:4,0]=(gamma-1)*k[0]*normal;defect[4,1:4]=(gamma-1)*k[0]*normal
            error=float(np.max(abs((J.sum(axis=0)-A-defect)*us[0][None,:]/fs[:,None])))
            assert error<2e-12
            uniform.append(dict(mode=r['mode'],case=c['name'],frozenK=k[0],predictedDefectError=error))
        if c['name'].startswith('actual_') and r['mode']==0:
            sizes=[]
            for side in range(2):
                correction=np.zeros((5,5));correction[1:4,0]=.5*(gamma-1)*k[side]*normal
                correction[4,1:4]=.5*(gamma-1)*k[side]*normal
                sizes.append(float(np.max(abs(correction*us[side][None,:]/fs[:,None]))))
            actual.append(dict(case=c['name'],pointIDs=c['pointIDs'],frozenK=k,centralMissingTermScaledSize=sizes,
                               actualGlobalRHSOrCausalAttribution=False))
    result=dict(passed=True,uniform=uniform,actual=actual,scope='algebraic central-matrix consistency only',
                decision='Do not enable untested accurate option or rebuild full solver from this local result alone.',
                newEngineExecutions=0,physicsQualified=False,newQualifiedArea=0,
                sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),folder/'cases.json',folder/'requests.json',folder/'native_result.csv']})
    write(target,result);print(json.dumps(dict(passed=True,maximumPredictionError=max(r['predictedDefectError'] for r in uniform),actual=actual)))


if __name__=='__main__':main()
