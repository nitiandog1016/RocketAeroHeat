"""独立核对 Git 内的空间归档、残差门槛和收支单位，并检查积分器的解析尺寸。"""
import hashlib
import json
from pathlib import Path
import zipfile
import numpy as np
from energy_budget import budget

DOC=Path(__file__).resolve().parent
DATA=DOC/"research"


def check_integrals():
    # 两个绕轴的双锥面。质量变化、能量变化、常热流、均匀穿流均可独立解析积分。
    xy=np.array([[1,0],[0,1],[-1,0],[2,0],[0,2],[-2,0]],float)
    a=dict(PointID=np.arange(6),x=xy[:,0],y=xy[:,1],Density=np.ones(6),Energy=np.ones(6)*2.5,
           Pressure=np.ones(6),Temperature=np.ones(6),Velocity_x=np.ones(6)*.25,Velocity_y=np.zeros(6),
           Heat_Flux=np.ones(6)*3)
    seed=dict(a,Density=a["Density"]-.1,Energy=a["Energy"]-.125)
    c=dict(shape=[2,3],freePressurePa=10,densityKgM3=2,freeTemperatureK=300,speedMps=.25*np.sqrt(5))
    b=budget(seed,a,.5,c);expected_volume=14*np.pi/3
    np.testing.assert_allclose([b["volumeM3"],b["massRateKgS"],b["storedTotalEnergyRateW"],b["wallIntoSolidW"]],
        [expected_volume,.4*expected_volume,2.5*expected_volume,6*np.sqrt(2)*np.pi],rtol=2e-14)
    np.testing.assert_allclose([b["farOutMassKgS"],b["farOutConductionW"],b["farOutViscousWorkW"],
                                b["farOutTotalEnergyW"]],0,atol=1e-8)


def verify():
    check_integrals()
    # 用一个具备共同终点的已知数组测试比较器，和本工况缺少共同终点的状态分开。
    from read_study import compare_pair
    areas=np.array([1.,2.,4.])
    fine=np.array([[10.,0.,0.],[10.,-5.,1.],[10.,-10.,2.]])
    coarse=fine.copy();coarse[2]+=[1.,2.,-3.]
    pair=compare_pair((coarse,np.array([False,False,True])),(fine,np.zeros(3,bool)),
                      areas,{"all":slice(None),"zero":slice(0,1)},2.)
    np.testing.assert_allclose(pair["all"]["heatRmsFraction"],np.sqrt(16/450))
    assert pair["all"]["worstOriginalFaceIDs"]==[3,3,3] and pair["all"]["reverseStatusDifferentFaces"]==1
    assert pair["zero"]["heatRmsFraction"] is None and pair["zero"]["powerDifferenceFraction"] is None
    status=json.loads((DATA/"status.json").read_bytes())
    records=json.loads((DATA/"runs.json").read_bytes())
    assert len(records)==6 and len({v["seedSHA256"] for v in records.values()})==1
    assert all(v["seedSHA256"]==status["originSHA256"] for v in records.values())
    for key,r in records.items():
        dt=25e-6 if key.endswith("half1") else 50e-6
        assert r["dtSeconds"]==dt==r["nativeCurTimeSeconds"]
        worst=max(r["last20Max"].values())
        assert worst==status["stepReceipts"][key]["maxLast20"] and worst>=-8 and not r["innerPassed"]
        b=r["savedFieldBudget"]
        mass=b["massRateKgS"]+b["farOutMassKgS"]
        energy=b["storedTotalEnergyRateW"]+b["farOutTotalEnergyW"]+b["wallIntoSolidW"]
        np.testing.assert_allclose([mass,energy,energy-b["freeTotalEnthalpyJkg"]*mass],
            [b["massImbalanceKgS"],b["totalEnergyImbalanceW"],b["excessEnergyImbalanceW"]],atol=1e-6)
        assert not b["nativeFluxClosureProven"] and r["rejectedFrameRegions"]["all"]["faces"]==9216
    for kind in ("full","half1"):
        cols=np.loadtxt(DATA/f"{kind}_columns.csv",delimiter=",",skiprows=1)
        layers=np.loadtxt(DATA/f"{kind}_layers.csv",delimiter=",",skiprows=1)
        np.testing.assert_allclose(cols.sum(axis=0),layers.sum(axis=0),rtol=2e-14,atol=0)
        assert cols.shape==(109,4) and layers.shape==(89,4)
        r=records["newton_"+kind]["spatialResidual"]
        for i,key in enumerate(("Residual_Density","Residual_Momentum_x","Residual_Momentum_y","Residual_Energy")):
            v=r[key]
            np.testing.assert_allclose(.5*np.log10(cols[:,i].sum()/9701),v["historyLogRMS"],atol=.001,rtol=0)
            np.testing.assert_allclose(cols[12:14,i].sum()/cols[:,i].sum(),v["baseRimAndFirstTailColumnShare"])
            np.testing.assert_allclose(layers[:4,i].sum()/layers[:,i].sum(),v["firstFourLayerShare"])
        share_key="fullEnergyRimShare" if kind=="full" else "halfEnergyRimShare"
        np.testing.assert_allclose(cols[12:14,-1].sum()/cols[:,-1].sum(),status[share_key])
    with zipfile.ZipFile(DATA/"executed.zip") as z:
        hashes=json.loads((DATA/"executed.json").read_bytes())
        assert {name:hashlib.sha256(z.read(name)).hexdigest() for name in z.namelist()}==hashes
    assert status["acceptedPhysicalSteps"]==0 and status["sameTimeHeatRmsFraction"] is None
    assert not status["sameTimePairAvailable"] and status["wholeBodyPowerW"] is None
    print("6 native records; same-origin checks, 8 spatial RMS values, 6 balances and analytic integration checks passed.")


if __name__=="__main__":
    verify()
