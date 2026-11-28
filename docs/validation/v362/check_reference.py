"""独立读取发布的球壁表，对比不可压解析解与 NACA TN 762 的 M0.30 球面解。"""
import argparse
import json
from pathlib import Path

import numpy as np

from read_sphere import common_nodes, metrics


HERE = Path(__file__).resolve().parent


def review(plot=None):
    wall = np.genfromtxt(HERE/"sphere_wall.csv", delimiter=",",names=True)
    theta = np.deg2rad(wall["angleDeg"])
    nodes = np.linspace(0,np.pi,129)
    assert len(wall)==128 and np.all((theta>nodes[:-1]) & (theta<nodes[1:]))
    weight = (np.r_[0,wall["ringWeight"]]+np.r_[wall["ringWeight"],0])/2
    v,cp,cp_poisson = [common_nodes(theta,wall[name],power) for name,power in
                         (("utOverU",1),("cpBernoulli",2),("cpNativePoisson",2))]
    computed = metrics(nodes,weight,v,cp,cp_poisson)
    raw = json.loads((HERE/"sphere.json").read_bytes())["commonNodes"]
    for method in ("bernoulli","nativePoisson"):
        for key in ("rmsCp","maxAbsCp"):
            np.testing.assert_allclose(computed[method][key],raw[method][key],atol=2e-14,rtol=0)
    np.testing.assert_allclose(computed["tangentialRMS"],raw["tangentialRMS"],atol=2e-14,rtol=0)

    # 作者公式 (17) 给定空气 γ=1.408；M² 和 M⁴ 两项分别保留。
    mach,gamma = .30,1.408
    sin1,sin3,sin5 = (np.sin(k*nodes) for k in (1,3,5))
    v0 = 1.5*sin1
    v2 = v0 + mach**2*(989*sin1-1215*sin3)/7040
    v4 = v2 + mach**4*(.11048*sin1-.17018*sin3+.07489*sin5)
    np.testing.assert_allclose([v2[64],v4[64]],[1.52818,1.53106],rtol=0,atol=5e-6)
    def pressure(speed):
        return 2/(gamma*mach**2)*np.expm1(gamma/(gamma-1)*np.log1p((gamma-1)*mach**2/2*(1-speed**2)))
    cp4,cp2,cp0 = pressure(v4),pressure(v2),1-v0**2
    def errors(delta):
        return {"rms":float(np.sqrt(np.average(delta**2,weights=weight))),
                "maxAbs":float(np.max(abs(delta)))}
    if plot is not None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,1,figsize=(7.2,6.2),sharex=True,constrained_layout=True)
        degrees=np.degrees(nodes)
        axes[0].plot(degrees,cp_poisson,color="#99775d",label="Native filtered Poisson pressure")
        axes[0].plot(degrees,cp,color="#226e87",linewidth=2,label="Bernoulli from computed wall speed")
        axes[0].plot(degrees,cp0,"--",color="#222222",label="Analytical incompressible sphere")
        axes[0].set(ylabel="Pressure coefficient",title="Slip sphere in a uniform stream")
        axes[1].plot(degrees,cp-cp0,color="#226e87",label="Relative to incompressible theory")
        axes[1].plot(degrees,cp-cp4,color="#ae523f",label="Relative to M0.30 theory (NACA TN 762)")
        axes[1].set(xlabel="Angle from the front stagnation point (degrees)",ylabel="Bernoulli Cp minus reference")
        for ax in axes:
            ax.grid(axis="y",alpha=.2)
            ax.legend(frameon=False,fontsize=8)
            ax.set_xticks([0,30,60,90,120,150,180])
        fig.savefig(plot,dpi=150)
        plt.close(fig)
    return {"source":"NACA TN 762, Equations (17), (23), Table I; printed pages 8, 9, 12",
            "url":"https://ntrs.nasa.gov/citations/19930081579",
            "pdfSHA256":"215f531320c88fe0cb20e5c242df12ec6e6f5555d871ccf8e9ffc75ed30ab0e0",
            "mach":mach,"gamma":gamma,"tableEquatorialSpeed":[1.52818,1.53106],
            "calculatedEquatorialSpeed":[float(v2[64]),float(v4[64])],
            "equatorialCp":{"numericalIncompressible":float(cp[64]),"analyticalIncompressible":-1.25,
                             "compressibleTheory":float(cp4[64])},
            "incompressibleTheoryVsMach030":{"speed":errors(v0-v4),"cp":errors(cp0-cp4)},
            "computedIncompressibleVsMach030":{"speed":errors(v-v4),"cp":errors(cp-cp4)},
            "secondVsFourthPowerTerms":{"speed":errors(v2-v4),"cp":errors(cp2-cp4)}}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure",type=Path)
    print(json.dumps(review(parser.parse_args().figure),indent=2,allow_nan=False))
