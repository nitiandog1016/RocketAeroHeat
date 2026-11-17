"""模块职责：绑定原46入口/有限球帽质量供体和原始公式，生成MATLAB联合积分输入。"""
from pathlib import Path
import numpy as np,hashlib,json
from phaseAwareEntropyInlet import ShockMassModel
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v254_validation'


def prepare():
    out=BASE/'prepare01';out.mkdir(exist_ok=False)
    e=read(ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json');c=e['condition']
    table=read(ROOT/'outputs/v253_validation/inlet01/shock_mass_source.json')
    old=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'))
    model=ShockMassModel(c['pressurePa'],c['temperatureK'],c['speedMps'])
    entropy=np.array([model.shock(u)['entropyJkgK'] for u in table['sinShockAngle']])
    probe=np.array([1e-3,2.5e-4,6.25e-5]);values=np.array([(entropy[0]-model.shock(np.sqrt(1-w))['entropyJkgK'])/w for w in probe])
    first=(4*values[1]-values[0])/3;D=(4*values[2]-values[1])/3;indicator=abs(D-first)/abs(D);assert indicator<1e-4
    dS_dPsi=-4*table['A']*D/(model.rho*model.speed*table['radiusM']**2)
    paths=[ROOT/'outputs/v253_validation/inlet01/shock_mass_source.json',ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz',
           ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    provenance={str(p.relative_to(ROOT)):sha(p) for p in paths}
    sourceText=json.dumps(provenance,sort_keys=True,separators=(',',':'))
    source=dict(pressurePa=c['pressurePa'],temperatureK=c['temperatureK'],speedMps=c['speedMps'],
        totalEnthalpyJkg=model.H,stagnationPressurePa=model.p0,stagnationTemperatureK=model.t0,
        stagnationEntropyJkgK=float(entropy[0]),radiusM=table['radiusM'],
        sourceMassCoordinate=table['massCoordinate'],sourceEntropyJkgK=entropy.tolist(),entropyDerivativeByMassAtZero=float(dS_dPsi),
        targetAnglesRad=np.r_[old['thetaRad'],np.pi/3].tolist(),sourceRole='engineering-candidate',
        sourceSHA256=hashlib.sha256(sourceText.encode()).hexdigest(),sourceText=sourceText,
        targetOriginalUpstreamFaceIDs=old['upstreamFaceIDs'].tolist(),targetOriginalDownstreamFaceIDs=old['downstreamFaceIDs'].tolist(),
        targetPhysicalWallXYZ=old['wallXYZ'].tolist(),targetActualMetricM=old['metricM'].tolist(),targetLabelJacobian=old['labelJacobian'].tolist(),
        sourcePressureModel='existing modified-Newtonian rotated fitted sphere; not downstream original-body pressure',
        actualProfileProvided=False,physicalQualification=False)
    save(out/'source.json',source)
    # 只绑定真正采用的原页及元数据；搜索失败与中断另留研究目录，不冒称取到D625。
    research=BASE/'research';selected=['NASA_CR_111921.pdf','NASA_CR_111922.pdf','cr111921-51.png','cr111921-52.png',
        'cr111921-53.png','cr111921-54.png','cr112180-32.png','cr112180-33.png','program_upright-1.png','program_upright-4.png']
    save(out/'method.json',dict(decision='adopt original wall-property integral mass closure plus implicit entropy feedback',
        formula='I=integral((p/p0)*(U/Uinf)*h^2 ds); theta=.4696*sqrt(2*rho_ws*mu_w*Uinf*I)/(rho_e*U*h); delta_m=5.27*theta',
        notThwaitesReferenceThickness=True,sourceProgramPages=[53,57],sourceProgramPDFPages=[58,62],sourceHeatEquations=[146,149,150,151,152,153,154,155,156,157,158,159],
        physicalMomentumDerivativeNotSolved=True,pressureGradientProxyExplicit=True,noActualVelocityTemperatureProfile=True,
        priorSearchNotFound='TN D-625 / TR R-118 full originals not obtained; exact needed CTHM/THM implementation verified in CR-111922',
        externalSourceLicense='NASA public-use metadata; no third-party runtime source redistributed',
        bindings={name:sha(research/name) for name in selected},derivativeRelativeIndicator=float(indicator),
        entropyDerivativeByCosSquaredGamma=float(D),entropyDerivativeByMassAtZero=float(dS_dPsi)))
    save(out/'report.json',dict(passed=True,targetRingPositions=46,extraTerminalAngleDeg=60,
        sourceSHA256=source['sourceSHA256'],sourceInputs=provenance,shockODERecomputed=False,oldBVPRecomputed=False,
        sourceEntropyNodes=len(entropy),derivativeRelativeIndicator=float(indicator),newQualifiedAreaM2=0))


if __name__=='__main__':prepare()
