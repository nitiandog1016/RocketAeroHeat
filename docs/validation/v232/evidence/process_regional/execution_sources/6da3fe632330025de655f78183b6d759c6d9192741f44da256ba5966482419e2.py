"""新具名输运组合：成熟线性mu/kappa表，复用v206量热及相态合同。"""
import numpy as np
from phaseAwareDryAir import PhaseAwareDryAir,MODEL,DOMAIN

TRANSPORT='rah-lemmon-dilute-icoTabulated-of13-v1'
class PhaseAwareViscousAir(PhaseAwareDryAir):
    def __init__(self):
        super().__init__()
        # 输运与压力无关，此合法低压仅满足物性接口范围，不授连续流资格。
        t=super().transport(np.ones_like(self.T),self.T)
        self.mu=np.asarray(t['dynamicViscosityPas']);self.kappa=np.asarray(t['conductivityWmK'])
    def interpolated_transport(self,p,t):
        state=self.evaluate(p,t);t=state['temperatureK']
        return dict(dynamicViscosityPas=np.interp(t,self.T,self.mu),conductivityWmK=np.interp(t,self.T,self.kappa))
    def dictionary(self):
        text=super().dictionary();old='transport { mu 0; Pr 0.72; }'
        rows=lambda values:'\n'.join('(%.17g %.17g)'%v for v in zip(self.T,values))
        return text.replace(old,'transport { mu { values (\n'+rows(self.mu)+'\n); } kappa { values (\n'+rows(self.kappa)+'\n); } }')
    def viscous_contract(self):
        return dict(caloric=self.contract(),transport=TRANSPORT,interpolation='piecewise linear',
            solverTransportType='icoTabulated',nativeWallHeatFluxToIntoWallSign=-1,
            turbulentPrandtlNumber=None,physicalQualification=False)
