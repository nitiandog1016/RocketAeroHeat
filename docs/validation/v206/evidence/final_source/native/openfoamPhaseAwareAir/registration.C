// 模块：独立注册新量热模型，不覆盖旧NIST/标准干空气/工程默认。
#include "psiThermo.H"
#include "pureMixture.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "rahPhaseAwareInternalEnergy.H"
#include "constTransport.H"
#include "thermo.H"
#include "forThermo.H"
#include "makeFluidThermo.H"
namespace Foam
{
    forThermo(constTransport,rahPhaseAwareInternalEnergy,hIcoTabulatedThermo,
        perfectGas,specie,makeFluidThermo,psiThermo,pureMixture);
}
