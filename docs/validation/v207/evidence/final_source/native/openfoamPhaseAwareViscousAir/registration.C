// 模块：将已验证量热反算与成熟温变mu/kappa表组合；不修改原生流动算子。
#include "psiThermo.H"
#include "pureMixture.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "rahPhaseAwareInternalEnergy.H"
#include "icoTabulatedTransport.H"
#include "thermo.H"
#include "forThermo.H"
#include "makeFluidThermo.H"
namespace Foam
{
    forThermo(icoTabulatedTransport,rahPhaseAwareInternalEnergy,hIcoTabulatedThermo,
        perfectGas,specie,makeFluidThermo,psiThermo,pureMixture);
}
