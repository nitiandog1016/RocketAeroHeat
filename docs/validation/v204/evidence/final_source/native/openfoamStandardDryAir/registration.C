// 模块：新具名标准干空气注册，不覆盖原NIST或原生默认气体。
#include "psiThermo.H"
#include "pureMixture.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "rahStandardDryAirInternalEnergy.H"
#include "constTransport.H"
#include "thermo.H"
#include "forThermo.H"
#include "makeFluidThermo.H"
namespace Foam
{
    forThermo(constTransport,rahStandardDryAirInternalEnergy,hIcoTabulatedThermo,
        perfectGas,specie,makeFluidThermo,psiThermo,pureMixture);
}
