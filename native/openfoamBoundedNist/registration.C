// 模块一：具名能量映射通过原生宏注册；不覆盖v202旧组合。
#include "psiThermo.H"
#include "pureMixture.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "rahBoundedNistInternalEnergy.H"
#include "constTransport.H"
#include "thermo.H"
#include "forThermo.H"
#include "makeFluidThermo.H"
namespace Foam
{
    forThermo
    (
        constTransport, rahBoundedNistInternalEnergy, hIcoTabulatedThermo,
        perfectGas, specie, makeFluidThermo, psiThermo, pureMixture
    );
}
