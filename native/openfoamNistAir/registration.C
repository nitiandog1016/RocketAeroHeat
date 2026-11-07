// 模块一：仅注册固定发行包的成熟模板组合，不实现或改变流动/热物性算法。
#include "psiThermo.H"
#include "pureMixture.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "sensibleInternalEnergy.H"
#include "constTransport.H"
#include "thermo.H"
#include "forThermo.H"
#include "makeFluidThermo.H"

namespace Foam
{
    forThermo
    (
        constTransport, sensibleInternalEnergy, hIcoTabulatedThermo,
        perfectGas, specie, makeFluidThermo, psiThermo, pureMixture
    );
}
