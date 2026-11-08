// 模块：固定来源60—2000K量热反算；最终新压力由独立每步状态守卫核验。
#ifndef rahPhaseAwareInternalEnergy_H
#define rahPhaseAwareInternalEnergy_H
#include "sensibleInternalEnergy.H"
#include "error.H"
#include <boost/math/tools/toms748_solve.hpp>
#include <cmath>
#include <algorithm>
namespace Foam
{
template<class Thermo>
class rahPhaseAwareInternalEnergy : public sensibleInternalEnergy<Thermo>
{
public:
    static word typeName(){return "rahPhaseAwareInternalEnergy";}
    scalar The(const Thermo& thermo,const scalar e,const scalar pOld,const scalar seed) const
    {
        // pOld属于压力修正之前的输入，不能冒充最终新压力的气相判定。
        if (!std::isfinite(e)||!std::isfinite(pOld)||pOld<=0||!std::isfinite(seed)||seed<60||seed>2000)
            FatalErrorInFunction<<"RAH_PHASE_INVERSE_INPUT e="<<e<<" pOld="<<pOld<<" seed="<<seed<<exit(FatalError);
        const scalar lo=60,hi=2000,elo=thermo.es(pOld,lo),ehi=thermo.es(pOld,hi);
        if (!std::isfinite(elo)||!std::isfinite(ehi)||elo>=ehi||e<elo||e>ehi)
            FatalErrorInFunction<<"RAH_PHASE_ENERGY_OUTSIDE_TABLE e="<<e<<" elo="<<elo<<" ehi="<<ehi
                <<" pOld="<<pOld<<" seed="<<seed<<exit(FatalError);
        if (e==elo) return lo;
        if (e==ehi) return hi;
        auto residual=[&](scalar T)
        {
            if (!std::isfinite(T)||T<lo||T>hi)
                FatalErrorInFunction<<"RAH_PHASE_ROOT_OUTSIDE_TABLE"<<exit(FatalError);
            return thermo.es(pOld,T)-e;
        };
        auto tolerance=[](scalar a,scalar b){return std::abs(b-a)<=1e-9;};
        std::uintmax_t count=60;
        const auto bracket=boost::math::tools::toms748_solve(residual,lo,hi,elo-e,ehi-e,tolerance,count);
        const scalar T=bracket.first+(bracket.second-bracket.first)/2;
        if (!tolerance(bracket.first,bracket.second)||std::abs(residual(T))>2e-10*std::max(scalar(1),std::abs(e)))
            FatalErrorInFunction<<"RAH_PHASE_INVERSE_NOT_CONVERGED"<<exit(FatalError);
        return T;
    }
};
}
#endif
