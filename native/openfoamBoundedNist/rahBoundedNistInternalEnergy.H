// 仅替换具名NIST气体的温度反算映射；流动内核及物性表保持原生实现。
#ifndef rahBoundedNistInternalEnergy_H
#define rahBoundedNistInternalEnergy_H
#include "sensibleInternalEnergy.H"
#include "error.H"
#include <boost/math/tools/toms748_solve.hpp>
#include <cmath>
#include <algorithm>

namespace Foam
{
template<class Thermo>
class rahBoundedNistInternalEnergy : public sensibleInternalEnergy<Thermo>
{
public:
    static word typeName() { return "rahBoundedNistInternalEnergy"; }

    scalar The(const Thermo& thermo, const scalar e, const scalar p, const scalar T0) const
    {
        // 模块一：只接受原合同输入；种子只核验，不参与区间扩大或外推。
        if (!std::isfinite(e) || !std::isfinite(p) || !std::isfinite(T0)
            || p <= 0 || p > 100000 || T0 < 100 || T0 > 2500)
        {
            FatalErrorInFunction << "RAH_INVALID_THERMO_INPUT e=" << e
                << " p=" << p << " seed=" << T0 << exit(FatalError);
        }
        const scalar lo=100, hi=2500;
        const scalar elo=thermo.es(p,lo), ehi=thermo.es(p,hi);
        if (!std::isfinite(elo) || !std::isfinite(ehi) || elo >= ehi)
        {
            FatalErrorInFunction << "RAH_INVALID_ENERGY_BRACKET" << exit(FatalError);
        }
        if (e < elo || e > ehi)
        {
            FatalErrorInFunction << "RAH_ES_OUTSIDE_SUPPORTED_RANGE e=" << e
                << " elo=" << elo << " ehi=" << ehi << " p=" << p
                << " seed=" << T0 << exit(FatalError);
        }
        if (e == elo) return lo;
        if (e == ehi) return hi;

        // 模块二：成熟Boost TOMS748保留异号区间；绝不改变物理输入能量。
        auto residual = [&](scalar T)
        {
            if (!std::isfinite(T) || T < lo || T > hi)
            {
                FatalErrorInFunction << "RAH_ROOT_TRIAL_OUTSIDE_BRACKET" << exit(FatalError);
            }
            return thermo.es(p,T)-e;
        };
        auto tolerance = [](scalar a, scalar b) { return std::abs(b-a) <= 1e-9; };
        // 2次端点 + 至多61次内部求值 + 1次末态检查，总计不超过64次。
        std::uintmax_t count=61;
        const auto bracket=boost::math::tools::toms748_solve
        (
            residual,lo,hi,elo-e,ehi-e,tolerance,count
        );
        const scalar T=bracket.first+(bracket.second-bracket.first)/2;

        // 模块三：区间宽度及真实能量残差都通过才返回，不以迭代次数冒充成功。
        if (!tolerance(bracket.first,bracket.second)
            || std::abs(residual(T)) > 2e-10*std::max(scalar(1),std::abs(e)))
        {
            FatalErrorInFunction << "RAH_INVERSE_NOT_CONVERGED evaluations=" << count
                << " e=" << e << " p=" << p << " seed=" << T0 << exit(FatalError);
        }
        return T;
    }
};
}
#endif
