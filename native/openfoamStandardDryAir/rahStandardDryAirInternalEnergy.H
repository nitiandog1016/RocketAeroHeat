// 标准干空气独立能量映射：保留成熟有界反算，只补充明确的联合物性域。
#ifndef rahStandardDryAirInternalEnergy_H
#define rahStandardDryAirInternalEnergy_H
#include "sensibleInternalEnergy.H"
#include "error.H"
#include <boost/math/tools/toms748_solve.hpp>
#include <cmath>
#include <algorithm>
#include <limits>
namespace Foam
{
template<class Thermo>
class rahStandardDryAirInternalEnergy : public sensibleInternalEnergy<Thermo>
{
public:
    static word typeName() { return "rahStandardDryAirInternalEnergy"; }
    static scalar pressureCeiling(const scalar T)
    {
        const scalar ratio=T/100;
        return std::min(scalar(1000000),scalar(10000)*((ratio*ratio)*ratio));
    }
    scalar The(const Thermo& thermo,const scalar e,const scalar p,const scalar seed) const
    {
        // 模块一：合法输入与真实能量区间；不修改es/p/种子或物理温度。
        if (!std::isfinite(e)||!std::isfinite(p)||!std::isfinite(seed)
            ||p<=0||p>1000000||seed<100||seed>2000)
            FatalErrorInFunction << "RAH_STANDARD_INVALID_INPUT e=" << e
                << " p=" << p << " seed=" << seed << exit(FatalError);
        const scalar hi=2000;
        const scalar tableLow=thermo.es(p,100),ehi=thermo.es(p,hi);
        if (!std::isfinite(tableLow)||!std::isfinite(ehi)||tableLow>=ehi)
            FatalErrorInFunction << "RAH_STANDARD_INVALID_BRACKET" << exit(FatalError);
        if (e<tableLow||e>ehi)
            FatalErrorInFunction << "RAH_STANDARD_ENERGY_OUTSIDE_TABLE e=" << e
                << " elo=" << tableLow << " ehi=" << ehi << " p=" << p << exit(FatalError);
        // 联合域本身确定求根下界；至多8个浮点邻值定位，不裁输入或扩大域。
        scalar lo=std::max(scalar(100),scalar(100)*std::cbrt(p/10000));
        for (unsigned i=0;i<8;++i)
        {
            const scalar previous=std::nextafter(lo,-std::numeric_limits<scalar>::infinity());
            if (p>pressureCeiling(lo)) lo=std::nextafter(lo,std::numeric_limits<scalar>::infinity());
            else if (previous>=100&&p<=pressureCeiling(previous)) lo=previous;
            else break;
        }
        const scalar previous=std::nextafter(lo,-std::numeric_limits<scalar>::infinity());
        if (p>pressureCeiling(lo)||(previous>=100&&p<=pressureCeiling(previous)))
            FatalErrorInFunction << "RAH_STANDARD_BOUNDARY_NOT_LOCATED" << exit(FatalError);
        const scalar elo=thermo.es(p,lo);
        if (e<elo)
            FatalErrorInFunction << "RAH_STANDARD_PT_OUTSIDE_DOMAIN e=" << e
                << " p=" << p << " minimumT=" << lo << " minimumE=" << elo << exit(FatalError);
        // 模块二：所有返回路径都核验联合域，反算试探值不当作已接受物理状态。
        auto accepted=[&](scalar T)->scalar
        {
            if (!std::isfinite(T)||T<lo||T>hi||p>pressureCeiling(T))
                FatalErrorInFunction << "RAH_STANDARD_PT_OUTSIDE_DOMAIN e=" << e
                    << " p=" << p << " T=" << T << " pMax=" << pressureCeiling(T) << exit(FatalError);
            const scalar rho=thermo.rho(p,T);
            if (!std::isfinite(rho)||rho<=0)
                FatalErrorInFunction << "RAH_STANDARD_NUMERICAL_RANGE" << exit(FatalError);
            return T;
        };
        if (e==elo) return accepted(lo);
        if (e==ehi) return accepted(hi);
        auto residual=[&](scalar T)
        {
            if (!std::isfinite(T)||T<lo||T>hi)
                FatalErrorInFunction << "RAH_STANDARD_ROOT_OUTSIDE_TABLE" << exit(FatalError);
            return thermo.es(p,T)-e;
        };
        auto tolerance=[](scalar a,scalar b){return std::abs(b-a)<=1e-9;};
        std::uintmax_t count=60;
        const auto bracket=boost::math::tools::toms748_solve(residual,lo,hi,elo-e,ehi-e,tolerance,count);
        const scalar T=bracket.first+(bracket.second-bracket.first)/2;
        // 模块三：总求值上限64次及双重数值门槛，旧映射与求解器均保持。
        if (!tolerance(bracket.first,bracket.second)||std::abs(residual(T))>2e-10*std::max(scalar(1),std::abs(e)))
            FatalErrorInFunction << "RAH_STANDARD_INVERSE_NOT_CONVERGED" << exit(FatalError);
        return accepted(T);
    }
};
}
#endif
