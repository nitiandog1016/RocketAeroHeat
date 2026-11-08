// 模块：Lemmon标准Air的最终气相/低密度范围；只检查，不改变任何状态。
#ifndef rahPhaseDomain_H
#define rahPhaseDomain_H
#include <cmath>
#include <algorithm>
namespace rahPhaseDomain
{
inline double ceiling(const double T)
{
    if (!std::isfinite(T)||T<60||T>2000) return -1;
    const double r=T/100;
    double upper=std::min(1000000.,10000.*((r*r)*r));
    if (T<=132.6312)
    {
        const double theta=1-T/132.6312;
        const double sum=-.1567266*std::sqrt(theta)-5.539635*theta
            +.7567212*std::pow(theta,2.5)-3.514322*std::pow(theta,4);
        const double dew=3785020.*std::exp(132.6312/T*sum);
        upper=std::min(upper,.1*dew);
    }
    return upper;
}
inline bool accepted(const double p,const double T)
{
    return std::isfinite(p)&&p>0&&p<=ceiling(T);
}
}
#endif
