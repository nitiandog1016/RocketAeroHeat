// 模块：调用与实际注册相同的模板，检查新表、联合域和反算真实能量残差。
#include "IFstream.H"
#include "dictionary.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "rahStandardDryAirInternalEnergy.H"
#include "thermo.H"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>
int main(int argc,char** argv)
{
    if (argc!=3) return 2;
    Foam::IFstream input(argv[1]);Foam::dictionary dict(input);
    typedef Foam::species::thermo<Foam::hIcoTabulatedThermo<Foam::perfectGas<Foam::specie>>,Foam::rahStandardDryAirInternalEnergy> Gas;
    Gas gas("air",dict);std::ifstream samples(argv[2]);std::cout<<std::setprecision(17);
    std::string mode,ps,vs,ss;
    while(samples>>mode>>ps>>vs>>ss)
    {
        // strtod保留可表示的subnormal，用于核验真实密度下溢拒绝路径。
        const double p=std::strtod(ps.c_str(),nullptr),v=std::strtod(vs.c_str(),nullptr),seed=std::strtod(ss.c_str(),nullptr);
        if (mode!="T"&&mode!="E") return 3;
        const double e=mode=="T"?gas.es(p,v):v;
        const double T=gas.The(e,p,seed);
        std::cout<<"ROW "<<p<<' '<<e<<' '<<seed<<' '<<T<<' '<<gas.es(p,T)<<' '
            <<gas.Cp(p,T)<<' '<<gas.Cv(p,T)<<' '<<gas.hs(p,T)<<' '<<gas.s(p,T)<<' '
            <<gas.rho(p,T)<<' '<<std::sqrt(gas.gamma(p,T)*gas.R()*T)<<' '
            <<Foam::rahStandardDryAirInternalEnergy<Gas>::pressureCeiling(T)<<'\n';
    }
    return samples.eof()?0:4;
}
