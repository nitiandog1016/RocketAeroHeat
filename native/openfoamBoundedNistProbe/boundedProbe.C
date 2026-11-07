// 已知真解/直接能量两类标量输入，核验实际注册映射使用的同一模板。
#include "IFstream.H"
#include "dictionary.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "rahBoundedNistInternalEnergy.H"
#include "thermo.H"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

int main(int argc,char** argv)
{
    if (argc!=3) return 2;
    Foam::IFstream input(argv[1]); Foam::dictionary dict(input);
    typedef Foam::species::thermo<Foam::hIcoTabulatedThermo<Foam::perfectGas<Foam::specie>>,Foam::rahBoundedNistInternalEnergy> Gas;
    Gas gas("air",dict);
    std::ifstream samples(argv[2]);
    std::cout << std::setprecision(17);
    std::string mode,ps,vs,ss;
    while(samples>>mode>>ps>>vs>>ss)
    {
        const double p=std::stod(ps),v=std::stod(vs),seed=std::stod(ss);
        if (mode!="T" && mode!="E") return 3;
        const double e=(mode=="T")?gas.es(p,v):v;
        const double T=gas.The(e,p,seed);
        std::cout << "ROW " << p << ' ' << e << ' ' << seed << ' ' << T << ' '
            << gas.es(p,T) << ' ' << gas.Cp(p,T) << ' ' << gas.Cv(p,T) << ' '
            << gas.hs(p,T) << ' ' << gas.s(p,T) << ' ' << gas.rho(p,T) << ' '
            << std::sqrt(gas.gamma(p,T)*gas.R()*T) << '\n';
    }
    return samples.eof()?0:4;
}
