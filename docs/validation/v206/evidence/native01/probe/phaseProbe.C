// 模块：独立检查旧压力反算与新压力接受两个不同阶段，调用实际注册模板。
#include "IFstream.H"
#include "dictionary.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "rahPhaseAwareInternalEnergy.H"
#include "rahPhaseDomain.H"
#include "thermo.H"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
int main(int argc,char** argv)
{
    if (argc!=3) return 2;
    Foam::IFstream input(argv[1]);Foam::dictionary dict(input);
    typedef Foam::species::thermo<Foam::hIcoTabulatedThermo<Foam::perfectGas<Foam::specie>>,Foam::rahPhaseAwareInternalEnergy> Gas;
    Gas gas("air",dict);std::ifstream samples(argv[2]);std::cout<<std::setprecision(17);
    std::string mode,ps,vs,ss,fs;
    while(samples>>mode>>ps>>vs>>ss>>fs)
    {
        const double pOld=std::strtod(ps.c_str(),nullptr),v=std::strtod(vs.c_str(),nullptr);
        const double seed=std::strtod(ss.c_str(),nullptr),pFinal=std::strtod(fs.c_str(),nullptr);
        if (mode!="T"&&mode!="E") return 3;
        const double e=mode=="T"?gas.es(pOld,v):v,T=gas.The(e,pOld,seed);
        if (!rahPhaseDomain::accepted(pFinal,T))
        {
            std::cerr<<"RAH_PHASE_PROBE_FINAL_REJECT p="<<pFinal<<" T="<<T<<'\n';return 10;
        }
        std::cout<<"ROW "<<pOld<<' '<<e<<' '<<seed<<' '<<T<<' '<<gas.es(pFinal,T)<<' '<<pFinal
            <<' '<<gas.Cp(pFinal,T)<<' '<<gas.Cv(pFinal,T)<<' '<<gas.hs(pFinal,T)<<' '<<gas.s(pFinal,T)
            <<' '<<gas.rho(pFinal,T)<<' '<<std::sqrt(gas.gamma(pFinal,T)*gas.R()*T)
            <<' '<<rahPhaseDomain::ceiling(T)<<'\n';
    }
    return samples.eof()?0:4;
}
