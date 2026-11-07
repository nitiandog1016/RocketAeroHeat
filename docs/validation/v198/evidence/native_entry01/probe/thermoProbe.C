// 模块一：只调用发行包原生物性类，核对气体常数、能量零点与温度逆算；不是CFD内核。
#include "IFstream.H"
#include "dictionary.H"
#include "specie.H"
#include "perfectGas.H"
#include "janafThermo.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
int main(int argc, char** argv)
{
    if (argc != 2) return 2;
    Foam::IFstream input(argv[1]);
    Foam::dictionary dict(input);
    typedef Foam::species::thermo<Foam::janafThermo<Foam::perfectGas<Foam::specie>>,Foam::sensibleInternalEnergy> Gas;
    Gas gas("air",dict);
    std::cout << std::setprecision(17);
    std::cout << "R " << gas.R() << " hf " << gas.hf() << " Tstd " << Foam::Tstd << '\n';
    // 模块二：不跨分段点求导，逆算种子是邻近温度；端点直接用端点作为种子。
    double samples[]={200,226.50908361133003,298.15,300,999,1000,1001,1500,2500};
    for(double T:samples)
    {
        double p=5000, seed=(T==200 || T==2500 || T==1000)?T:T*1.000001;
        double a=std::sqrt(gas.gamma(p,T)*gas.R()*T);
        std::cout << "ROW " << T << ' ' << gas.Cp(p,T) << ' ' << gas.Cv(p,T)
            << ' ' << gas.ha(p,T) << ' ' << gas.hs(p,T) << ' ' << gas.es(p,T)
            << ' ' << gas.rho(p,T) << ' ' << a << ' ' << gas.Tes(gas.es(p,T),p,seed) << '\n';
    }
    return 0;
}
