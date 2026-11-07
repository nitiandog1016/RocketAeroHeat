// 模块一：读取同一原生Cp表，查询实际成熟类；只检验物性，不推进流场。
#include "IFstream.H"
#include "dictionary.H"
#include "specie.H"
#include "perfectGas.H"
#include "hIcoTabulatedThermo.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

int main(int argc, char** argv)
{
    if (argc != 3) return 2;
    Foam::IFstream input(argv[1]);
    Foam::dictionary dict(input);
    typedef Foam::species::thermo<Foam::hIcoTabulatedThermo<Foam::perfectGas<Foam::specie>>,Foam::sensibleInternalEnergy> Gas;
    Gas gas("air",dict);
    std::cout << std::setprecision(17);
    std::cout << "META " << gas.R() << ' ' << gas.hf() << ' '
        << Foam::constant::thermodynamic::Tstd << ' ' << Foam::constant::thermodynamic::Pstd << '\n';
    // 模块二：样本由输入文件声明，含端点/分段/非网格节点；逆算种子也明确保存。
    std::ifstream samples(argv[2]);
    double p,T,seed;
    while(samples>>p>>T>>seed)
    {
        std::cout << "ROW " << p << ' ' << T << ' ' << gas.Cp(p,T) << ' ' << gas.Cv(p,T)
            << ' ' << gas.ha(p,T) << ' ' << gas.hs(p,T) << ' ' << gas.es(p,T) << ' ' << gas.s(p,T)
            << ' ' << gas.rho(p,T) << ' ' << std::sqrt(gas.gamma(p,T)*gas.R()*T)
            << ' ' << gas.Tes(gas.es(p,T),p,seed) << '\n';
    }
    return samples.eof()?0:3;
}
