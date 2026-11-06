// 模块用途。实际CVariable原生更新与CNSVariable热状态验证，输入为独立构造及v166状态。
#include "variables/CVariable.hpp"
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2IdealAir.hpp"
#include "fluid/SSTEnergyUpdate.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <omp.h>

void check(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
int main(int argc,char** argv) {
  try {
    check(argc==4,"usage: fixture config inputs output");
    CConfig cfg(argv[1],SU2_COMPONENT::SU2_CFD);
    std::ifstream input(argv[2]);std::string line;std::getline(input,line);
    std::ofstream output(argv[3]);output<<std::setprecision(17);
    output<<"case,accepted,limited,alpha,native_alpha,k,omega,e,T,refinements\n";
    unsigned accepted=0,rejected=0,limited=0,zero=0;
    std::vector<rah_sst_energy::Input> valid_inputs;
    while(std::getline(input,line)) {
      std::stringstream stream(line);std::string part;std::vector<double> row;
      while(std::getline(stream,part,','))row.push_back(std::stod(part));
      check(row.size()==18,"input columns");
      const unsigned id=static_cast<unsigned>(row[0]);const bool expect_reject=row[1]!=0;
      rah_sst_energy::Input in;in.dimension=static_cast<unsigned>(row[2]);
      for(unsigned j=0;j<5;++j)in.flow[j]=row[3+j];
      for(unsigned j=0;j<2;++j) {in.old[j]=row[8+j];in.increment[j]=row[10+j];in.lower[j]=row[12+j];in.upper[j]=row[14+j];}
      in.native_alpha=row[16];const int expected_kind=static_cast<int>(row[17]);
      rah_sst_energy::Result result;
      try {result=rah_sst_energy::constrain(in);}
      catch(const std::domain_error&) {check(expect_reject,"unexpected kernel rejection");++rejected;output<<id<<",0,0,0,0,0,0,0,0,0\n";continue;}
      check(!expect_reject,"invalid input accepted");++accepted;
      check(static_cast<int>(result.limited)==(expected_kind!=0),"wrong limiting branch");
      if(result.limited)++limited;
      if(expected_kind==2){check(result.alpha==0.,"outward boundary must use zero step");++zero;}
      if(!result.limited)check(result.alpha==in.native_alpha,"legal native alpha changed");
      CVariable scalar(1,in.dimension,2,&cfg);
      // 刻意令current与old不同，真实AddClippedSolution必须只用old量。
      const double irrelevant[2]={123456,234567};scalar.SetSolution(0,irrelevant);scalar.SetSolution_Old(0,in.old.data());
      for(unsigned j=0;j<2;++j)scalar.AddClippedSolution(0,j,result.alpha*in.increment[j],in.lower[j],in.upper[j],in.flow[0],in.flow[0]);
      check(scalar.GetSolution(0,0)==result.k && scalar.GetSolution(0,1)==result.omega,"native clipped update differs from prediction");
      double velocity[3]={in.flow[1]/in.flow[0],in.flow[2]/in.flow[0],in.flow[3]/in.flow[0]};
      CNSVariable flow(in.flow[0],velocity,in.flow[4]/in.flow[0],1,in.dimension,in.dimension+2,&cfg);
      double conserved[5]={in.flow[0],in.flow[1],in.flow[2],in.flow[3],in.flow[4]};
      if(in.dimension==2)conserved[3]=in.flow[4];
      flow.SetSolution(0,conserved);flow.SetSolution_Old(0,conserved);CSU2IdealAir air;
      check(flow.SetPrimVar(0,0.,result.k,&air),"accepted turbulence makes actual CNSVariable reject flow");
      for(unsigned j=0;j<in.dimension+2;++j)check(flow.GetSolution(0,j)==conserved[j],"flow conserved component changed");
      check(flow.GetTemperature(0)>=60 && flow.GetTemperature(0)<=2000,"temperature outside original domain");
      valid_inputs.push_back(in);
      output<<id<<",1,"<<result.limited<<','<<result.alpha<<','<<in.native_alpha<<','<<result.k<<','<<result.omega<<','<<result.energy<<','<<flow.GetTemperature(0)<<','<<result.refinements<<'\n';
    }
    unsigned parallel_ok=0;
    #pragma omp parallel for num_threads(2) reduction(+:parallel_ok)
    for(int i=0;i<32;++i) {
      const auto& in=valid_inputs[i%valid_inputs.size()];const auto a=rah_sst_energy::constrain(in);
      const auto b=rah_sst_energy::constrain(in);
      if(a.alpha==b.alpha && a.k==b.k && a.omega==b.omega)++parallel_ok;
    }
    check(parallel_ok==32,"parallel result changed");
    std::cout<<"{\"passed\":true,\"accepted\":"<<accepted<<",\"rejected\":"<<rejected<<",\"limited\":"<<limited<<",\"zero\":"<<zero<<",\"parallel\":"<<parallel_ok<<"}\n";
    return 0;
  } catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
