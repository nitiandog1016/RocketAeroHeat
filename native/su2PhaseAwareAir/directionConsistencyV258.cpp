// 模块职责：一次失败方向差分/块预条件诊断；不执行KINSOL、时间步或新的非线性更新。
#include "FullResidualProviderV257.hpp"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>

namespace budget258 { double entered=0; }
class Direction258 final : public FullResidualProviderV257 {
  CSysVector<su2mixedfloat> fr,fz,fp,tr,tz,tp;
  std::array<double,6> scales{},rscales{};
  std::ofstream events;
  void timeCheck(){if(omp_get_wtime()-budget258::entered>150)throw std::runtime_error("SOFT_BUDGET_BEFORE_OPERATION");}
  std::vector<double> readVector(const std::string& path){
    std::ifstream f(path,std::ios::binary);char magic[8];uint64_t n=0;f.read(magic,8);f.read(reinterpret_cast<char*>(&n),8);
    if(!f||std::memcmp(magic,"RAH256V1",8)||n!=entries.size())throw std::runtime_error("SOURCE_VECTOR_IDENTITY");
    std::vector<double> x(n);f.read(reinterpret_cast<char*>(x.data()),n*8);if(!f||f.peek()!=EOF)throw std::runtime_error("SOURCE_VECTOR_LENGTH");return x;
  }
  bool sameFile(const std::string& a,const std::string& b){
    std::ifstream x(a,std::ios::binary),y(b,std::ios::binary);if(!x||!y)return false;std::array<char,65536> xx{},yy{};
    do{x.read(xx.data(),xx.size());y.read(yy.data(),yy.size());if(x.gcount()!=y.gcount()||std::memcmp(xx.data(),yy.data(),x.gcount()))return false;}while(x.gcount());
    return x.eof()&&y.eof();
  }
  // 模块一：完整状态/原六式残差逐阶段保存；拒绝时保留真实输入，不造残差。
  void evaluateSaved(const std::vector<double>& x,std::vector<double>& f,const std::string& label){
    timeCheck();int status=evaluate(x,f,label);
    events<<label<<','<<status<<','<<c->GetNonphysical_Reconstr()<<','<<omp_get_wtime()-budget258::entered<<'\n';events.flush();
    if(status){dump(label+"_refused.bin",x);throw std::runtime_error("PROBE_STATE_REFUSED_"+label);}
    state(label+"_state.bin");full_residual(label);
  }
  // 模块二：自由向量→完整原生能量坐标；强约束行填其真实固定零方向。
  void expand(const std::vector<double>& v){
    fr.SetValZero();tr.SetValZero();for(size_t k=0;k<entries.size();++k){auto e=entries[k];
      if(e.second<5)fr(e.first,e.second)=v[k];else tr(e.first,0)=v[k];}
    for(unsigned long i=0;i<g->GetnPoint();++i)fr(i,4)+=300000*fr(i,0);
  }
  void pairFile(const std::string& name,const CSysVector<su2mixedfloat>& af,const CSysVector<su2mixedfloat>& at,
                const CSysVector<su2mixedfloat>& bf,const CSysVector<su2mixedfloat>& bt){
    std::ofstream f(destination+"/"+name,std::ios::binary);uint64_t n=g->GetnPoint(),cols=13;
    f.write("RAH258B1",8);f.write(reinterpret_cast<char*>(&n),8);f.write(reinterpret_cast<char*>(&cols),8);
    for(unsigned long i=0;i<n;++i){std::array<double,13>a{};a[0]=g->nodes->GetGlobalIndex(i);
      for(unsigned j=0;j<6;++j){a[1+j]=j<5?af(i,j):at(i,0);a[7+j]=j<5?bf(i,j):bt(i,0);
        if(!std::isfinite(a[1+j])||!std::isfinite(a[7+j]))throw std::runtime_error("NATIVE_MATRIX_NONFINITE");}
      f.write(reinterpret_cast<char*>(a.data()),sizeof(a));
    }f.close();if(!f)throw std::runtime_error("MATRIX_PAIR_WRITE");
  }
public:
  using FullResidualProviderV257::FullResidualProviderV257;
  void run(const std::string& folder){
    initialize(folder);timeCheck();events.open(folder+"/events.csv");events<<std::setprecision(17)<<"frame,status,nonphysicalRecon,elapsedFromProcessEntry\n";
    const auto N=g->GetnPoint();fr.Initialize(N,N,5);fz.Initialize(N,N,5);fp.Initialize(N,N,5);tr.Initialize(N,N,1);tz.Initialize(N,N,1);tp.Initialize(N,N,1);
    auto base=readVector("../../v257_validation/rocket01/initial_unknowns.bin"),direction=readVector("../../v257_validation/rocket01/linear_refused_trial.bin");
    auto actual=pack();if(actual!=base)throw std::runtime_error("INITIAL_UNKNOWN_NOT_V257");actual.clear();actual.shrink_to_fit();
    for(size_t k=0;k<base.size();++k)direction[k]-=base[k];dump("saved_direction.bin",direction);
    std::vector<double> f0,probe;evaluateSaved(base,f0,"initial");
    if(!sameFile(destination+"/initial_state.bin","../../v257_validation/rocket01/initial_state.bin")||
       !sameFile(destination+"/initial_full_residual.bin","../../v257_validation/rocket01/initial_full_residual.bin"))throw std::runtime_error("INITIAL_NOT_BITWISE_V257");
    expand(direction);s[FLOW_SOL]->Jacobian.MatrixVectorProduct(fr,fp,g,c);s[TURB_SOL]->Jacobian.MatrixVectorProduct(tr,tp,g,c);
    pairFile("direction_matrix.bin",fr,tr,fp,tp);
    // 唯一右端为原完整F的负值；仅测试已采用的纯空间块ILU，不额外加时间对角。
    probe=f0;for(auto& v:probe)v=-v;expand(probe);timeCheck();
    s[FLOW_SOL]->Jacobian.BuildILUPreconditioner();s[TURB_SOL]->Jacobian.BuildILUPreconditioner();
    s[FLOW_SOL]->Jacobian.ComputeILUPreconditioner(fr,fz,g,c);s[TURB_SOL]->Jacobian.ComputeILUPreconditioner(tr,tz,g,c);
    s[FLOW_SOL]->Jacobian.MatrixVectorProduct(fz,fp,g,c);s[TURB_SOL]->Jacobian.MatrixVectorProduct(tz,tp,g,c);
    pairFile("ilu_matrix.bin",fz,tz,fp,tp);
    // 模块三：固定原生DQ步长及half。所有向量积/加法使用同一serial NVector实现。
    std::ifstream sf("../../v257_validation/rocket01/scaling.csv");std::string line;std::getline(sf,line);
    for(unsigned j=0;j<6;++j){if(!std::getline(sf,line))throw std::runtime_error("SCALING_ROWS");std::replace(line.begin(),line.end(),',',' ');
      std::istringstream stream(line);unsigned index;stream>>index>>scales[j]>>rscales[j];if(!stream||index!=j||!(scales[j]>0&&rscales[j]>0))throw std::runtime_error("SCALING_VALUES");}
    SUNContext context=nullptr;if(SUNContext_Create(SUN_COMM_NULL,&context))throw std::runtime_error("SUN_CONTEXT");
    N_Vector u=N_VNew_Serial(base.size(),context),d=N_VClone(u),us=N_VClone(u),work=N_VClone(u),trial=N_VClone(u);
    if(!u||!d||!us||!work||!trial)throw std::runtime_error("NVECTOR_ALLOC");
    for(size_t k=0;k<base.size();++k){NV_Ith_S(u,k)=base[k];NV_Ith_S(d,k)=direction[k];NV_Ith_S(us,k)=1/scales[entries[k].second];}
    N_VProd(d,us,work);N_VProd(u,us,trial);
    const double dot=N_VDotProd(trial,work),vv=N_VDotProd(work,work),l1=N_VL1Norm(work);
    const double sigma=(dot>=0?1:-1)*std::sqrt(std::numeric_limits<double>::epsilon())*std::max(std::abs(dot),l1)/vv;
    if(!(vv>0)||!std::isfinite(sigma)||sigma==0)throw std::runtime_error("INVALID_FIXED_DIRECTION");
    std::ofstream step(folder+"/step.json");step<<std::setprecision(17)<<"{\"sigma\":"<<sigma<<",\"weightedDotUD\":"<<dot
      <<",\"weightedDirectionSquaredNorm\":"<<vv<<",\"weightedDirectionL1\":"<<l1<<",\"machineEpsilon\":"<<std::numeric_limits<double>::epsilon()
      <<",\"directionIsConstrainedRejectedTrialDifference\":true,\"nonlinearUpdates\":0}\n";step.close();if(!step)throw std::runtime_error("STEP_WRITE");
    std::vector<double> x(base.size());
    N_VLinearSum(1,u,sigma,d,trial);std::copy(N_VGetArrayPointer(trial),N_VGetArrayPointer(trial)+x.size(),x.begin());evaluateSaved(x,probe,"sigma");
    N_VLinearSum(1,u,sigma/2,d,trial);std::copy(N_VGetArrayPointer(trial),N_VGetArrayPointer(trial)+x.size(),x.begin());evaluateSaved(x,probe,"half");
    evaluateSaved(base,probe,"terminal_restored");
    if(probe!=f0)throw std::runtime_error("RESTORE_F_CHANGED");
    std::ofstream finish(folder+"/native_completion.json");finish<<std::setprecision(17)
      <<"{\"normalReturn\":true,\"fullPhysicalEvaluations\":"<<physicalCalls<<",\"recoverableRejects\":"<<rejected
      <<",\"initialBitwiseV257\":true,\"terminalIsRestoredInitial\":true,\"nonlinearUpdates\":0,\"timeUpdates\":0,\"ILUSetups\":1,\"ILUSolves\":1,\"matrixProducts\":2,\"KINSOLExecuted\":false,\"physicalQualification\":false,\"newQualifiedAreaM2\":0,\"elapsedFromProcessEntry\":"<<omp_get_wtime()-budget258::entered<<"}\n";
    finish.close();if(!finish)throw std::runtime_error("COMPLETION_WRITE");
    N_VDestroy(u);N_VDestroy(d);N_VDestroy(us);N_VDestroy(work);N_VDestroy(trial);SUNContext_Free(&context);rah_single_observer::rejection=nullptr;active=nullptr;
  }
};
int main(int argc,char** argv){
  budget258::entered=omp_get_wtime();std::string folder=argc>2?argv[2]:".";
  try{if(argc!=3)throw std::runtime_error("CONFIG_OUTPUT_REQUIRED");SU2_MPI::Init(&argc,&argv);
    Direction258 driver(argv[1],1,SU2_MPI::GetComm());driver.run(folder);driver.Finalize();SU2_MPI::Finalize();return 0;
  }catch(const std::exception& e){std::ofstream f(folder+"/native_failure.txt");f<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}
}
