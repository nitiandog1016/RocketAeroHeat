// 模块职责：固定KINSOL全局线搜索接原完整flow-SA空间算子；不自研推进/差分内核。
#include "FullResidualProviderV257.hpp"
#include <kinsol/kinsol.h>
#include <kinsol/kinsol_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_logger.h>

class Coupled257 final : public FullResidualProviderV257 {
  SUNContext ctx=nullptr;SUNLinearSolver ls=nullptr;void* kin=nullptr;
  N_Vector u=nullptr,us=nullptr,fs=nullptr,constraints=nullptr;
  CSysVector<su2mixedfloat> fr,fd,tr,td;
  std::vector<double> current,input,output;
  std::ofstream events,frames,linear;
  double start=omp_get_wtime(),linearTarget=0;
  long setups=0,jvCalls=0,sysCalls=0,precSolves=0,observedChanges=0,rejectFiles=0;
  bool first=true,linearRefused=false,budgetStop=false;
  std::string error;
  const std::array<double,6> residualScale{{0.00021287581062475784,0.2749918398779103,
    0.0888042842240357,0.10889617958483502,237.73876378982774,1.018701759607621e-5}};
  static void check(int code){if(code!=0)throw std::runtime_error("SUNDIALS_API_"+std::to_string(code));}
  void readVector(N_Vector x){auto p=N_VGetArrayPointer(x);input.assign(p,p+entries.size());}
  long iteration()const{long n=0;if(kin)KINGetNumNonlinSolvIters(kin,&n);return n;}
  bool expired(){budgetStop=budgetStop||omp_get_wtime()-start>840;return budgetStop;}
  static bool sameFile(const std::string& a,const std::string& b){
    std::ifstream x(a,std::ios::binary),y(b,std::ios::binary);if(!x||!y)return false;
    std::array<char,65536> xx{},yy{};
    do{x.read(xx.data(),xx.size());y.read(yy.data(),yy.size());
      if(x.gcount()!=y.gcount()||std::memcmp(xx.data(),yy.data(),x.gcount()))return false;
    }while(x.gcount());return x.eof()&&y.eof();
  }
  // 模块一：保存同状态的完整六式残差、弱热壁导热量及真实当前/返回状态。
  void wallFrame(const std::string& label){
    auto* n=s[FLOW_SOL]->GetNodes();std::ofstream f(destination+"/"+label+"_wall.csv");
    f<<std::setprecision(17)<<"globalID,neighborID,pressurePa,wallFluidK,qIntoWallWm2,distanceM,neighborK,k,Cp,mut,Prt\n";
    unsigned count=0;
    for(unsigned m=0;m<c->GetnMarker_All();++m)if(c->GetMarker_All_KindBC(m)==ISOTHERMAL)
      for(unsigned v=0;v<g->nVertex[m];++v){auto i=g->vertex[m][v]->GetNode(),j=g->vertex[m][v]->GetNormal_Neighbor();
        auto normal=g->vertex[m][v]->GetNormal();double aa=0,projection=0;
        for(unsigned d=0;d<3;++d){aa+=normal[d]*normal[d];projection+=normal[d]*(g->nodes->GetCoord(j,d)-g->nodes->GetCoord(i,d));}
        double dist=std::abs(projection)/std::sqrt(aa),k=n->GetThermalConductivity(i),cp=n->GetSpecificHeatCp(i),
          mut=n->GetEddyViscosity(i),prt=c->GetPrandtl_Turb(),Tj=n->GetTemperature(j);
        double q=(k+cp*mut/prt)*(Tj-300)/dist;
        if(!(dist>0)||!std::isfinite(q))throw std::runtime_error("WALL_DIAGNOSTIC_INVALID");
        f<<g->nodes->GetGlobalIndex(i)<<','<<g->nodes->GetGlobalIndex(j)<<','<<n->GetPressure(i)<<','<<n->GetTemperature(i)<<','
          <<q<<','<<dist<<','<<Tj<<','<<k<<','<<cp<<','<<mut<<','<<prt<<'\n';++count;
      }
    f.close();if(!f||count!=9878)throw std::runtime_error("WALL_WRITE_OR_COUNT");
  }
  void snapshot(const std::string& label,const std::vector<double>& x){
    state(label+"_state.bin");full_residual(label);dump(label+"_unknowns.bin",x);wallFrame(label);
    std::array<long double,6> sums{};for(unsigned long i=0;i<g->GetnPoint();++i)for(unsigned j=0;j<6;++j){
      double r=j<5?s[FLOW_SOL]->LinSysRes(i,j):s[TURB_SOL]->LinSysRes(i,0);
      if(j==4)r-=300000*s[FLOW_SOL]->LinSysRes(i,0);sums[j]+=r*r;
    }
    frames<<label<<','<<iteration()<<','<<observedChanges<<','<<c->GetNonphysical_Reconstr();
    for(auto r:sums)frames<<','<<std::sqrt(double(r/g->GetnPoint()));frames<<'\n';frames.flush();
    if(!frames)throw std::runtime_error("FRAME_LOG_WRITE");
  }
  int failure(const std::exception& e){error=e.what();events<<"exception,"<<iteration()<<",-1,"<<error<<'\n';events.flush();return -1;}
  int evaluateVector(N_Vector x,N_Vector f,const std::string& role){
    readVector(x);int flag=evaluate(input,output,role);
    events<<role<<','<<iteration()<<','<<flag<<','<<c->GetNonphysical_Reconstr()<<'\n';events.flush();
    if(flag){
      // 拒绝态保留原向量，不发布伪残差。16份为输出预算上限，超限终止本控制。
      if(rejectFiles>=16){budgetStop=true;return -1;}
      dump("rejected_"+std::to_string(++rejectFiles)+"_"+role+".bin",input);return flag;
    }
    if(f)std::copy(output.begin(),output.end(),N_VGetArrayPointer(f));return 0;
  }
  int system(N_Vector x,N_Vector f,bool jv){
    try{
      if(expired())return -1;
      if(jv){++jvCalls;return evaluateVector(x,f,"native_DQ_Jv");}
      ++sysCalls;
      if(!first){
        auto flag=SUNLinSolLastFlag(ls);double norm=SUNLinSolResNorm(ls);
        bool pass=setups>0&&flag==SUN_SUCCESS&&std::isfinite(norm)&&norm>=0&&norm<=linearTarget;
        linear<<iteration()<<','<<flag<<','<<norm<<','<<linearTarget<<','<<pass<<'\n';linear.flush();
        if(!pass){linearRefused=true;readVector(x);dump("linear_refused_trial.bin",input);return -1;}
      }
      int code=evaluateVector(x,f,first?"initial":"nonlinear_trial");if(code)return code;
      if(first){
        current=input;snapshot("initial",input);
        if(!sameFile(destination+"/initial_state.bin","../../v256_validation/rocket02/initial_state.bin")||
           !sameFile(destination+"/initial_full_residual.bin","../../v256_validation/rocket02/base_full_residual.bin"))
          throw std::runtime_error("INITIAL_NOT_BITWISE_V256");
        first=false;
      }
      return 0;
    }catch(const std::exception& e){return failure(e);}
  }
  static int sysfn(N_Vector x,N_Vector f,void* p){return static_cast<Coupled257*>(p)->system(x,f,false);}
  static int jvfn(N_Vector x,N_Vector f,void* p){return static_cast<Coupled257*>(p)->system(x,f,true);}
  // 模块二：仅同状态空间Jacobian的原生块ILU0；无时间对角，无内部Krylov。
  int setup(N_Vector x,N_Vector fval,N_Vector fscale){
    try{
      if(expired())return -1;
      int code=evaluateVector(x,nullptr,"preconditioner_current_state");if(code)return -1;
      auto expected=N_VGetArrayPointer(fval);
      if(std::memcmp(expected,output.data(),output.size()*sizeof(double)))throw std::runtime_error("SETUP_RESIDUAL_NOT_SAME_STATE");
      if(input!=current){++observedChanges;current=input;snapshot("current_"+std::to_string(observedChanges),current);}
      linearTarget=(1e-5+std::numeric_limits<double>::epsilon())*N_VWL2Norm(fval,fscale);
      s[FLOW_SOL]->Jacobian.BuildILUPreconditioner();s[TURB_SOL]->Jacobian.BuildILUPreconditioner();++setups;
      events<<"ILU_setup,"<<iteration()<<",0,"<<linearTarget<<'\n';events.flush();return 0;
    }catch(const std::exception& e){return failure(e);}
  }
  static int psetup(N_Vector x,N_Vector,N_Vector f,N_Vector fs,void* p){return static_cast<Coupled257*>(p)->setup(x,f,fs);}
  int solve(N_Vector v){
    try{
      if(expired())return -1;
      ++precSolves;fr.SetValZero();tr.SetValZero();auto a=N_VGetArrayPointer(v);
      for(size_t k=0;k<entries.size();++k){auto e=entries[k];if(e.second<5)fr(e.first,e.second)=a[k];else tr(e.first,0)=a[k];}
      for(unsigned long i=0;i<g->GetnPoint();++i)fr(i,4)+=300000*fr(i,0);
      s[FLOW_SOL]->Jacobian.ComputeILUPreconditioner(fr,fd,g,c);s[TURB_SOL]->Jacobian.ComputeILUPreconditioner(tr,td,g,c);
      for(size_t k=0;k<entries.size();++k){auto e=entries[k];double val=e.second<5?fd(e.first,e.second):td(e.first,0);
        if(e.second==4)val-=300000*fd(e.first,0);
        if(!std::isfinite(val))throw std::runtime_error("ILU_NONFINITE_DIRECTION");a[k]=val;
      }return 0;
    }catch(const std::exception& e){return failure(e);}
  }
  static int psolve(N_Vector,N_Vector,N_Vector,N_Vector,N_Vector v,void* p){return static_cast<Coupled257*>(p)->solve(v);}
public:
  using FullResidualProviderV257::FullResidualProviderV257;
  // 模块三：固定原生求解器设置；统计尝试次数与状态实际变化分别保存。
  void run(const std::string& folder){
    initialize(folder);events.open(folder+"/events.csv");events<<std::setprecision(17)<<"role,iteration,status,detail\n";
    frames.open(folder+"/frames.csv");frames<<std::setprecision(17)<<"frame,nativeAttempt,observedCurrentChanges,nonphysicalRecon,mass,momX,momY,momZ,physicalEnergy,SA\n";
    linear.open(folder+"/linear_gate.csv");linear<<std::setprecision(17)<<"attempt,flag,reportedResidual,target,passed\n";
    const auto N=g->GetnPoint();fr.Initialize(N,N,5);fd.Initialize(N,N,5);tr.Initialize(N,N,1);td.Initialize(N,N,1);
    check(SUNContext_Create(SUN_COMM_NULL,&ctx));SUNLogger logger=nullptr;check(SUNContext_GetLogger(ctx,&logger));
    check(SUNLogger_SetInfoFilename(logger,(folder+"/kinsol_info.log").c_str()));
    check(SUNLogger_SetErrorFilename(logger,(folder+"/kinsol_error.log").c_str()));
    check(SUNLogger_SetWarningFilename(logger,(folder+"/kinsol_warning.log").c_str()));
    u=N_VNew_Serial(entries.size(),ctx);us=N_VClone(u);fs=N_VClone(u);constraints=N_VClone(u);
    if(!u||!us||!fs||!constraints)throw std::runtime_error("NVECTOR_ALLOCATION");
    const double rho=c->GetDensity_FreeStreamND(),mu=c->GetViscosity_FreeStreamND();
    const auto vel=c->GetVelocity_FreeStreamND();double speed=0;for(int j=0;j<3;++j)speed+=vel[j]*vel[j];speed=std::sqrt(speed);
    std::array<double,6> scales{{rho,rho*speed,rho*speed,rho*speed,rho*speed*speed,mu/rho}};
    auto initial=pack();for(size_t k=0;k<entries.size();++k){unsigned j=entries[k].second;
      if(!(scales[j]>0))throw std::runtime_error("NONPOSITIVE_SCALE");
      NV_Ith_S(u,k)=initial[k];NV_Ith_S(us,k)=1/scales[j];NV_Ith_S(fs,k)=1/residualScale[j];NV_Ith_S(constraints,k)=j==0?2:0;
    }
    std::ofstream scaleFile(folder+"/scaling.csv");scaleFile<<std::setprecision(17)<<"component,variableCharacteristic,residualCharacteristic\n";
    for(unsigned j=0;j<6;++j)scaleFile<<j<<','<<scales[j]<<','<<residualScale[j]<<'\n';scaleFile.close();
    initial.clear();initial.shrink_to_fit();
    kin=KINCreate(ctx);if(!kin)throw std::runtime_error("KIN_ALLOCATION");
    check(KINInit(kin,sysfn,u));check(KINSetUserData(kin,this));
    ls=SUNLinSol_SPGMR(u,SUN_PREC_RIGHT,20,ctx);if(!ls)throw std::runtime_error("SPGMR_ALLOCATION");
    check(SUNLinSol_SPGMRSetMaxRestarts(ls,0));check(KINSetLinearSolver(kin,ls,nullptr));
    check(KINSetPreconditioner(kin,psetup,psolve));check(KINSetJacTimesVecSysFn(kin,jvfn));
    check(KINSetEtaForm(kin,KIN_ETACONSTANT));check(KINSetEtaConstValue(kin,1e-5));check(KINSetNoMinEps(kin,SUNTRUE));
    check(KINSetMaxSetupCalls(kin,1));check(KINSetNumMaxIters(kin,3));check(KINSetFuncNormTol(kin,1e-8));
    check(KINSetScaledStepTol(kin,1e-12));check(KINSetConstraints(kin,constraints));
    const int flag=KINSol(kin,u,KIN_LINESEARCH,us,fs);
    long nni=0,nfe=0,back=0,nli=0,linfail=0,lastlin=0,jtimes=0;
    double norm=0,step=0;KINGetNumNonlinSolvIters(kin,&nni);KINGetNumFuncEvals(kin,&nfe);KINGetNumBacktrackOps(kin,&back);
    KINGetNumLinIters(kin,&nli);KINGetNumLinConvFails(kin,&linfail);KINGetLastLinFlag(kin,&lastlin);
    KINGetNumJtimesEvals(kin,&jtimes);KINGetFuncNorm(kin,&norm);KINGetStepLength(kin,&step);
    // 返回uu先保存，再刷新同状态诊断；与最后的Jv试探或被拒状态严格分开。
    readVector(u);dump("returned_unknowns.bin",input);const bool changed=input!=current;
    if(changed)++observedChanges;
    int terminal=evaluate(input,output,"returned_state_refresh");
    if(terminal)throw std::runtime_error("RETURNED_STATE_REFUSED");snapshot("terminal",input);
    char* flagName=KINGetReturnFlagName(flag);
    std::ofstream report(folder+"/native_completion.json");report<<std::boolalpha<<std::setprecision(17)
      <<"{\"normalReturn\":true,\"KINSOLExecuted\":true,\"nativeFlag\":"<<flag<<",\"nativeFlagName\":\""<<flagName
      <<"\",\"unknowns\":"<<entries.size()<<",\"nativeNonlinearAttempts\":"<<nni<<",\"observedCurrentStateChanges\":"<<observedChanges
      <<",\"nativeFunctionCalls\":"<<nfe<<",\"nativeLinearIterations\":"<<nli<<",\"nativeLinearFailures\":"<<linfail<<",\"nativeLastLinearFlag\":"<<lastlin
      <<",\"nativeJvCalls\":"<<jtimes<<",\"nativeBacktracks\":"<<back<<",\"nativeWeightedL2Norm\":"<<norm<<",\"nativeStepLength\":"<<step
      <<",\"setups\":"<<setups<<",\"preconditionerSolves\":"<<precSolves<<",\"jvSystemCalls\":"<<jvCalls<<",\"systemCalls\":"<<sysCalls
      <<",\"fullPhysicalEvaluations\":"<<physicalCalls<<",\"recoverableRejectedStates\":"<<rejected<<",\"linearGateRefused\":"<<linearRefused
      <<",\"budgetStop\":"<<budgetStop<<",\"interfaceError\":\""<<error<<"\",\"initialStateAndResidualBitwiseV256\":"<<!first
      <<",\"timeUpdates\":0,\"physicalQualification\":false,\"newQualifiedAreaM2\":0,\"profileReadyForCFD\":false,\"elapsedSeconds\":"<<omp_get_wtime()-start<<"}\n";
    report.close();free(flagName);if(!report)throw std::runtime_error("COMPLETION_WRITE");
    rah_single_observer::rejection=nullptr;active=nullptr;
    KINFree(&kin);SUNLinSolFree(ls);ls=nullptr;N_VDestroy(u);N_VDestroy(us);N_VDestroy(fs);N_VDestroy(constraints);SUNContext_Free(&ctx);
  }
};
int main(int argc,char** argv){
  std::string destination=argc>2?argv[2]:".";
  try{if(argc!=3)throw std::runtime_error("CONFIG_OUTPUT_REQUIRED");SU2_MPI::Init(&argc,&argv);
    Coupled257 driver(argv[1],1,SU2_MPI::GetComm());driver.run(destination);driver.Finalize();SU2_MPI::Finalize();return 0;
  }catch(const std::exception& e){std::ofstream f(destination+"/native_failure.txt");f<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}
}
