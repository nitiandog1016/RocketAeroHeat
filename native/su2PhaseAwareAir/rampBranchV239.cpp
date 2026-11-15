// 模块职责：成熟阶次日程只读控制与真实flow-SA入口；不改通量、物性、松弛或原迭代顺序。
#include "hllcCoreV239.hpp"
namespace ramp239 {
bool control=false;std::ofstream trace;unsigned calls=0,monitors=0;double used=0;
double expected(unsigned step){return .5*(1-std::cos(PI_NUMBER*std::min(100u,step>1?step-2:0u)/100.));}
void check_identity(CConfig* c){
 const auto p=c->GetMUSCLRampParam();
 if(!c->GetMUSCL_Flow()||!c->GetMUSCLRamp()||p.Kind_MUSCLRamp!=MUSCL_RAMP_TYPE::SMOOTH_FUNCTION||
    p.rampMUSCLCoeff[0]!=0||p.rampMUSCLCoeff[1]!=1||p.rampMUSCLCoeff[2]!=100||p.RampMUSCLPower!=1||
    c->GetNewtonKrylovRelaxation()!=1||c->GetMUSCL_Kappa_Flow()!=0)throw std::runtime_error("RAMP239_IDENTITY");
}
void record(const char* phase,CConfig* c,CSolver* solver){
 const unsigned step=c->GetInnerIter()+1;const double r=c->GetMUSCLRampValue(),e=std::string(phase)=="after_monitor"?.5*(1-std::cos(PI_NUMBER*std::min(100u,step-1)/100.)):expected(step);
 if(std::abs(r-e)>2e-13)throw std::runtime_error("RAMP239_STAGE_MISMATCH");
 double grad[3]={1,2,-.5},edge[3]={.3,-.2,.4};
 const double projection=solver->MUSCL_Reconstruction(grad,edge,.7,c->GetMUSCL_Kappa_Flow(),r*c->GetNewtonKrylovRelaxation());
 trace<<step<<','<<phase<<','<<r<<','<<c->GetNewtonKrylovRelaxation()<<','<<projection<<'\n';trace.flush();
 if(!trace)throw std::runtime_error("RAMP239_WRITE");
}
void stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short system){
 if(c->GetContainerPosition(system)==FLOW_SOL&&std::string(label)=="after_sources")record("flow_after_sources",c,s[FLOW_SOL]);
 events238::stage(label,g,s,c,system);
}
}

// 模块一：控制只走成熟调度；生产Iterate仍调用完整父类，Monitor后才更新下一次使用的阶次。
class RampIterationV239 final:public CFluidIteration {
 public:using CFluidIteration::CFluidIteration;
 void Iterate(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,
   CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   ramp239::check_identity(c[z]);ramp239::record("before_iterate",c[z],s[z][inst][0][FLOW_SOL]);++ramp239::calls;
   if(ramp239::control)return;
   observe::phase="native-iterate/preprocessing";CFluidIteration::Iterate(o,in,g,s,n,c,sm,gm,b,z,inst);
   observe::phase="accepted-flow-SA/end-iteration";observe::end_iteration();branch236::benefit();
 }
 bool Monitor(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,
   CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   const bool native=CFluidIteration::Monitor(o,in,g,s,n,c,sm,gm,b,z,inst);
   ramp239::record("after_monitor",c[z],s[z][inst][0][FLOW_SOL]);++ramp239::monitors;
   if(ramp239::control)return false; // 此控制没有流动残差，不让初始化默认值充当收敛证据。
   if(native&&!branch236::stop)branch236::stopReason="native-convergence-stop-unqualified";
   return native||branch236::stop;
 }
};

class RampDriverV239 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void launch(const std::string& expected,const std::string& facefile){using namespace observe;
   cfg=config_container[0];geo=geometry_container[0][0][0];sol=solver_container[0][0][0];auto*** num=numerics_container[0][0][0];
   ramp239::check_identity(cfg);
   if(cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||cfg->GetKind_Upwind_Flow()!=UPWIND::HLLC||cfg->GetKind_Solver()!=MAIN_SOLVER::RANS||
      cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!cfg->GetSAParsedOptions().ft2||
      cfg->GetnMGLevels()!=0||cfg->GetTime_Domain()||cfg->GetMultizone_Problem()||omp_get_max_threads()!=1||
      !dynamic_cast<CSU2PhaseAwareAir*>(sol[FLOW_SOL]->GetFluidModel())||sol[FLOW_SOL]->GetnVar()!=5||sol[TURB_SOL]->GetnVar()!=1||
      !dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_TERM])||!dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_BOUND_TERM]))throw std::runtime_error("RAMP239_NATIVE_IDENTITY");
   ramp239::trace.open(folder+"/ramp.csv");ramp239::trace<<std::setprecision(17)<<"update,stage,ramp,newtonRelaxation,nativeSyntheticProjection\n";
   delete iteration_container[0][0];iteration_container[0][0]=new RampIterationV239(cfg);
   if(ramp239::control){
     if(geo->GetnPoint()!=24||cfg->GetnInner_Iter()!=103)throw std::runtime_error("RAMP239_CONTROL_IDENTITY");
     std::vector<double> initial;std::ofstream f(folder+"/state_pairs.csv");f<<std::setprecision(17)<<"point,component,before,after\n";
     for(unsigned i=0;i<24;++i)for(unsigned j=0;j<6;++j)initial.push_back(j<5?sol[FLOW_SOL]->GetNodes()->GetSolution(i,j):sol[TURB_SOL]->GetNodes()->GetSolution(i,0));
     StartSolver();unsigned mismatch=0;
     for(unsigned i=0;i<24;++i)for(unsigned j=0;j<6;++j){double value=j<5?sol[FLOW_SOL]->GetNodes()->GetSolution(i,j):sol[TURB_SOL]->GetNodes()->GetSolution(i,0);const auto old=initial[i*6+j];mismatch+=std::memcmp(&old,&value,8)!=0;f<<i<<','<<j<<','<<old<<','<<value<<'\n';}
     if(mismatch||ramp239::calls!=103||ramp239::monitors!=103||rah_phase_step_v237::checked!=0)throw std::runtime_error("RAMP239_READONLY_CONTROL");
     std::ofstream report(folder+"/control.json");report<<"{\"passed\":true,\"schedulerCalls\":103,\"flowUpdates\":0,\"SAUpdates\":0,\"sixStateBitMismatches\":0,\"phaseGuardCalls\":0,\"physicalQualification\":false}\n";return;
   }
   if(cfg->GetnInner_Iter()!=400)throw std::runtime_error("RAMP239_UPDATE_BUDGET");
   phase="uniform-initial-check";uniform_initial(expected,facefile);events238::open(folder);
   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");steps.open(folder+"/steps.csv");walls.open(folder+"/wall_stages.csv");
   branch236::trends.open(folder+"/surface_trend.csv");branch236::decisions.open(folder+"/benefit_decisions.csv");
   for(auto* f:{&ledger,&events,&steps,&walls,&branch236::trends,&branch236::decisions})*f<<std::setprecision(17);
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   events<<"update,phase,localID,globalID,recovered,reason";for(auto word:{"trial","old","returned"})for(int j=0;j<5;++j)events<<','<<word<<j;events<<",SAcurrent,SAold\n";
   steps<<"update,lastFlowRejectionCount,totalObservedRejections,maxWallSpeedComponent,maxWallSA,minT,maxT,maxFlowVsUpdatedSAMuTDifference\n";
   branch236::trends<<"update,pressureAreaL1Change,heatAreaL1Change\n";
   branch236::decisions<<"check,windowEnd,ratio0,ratio1,ratio2,ratio3,ratio4,ratio5,geometricMeanRatio,oldPressureChange,oldHeatChange,newPressureChange,newHeatChange,noSignificantDeterioration,benefit\n";
   std::ofstream schema(folder+"/schema.json");schema<<"{\"format\":\"RAH235LE\",\"columns\":\""<<columns<<",reserved0\",\"internalEnergyOffsetJkg\":300000,\"physicalQualification\":false,\"branch\":\"phase-air-hllc-native-muscl-ramp-v239\"}\n";schema.close();
   rah_single_observer::stage=ramp239::stage;rah_single_observer::rejection=rejection;enabled=true;interval=100;
   state_frame(0);branch236::wall_frame(0);StartSolver();if(saved!=completed)state_frame(completed);if(completed%50)branch236::wall_frame(completed);
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);
   if(!forced||!events238::pending.empty())throw std::runtime_error("RAMP239_FINAL_SAVE_OR_ACCEPTANCE");
   std::ofstream accept(folder+"/acceptance_summary.json");accept<<"{\"checked\":"<<rah_phase_step_v237::checked<<",\"limited\":"<<rah_phase_step_v237::limited<<",\"forcedFinalNativeOutput\":true}\n";
   rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;
   std::ofstream f(folder+"/native_completion.json");f<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"stopReason\":\""<<branch236::stopReason<<"\",\"initialSaved\":true,\"uniformFreestreamInitial\":true,\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{
 if(argc!=5&&argc!=4)throw std::runtime_error("RAMP239_ARGUMENTS");
 ramp239::control=argc==4&&std::string(argv[3])=="--schedule-control";
 if(argc==4&&!ramp239::control)throw std::runtime_error("RAMP239_MODE");
 observe::folder=argv[2];std::set_terminate(observe::failed);SU2_MPI::Init(&argc,&argv);
 RampDriverV239 driver(argv[1],1,SU2_MPI::GetComm());driver.launch(ramp239::control?"":argv[3],ramp239::control?"":argv[4]);driver.Finalize();SU2_MPI::Finalize();return 0;
 }catch(const std::exception& e){std::ofstream f(observe::folder+"/caught_failure.txt");f<<"phase="<<observe::phase<<"\ncompleted="<<observe::completed<<"\nwhat="<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}}
