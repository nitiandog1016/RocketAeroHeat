bool CNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {
  // 模块一。使用本次实际SST k计算热状态；合法候选不改任何守恒量。
  FlowPrimitiveK(iPoint)=turb_ke;
  bool RightVol=true;
  su2double density=0, staticEnergy=0;
  auto query=[&]() {
    SetVelocity(iPoint);
    density=GetDensity(iPoint);
    staticEnergy=GetEnergy(iPoint)-0.5*Velocity2(iPoint)-turb_ke;
    FluidModel->SetTDState_rhoe(density,staticEnergy);
  };
  try {
    query();
  } catch(const std::domain_error& error) {
    if(dynamic_cast<CSU2IdealAir*>(FluidModel)==nullptr) throw;
    const std::string reason=error.what();
    const double trial_energy=staticEnergy;
    rah_step::State trial{},old{},accepted{};
    for(unsigned long j=0;j<nVar;++j) {
      const unsigned long column=(j==nVar-1)?4:j;
      trial[column]=Solution(iPoint,j);
      old[column]=Solution_Old(iPoint,j);
      Solution(iPoint,j)=Solution_Old(iPoint,j);
    }
    // 模块二。旧态也必须在当前k下合法；不得裁温、清k或使用陈旧缓存。
    try {
      query();
    } catch(const std::domain_error& old_error) {
      rah_step::record(iPoint,nDim,trial,old,old,trial_energy,staticEnergy,staticEnergy,
                       turb_ke,0,0,false,old_error.what());
      throw;
    }
    const double old_energy=staticEnergy;
    double theta=1.0;
    unsigned attempts=0;
    bool found=false;
    // 模块三。整组守恒分量共用一个有限半步，真实EOS同时审核上下温界。
    for(unsigned attempt=1;attempt<=32;++attempt) {
      theta*=0.5;
      attempts=attempt;
      for(unsigned long j=0;j<nVar;++j) {
        const unsigned long column=(j==nVar-1)?4:j;
        Solution(iPoint,j)=old[column]+theta*(trial[column]-old[column]);
      }
      try {
        query();
        found=true;
        break;
      } catch(const std::domain_error&) {}
    }
    if(!found) {
      theta=0;
      for(unsigned long j=0;j<nVar;++j) Solution(iPoint,j)=Solution_Old(iPoint,j);
      query();
    }
    for(unsigned long j=0;j<nVar;++j) accepted[(j==nVar-1)?4:j]=Solution(iPoint,j);
    rah_step::record(iPoint,nDim,trial,old,accepted,trial_energy,old_energy,staticEnergy,
                     turb_ke,theta,attempts,true,reason.c_str());
    RightVol=false;
  }

  // 模块四。保留其他模型原有的非物理状态检查与恢复路径。
  bool check_dens=SetDensity(iPoint);
  bool check_press=SetPressure(iPoint,FluidModel->GetPressure());
  bool check_sos=SetSoundSpeed(iPoint,FluidModel->GetSoundSpeed2());
  bool check_temp=SetTemperature(iPoint,FluidModel->GetTemperature());
  if(check_dens || check_press || check_sos || check_temp) {
    for(unsigned long j=0;j<nVar;++j) Solution(iPoint,j)=Solution_Old(iPoint,j);
    query();
    SetDensity(iPoint);
    SetPressure(iPoint,FluidModel->GetPressure());
    SetSoundSpeed(iPoint,FluidModel->GetSoundSpeed2());
    SetTemperature(iPoint,FluidModel->GetTemperature());
    RightVol=false;
  }

