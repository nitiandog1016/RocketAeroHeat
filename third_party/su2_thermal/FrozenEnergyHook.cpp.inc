    // 原生alpha为0时才评估恢复；诊断模式只记录，默认环境未指定时完全不执行。
    ++rah_checked;
    if (rah_enabled && localUnderRelaxation == 0.0) {
      ++rah_frozen;
      std::array<double,5> q{}, delta{};
      for (unsigned j=0;j<nVar;++j) {
        q[j]=SU2_TYPE::GetValue(nodes->GetSolution(iPoint,j));
        delta[j]=SU2_TYPE::GetValue(LinSysSol(iPoint,j));
      }
      const auto d=rah::frozenEnergyRecovery(q,delta,nDim,0.0,SU2_TYPE::GetValue(allowableRatio));
      if (d.eligible) ++rah_eligible;
      if (rah_apply && d.eligible) {localUnderRelaxation=d.relaxation;++rah_applied;}
      SU2_OMP_CRITICAL {
        static rah::FrozenEnergyLog& log=rah_frozen_log();
        log.record(config->GetInnerIter(),iPoint,nDim,0.0,SU2_TYPE::GetValue(allowableRatio),q,delta,d,
                   SU2_TYPE::GetValue(localUnderRelaxation));
      } END_SU2_OMP_CRITICAL
    }
