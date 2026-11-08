# 相态有界标准Air的隔离Eilmer接口

`gas/rah_phase_aware_air.d`是独立GasModel类，具名`rah-phase-aware-dry-air-eilmer-v1`。它嵌入项目已有`data/thermo/phase_aware_dry_air_v1.json`并核对SHA，复用60—2000K/相态压力包络、Cp积分、hs298.15=0及u=h−RT。输运采用v207相同节点的线性mu/k插值。

模块职责：数据绑定与积分；联合域接受；原生Brent实根与复数隐函数切向；EOS/声速/热力学/输运。温压能量越界抛GasModelException，合法负参考内能不会被拒绝。Cp/输运节点使用实际返回T所在区间，精确节点取右侧（2000K取最后区间）；不能把表节点称全阶光滑。

`probe.d`经隔离的真实`gas.init_gas_model`登记分派调用，不含流动方程或替代物性。`scripts/prepareEilmerPhaseAwareAir.py`只生成隔离factory覆盖，原共享源和旧模型不改。准确执行版本、失败及补行证据见v210。未知模型由上游factory抛Error而退出，应作独立进程拒绝检查，不能假定能被普通Exception批处理捕获。

编译需要现有固定Eilmer源码、LDC1.41、Lua5.4.3/EQC静态库，`-J`指向项目`data/thermo`；实数和复数不能混用ABI。v210只验证独立helper，不分发运行库，不代表实际lmr主程序已登记或区域候选准入。不能直接运行历史脚本重新构建：对应CLI交付后关闭。下一步按目标文件核对准确NK/SA资源、真实棱柱梯度及出口/截断条件。
