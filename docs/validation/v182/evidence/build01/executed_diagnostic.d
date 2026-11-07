/**
 * RocketAeroHeat：固定Eilmer同状态导数检查入口。
 * 复用官方checkjacobian的原生初始化/稀疏矩阵/Jv函数，仅在独立程序中观测。
 * 不推进时间，不修改生产内核；实际sigma、确定方向、同阶基态和恢复账分开记录。
 */
module rah.eilmer_operator_diagnostic;

import std.stdio;
import std.conv : to;
import std.math;
import std.json : JSONValue;
import std.exception : enforce;
import nm.number : number;
import nm.smla;
import ntypes.complex;
import util.json_helper;
import lmr.globalconfig;
import lmr.globaldata : localFluidBlocks;
import lmr.lmrconfig;
import lmr.newtonkrylovsolver;

double[] values(int field)
{
    // 模块一：统一按cell/守恒变量序读取，避免近零分量被总能量尺度掩盖。
    double[] result;
    foreach (blk; localFluidBlocks) foreach (cell; blk.cells) {
        foreach (ivar; 0 .. GlobalConfig.cqi.n) {
            result ~= field == -1 ? cell.U[0][ivar].re : cell.dUdt[field][ivar].re;
        }
    }
    return result;
}

void restore(double[] original)
{
    // 恢复准确守恒基态；梯度/边界由原生残差算子重算，不代填缓存。
    size_t i=0;
    foreach (blk; localFluidBlocks) foreach (cell; blk.cells) {
        foreach (ivar; 0 .. GlobalConfig.cqi.n) {
            cell.U[0][ivar]=original[i++];
            cell.U[1][ivar]=cell.U[0][ivar];
        }
        cell.decode_conserved(0,0,blk.omegaz);
    }
}

double[] directionalResidual(double[] original, double[] direction, double step)
{
    // 模块二：独立中心差分直接调用残差，与原生前向Jv实现分开。
    restore(original);
    size_t i=0;
    foreach (blk; localFluidBlocks) foreach (cell; blk.cells) {
        foreach (ivar; 0 .. GlobalConfig.cqi.n) {
            cell.U[1][ivar]=original[i]+step*direction[i];i++;
        }
        cell.decode_conserved(0,1,blk.omegaz);
    }
    evalResidual(1,0,0);
    return values(1);
}

int main(string[] args)
{
    writeln("RAH_NATIVE_OPERATOR_V182_ENTRY");
    if (args.length != 3) { writeln("usage: diagnostic snapshot output.csv");return 2; }
    int snapshot=to!int(args[1]);
    initNewtonKrylovSimulation(snapshot,1,1,"30");
    JSONValue settings=readJSONfile(lmrCfg.nkCfgFile);
    nkCfg.readValuesFromJSON(settings);
    nkPhases.length=nkCfg.numberOfPhases;
    foreach (i,ref phase; nkPhases) phase.readValuesFromJSON(settings["NewtonKrylovPhase_"~to!string(i)]);
    activePhase=nkPhases[$-1];
    enforce(localFluidBlocks.length==1,"诊断只支持一个已审核原生块");
    auto blk=localFluidBlocks[0];
    const realStep=0x1p-26;
    version(complex_numbers) {
        writeln("RAH_NUMBER_TYPE complex");
        nkCfg.useRealValuedFrechetDerivative=false;
        nkCfg.frechetDerivativePerturbation=1e-30;
        nkCfg.preconditionerPerturbation=1e-30;
    } else {
        writeln("RAH_NUMBER_TYPE real");
        nkCfg.useRealValuedFrechetDerivative=true;
        nkCfg.frechetDerivativePerturbation=realStep;
        nkCfg.preconditionerPerturbation=realStep;
    }
    allocateGlobalGMRESWorkspace();
    blk.allocate_GMRES_workspace(1,true);
    auto original=values(-1);
    double[] direction;direction.length=original.length;
    double sumMagnitude=0, normDirection=0;
    foreach (i,u;original) {
        direction[i]=fmax(fabs(u),1e-3)*sin((i+1)*0.7);
        sumMagnitude+=fabs(u)+1;
        normDirection+=direction[i]*direction[i];
    }
    normDirection=sqrt(normDirection);
    auto centralStep=pow(2.0,-52.0/3.0)*sumMagnitude/(original.length*fmax(normDirection,sqrt(1e-7)));

    // 模块三：先观察生产二阶残差/一阶Jv的R0构造及恢复，不由相似历史猜因果。
    restore(original);evalResidual(0,0,0);setResiduals();
    auto highBefore=values(0);
    if (nkCfg.useRealValuedFrechetDerivative) setR0(0,0);
    auto highAfterR0=values(0);
    auto stateAfterR0=values(-1);
    blk.zed[]=direction[];
    double sigma;
    version(complex_numbers) { sigma=1e-30; }
    else { sigma=computePerturbationSize(); }
    evalJacobianVectorProduct(sigma,0,0);
    auto mixedJv=blk.zed.dup;

    // 模块四：矩阵和Jv都用生产Jv的一阶；固定通量/方向，形成同算子比较。
    activePhase.residualInterpolationOrder=activePhase.jacobianInterpolationOrder;
    foreach (b;localFluidBlocks) b.set_interpolation_order(activePhase.jacobianInterpolationOrder);
    restore(original);evalResidual(0,0,0);setResiduals();setR0(0,0);
    auto baseline=values(0);
    blk.initialize_jacobian(activePhase.jacobianInterpolationOrder,nkCfg.preconditionerPerturbation,0);
    blk.evaluate_jacobian(GlobalConfig.flux_calculator);
    double[] sparseJv;sparseJv.length=original.length;
    multiply(blk.flowJacobian.local,direction,sparseJv);
    auto stateAfterMatrix=values(-1);
    evalResidual(0,0,0);auto afterMatrix=values(0);

    restore(original);evalResidual(0,0,0);setResiduals();setR0(0,0);
    blk.zed[]=direction[];
    version(complex_numbers) { sigma=1e-30; }
    else { sigma=computePerturbationSize(); }
    evalJacobianVectorProduct(sigma,0,0);
    auto nativeJv=blk.zed.dup;
    auto stateAfterJv=values(-1);
    evalResidual(0,0,0);auto afterJv=values(0);

    auto plus=directionalResidual(original,direction,centralStep);
    auto minus=directionalResidual(original,direction,-centralStep);
    restore(original);evalResidual(0,0,0);auto restored=values(0);
    auto finalState=values(-1);
    auto file=File(args[2],"w");
    file.writeln("index,cell,variable,U,z,R,J_sparse,J_native,J_central,J_mixed,R_after_matrix,R_after_jv,R_restored,U_after_matrix,U_after_jv,U_after_R0,U_final,R_high_before,R_high_after_R0");
    foreach (i,u;original) {
        file.writefln("%d,%d,%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e",
            i,i/GlobalConfig.cqi.n,i%GlobalConfig.cqi.n,u,direction[i],baseline[i],sparseJv[i],nativeJv[i],
            (plus[i]-minus[i])/(2*centralStep),mixedJv[i],afterMatrix[i],afterJv[i],restored[i],
            stateAfterMatrix[i],stateAfterJv[i],stateAfterR0[i],finalState[i],highBefore[i],highAfterR0[i]);
    }
    file.close();
    writefln("RAH_SIGMA %.17e CENTRAL %.17e VARIABLES %d",sigma,centralStep,original.length);
    writeln("RAH_OPERATOR_COMPLETED");
    return 0;
}
