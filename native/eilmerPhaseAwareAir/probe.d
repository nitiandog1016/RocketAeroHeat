// 模块职责：经真实factory/GasModel调用测试状态，不包含CFD或替代物性公式。
module rah_phase_air_probe;
import std.stdio;
import std.file : readText;
import std.string;
import std.conv : to;
import nm.number;
import ntypes.complex;
import gas.gas_model;
import gas.gas_state;
import gas.init_gas_model;

number scalar(double x, double tangent) {
    version(complex_numbers) return complex(x,tangent);
    else { if (tangent!=0.0) throw new Exception("real helper received tangent"); return x; }
}
double imag(number x) { version(complex_numbers) return x.im; else return 0.0; }

int main(string[] args) {
    if (args.length!=4) { stderr.writeln("gas.lua cases.tsv invalid-gas.lua required"); return 2; }
    GasModel gm=init_gas_model(args[1]);
    if (gm.type_str!="rah-phase-aware-dry-air-eilmer-v1") throw new Exception("factory identity mismatch");
    writeln("MODEL ",gm.type_str);
    foreach (line; readText(args[2]).splitLines()) {
        if (line.length==0 || line[0]=='#') continue;
        auto parts=line.split(); string id=parts[0], mode=parts[1];
        try {
            if (mode=="BAD_MODEL") { init_gas_model(args[3]); throw new Exception("bad model accepted"); }
            if (parts.length!=6) throw new Exception("invalid query row");
            number x=scalar(to!double(parts[2]),to!double(parts[4]));
            number y=scalar(to!double(parts[3]),to!double(parts[5]));
            GasState q=GasState(gm); q.massf[0]=to!number(1.0);
            q.T=to!number(180.0); q.p=to!number(1.0);
            if (mode.startsWith("RT_")) {
                q.p=x; q.T=y; gm.update_thermo_from_pT(q);
                auto targetH=gm.enthalpy(q), targetS=gm.entropy(q);
                // 反算故意使用不同种子；原方法不依赖种子，不能碰巧沿用目标T。
                if (mode=="RT_RU") { q.T=to!number(900.0); gm.update_thermo_from_rhou(q); }
                else if (mode=="RT_RP") { q.T=to!number(900.0); gm.update_thermo_from_rhop(q); }
                else if (mode=="RT_RHOT") gm.update_thermo_from_rhoT(q);
                else if (mode=="RT_PU") { q.T=to!number(900.0); gm.update_thermo_from_pu(q); }
                else if (mode=="RT_PS") { q.T=to!number(900.0); gm.update_thermo_from_ps(q,targetS); }
                else if (mode=="RT_HS") { q.T=to!number(900.0); gm.update_thermo_from_hs(q,targetH,targetS); }
                else throw new Exception("unknown roundtrip");
            } else if (mode=="PT") { q.p=x; q.T=y; gm.update_thermo_from_pT(q); }
            else if (mode=="RU") { q.rho=x; q.u=y; gm.update_thermo_from_rhou(q); }
            else if (mode=="RP") { q.rho=x; q.p=y; gm.update_thermo_from_rhop(q); }
            else if (mode=="RHOT") { q.rho=x; q.T=y; gm.update_thermo_from_rhoT(q); }
            else if (mode=="PU") { q.p=x; q.u=y; gm.update_thermo_from_pu(q); }
            else if (mode=="PS") { q.p=x; gm.update_thermo_from_ps(q,y); }
            else if (mode=="HS") { gm.update_thermo_from_hs(q,x,y); }
            else if (mode=="MIN") gm.minimum_mixture_energy(q);
            else throw new Exception("unknown query mode");
            gm.update_sound_speed(q); gm.update_trans_coeffs(q);
            number[] values=[q.p,q.T,q.rho,q.u,gm.Cp(q),gm.Cv(q),gm.enthalpy(q),gm.entropy(q),q.a,q.mu,q.k];
            write(id," OK"); foreach (v; values) writef(" %.17g %.17g",v.re,imag(v)); writeln();
        } catch (GasModelException err) { writeln(id," REJECT ",err.msg.replace("\n"," ")); }
          catch (Exception err) { writeln(id," INTERNAL ",err.msg.replace("\n"," ")); }
    }
    return 0;
}
