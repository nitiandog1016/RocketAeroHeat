// RocketAeroHeat：隔离的标准干空气适配；不修改Eilmer流动求解内核。
// 模块职责：相容量热/输运、严格状态接受、成熟实数反算与NK复数切向。
module gas.rah_phase_aware_air;

import std.math;
import std.json;
import std.conv : to;
import std.digest.sha : sha256Of;
import std.digest : toHexString;
import std.string : toLower;
import nm.number;
import nm.brent;
import ntypes.complex;
import gas.gas_model;
import gas.gas_state;

class RAHPhaseAwareAir : GasModel {
    enum identity = "rah-phase-aware-dry-air-eilmer-v1";
    enum dataSHA = "7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985";
    enum double Rair = 287.0491267875601;
    enum double Tmin = 60.0, Tmax = 2000.0, Tref = 298.15, Pref = 100000.0;
    private double[] ts, cps, hs, ss, slopes, mus, ks;

    this() {
        // 模块一：嵌入并绑定既有数据；初始化只做积分和原输运节点计算。
        enum raw = import("phase_aware_dry_air_v1.json");
        if (toLower(toHexString(sha256Of(raw))) != dataSHA)
            throw new GasModelException("RAH_DATA_SHA");
        auto d = parseJSON(raw);
        if (d["model"].str != "rah-phase-aware-dry-air-of13-v1" ||
            d["domain"].str != "rah-lemmon-gas-low-density-60-2000-v1" ||
            d["gasConstantJkgK"].floating != Rair)
            throw new GasModelException("RAH_DATA_IDENTITY");
        size_t n = d["temperatureK"].array.length;
        if (n < 2 || d["cpJkgK"].array.length != n) throw new GasModelException("RAH_TABLE");
        ts.length=n; cps.length=n; hs.length=n; ss.length=n;
        slopes.length=n-1; mus.length=n; ks.length=n;
        size_t refIndex=n;
        foreach (i; 0..n) {
            ts[i]=d["temperatureK"][i].floating; cps[i]=d["cpJkgK"][i].floating;
            if (!isFinite(ts[i]) || !isFinite(cps[i]) || cps[i]<=Rair || (i && ts[i]<=ts[i-1]))
                throw new GasModelException("RAH_TABLE");
            if (ts[i]==Tref) refIndex=i;
            double x=log(ts[i]/103.3);
            double omega=exp(.431-.4623*x+.08406*x*x+.005341*x*x*x-.00331*x*x*x*x);
            mus[i]=2.66958e-8*sqrt(28.9586*ts[i])/(.36*.36*omega);
            double tau=132.6312/ts[i];
            ks[i]=.001308*(mus[i]*1e6)+.001405*pow(tau,-1.1)-.001036*pow(tau,-.3);
        }
        if (ts[0]!=Tmin || ts[$-1]!=Tmax || refIndex==n) throw new GasModelException("RAH_TABLE_DOMAIN");
        hs[0]=0; ss[0]=0;
        foreach (i; 0..n-1) {
            double dt=ts[i+1]-ts[i]; slopes[i]=(cps[i+1]-cps[i])/dt;
            hs[i+1]=hs[i]+(cps[i]+cps[i+1])*dt/2;
            ss[i+1]=ss[i]+(cps[i]-slopes[i]*ts[i])*log(ts[i+1]/ts[i])+slopes[i]*dt;
        }
        double h0=hs[refIndex], s0=ss[refIndex];
        foreach (i; 0..n) { hs[i]-=h0; ss[i]-=s0; }
        type_str=identity; _n_species=1; _n_modes=0;
        _species_names=["rah-standard-dry-air"]; _mol_masses=[.02896546];
        create_species_reverse_lookup();
    }

    // 模块二：最终联合域在每个公开状态入口验收，复数扰动仅用于局部切向。
    @nogc static void finite(number x) {
        if (!isFinite(x.re)) throw new GasModelException("RAH_NONFINITE");
        version(complex_numbers) { if (!isFinite(x.im)) throw new GasModelException("RAH_NONFINITE_TANGENT"); }
    }
    @nogc static double pressureCeiling(double t) {
        if (!isFinite(t) || t<Tmin || t>Tmax) throw new GasModelException("RAH_T_DOMAIN");
        double ratio=t/100.; double p=fmin(1e6,1e4*((ratio*ratio)*ratio));
        if (t<=132.6312) {
            double th=1-t/132.6312;
            double val=-.1567266*sqrt(th)-5.539635*th+.7567212*pow(th,2.5)-3.514322*pow(th,4);
            p=fmin(p,.1*(3785020.*exp(132.6312/t*val)));
        }
        return p;
    }
    @nogc static void accept(number p, number t) {
        finite(p); finite(t);
        if (p.re<=0 || p.re>pressureCeiling(t.re)) throw new GasModelException("RAH_PT_DOMAIN");
    }
    @nogc private size_t interval(number t) const {
        finite(t);
        if (t.re<Tmin || t.re>Tmax) throw new GasModelException("RAH_T_DOMAIN");
        size_t lo=0, hi=ts.length;
        while (lo<hi) { size_t mid=(lo+hi)/2; if (ts[mid]<=t.re) lo=mid+1; else hi=mid; }
        // 仅选取已接受端点的末区间；这里不裁温度或外推距离。
        return lo==ts.length ? ts.length-2 : lo-1;
    }
    @nogc number cp(number t) const { auto i=interval(t); return cps[i]+slopes[i]*(t-ts[i]); }
    @nogc number h(number t) const {
        auto i=interval(t); number dt=t-ts[i]; return hs[i]+cps[i]*dt+.5*slopes[i]*dt*dt;
    }
    @nogc number s0(number t) const {
        auto i=interval(t); return ss[i]+(cps[i]-slopes[i]*ts[i])*log(t/ts[i])+slopes[i]*(t-ts[i]);
    }
    @nogc number e(number t) const { return h(t)-Rair*t; }

    // 模块三：成熟Brent求实根；显式隐函数导数保留复数NK所需一阶切向。
    // kind: 0=h，1=u，2=相对热熵s0。实根允许负能量，禁止超表目标。
    @nogc number inverse(number target, int kind, number shift=to!number(0.0)) const {
        finite(target); finite(shift);
        if (kind<0 || kind>2) throw new GasModelException("RAH_INVERSE_KIND");
        double value(double t) { return (kind==0 ? h(to!number(t)) : kind==1 ? e(to!number(t)) : s0(to!number(t))).re-shift.re; }
        double left=value(Tmin), right=value(Tmax), v=target.re;
        if (v<left || v>right) throw new GasModelException("RAH_CALORIC_DOMAIN");
        double root;
        if (v==left) root=Tmin;
        else if (v==right) root=Tmax;
        else {
            double f(double t) { return value(t)-v; }
            try { root=nm.brent.solve!(f,double)(Tmin,Tmax,2e-11); }
            catch (Exception err) { throw new GasModelException("RAH_BRENT_FAILURE"); }
        }
        number t=to!number(root);
        version(complex_numbers) {
            double deriv=cp(t).re;
            if (kind==1) deriv-=Rair;
            else if (kind==2) deriv/=root;
            t=complex(root,(target.im+shift.im)/deriv);
        }
        return t;
    }
    @nogc private void fill(ref GasState q) const {
        accept(q.p,q.T);
        number rho=q.p/(Rair*q.T), u=e(q.T);
        finite(rho); finite(u);
        if (rho.re<=0) throw new GasModelException("RAH_DENSITY_RANGE");
        q.rho=rho; q.u=u;
    }
    @nogc override void update_thermo_from_pT(ref GasState q) const { fill(q); }
    @nogc override void update_thermo_from_pu(ref GasState q) const {
        number t=inverse(q.u,1); accept(q.p,t); q.T=t; fill(q);
    }
    @nogc override void update_thermo_from_rhou(ref GasState q) const {
        finite(q.rho); if (q.rho.re<=0) throw new GasModelException("RAH_RHO_DOMAIN");
        number t=inverse(q.u,1), p=q.rho*Rair*t; accept(p,t); q.T=t; q.p=p;
        // rho/u为实际输入，不用反算舍入后的值覆盖守恒输入。
    }
    @nogc override void update_thermo_from_rhoT(ref GasState q) const {
        finite(q.rho); finite(q.T); if (q.rho.re<=0) throw new GasModelException("RAH_RHO_DOMAIN");
        number p=q.rho*Rair*q.T; accept(p,q.T); q.p=p; q.u=e(q.T);
    }
    @nogc override void update_thermo_from_rhop(ref GasState q) const {
        finite(q.rho); finite(q.p); if (q.rho.re<=0) throw new GasModelException("RAH_RHO_DOMAIN");
        number t=q.p/(Rair*q.rho); accept(q.p,t); q.T=t; q.u=e(t);
    }
    @nogc override void update_thermo_from_ps(ref GasState q, number s) const {
        finite(q.p); finite(s); if (q.p.re<=0) throw new GasModelException("RAH_PT_DOMAIN");
        number t=inverse(s,2,Rair*(log(q.p)-log(Pref))); accept(q.p,t); q.T=t; fill(q);
    }
    @nogc override void update_thermo_from_hs(ref GasState q, number hh, number s) const {
        finite(s); number t=inverse(hh,0), p=Pref*exp((s0(t)-s)/Rair);
        accept(p,t); q.T=t; q.p=p; fill(q);
    }

    // 模块四：量热/声速/输运共用同一状态与同一节点的右侧导数。
    @nogc override void update_sound_speed(ref GasState q) const {
        accept(q.p,q.T); number c=cp(q.T); q.a=sqrt(c/(c-Rair)*Rair*q.T); finite(q.a);
    }
    @nogc override void update_trans_coeffs(ref GasState q) const {
        accept(q.p,q.T); auto i=interval(q.T); number frac=(q.T-ts[i])/(ts[i+1]-ts[i]);
        q.mu=mus[i]+frac*(mus[i+1]-mus[i]); q.k=ks[i]+frac*(ks[i+1]-ks[i]);
        finite(q.mu); finite(q.k);
        if (q.mu.re<=0 || q.k.re<=0) throw new GasModelException("RAH_TRANSPORT_RANGE");
    }
    @nogc override number dudT_const_v(in GasState q) const { accept(q.p,q.T); return cp(q.T)-Rair; }
    @nogc override number dhdT_const_p(in GasState q) const { accept(q.p,q.T); return cp(q.T); }
    @nogc override number dpdrho_const_T(in GasState q) const { accept(q.p,q.T); return Rair*q.T; }
    @nogc override number gas_constant(in GasState q) const { accept(q.p,q.T); return to!number(Rair); }
    @nogc override number internal_energy(in GasState q) const { accept(q.p,q.T); return e(q.T); }
    @nogc override number enthalpy(in GasState q) const { accept(q.p,q.T); return h(q.T); }
    @nogc override number entropy(in GasState q) const { accept(q.p,q.T); return s0(q.T)-Rair*(log(q.p)-log(Pref)); }
    @nogc override void minimum_mixture_energy(ref GasState q) const {
        // 原接口本来会重置状态；这里只返回合法参考压力下的量热下界。
        // 不是指定rho或p下的可接受性证明，真实输入仍经accept严格验证。
        q.p=to!number(.5*pressureCeiling(Tmin)); q.T=to!number(Tmin); fill(q);
    }
}
