
#include <cmath>
#include <cstdio>
#include <algorithm>
using su2double = double;
constexpr double EPS = 1e-16, TWO3 = 2.0/3.0;
struct CConfig {};
struct CNumerics {
  template<class T=double> struct ResidualType {
    const double* residual;
    ResidualType(const double* r, double**, void*) : residual(r) {}
  };
};
// 最小存储壳，仅提供逐字原函数读取的成员；不是整个 SU2 类/离散求解器。
struct CSourceAxisymmetric_Flow : CNumerics {
  unsigned short nVar=4, nDim=2;
  bool implicit=true, viscous=true, rans=false;
  double yinv=0.0, Volume=.7, turb_ke_i=0, Prandtl_Turb=.9;
  double Coord_i[2]{}, U_i[4]{}, V_i[12]{}, V_j[12]{}, S_i[2]{};
  double PrimVar_Grad_i[3][2]{}, AuxVar_Grad_i[3][2]{};
  double residual[4]{}, matrix[4][4]{};
  double* jacobian[4]{matrix[0], matrix[1], matrix[2], matrix[3]};
  void ResidualDiffusion();
};
struct CSourceGeneralAxisymmetric_Flow : CSourceAxisymmetric_Flow {
  ResidualType<> ComputeResidual(const CConfig*);
};

void CSourceAxisymmetric_Flow::ResidualDiffusion(){

  if (!rans){ turb_ke_i = 0.0; }

  su2double laminar_viscosity_i    = V_i[nDim+5];
  su2double eddy_viscosity_i       = V_i[nDim+6];
  su2double thermal_conductivity_i = V_i[nDim+7];
  su2double heat_capacity_cp_i     = V_i[nDim+8];

  su2double total_viscosity_i = laminar_viscosity_i + eddy_viscosity_i;
  su2double total_conductivity_i = thermal_conductivity_i + heat_capacity_cp_i*eddy_viscosity_i/Prandtl_Turb;

  su2double u = U_i[1]/U_i[0];
  su2double v = U_i[2]/U_i[0];

  residual[0] -= 0.0;
  residual[1] -= Volume*(yinv*total_viscosity_i*(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0])
                         -TWO3*AuxVar_Grad_i[0][0]);
  residual[2] -= Volume*(yinv*total_viscosity_i*2*(PrimVar_Grad_i[2][1]-v*yinv)
                         -TWO3*AuxVar_Grad_i[0][1]);
  residual[3] -= Volume*(yinv*(total_viscosity_i*(u*(PrimVar_Grad_i[2][0]+PrimVar_Grad_i[1][1])
                                                 +v*TWO3*(2*PrimVar_Grad_i[2][1]-PrimVar_Grad_i[1][0]
                                                 -v*yinv+U_i[0]*turb_ke_i))
                                                 +total_conductivity_i*PrimVar_Grad_i[0][1])
                                                 -TWO3*(AuxVar_Grad_i[1][1]+AuxVar_Grad_i[2][0]));
}
CNumerics::ResidualType<> CSourceGeneralAxisymmetric_Flow::ComputeResidual(const CConfig* config) {
  unsigned short iVar, jVar;

  if (Coord_i[1] > EPS) {

    yinv = 1.0/Coord_i[1];

    su2double Density_i = U_i[0];
    su2double Velocity1_i = U_i[1]/U_i[0];
    su2double Velocity2_i = U_i[2]/U_i[0];
    su2double Energy_i = U_i[3]/U_i[0];

    su2double Pressure_i = V_j[3];
    su2double Enthalpy_i = Energy_i + Pressure_i/Density_i;

    /*--- Inviscid component of the source term. ---*/

    residual[0] = yinv*Volume*U_i[2];
    residual[1] = yinv*Volume*U_i[1]*Velocity2_i;
    residual[2] = yinv*Volume*U_i[2]*Velocity2_i;
    residual[3] = yinv*Volume*U_i[2]*Enthalpy_i;

    if (implicit) {

      su2double dPdrho_e_i = S_i[0];
      su2double dPde_rho_i = S_i[1];

      jacobian[0][0] = 0.0;
      jacobian[0][1] = 0.0;
      jacobian[0][2] = 1.0;
      jacobian[0][3] = 0.0;

      jacobian[1][0] = -Velocity1_i*Velocity2_i;
      jacobian[1][1] = Velocity2_i;
      jacobian[1][2] = Velocity1_i;
      jacobian[1][3] = 0.0;

      jacobian[2][0] = -Velocity2_i*Velocity2_i;
      jacobian[2][1] = 0.0;
      jacobian[2][2] = 2*Velocity2_i;
      jacobian[2][3] = 0.0;

      jacobian[3][0] = Velocity2_i*(dPdrho_e_i + dPde_rho_i/Density_i*(Velocity1_i*Velocity1_i
                                                                        + Velocity2_i*Velocity2_i
                                                                        - Energy_i) - Enthalpy_i);
      jacobian[3][1] = -Velocity1_i*Velocity2_i/Density_i *dPde_rho_i;
      jacobian[3][2] = Enthalpy_i - Velocity2_i*Velocity2_i/Density_i *dPde_rho_i;
      jacobian[3][3] = Velocity2_i + Velocity2_i/Density_i *dPde_rho_i;

      for (iVar=0; iVar < nVar; iVar++)
        for (jVar=0; jVar < nVar; jVar++)
          jacobian[iVar][jVar] *= yinv*Volume;

    }

    /*--- Add the viscous terms if necessary. ---*/

    if (viscous) ResidualDiffusion();

  }

  else {
    /*--- At the axis of symmetry, use L'Hôpital's rule: lim(v/r) = dv/dr ---*/
    const su2double dv_dr = PrimVar_Grad_i[2][1];  // ∂v/∂r (radial velocity gradient)
    const su2double u = U_i[1]/U_i[0];            // axial velocity u
    const su2double rho = U_i[0];                 // density

    /* Compute pressure and enthalpy consistently with the general-gas formulation. */
    const su2double Density_i = rho;
    const su2double Energy_i  = U_i[3]/U_i[0];
    const su2double Pressure_i = V_j[3];
    const su2double Enthalpy_i = Energy_i + Pressure_i/Density_i;

    /*--- Apply L'Hôpital's rule to axisymmetric source terms ---*/
    residual[0] = Volume * rho * dv_dr;               // ρ(∂v/∂r)
    residual[1] = Volume * rho * u * dv_dr;           // ρu(∂v/∂r)
    residual[2] = 0.0;                                // ρv(∂v/∂r) = 0 since v=0 at axis
    residual[3] = Volume * rho * Enthalpy_i * dv_dr;  // ρH(∂v/∂r)

    if (implicit) {
      /* For now, set Jacobian to zero at axis (can be improved later to help with convergence). */
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++)
          jacobian[iVar][jVar] = 0.0;
      }
    }

    /*--- Add the viscous terms if necessary. ---*/
    if (viscous) ResidualDiffusion();

  }

  return ResidualType<>(residual, jacobian, nullptr);
}

void fill(CSourceGeneralAxisymmetric_Flow& x, double r, int mode) {
  // 光滑偶/奇场：u=10+a*r²，v=c*r，T=300+b*r²；轴线极限有解析值。
  const double a=mode==2 ? 2.0 : 0.0;
  const double b=mode==2 ? 5.0 : 0.0;
  const double c=mode==1 ? .3 : 0.0;
  const double rho=1.2, u=10+a*r*r, v=c*r, mu=2e-5, k=.03;
  const double p=100000, e=250000;
  x.Coord_i[1]=r;
  x.U_i[0]=rho; x.U_i[1]=rho*u; x.U_i[2]=rho*v;
  x.U_i[3]=rho*(e+.5*(u*u+v*v));
  x.V_j[3]=p; x.V_i[7]=mu; x.V_i[8]=0; x.V_i[9]=k; x.V_i[10]=1005;
  x.S_i[0]=p/rho; x.S_i[1]=.4*rho;
  x.PrimVar_Grad_i[0][1]=2*b*r;
  x.PrimVar_Grad_i[1][1]=2*a*r;
  x.PrimVar_Grad_i[2][1]=c;
  // 原辅助变量为 mu*v/r、mu*v²/r、mu*u*v/r。第二项径向导数不能遗漏。
  for(auto& row : x.AuxVar_Grad_i) for(auto& v : row) v=0;
  x.AuxVar_Grad_i[1][1]=mu*c*c;
}
int main() {
  std::puts("mode,viscous,preceding_radius,radius,R0,R1,R2,R3");
  for(int mode : {1,2}) for(bool viscous : {false,true}) {
    for(double before : {0.0,.1,.2}) {
      CSourceGeneralAxisymmetric_Flow x; x.viscous=viscous;
      if(before>0) {fill(x,before,mode); x.ComputeResidual(nullptr);}
      fill(x,0,mode); x.ComputeResidual(nullptr);
      std::printf("%d,%d,%.17g,0,%.17g,%.17g,%.17g,%.17g\n",mode,viscous,before,
                  x.residual[0],x.residual[1],x.residual[2],x.residual[3]);
    }
    CSourceGeneralAxisymmetric_Flow x; x.viscous=viscous;
    for(double radius : {.01,.001,.0001}) {
      fill(x,radius,mode); x.ComputeResidual(nullptr);
      std::printf("%d,%d,0,%.17g,%.17g,%.17g,%.17g,%.17g\n",mode,viscous,radius,
                  x.residual[0],x.residual[1],x.residual[2],x.residual[3]);
    }
  }
}
