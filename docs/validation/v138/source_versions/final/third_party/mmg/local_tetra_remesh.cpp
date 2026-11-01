// SPDX-License-Identifier: LGPL-3.0-or-later
// 将已有带孔四面体域传入MMG；全部边界和未选单元冻结，原生输出另行审核。
#include <mmg/mmg3d/libmmg3d.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

void ok(int v) { if(v!=1) throw std::runtime_error("MMG API operation failed"); }
void save(MMG5_pMesh mesh,const std::string& path) {
  MMG5_int nv,nt,np,nf,nq,ne;
  ok(MMG3D_Get_meshSize(mesh,&nv,&nt,&np,&nf,&nq,&ne));
  if(np || nq) throw std::runtime_error("unexpected non-tetrahedral element");
  std::ofstream out(path);if(!out)throw std::runtime_error("output open failed");
  out<<std::setprecision(17)<<nv<<' '<<nt<<' '<<nf<<'\n';
  for(MMG5_int i=0;i<nv;++i) {
    double x,y,z;MMG5_int ref;int corner,req;
    ok(MMG3D_Get_vertex(mesh,&x,&y,&z,&ref,&corner,&req));
    out<<x<<' '<<y<<' '<<z<<' '<<ref<<' '<<req<<'\n';
  }
  for(MMG5_int i=0;i<nt;++i) {
    MMG5_int a,b,c,d,ref;int req;ok(MMG3D_Get_tetrahedron(mesh,&a,&b,&c,&d,&ref,&req));
    out<<a<<' '<<b<<' '<<c<<' '<<d<<' '<<req<<'\n';
  }
  for(MMG5_int i=0;i<nf;++i) {
    MMG5_int a,b,c,ref;int req;ok(MMG3D_Get_triangle(mesh,&a,&b,&c,&ref,&req));
    out<<a<<' '<<b<<' '<<c<<' '<<ref<<' '<<req<<'\n';
  }
}

int main(int argc,char** argv) {
  MMG5_pMesh mesh=nullptr;MMG5_pSol met=nullptr;
  try {
    if(argc!=3)throw std::runtime_error("input.dat output_prefix required");
    std::ifstream in(argv[1]);MMG5_int nv,nt,nf;
    if(!(in>>nv>>nt>>nf)||nv<4||nt<1||nf<4)throw std::runtime_error("invalid input header");
    ok(MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,MMG5_ARG_end));
    ok(MMG3D_Set_meshSize(mesh,nv,nt,0,nf,0,0));
    for(MMG5_int i=1;i<=nv;++i) {
      double x,y,z;MMG5_int ref;int req;
      if(!(in>>x>>y>>z>>ref>>req))throw std::runtime_error("truncated vertices");
      ok(MMG3D_Set_vertex(mesh,x,y,z,ref,i));if(req)ok(MMG3D_Set_requiredVertex(mesh,i));
    }
    for(MMG5_int i=1;i<=nt;++i) {
      MMG5_int a,b,c,d;int req;
      if(!(in>>a>>b>>c>>d>>req))throw std::runtime_error("truncated tetrahedra");
      ok(MMG3D_Set_tetrahedron(mesh,a,b,c,d,1,i));if(req)ok(MMG3D_Set_requiredTetrahedron(mesh,i));
    }
    for(MMG5_int i=1;i<=nf;++i) {
      MMG5_int a,b,c;
      if(!(in>>a>>b>>c))throw std::runtime_error("truncated boundary");
      ok(MMG3D_Set_triangle(mesh,a,b,c,1,i));ok(MMG3D_Set_requiredTriangle(mesh,i));
    }
    std::string extra;if(in>>extra)throw std::runtime_error("extra input data");
    ok(MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,1));
    ok(MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_mem,4096));
    ok(MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nosurf,1));
    // 本运行时未编入SCOTCH；其条件编译下的renum参数不存在，不能调用该分支。
    ok(MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_angle,0));
    ok(MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nofem,1));
    ok(MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_optim,1));
    save(mesh,std::string(argv[2])+"_import.dat");
    int result=MMG3D_mmg3dlib(mesh,met);
    std::cout<<"MMG_NATIVE_STATUS "<<result<<'\n';
    if(result!=MMG5_STRONGFAILURE)save(mesh,std::string(argv[2])+"_result.dat");
    MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,MMG5_ARG_end);
    return result==MMG5_SUCCESS?0:result;
  }catch(const std::exception& e) {
    std::cerr<<e.what()<<'\n';
    if(mesh)MMG3D_Free_all(MMG5_ARG_start,MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,MMG5_ARG_end);
    return 3;
  }
}
