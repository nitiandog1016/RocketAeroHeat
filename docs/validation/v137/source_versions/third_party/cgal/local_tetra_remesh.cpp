// SPDX-License-Identifier: GPL-3.0-or-later
// 局部四面体研究入口：保留输入连接，冻结边界，仅处理指定子域。
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/version.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Tr = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;

struct Selected {
  using key_type = Tr::Cell_handle;
  using value_type = bool;
  using reference = bool;
  using category = boost::read_write_property_map_tag;
  friend bool get(const Selected&, const key_type& c) { return c->subdomain_index()==1; }
  // CGAL在局部操作中传播子域编号；选择条件由编号重新计算。
  friend void put(Selected&, const key_type&, bool) {}
};

void save(const Tr& tr, const std::string& path) {
  std::ofstream out(path);
  if(!out) throw std::runtime_error("cannot write output");
  std::map<Tr::Vertex_handle,std::size_t> ids;
  out << std::setprecision(17) << tr.number_of_vertices() << ' '
      << tr.number_of_finite_cells() << '\n';
  for(auto v:tr.finite_vertex_handles()) {
    ids.emplace(v,ids.size());
    out << v->point().x() << ' ' << v->point().y() << ' ' << v->point().z() << '\n';
  }
  for(auto c:tr.finite_cell_handles()) {
    for(int i=0;i<4;++i) out << ids.at(c->vertex(i)) << ' ';
    out << c->subdomain_index() << '\n';
  }
}

int main(int argc,char** argv) {
  try {
    if(argc!=5) throw std::runtime_error("input.mesh output_prefix target_length iterations required");
    double h=std::stod(argv[3]);int iterations=std::stoi(argv[4]);
    if(!(h>0) || iterations!=1) throw std::runtime_error("one fixed iteration and positive size required");
    std::cout << "CGAL " << CGAL_VERSION_STR << " fixed boundaries, selected subdomain 1\n" << std::flush;
    Tr tr;std::ifstream input(argv[1]);
    if(!input || !CGAL::IO::read_MEDIT(input,tr)) throw std::runtime_error("connectivity import failed");
    if(!tr.is_valid(false)) throw std::runtime_error("input violates Triangulation_3 contract");
    save(tr,std::string(argv[2])+"_import.txt");
    CGAL::tetrahedral_isotropic_remeshing(tr,h,
      CGAL::parameters::number_of_iterations(1).remesh_boundaries(false).cell_is_selected_map(Selected{}));
    if(!tr.is_valid(false)) throw std::runtime_error("output triangulation invalid");
    save(tr,std::string(argv[2])+"_result.txt");
    std::cout << "valid output " << tr.number_of_vertices() << ' ' << tr.number_of_finite_cells() << '\n';
    return 0;
  } catch(const std::exception& e) {
    std::cerr << e.what() << '\n';return 2;
  }
}
