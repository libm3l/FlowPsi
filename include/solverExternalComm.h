#ifndef SOLVEREXTERNALCOMM_H
#define SOLVEREXTERNALCOMM_H

#include <Loci.h>


namespace flowPsi {
  void get_interface_info(fact_db &facts);

  void find_mind_noslip(fact_db &facts) ;
  void find_mind_surf(fact_db &facts) ;
  void find_mind_surf_node(fact_db &facts) ;

  void create_cell_stencil(fact_db &facts) ;

  bool check_boundary_conditions(fact_db &facts) ;
  bool check_scalar_units(const options_list &o, std::string option,
                          std::string unit) ;
  bool check_vector_units(const options_list &ol,std::string vname,
                          std::string units) ;
}

#endif
