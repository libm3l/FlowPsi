#ifndef GRID_COMPONENT_H
#define GRID_COMPONENT_H
#include "Loci.h"
#include "flowTypes.h"
#include <vector>
namespace flowPsi {
  
/*
 * defines prescribed gust variables
 */
  struct gustSplines {
    std::vector<real> xg,vg;
    void gust_initialize(std::vector<real> xi,std::vector<real> vi) {
      xg = xi ;
      vg = vi ;
    }
  } ;
}
#endif
