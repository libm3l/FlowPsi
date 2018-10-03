#ifndef SOLVEREXTERNALCOMM_H
#define SOLVEREXTERNALCOMM_H

#include <Loci.h>

namespace flowPsi {
  bool get_interface_info(fact_db &facts);
}


typedef struct InterfStr{
  int INTF_comm;
  int INTF_Root_ind;
}InterfStr_t;

 
namespace Loci {
  template <> struct data_schema_traits <InterfStr_t > {
    typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        CompoundDatatypeP cmpd = CompoundFactory(InterfStr_t());
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_comm);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_Root_ind);
        return DatatypeP(cmpd);
      }
    };
}



#endif
