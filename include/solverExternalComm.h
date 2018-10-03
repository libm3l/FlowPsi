#ifndef SOLVEREXTERNALCOMM_H
#define SOLVEREXTERNALCOMM_H

#include <Loci.h>

namespace flowPsi {
  bool get_interface_info(fact_db &facts);
}


typedef struct InterfStr{
  int INTF_comm;
  int INTF_Root_ind;
  int INTF_Number;
  int INTF_Activated;
  int INTF_domsize;
  Array<char,24> INTF_Name;
}InterfStr_t;

 
namespace Loci {
  template <> struct data_schema_traits <InterfStr_t > {
    typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        CompoundDatatypeP cmpd = CompoundFactory(InterfStr_t());
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_comm);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_Root_ind);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_Number);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_Activated);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_Name);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_domsize);
        return DatatypeP(cmpd);
      }
    };
}



#endif
