#ifndef SOLVEREXTERNALCOMM_H
#define SOLVEREXTERNALCOMM_H

#include <Loci.h>
#include "libm3l.h"
#include "lsipdx.h"

namespace flowPsi {
  bool get_interface_info(fact_db &facts);
}

// node_t *Gnode=NULL;

typedef struct InterfStr{
  int INTF_comm;        /*   communicator */
  int INTF_Root_ind;    /* =0 interface data not on root partition, =1 interface data on root partition */
  int INTF_Number;      /* number of interface in order as they are processed, start from 1 */
  int INTF_Activated;   /* =0 interface not active   =1 interface active */
  int INTF_domsize;     /* size of interface data set on partition */
  int INTF_totdomsize;     /* total size of interface data*/
  int INTF_type     ;     /* integer defining interface type */
  Array<char,24> INTF_Name;  /* interface name */
  bool INTF_docommunicate;
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
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_totdomsize);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_docommunicate);
        LOCI_INSERT_TYPE(cmpd, InterfStr_t, INTF_type);
        return DatatypeP(cmpd);
      }
    };
}



#endif
