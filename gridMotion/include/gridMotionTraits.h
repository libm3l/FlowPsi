//#############################################################################
//#
//# Copyright 2017, Mississippi State University
//#
//# This file is part of the flowPsi computational fluid dynamics solver.
//#
//# The flowPsi solver is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The flowPsi solver is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with the flowPsi solver.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef GRIDMOTIONTRAITS_H
#define GRIDMOTIONTRAITS_H

#include "gridMotion/gridTypes.h"
#include "gridMotion/treeApprox.h"
#include "gridMotion/rotor.h"
#include <Loci>

namespace Loci { 

  // This hack seems to be necessary due to template instantiation problems 
  // on some compilers
  template<> inline std::ostream &
    operator<<(std::ostream &s, const std::pair<gridMotion::vect3d,gridMotion::vect3d> &v) {
    s<<"["<<v.first<<","<<v.second<<"]";
    return s;
  }

  template<> inline std::istream &
    operator>>(std::istream &s, std::pair<gridMotion::vect3d,gridMotion::vect3d> &i) {
    char ch ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!='[') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a '[' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.first ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=',') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ',' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    s >> i.second ;
    
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=']') {
      std::cerr << "Incorrect format when reading interval" << std::endl ;
      std::cerr << "expected a ']' but got a '" << ch << "'" << std::endl ;
      s.putback(ch) ;
      return s ;
    }
    return s;
  }

  template<> struct data_schema_traits<gridMotion::Quaternion> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(gridMotion::Quaternion()) ;
      LOCI_INSERT_TYPE(ct,gridMotion::Quaternion,x) ;
      LOCI_INSERT_TYPE(ct,gridMotion::Quaternion,y) ;
      LOCI_INSERT_TYPE(ct,gridMotion::Quaternion,z) ;
      LOCI_INSERT_TYPE(ct,gridMotion::Quaternion,w) ;
      return DatatypeP(ct) ;
    }
  } ;

  template <>

    struct data_schema_traits<gridMotion::tree_info> {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(gridMotion::tree_info()) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, left) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, right) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, start) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, end) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, coord) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, split) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, q) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, centroid) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, weight) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, displacement) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, rotation) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, radius) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, err) ;
      LOCI_INSERT_TYPE(ct, gridMotion::tree_info, drot) ;
      return DatatypeP(ct) ;
    }
  };
  template<>
    struct data_schema_traits<gridMotion::NodeData> {
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(gridMotion::NodeData());
      LOCI_INSERT_TYPE(ct, gridMotion::NodeData, pos);
      LOCI_INSERT_TYPE(ct, gridMotion::NodeData, disp);
      LOCI_INSERT_TYPE(ct, gridMotion::NodeData, rot);
      LOCI_INSERT_TYPE(ct, gridMotion::NodeData, weight);
      return DatatypeP(ct);
    }
  };

  template<>
  struct data_schema_traits<gridMotion::Rotor> {
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(gridMotion::Rotor());
      LOCI_INSERT_TYPE(ct, gridMotion::Rotor, alpha);
      LOCI_INSERT_TYPE(ct, gridMotion::Rotor, beta);
      return DatatypeP(ct);
    }
  };
  template<>
  struct data_schema_traits<gridMotion::AffineMatrix> {
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      CompoundDatatypeP ct = CompoundFactory(gridMotion::AffineMatrix());
      LOCI_INSERT_TYPE(ct, gridMotion::AffineMatrix, M);
      return DatatypeP(ct);
    }
  };

}

#endif
