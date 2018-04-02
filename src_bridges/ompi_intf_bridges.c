//#############################################################################
//#
//# Copyright 2015, Mississippi State University
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


/*
 * 
 * Description: This function calculates forces/moments
 *    It uses some of the functions and types from IOintegrated.loci
 *    It is executed conditionally based on comm_freq
 *
 *
 * History:
 * Version   Author:               Date       Patch number  CLA     Comment
 * -------   -------               --------   --------      ---     -------
 * 1.1       Adam Jirasek    2018-03-21                         Initial implementation
 *
 *
 *
 * 
 */

#include <Loci.h>
#include <string>
#include <stdio.h>
#include <omp.h>
#include "src_bridges.h"

using std::list ;
using std::vector ;
using std::endl ;
using std::cerr ;


int ompi_intf_bridges(fact_db &facts){
//int ompi_intf_bridges(){
/*
 * example of function which opens OpenMPI processses 
 * and access iterface data in paralell
 */

/*
 * get access to current fact database
 */
      int ithread, tid;
      int nProcessors,nthreads;
/*
 * get poonter on interfaces in facts 
 */
      param<options_list> int_info ;
      int_info = facts.get_variable("interfaces") ;
      options_list::option_namelist nl = int_info->getOptionNameList() ;
      options_list::option_namelist::iterator li;
/*
 * set li to the beginning of list of interfaces
 */
      li = nl.begin();
      string bname;

      ithread = 1;
      nProcessors = 3;
/*
 * set number of processors to required number 
 * which is equal to number of interfaces
 * each process will then handle its own interface
 */
      omp_set_num_threads(nProcessors);  
/*
 * set li to the beginning of interface list
 */
      li=nl.begin();
/*
 * spawn openmpi processes
 */
#pragma omp parallel private(tid, bname)
{
 
#pragma omp critical
{   
/*
 * if not the first process, increment li and set bname
 * this is a critical section, protect it
 */
       if(omp_get_thread_num() > 0)
         li++;

       bname = *li ;
}   

      Loci::option_value_type vt = int_info->getOptionValueType(bname);
      Loci::option_values ov = int_info->getOption(bname) ;
      options_list::arg_list value_list ;
      string name ;
      std::cout << "-------------  OPENMPI :::      " << omp_get_thread_num() << "  " << bname << endl ;

/*
 * get options in intf
 */


//#pragma omp critical
//{
      param<options_list> bc_info;
      switch(vt) {
      case Loci::NAME :
        ov.get_value(name) ;
        bc_info->setOption(bname,name) ;
        break ;
      case Loci::FUNCTION:
        ov.get_value(name) ;
        ov.get_value(value_list) ;
        bc_info->setOption(bname,name,value_list) ;
  //      std::cout << "-------------  OPENMPI IN FUNCTION  " <<  omp_get_thread_num() << "  " << value_list << endl ;
        break ;
      default:
        cerr << "setup_interface can not interpret value assigned to " << bname 
             << " in interfaces" << endl ;
        exit(-1) ;
      }
/*
 * parse options
 */
      options_list ol ;
      ol.Input(value_list) ;
/*
 * at the moment bcv is going to be intfr_local or intrf_global
 */
      Loci::variable bcv(name) ;
/*
 * loop through options and add them to bvars
 */
      options_list::option_namelist nlb = ol.getOptionNameList() ;
      Loci::variableSet bvars ;

      options_list::option_namelist::iterator lii;
      for(lii=nlb.begin();lii!=nlb.end();++lii){

        bvars += Loci::variable(*lii) ;
   //   std::cout << "-------------  OPENMPI   OPTION LIST  " <<  omp_get_thread_num() << "  "  << bvars  << endl ;

        if( *lii == "boundary_conditions" ){
            std::cout << "-------------  OPENMPI BOUNDARY CONDITIONS FOUND ------------------- "<< endl ;
            Loci::option_values oss = ol.getOption(*lii) ;
            std::cout << "-------------  OPENMPI BC list in intf is  " <<  omp_get_thread_num() << "  " << oss << endl ;
        }
      }
/*
 * get value of comm_freq from interface
 */
      int comm_freq = 0 ;
      double param = 0;
      ol.getOptionUnits("comm_freq","",param);
      comm_freq = (int)param;
      std::cout << "-------------  OPENMPI value of comm_freq is  " <<  omp_get_thread_num() << "  " << comm_freq  << endl ;
//}
/*
 * close MPI processes
 */
# pragma omp barrier            
}
//omp_set_num_threads(1);   

}
