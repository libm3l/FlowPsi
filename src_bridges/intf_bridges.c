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
#include "src_bridges_types.h"
#include "src_bridges.h"

using std::list ;
using std::vector ;
using std::endl ;
using std::cerr ;

/*
 * this is a function which loops over interfaces and communicates as required
 * it was used in rule in solverExternalComm.loci 
 * this is a serial communication which is not ideal as a running process
 * block all other processes
 */

int intf_bridges(fact_db &facts){
/*
 * get pointer on interfaces in facts 
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
/*
 * set li to the beginning of interface list
 */
      for(li=nl.begin();li!=nl.end();++li) {

        bname = *li ;

        Loci::option_value_type vt = int_info->getOptionValueType(bname);
        Loci::option_values ov = int_info->getOption(bname) ;
        options_list::arg_list value_list ;
        string name ;
/*
 * get options in intf
 */
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
      Loci::option_values oss;
      int portno,comm_freq;
      double param;
      string type,tag,I_channel, O_channel, intf_name,IP, BCs, I_channel_name, O_channel_name;
      
      comm_struct_t *pcomm_str, comm_str;
      pcomm_str = &comm_str;

      options_list::option_namelist::iterator lii;
      for(lii=nlb.begin();lii!=nlb.end();++lii){

        bvars += Loci::variable(*lii) ;

        if( *lii == "boundary_conditions" ){
            oss = ol.getOption(*lii) ;
            ol.getOption(*lii,BCs);
        }
        else if(*lii == "I_channel" ){
            ol.getOption(*lii,I_channel_name) ;
            pcomm_str->I_channel = I_channel_name.c_str();
        }
        else if(*lii == "O_channel" ){
	         ol.getOption(*lii,O_channel_name) ;
            pcomm_str->O_channel = O_channel_name.c_str();
            std::cout << "-------------  OPENMPI O_channel name is  " <<  omp_get_thread_num() << "  " << pcomm_str->O_channel<< endl ;
        }
        else if(*lii == "tag" ){
	        ol.getOption(*lii,tag) ;
            pcomm_str->tag = tag.c_str();
        }
        else if(*lii == "name" ){
	        ol.getOption(*lii,intf_name) ;
            pcomm_str->intf_name = intf_name.c_str();
        }
        else if(*lii == "type" ){
	        ol.getOption(*lii,type) ;
            pcomm_str->type = type.c_str();
        }
        else if(*lii == "IP" ){
            ol.getOption(*lii,IP) ;
            pcomm_str->IP = IP.c_str();
        }
      }
/*
 * get value of comm_freq from interface
 */
      ol.getOptionUnits("comm_freq","",param);
      pcomm_str->comm_freq = (int)param;
/*
 * get port number
 */
      ol.getOptionUnits("portno","",param);
      pcomm_str->portno = (int)param;
/*
 * get comm channel name
 */
      test_bridge1(pcomm_str);
  }
  return 1;
}
