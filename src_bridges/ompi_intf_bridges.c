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

//#include <Loci.h>
//#include "flowTypes.h"
//#include <string>
#include <stdio.h>
#include <omp.h>
#include "src_bridges.h"


//int ompi_intf_bridges(fact_db &facts){
int ompi_intf_bridges(){
/*
 * example of function which opens OpenMPI processses 
 * and access iterface data in paralell
 */

/*
 * get access to current fact database
 */
      int ithread, tid;
      int threads = 3;
      int nProcessors,nthreads;
      char* bname;

      ithread = 1;
      nProcessors = 3;
      omp_set_num_threads(nProcessors);  
 
#pragma omp parallel private(tid, bname)
{
 
#pragma omp critical
{   
       ithread++;
}   
      tid = omp_get_thread_num(); 
      if (omp_get_thread_num() == 0){
         nthreads = omp_get_num_threads();
         printf("Number of threads = %d\n", nthreads);
      }

      printf(" TID is %d \n", tid);
/*
 * close MPI processes
 */
# pragma omp barrier            
}
//omp_set_num_threads(1);   

}
