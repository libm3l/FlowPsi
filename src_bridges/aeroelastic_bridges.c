 /*#############################################################################
 *     Copyright (C) 2018  Adam Jirasek
 * 
 *     This file is part of the flowPsi computational fluid dynamics solver.
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *     
 *     contact: libm3l@gmail.com
 * 
 *#############################################################################
 */



/*
 *     Description: - contains bridge function to lsipdx engine
 * 
 * History:
 * Version   Author:               Date       Patch number  CLA     Comment
 * -------   -------               --------   --------      ---     -------
 * 1-beta-6  Adam Jirasek         2018-03-21                        Initial Implementation
 *
 *
 *     Description
 * 
 */

/*
 * check that there is support for libm3l and lsipdx
 */
#include "extcomm_def.h"

#ifdef LIBM3LSIPDX

#include "libm3l.h"
#include "lsipdx.h"
#include "src_bridges_types.h"
#include "src_bridges.h"

/*
 *     Description: - this is a bridge_function for bridge_prescribed_quaternion
 *       this function creates data which are located in head node
 *       CFD_2_SIM and consists of double array of 6 called ForcesMoments and single double value called Time
 *
 *  H-DIR           CFD_2_SIM               2
 *     -D               Forces           1      6   
 *     -D               Time            1      1   
 *
 *      this data set is recevied through channel CFD2SIM
 *
 *      as a respons, it receives a data set in head node SIM_2_CFD
 *
 *  H-DIR           SIM_2_CFD               3
 *   -D                     Quaternion              1      4
 *   -D                     RotCenter               1      3   
 *   -D                     TransVec                1      3   
 *
 *      which contains double array of 4 called Quaternion, double array of 3 called RotCenter and double array of 
 *      3 called TransVec
 *
 *      to create and manipulate with this type of data set, you need to link two libraries - libm3l and lsipdx 
 *      for more details go to www.github.com/libm3l/libm3l and www.github.com/libm3l/lsipdx
 *
 *     Input parameters:
 *     t_time  - simulation time
 *     ForceX, ForceY, ForceZ - value of aerodynamic forces
 *     comm_str - structure containing details for lsipdx engine
 *
 *     Return value:
 *      Alpha, Qx, Qy, Qz - values of quaternion
 *      TransX, TransY, TransZ vector of translation
 *      RotCX,  RotCY,  RotCZ - vector of rotation center
 * 
 * History:
 * Version   Author:               Date       Patch number  CLA     Comment
 * -------   -------               --------   --------      ---     -------
 * 1-beta-6  Adam Jirasek         2018-03-21                        Initial Implementation
 *
 *
 *     Description
 * 
 */
CPP_C node_t *bridge_aeroelastic(node_t *Gnode, comm_struct_t *comm_str){

	lmint_t sockfd;
    node_t *Snode, *GlobDataNode;

	client_fce_struct_t InpPar, *PInpPar;
	opts_t *Popts_1, opts, opts_1, *Popts;
    find_t *Founds=NULL;
/*
 * Set parameters which are needed for opening the socket for sending data
 */
    PInpPar = &InpPar;
	PInpPar->channel_name = comm_str->O_channel;   /* name of channel where to send data*/
	PInpPar->SR_MODE = 'S';         /* set R or S, because we will write outgoing data to this 
                                                   socket, set it to S */
/*
 *  this parameter determines two modes.
 *  for this case it is set to D as direct, meaning we will open the socket and get the data only
 *  had we planned to send data back through the same socket, we would use A as alternate
 *
 *  the second parameter is either Y or N. It determines if we should close the socket afetr receiving data
 *  or keep it opened. In this scenario we close the socket as soon as we receive the data so set it N
 *
 *  NOTE: these nodes have to be the same for both processign accessing this channel (ie. sendign process and receviing process
 *        need to use the same settings - see ../External_Processes/Client1_FakedSimulink.c)
 *
 *  for more details on modes see:  
 *  Adam Jirasek and Arthur W. Rizzi: libm3l and lsipdx - Utilities for Inter-Process Data Transfer and Synchronization", 
 *  52nd Aerospace Sciences Meeting, AIAA SciTech Forum, (AIAA 2014-1045) https://doi.org/10.2514/6.2014-1045
 */
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("socket_Flow2Aeroel: wrong client mode");
/*
 * this is just a structure containing additional data
 */
	Popts   = &opts;
	Popts_1 = &opts_1;
/*
 *  use the set options to set parameters for connection
 */
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	m3l_set_Find(&Popts);
/*
 * find global_data set
 */
    if( (Founds = m3l_Locate(Gnode, "./Interface/global_data", "./*/*",  (lmchar_t *)NULL)) != NULL){
                    
       GlobDataNode = m3l_get_Found_node(Founds, 0);
       m3l_DestroyFound(&Founds);

    }
    else{
      Perror("Did not find glob data set");
    }
                
/*
 * print data on screen
 */
    if(m3l_Cat(GlobDataNode, "--detailed", "-P", "-L",  "*",   (char *)NULL) != 0)
 	   Error("CatData");
/*
 * open socket
 */
#pragma omp critical
{
	if( (sockfd = open_connection_to_server(comm_str->IP, comm_str->portno, PInpPar, Popts_1)) < 1)
		Error("socket_Flow2Aeroel: Error when opening socket");
}
/*
 * send data 
 */
	if ( client_sender(Gnode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL) !=1 )
		Error("socket_Flow2Aeroel: client_sender()");
/*
 * close socket
 */
#pragma omp critical
{
	if( close(sockfd) == -1)
		Perror("socket_Flow2Aeroel: close");
}
/*
 * free borrowed memory
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_Flow2Aeroel: m3l_Umount");
/*
 * receive data 
 */
	PInpPar = &InpPar;
	PInpPar->channel_name = comm_str->I_channel;       /* name of channel */
	PInpPar->SR_MODE = 'R';              /* set R or S, because we will read incoming data from this 
                                                   socket, set it to R */
/*
 *  this parameter determines two modes.
 *  for this case it is set to D as direct, meaning we will open the socket and get the data only
 *  had we planned to send data back through the same socket, we would use A as alternate
 *
 *  the second parameter is either Y or N. It determines if we should close the socket afetr receiving data
 *  or keep it opened. In this scenario we close the socket as soon as we receive the data so set it N
 *
 *  NOTE: these nodes have to be the same for both processign accessing this channel (ie. sendign process and receviing process
 *        need to use the same settings - see ../External_Processes/Client1_FakedSimulink.c)
 *
 *  for more details on modes see:  
 *  Adam Jirasek and Arthur W. Rizzi: libm3l and lsipdx - Utilities for Inter-Process Data Transfer and Synchronization", 
 *  52nd Aerospace Sciences Meeting, AIAA SciTech Forum, (AIAA 2014-1045) https://doi.org/10.2514/6.2014-1045
 */
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("socket_Flow2Aeroel: wrong client mode");

	Popts   = &opts;
	Popts_1 = &opts_1;
/*
 * set socket settings
 */
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	m3l_set_Find(&Popts);
/*
 * open socket for reading data 
 */
#pragma omp critical
{
	if( (sockfd = open_connection_to_server(comm_str->IP, comm_str->portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");
}
/*
 * receive data from socket
 */
	if ( (Snode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL)) == NULL)
		Error("socket_Flow2Aeroel: client_receiver()");
/*
 * close socket 
 */
#pragma omp critical
{
	if( close(sockfd) == -1)
		Perror("socket_Flow2Aeroel: close");
}


//=======================================================================

        return Snode;

}


CPP_C lmint_t bridge_aeroelastic_rigid(node_t *Gnode, lmdouble_t * Alpha, lmdouble_t * Qx, lmdouble_t * Qy, lmdouble_t * Qz,
                   lmdouble_t * TransX, lmdouble_t * TransY, lmdouble_t * TransZ, 
                   lmdouble_t * RotCX, lmdouble_t * RotCY, lmdouble_t * RotCZ,
                   comm_struct_t *comm_str){

	lmint_t sockfd;
    node_t *Snode, *GlobDataNode;

	client_fce_struct_t InpPar, *PInpPar;
	opts_t *Popts_1, opts, opts_1, *Popts;
    find_t *Founds=NULL;
/*
 * Set parameters which are needed for opening the socket for sending data
 */
    PInpPar = &InpPar;
	PInpPar->channel_name = comm_str->O_channel;   /* name of channel where to send data*/
	PInpPar->SR_MODE = 'S';         /* set R or S, because we will write outgoing data to this 
                                                   socket, set it to S */
/*
 *  this parameter determines two modes.
 *  for this case it is set to D as direct, meaning we will open the socket and get the data only
 *  had we planned to send data back through the same socket, we would use A as alternate
 *
 *  the second parameter is either Y or N. It determines if we should close the socket afetr receiving data
 *  or keep it opened. In this scenario we close the socket as soon as we receive the data so set it N
 *
 *  NOTE: these nodes have to be the same for both processign accessing this channel (ie. sendign process and receviing process
 *        need to use the same settings - see ../External_Processes/Client1_FakedSimulink.c)
 *
 *  for more details on modes see:  
 *  Adam Jirasek and Arthur W. Rizzi: libm3l and lsipdx - Utilities for Inter-Process Data Transfer and Synchronization", 
 *  52nd Aerospace Sciences Meeting, AIAA SciTech Forum, (AIAA 2014-1045) https://doi.org/10.2514/6.2014-1045
 */
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("socket_Flow2Aeroel: wrong client mode");
/*
 * this is just a structure containing additional data
 */
	Popts   = &opts;
	Popts_1 = &opts_1;
/*
 *  use the set options to set parameters for connection
 */
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	m3l_set_Find(&Popts);
/*
 * find global_data set
 */
    if( (Founds = m3l_Locate(Gnode, "./Interface/global_data", "./*/*",  (lmchar_t *)NULL)) != NULL){
                    
       GlobDataNode = m3l_get_Found_node(Founds, 0);
       m3l_DestroyFound(&Founds);

    }
    else{
      Perror("Did not find glob data set");
    }
                
/*
 * print data on screen
 */
    if(m3l_Cat(GlobDataNode, "--detailed", "-P", "-L",  "*",   (char *)NULL) != 0)
 	   Error("CatData");
/*
 * open socket
 */
#pragma omp critical
{
	if( (sockfd = open_connection_to_server(comm_str->IP, comm_str->portno, PInpPar, Popts_1)) < 1)
		Error("socket_Flow2Aeroel: Error when opening socket");
}
/*
 * send data 
 */
	if ( client_sender(Gnode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL) !=1 )
		Error("socket_Flow2Aeroel: client_sender()");
/*
 * close socket
 */
#pragma omp critical
{
	if( close(sockfd) == -1)
		Perror("socket_Flow2Aeroel: close");
}
/*
 * free borrowed memory
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_Flow2Aeroel: m3l_Umount");
/*
 * receive data 
 */
	PInpPar = &InpPar;
	PInpPar->channel_name = comm_str->I_channel;       /* name of channel */
	PInpPar->SR_MODE = 'R';              /* set R or S, because we will read incoming data from this 
                                                   socket, set it to R */
/*
 *  this parameter determines two modes.
 *  for this case it is set to D as direct, meaning we will open the socket and get the data only
 *  had we planned to send data back through the same socket, we would use A as alternate
 *
 *  the second parameter is either Y or N. It determines if we should close the socket afetr receiving data
 *  or keep it opened. In this scenario we close the socket as soon as we receive the data so set it N
 *
 *  NOTE: these nodes have to be the same for both processign accessing this channel (ie. sendign process and receviing process
 *        need to use the same settings - see ../External_Processes/Client1_FakedSimulink.c)
 *
 *  for more details on modes see:  
 *  Adam Jirasek and Arthur W. Rizzi: libm3l and lsipdx - Utilities for Inter-Process Data Transfer and Synchronization", 
 *  52nd Aerospace Sciences Meeting, AIAA SciTech Forum, (AIAA 2014-1045) https://doi.org/10.2514/6.2014-1045
 */
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("socket_Flow2Aeroel: wrong client mode");

	Popts   = &opts;
	Popts_1 = &opts_1;
/*
 * set socket settings
 */
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
	m3l_set_Find(&Popts);
/*
 * open socket for reading data 
 */
#pragma omp critical
{
	if( (sockfd = open_connection_to_server(comm_str->IP, comm_str->portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");
}
/*
 * receive data from socket
 */
	if ( (Snode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL)) == NULL)
		Error("socket_Flow2Aeroel: client_receiver()");
/*
 * close socket 
 */
#pragma omp critical
{
	if( close(sockfd) == -1)
		Perror("socket_Flow2Aeroel: close");
}


//=======================================================================

        return Snode;

}

#else

   #include "src_bridges_types.h"
   #include "src_bridges.h"

CPP_C int bridge_aeroelastic(double time){

      printf("This version of flowPsi does not support external communication\n");
      exit;

   }

#endif
