/*
 *     Copyright (C) 2018  Adam Jirasek
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
 */
/*
 * 
 * Description: contains example of module receving forces/moments and sending back
 *       values of quaternion, rotation center and translation
 *
 *       this function receives data from a separatelly run process. The data is contained in a head node
 *       CFD_2_SIM and consists of double array of 6 called ForcesMoments and single double value called Time
 *
 *  H-DIR           CFD_2_SIM               2
 *     -D               ForcesMoments           1      6   
 *     -D               Time            1      1   
 *
 *      this data set is recevied through channel CFD2SIM
 *
 *      as a respons, it creates a data set in head node SIM_2_CFD
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
#include "libm3l.h"
#include "lsipdx.h"

int main(int argc, char *argv[])
{
/*
 *  definition of variables
 */
	node_t *Gnode=NULL, *Snode=NULL, *FoundNode=NULL, *TmpNode=NULL;
	size_t i, niter, dim[1], tot_dim, narray;

	lmint_t sockfd, portno;

        socklen_t clilen;
        struct sockaddr_in cli_addr;
/*
 *   here define name of channel through which to receive the data and through which to send the data
 */
	lmchar_t *name ="CFD2SIM";
	lmchar_t *name1="SIM2CFD";

	lmdouble_t *P, dy, *tmpfloat, *x, *y, *z, *time, sign;
        lmdouble_t yaw,pitch,roll,cy,sy,cr,sr,cp,sp,qw,qx,qy,qz;
	lmdouble_t psi;
	
	find_t *SFounds;
	
	opts_t opts, *Popts_1;

	client_fce_struct_t InpPar, *PInpPar;

	PInpPar = &InpPar;
/*
 * get port number
 */
	if (argc < 3) {
		fprintf(stderr,"ERROR, no IPaddress and port number provided\n");
		exit(1);
	}
 	portno = atoi(argv[2]);
/*
 * open socket - hostname is in argv[1] and port number is in portno
 */
	niter = 0;
 	while(1){
/*
 * this is an infinite loop, in each cycle it receves and send data. The loop is infinite
 */

 		printf("\n\n--------------------------------    i = %ld\n\n", ++niter);
/*
 * open socket - set parameters which are needed for opening the socket for receving data
 *  
 */
		PInpPar->channel_name = name;   /* name of channel */
		PInpPar->SR_MODE = 'R';         /* set R or S, because we will read incoming data from this 
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
 *        need to use the same settings - see ../src_bridges/test_bridge.c)
 *
 *  for more details on modes see:  
 *  Adam Jirasek and Arthur W. Rizzi: libm3l and lsipdx - Utilities for Inter-Process Data Transfer and Synchronization", 
 *  52nd Aerospace Sciences Meeting, AIAA SciTech Forum, (AIAA 2014-1045) https://doi.org/10.2514/6.2014-1045
 */
		if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
			Error("wrong client mode");
/*
 * this is just a structure containing additional data
 */
		Popts_1 = &opts;
/*
 *  use the set options to set parameters for connection
 */
		m3l_set_Send_receive_tcpipsocket(&Popts_1);
/* 
 * open socket, use host name and port number as input parameters, use PInpPar which contains 
 * structure with mode for this socket
 */
		if( (sockfd = open_connection_to_server(argv[1], portno, PInpPar, Popts_1)) < 1)
			Error("client_sender: Error when opening socket");		
/*
 * receive data set
 */
		Gnode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL);
/* 
 * print data set on screen
 */	
 		if(m3l_Cat(Gnode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
 			Error("CatData");
/*
 * locate array ForcesMoments which is located in CFD_2_SIM
 */
	if( (SFounds = m3l_Locate(Gnode, "/CFD_2_SIM/ForcesMoments", "/*/*",  (lmchar_t *)NULL)) != NULL){
/*
 * this does not have to be done, but just for sake of programming, check that CFD_2_SIM contains only
 * one array ForcesMoments
 */
		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_FlowPsi2out: More then one ForcesMoments data set found");
/* 
 * get pointer to list of found nodes, get its total dimensions
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_FlowPsi2out: Did not find 1st data pointer");
		tot_dim = m3l_get_List_totdim(FoundNode);
		narray = tot_dim;
		
		printf(" Size of array is %ld\n", narray);
/*
 *  ... and get double data it contains
 */	
		if( (x = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_FlowPsi2out: Did not find DX data pointer");
/* 
 * free memory allocated in m3l_Locate - this has to be done awlays once we do not use Found data set any more
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_FlowPsi2out: ForcesMoments not found\n");
	}
/*
 * similarly, find data set containing Time and get the value of time
 */
	SFounds = m3l_Locate(Gnode, "/CFD_2_SIM/Time", "/*/*",  (lmchar_t *)NULL);
	FoundNode = m3l_get_Found_node(SFounds, 0);
	time = (lmdouble_t *)m3l_get_data_pointer(FoundNode);
	m3l_DestroyFound(&SFounds);

	printf("Time is %lf  tot-dim is %ld\n", *time, tot_dim);
/*
 *  close socket
 */
	if( close(sockfd) == -1)
		Perror("close");
/*
 * this section prepares data which will be sent back and send them back
 * open socket for sending data back
 */	
	PInpPar->channel_name = name1;   /* name of channel where to send data*/
	PInpPar->SR_MODE = 'S';          /* set R or S, because we will write outgoing data to this 
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
 *        need to use the same settings - see ../src_bridges/test_bridge.c)
 *
 *  for more details on modes see:  
 *  Adam Jirasek and Arthur W. Rizzi: libm3l and lsipdx - Utilities for Inter-Process Data Transfer and Synchronization", 
 *  52nd Aerospace Sciences Meeting, AIAA SciTech Forum, (AIAA 2014-1045) https://doi.org/10.2514/6.2014-1045
 */
	if ( (PInpPar->mode = get_exchange_channel_mode('D', 'N')) == -1)
		Error("wrong client mode");
	Popts_1 = &opts;
/*
 *  use the set options to set parameters for connection
 */
	m3l_set_Send_receive_tcpipsocket(&Popts_1);
/* 
 * open socket, use host name and port number as input parameters, use PInpPar which contains 
 * structure with mode for this socket
 */	
	if( (sockfd = open_connection_to_server(argv[1], portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");
/*
 * Create Snode containing outgoing data and call it SIM_2_CFD
 */
	if(  (Snode = m3l_Mklist("SIM_2_CFD", "DIR", 0, 0, (node_t **)NULL, (const char *)NULL, (const char *)NULL, (char *)NULL)) == 0)
		Perror("m3l_Mklist");
/*
 * add to it double array of 4 called Quaternion and store it in SIM_2_CFD
 */
	dim[0] = 4;
	if(  (TmpNode = m3l_Mklist("Quaternion", "D", 1, dim, &Snode, "/SIM_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
/*
 * get pointer to recently created array
 */
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
/*
 * calculate angle
 */
	pitch = -2.41*sin(*time*2*3.1415926*50.32);
//	pitch = 15*sin(*time*2*3.1415926*5);
        roll = 0;
        yaw = 0;
	printf("Pitch angle is %lf\n", pitch);
        pitch = pitch *3.1415926/180.;
/*
 * transofrm to quaternion
 */
	cy = cos(yaw * 0.5);
	sy = sin(yaw * 0.5);
	cr = cos(roll * 0.5);
	sr = sin(roll * 0.5);
	cp = cos(pitch * 0.5);
	sp = sin(pitch * 0.5);

	qw = cy * cr * cp + sy * sr * sp;
	qx = cy * sr * cp - sy * cr * sp;
	qy = cy * cr * sp + sy * sr * cp;
	qz = sy * cr * cp - cy * sr * sp;
/*
 * store values in quaternion array
 */
	tmpfloat[0] = qw;
	tmpfloat[1] = qx;
	tmpfloat[2] = qy;
	tmpfloat[3] = qz;
/*
 *  create double array of 3 called RotCenter containing point of rotation, store it in SIM_2_CFD and fill it with data
 */
	dim[0] = 3;
	if(  (TmpNode = m3l_Mklist("RotCenter", "D", 1, dim, &Snode, "/SIM_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
	tmpfloat[0] = 0.25;
	tmpfloat[1] = 0.;
	tmpfloat[2] = 0.;
/*
 *  create double array of 3 called TransVec containing translation vector, store it in SIM_2_CFD and fill it with data
 */
	dim[0] = 3;
	if(  (TmpNode = m3l_Mklist("TransVec", "D", 1, dim, &Snode, "/SIM_2_CFD", "./", (char *)NULL)) == 0)
		Error("m3l_Mklist");
	tmpfloat = (lmdouble_t *)m3l_get_data_pointer(TmpNode);
	tmpfloat[0] = 0;
	tmpfloat[1] = 0;
	tmpfloat[2] = 0;
/*
 * print data on screen
 */
 	if(m3l_Cat(Snode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
 		Error("CatData");
/*
 * send data through socket
 */
	client_sender(Snode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL);
/*
 * during this iteration, we allocated memory, either by getting data set from socket or 
 * creatign in by ourselves. Deallocate
 */	
	if(m3l_Umount(&Gnode) != 1)
		Perror("m3l_Umount");
	if(m3l_Umount(&Snode) != 1)
		Perror("m3l_Umount");
/* 
 * close socket
 */
	if( close(sockfd) == -1)
		Perror("close");
 	}


     return 0; 
}
