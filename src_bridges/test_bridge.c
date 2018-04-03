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
 *     Description: - this is a test of bridge function
 *       this function creates data which are located in head node
 *       CFD_2_SIM and consists of double array of 6 called ForcesMoments and single double value called Time
 *
 *  H-DIR           CFD_2_SIM               2
 *     -D               ForcesMoments           1      6   
 *     -D               Time            1      1   
 *
 *      this data set is recevied through channel CFD2SIM
 *
 *      as a respons, it receives a data set in head node SIM_2_CFD
 *
 *  H-DIR           SIM_2_CFD               3
 *   -D                     Angles                  1      3   
 *   -D                     RotCenter               1      3   
 *   -D                     TransVec                1      3   
 *
 *      which contains double array of 3 called Angles, double array of 3 called RotCenter and double array of 
 *      3 called TransVec
 *
 *      to create and manipulate with this type of data set, you need to link two libraries - libm3l and lsipdx 
 *      for more details go to www.github.com/libm3l/libm3l and www.github.com/libm3l/lsipdx
 *
 *     Input parameters:
 * 
 *
 *     Return value:
 * 
 * History:
 * Version   Author:               Date       Patch number  CLA     Comment
 * -------   -------               --------   --------      ---     -------
 * 1.1       Adam Jirasek    2018-03-21                       Initial implementation
 *
 *
 *     Description
 * 
 */

#include "libm3l.h"
#include "lsipdx.h"
#include "src_bridges_types.h"
#include "src_bridges.h"


lmint_t test_bridge(lmdouble_t ForceX, lmdouble_t ForceY, lmdouble_t ForceZ , 
                    lmdouble_t Alpha, lmdouble_t Beta, lmdouble_t Gamma,
                    lmdouble_t TransX, lmdouble_t TransY, lmdouble_t TransZ,
                    lmdouble_t RotCX, lmdouble_t RotCY, lmdouble_t RotCZ){


	node_t *Gnode=NULL, *TmpNode = NULL, *FoundNode = NULL;
	lmsize_t dim[1], i, tot_dim;

//	lmchar_t hostname[80], channel_name[80];

	lmint_t sockfd, portno;

	lmchar_t *name ="CFD2SIM";
	lmchar_t *name1="SIM2CFD";
        lmchar_t *hostname ="localhost";

	lmdouble_t *tmpfloat;
	client_fce_struct_t InpPar, *PInpPar;
	opts_t *Popts_1, opts, opts_1, *Popts;
	find_t *SFounds;
/*
 * Set parameters which are needed for opening the socket for sending data
 */
        PInpPar = &InpPar;
	PInpPar->channel_name = name;   /* name of channel where to send data*/
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
		Error("socket_FlowPsi2simulink: wrong client mode");
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

        portno = 31000;

//=======================================================================

/*
 * create data structure which will be sent, call it CFD_2_SIM
 */        
	if(  (Gnode = m3l_Mklist("CFD_2_SIM", "DIR", 0, 0, (node_t **)NULL, (const char *)NULL, (const char *)NULL, (char *)NULL)) == 0)
		Perror("socket_FlowPsi2simulink: m3l_Mklist");
	
	dim[0] = 6;
/*
 * create and store double array of 6 with name ForcesMoments in CFD_2_SIM
 */
	if(  (TmpNode = m3l_Mklist("ForcesMoments", "D", 1, dim, &Gnode, "/CFD_2_SIM", "./", (char *)NULL)) == 0)
		Error("socket_FlowPsi2simulink: m3l_Mklist");
        for (i=0; i<6; i++)
	  TmpNode->data.df[i] = 1.*i;

	TmpNode->data.df[0] = ForceX;
	TmpNode->data.df[1] = ForceY;
	TmpNode->data.df[2] = ForceZ;
/*
 * add time, store it in double variable Time in CFD_2_SIM
 */
	dim[0] = 1;
	if(  (TmpNode = m3l_Mklist("Time", "D", 1, dim, &Gnode, "/CFD_2_SIM", "./", (char *)NULL)) == 0)
		Error("socket_FlowPsi2simulink: m3l_Mklist");
	TmpNode->data.df[0] = 1.54;
/*
 * print data on screen
 */
        if(m3l_Cat(Gnode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
 	  Error("CatData");
/*
 * open socket
 */
	if( (sockfd = open_connection_to_server(hostname, portno, PInpPar, Popts_1)) < 1)
		Error("socket_FlowPsi2simulink: Error when opening socket");
/*
 * send data 
 */
	if ( client_sender(Gnode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL) !=1 )
		Error("socket_FlowPsi2simulink: client_sender()");
/*
 * close socket
 */
	if( close(sockfd) == -1)
		Perror("socket_FlowPsi2simulink: close");
/*
 * free borrowed memory
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_FlowPsi2simulink: m3l_Umount");
/*
 * receive data 
 */
	PInpPar = &InpPar;
	PInpPar->channel_name = name1;       /* name of channel */
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
		Error("socket_FlowPsi2simulink: wrong client mode");

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
	if( (sockfd = open_connection_to_server(hostname, portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");
/*
 * receive data from socket
 */
	if ( (Gnode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL)) == NULL)
		Error("socket_FlowPsi2simulink: client_receiver()");
/*
 * close socket 
 */
	if( close(sockfd) == -1)
		Perror("socket_FlowPsi2simulink: close");
/*
 * print data on screen
 */
        if(m3l_Cat(Gnode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
 	  Error("CatData");
/*
 * find Angles - rotation matrix and copy the values to FlowPsi allocated memory
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/Angles", "/*/*",  (lmchar_t *)NULL)) != NULL){
/*
 * check that only one data set called Angles is present
 */
		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_FlowPsi2simulink: More then one Angles data set found");
/* 
 * point to the first list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find 1st data pointer");
/*
 * get array dimensions
 */
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_FlowPsi2simulink: Wrong dimensions of Angles array");
/*
 * get pointer on actual, double array data
 */
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find Angles data pointer");

                Alpha = tmpfloat[0];
                Beta  = tmpfloat[1];
                Gamma = tmpfloat[2];

		//for (i=0; i<tot_dim; i++)
		//	Angles[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_FlowPsi2simulink: Angles not found\n");
	}
/*
 * find center of rotation, procedure identical as in case of Angles
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/RotCenter", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_FlowPsi2simulink: More then one RotCenter data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find 1st data pointer");
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_FlowPsi2simulink: Wrong dimensions of RotCenter array");
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find RotCenter data pointer");

                RotCX = tmpfloat[0];
                RotCY = tmpfloat[1];
                RotCZ = tmpfloat[2];
		//for (i=0; i<tot_dim; i++)
		//	RotCenter[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_FlowPsi2simulink: RotCenter not found\n");
	}
/*
 * find center of translation, procedure identical as in case of Angles
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/TransVec", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_FlowPsi2simulink: More then one TransVec data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find 1st data pointer");
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_FlowPsi2simulink: Wrong dimensions of TransVec array");
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find TransVec data pointer");

                TransX = tmpfloat[0];
                TransY = tmpfloat[1];
                TransZ = tmpfloat[2];
		//for (i=0; i<tot_dim; i++)
		//	TransVec[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_FlowPsi2simulink: TransVec not found\n");
	}
/*
 * free borrowed memory (it was allocated during process of receiving data)
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_FlowPsi2simulink: m3l_Umount");


//=======================================================================

        return 1;

}


lmint_t test_bridge1(comm_struct_t *pcomm_str){

        lmdouble_t ForceX; lmdouble_t ForceY; lmdouble_t ForceZ ; 
                    lmdouble_t Alpha; lmdouble_t Beta; lmdouble_t Gamma;
                    lmdouble_t TransX; lmdouble_t TransY; lmdouble_t TransZ;
                    lmdouble_t RotCX; lmdouble_t RotCY; lmdouble_t RotCZ;

	node_t *Gnode=NULL, *TmpNode = NULL, *FoundNode = NULL;
	lmsize_t dim[1], i, tot_dim;

	lmint_t sockfd;

	lmdouble_t *tmpfloat;
	client_fce_struct_t InpPar, *PInpPar;
	opts_t *Popts_1, opts, opts_1, *Popts;
	find_t *SFounds;

        ForceX = 0;
        ForceY = 0;
        ForceZ = 0;

  //      printf("INPUT NAME EXAMPLE IS %s\n",pcomm_str->I_channel);
  //      printf("OUTPUT NAME EXAMPLE IS %s\n",pcomm_str->O_channel);
  //      printf("PORTNO  EXAMPLE IS %d\n",pcomm_str->portno);

       // return 1;


/*
 * Set parameters which are needed for opening the socket for sending data
 */
        PInpPar = &InpPar;
	PInpPar->channel_name = pcomm_str->O_channel;   /* name of channel where to send data*/
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
		Error("socket_FlowPsi2simulink: wrong client mode");
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

//=======================================================================

/*
 * create data structure which will be sent, call it CFD_2_SIM
 */        
	if(  (Gnode = m3l_Mklist("CFD_2_SIM", "DIR", 0, 0, (node_t **)NULL, (const char *)NULL, (const char *)NULL, (char *)NULL)) == 0)
		Perror("socket_FlowPsi2simulink: m3l_Mklist");
	
	dim[0] = 6;
/*
 * create and store double array of 6 with name ForcesMoments in CFD_2_SIM
 */
	if(  (TmpNode = m3l_Mklist("ForcesMoments", "D", 1, dim, &Gnode, "/CFD_2_SIM", "./", (char *)NULL)) == 0)
		Error("socket_FlowPsi2simulink: m3l_Mklist");
        for (i=0; i<6; i++)
	  TmpNode->data.df[i] = 1.*i;

	TmpNode->data.df[0] = ForceX;
	TmpNode->data.df[1] = ForceY;
	TmpNode->data.df[2] = ForceZ;
/*
 * add time, store it in double variable Time in CFD_2_SIM
 */
	dim[0] = 1;
	if(  (TmpNode = m3l_Mklist("Time", "D", 1, dim, &Gnode, "/CFD_2_SIM", "./", (char *)NULL)) == 0)
		Error("socket_FlowPsi2simulink: m3l_Mklist");
	TmpNode->data.df[0] = 1.54;
/*
 * print data on screen
 */
        if(m3l_Cat(Gnode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
 	  Error("CatData");
/*
 * open socket
 */
#pragma omp critical
{
        printf(" Opening channel %s  %d   %s  \n", pcomm_str->IP, pcomm_str->portno, PInpPar->channel_name);
	if( (sockfd = open_connection_to_server(pcomm_str->IP, pcomm_str->portno, PInpPar, Popts_1)) < 1)
		Error("socket_FlowPsi2simulink: Error when opening socket");
}
/*
 * send data 
 */
	if ( client_sender(Gnode, sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL) !=1 )
		Error("socket_FlowPsi2simulink: client_sender()");
/*
 * close socket
 */
#pragma omp critical
{
	if( close(sockfd) == -1)
		Perror("socket_FlowPsi2simulink: close");
}
/*
 * free borrowed memory
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_FlowPsi2simulink: m3l_Umount");
/*
 * receive data 
 */
	PInpPar = &InpPar;
	PInpPar->channel_name = pcomm_str->I_channel;    /* name of channel */
	PInpPar->SR_MODE = 'R';                          /* set R or S, because we will read incoming  
                                                            data from this socket, set it to R */
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
		Error("socket_FlowPsi2simulink: wrong client mode");

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
	if( (sockfd = open_connection_to_server(pcomm_str->IP, pcomm_str->portno, PInpPar, Popts_1)) < 1)
		Error("client_sender: Error when opening socket");
}
/*
 * receive data from socket
 */
	if ( (Gnode = client_receiver(sockfd, PInpPar, (opts_t *)NULL, (opts_t *)NULL)) == NULL)
		Error("socket_FlowPsi2simulink: client_receiver()");
/*
 * close socket 
 */
#pragma omp critical
{
	if( close(sockfd) == -1)
		Perror("socket_FlowPsi2simulink: close");
}
/*
 * print data on screen
 */
        if(m3l_Cat(Gnode, "--all", "-P", "-L",  "*",   (char *)NULL) != 0)
 	  Error("CatData");
/*
 * find Angles - rotation matrix and copy the values to FlowPsi allocated memory
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/Angles", "/*/*",  (lmchar_t *)NULL)) != NULL){
/*
 * check that only one data set called Angles is present
 */
		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_FlowPsi2simulink: More then one Angles data set found");
/* 
 * point to the first list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find 1st data pointer");
/*
 * get array dimensions
 */
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_FlowPsi2simulink: Wrong dimensions of Angles array");
/*
 * get pointer on actual, double array data
 */
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find Angles data pointer");

                Alpha = tmpfloat[0];
                Beta  = tmpfloat[1];
                Gamma = tmpfloat[2];

		//for (i=0; i<tot_dim; i++)
		//	Angles[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_FlowPsi2simulink: Angles not found\n");
	}
/*
 * find center of rotation, procedure identical as in case of Angles
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/RotCenter", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_FlowPsi2simulink: More then one RotCenter data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find 1st data pointer");
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_FlowPsi2simulink: Wrong dimensions of RotCenter array");
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find RotCenter data pointer");

                RotCX = tmpfloat[0];
                RotCY = tmpfloat[1];
                RotCZ = tmpfloat[2];
		//for (i=0; i<tot_dim; i++)
		//	RotCenter[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_FlowPsi2simulink: RotCenter not found\n");
	}
/*
 * find center of translation, procedure identical as in case of Angles
 */
	if( (SFounds = m3l_Locate(Gnode, "/SIM_2_CFD/TransVec", "/*/*",  (lmchar_t *)NULL)) != NULL){

		if( m3l_get_Found_number(SFounds) != 1)
			Error("socket_FlowPsi2simulink: More then one TransVec data set found");
/* 
 * pointer to list of found nodes
 */
		if( (FoundNode = m3l_get_Found_node(SFounds, 0)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find 1st data pointer");
		if( (tot_dim = m3l_get_List_totdim(FoundNode)) != 3)
			Error("socket_FlowPsi2simulink: Wrong dimensions of TransVec array");
		if( (tmpfloat = (lmdouble_t *)m3l_get_data_pointer(FoundNode)) == NULL)
			Error("socket_FlowPsi2simulink: Did not find TransVec data pointer");

                TransX = tmpfloat[0];
                TransY = tmpfloat[1];
                TransZ = tmpfloat[2];
		//for (i=0; i<tot_dim; i++)
		//	TransVec[i]  = tmpfloat[i];
/* 
 * free memory allocated in m3l_Locate
 */
		m3l_DestroyFound(&SFounds);
	}
	else
	{
		Error("socket_FlowPsi2simulink: TransVec not found\n");
	}
/*
 * free borrowed memory (it was allocated during process of receiving data)
 */
	if(m3l_Umount(&Gnode) != 1)
		Perror("socket_FlowPsi2simulink: m3l_Umount");


//=======================================================================

        return 1;

}
