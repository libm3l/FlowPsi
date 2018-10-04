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
 *     Description: - this is a header file for test of bridge function
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
#ifdef LIBM3LSIPDX

#include <libm3l.h> 

#ifndef __SRC_BRIDGES_TYPES_H__
#define __SRC_BRIDGES_TYPES_H__

typedef struct comm_struct {
    lmchar_t type[80];         /* interface type */
    lmchar_t VolTag[80];    /* volume tag  */
    lmchar_t I_channel[80];    /* channel in  */
    lmchar_t O_channel[80];    /* channel out */
    lmchar_t IP[80];           /* IP address of server */
    lmchar_t intf_name[80]; /* interface name */
    lmint_t portno;            /* port number */
    lmint_t comm_freq;         /* communicatin frequnecy */
    
}comm_struct_t;

#endif

#else

    typedef struct comm_struct {
    char type[80];
    char VolTag;
    char I_channel[80];
    char O_channel[80];
    char IP[80];    
    char intf_name;
    
    int portno;
    int comm_freq;
    
}comm_struct_t;

#endif
