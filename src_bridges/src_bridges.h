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
 * 1-beta-6  Adam Jirasek         2018-03-21                        Initial Implementation
 *
 *
 *     Description
 * 
 */

#ifndef __SRC_BRIDGES_H__
#define __SRC_BRIDGES_H__

#include <Loci.h>
#include <libm3l.h>

#ifndef CPP_C
#ifdef __cplusplus 
#define CPP_C "C"
#else
#define CPP_C
#endif
#endif

extern CPP_C lmint_t test_bridge(lmdouble_t, lmdouble_t, lmdouble_t,
				  lmdouble_t, lmdouble_t, lmdouble_t,
				  lmdouble_t, lmdouble_t, lmdouble_t,
				  lmdouble_t, lmdouble_t, lmdouble_t);

extern CPP_C lmint_t test_bridge_quaternion(lmdouble_t, lmdouble_t, lmdouble_t, lmdouble_t,
				  lmdouble_t, lmdouble_t, lmdouble_t, lmdouble_t,
				  lmdouble_t, lmdouble_t, lmdouble_t,
				  lmdouble_t, lmdouble_t, lmdouble_t);

#endif

