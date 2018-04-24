##############################################################################
#     Copyright (C) 2018  Adam Jirasek
# 
#     This file is part of the flowPsi computational fluid dynamics solver.
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#     
#     contact: libm3l@gmail.com
# 
##############################################################################


#
#     Description: - This is a Makefile for src_bridges
#
#     Input parameters:
# 
#
#     Return value:
# 
# History:
# Version   Author:               Date       Patch number  CLA     Comment
# -------   -------               --------   --------      ---     -------
# 1-beta-6  Adam Jirasek         2018-03-21                        Initial Implementation
#
#
#     Description
# 
#

include ../flowpsi.conf

#FILES = $(shell find ./ -name "*.c" | sed 's/.\///' )
FILES = $(shell ls  *.c )
PATHL  = $(shell pwd )
HFILES = $(shell ls  *.h  |   sed '/lsipdx_header.h/d'  |   sed '/lsipdx.h/d' | sed '/libm3l.h/d' | sort -df)
EMPTYL = $(shell find $(LIBM3LPATH) \( -name libm3l.so -o -iname libm3l.h -o -iname libm3l.so.1.0 \) | wc -l )
EMPTYLS = $(shell find $(LSIPDXPATH) \( -name liblsipdx.so -o -iname lsipdx.h -o -iname liblsipdx.so.1.0 \) | wc -l )

OBJS=$(FILES:%.c=%.o)

CC = mpic++
CFLAGS = -fbounds-check -fstack-check -g 
LFLAGS = -fbounds-check -fstack-check -g 

INCLUDES = -I$(LIBM3LPATH) -I$(LSIPDXPATH) -I$(LOCI_BASE)/include -I../include

#
#Autodependencies with GNU make
#Scott McPeak, November 2001 
#
# link
#
main:
#
#  delete empty links
#
	@find -L -type l -delete

ifeq ($(EMPTYL),3)

	make prog
else
#
#  did not find libm3l files, posibly LIBM3LPATH#L wrong
#
	@echo ERROR: Did not find libml3 library and header file
	@echo
	@echo Check LIBM3LPATH variable setting  in flowpsi.conf file
	@echo
endif

ifeq ($(EMPTYLS),3)

	make prog
else
#
#  did not find libm3l files, posibly LIBM3LPATH#L wrong
#
	@echo ERROR: Did not find lsipdx library and header file
	@echo
	@echo Check LSIPDXPATH variable setting in flowpsi.conf file
	@echo
endif


prog: $(OBJS)

	$(CC) -c $(OBJS)
	ar  rcs libbridges.a  *.o    

-include $(OBJS:.o=.d)

%.o: %.c
	$(CC) -c -O  -Wall -I$(LIBM3LPATH) -I$(LSIPDXPATH) -I$(LOCI_BASE)/include -I../include -$(LIBOPENMPI) $*.c -o $*.o
	$(CC) -MM -O  -Wall   -I$(LIBM3LPATH) -I$(LSIPDXPATH) -I$(LOCI_BASE)/include -I../include $*.c > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -f *.d  *.o libbridges.a  libbridges.so  libbridges.so.1.0
