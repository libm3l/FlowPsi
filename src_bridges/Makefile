include ../flowpsi.conf

#FILES = $(shell find ./ -name "*.c" | sed 's/.\///' )
FILES = $(shell ls  *.c )
PATHL  = $(shell pwd )
HFILES = $(shell ls  *.h  |   sed '/lsipdx_header.h/d'  |   sed '/lsipdx.h/d' | sed '/libm3l.h/d' | sort -df)
EMPTYL = $(shell find $(LIBM3LPATH) \( -name libm3l.so -o -iname libm3l.h -o -iname libm3l.so.1.0 \) | wc -l )
EMPTYLS = $(shell find $(LSIPDXPATH) \( -name liblsipdx.so -o -iname lsipdx.h -o -iname liblsipdx.so.1.0 \) | wc -l )

# OBJS=$(FILES:%.c=%.o)
OBJS = $(shell ls *.c | sed 's/\.c/\.o/'  |   sed '/Server_Main.o/d' | sort -df)

CC = mpic++
CFLAGS = -fbounds-check -fstack-check -g 
LFLAGS = -fbounds-check -fstack-check -g 

INCLUDES = -I$(LOCI_BASE)

#CFLAGS =
#LFLAGS =
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
	$(CC) -shared -Wl,-soname,libbridges.so.1.0 -o libbridges.so.1.0   $(OBJS)
	ln -sf libbridges.so.1.0 libbridges.so	
	ar  rcs libbridges.a  *.o    

-include $(OBJS:.o=.d)

%.o: %.c
	$(CC) -c -O  -Wall -fPIC -I$(LIBM3LPATH) -I$(LSIPDXPATH) -I$(LOCI_BASE)/include -$(LIBOPENMPI) $*.c -o $*.o
	$(CC) -MM -O  -Wall -fPIC  -I$(LIBM3LPATH) -I$(LSIPDXPATH) -I$(LOCI_BASE)/include $*.c > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -f *.d  *.o libbridges.a  libbridges.so  libbridges.so.1.0
