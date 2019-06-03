include ../flowpsi.conf
include ../revision.conf
include $(LOCI_BASE)/Loci.conf

#COPTLESS = -O0 -g
#COPT=$(COPTLESS)


INCLUDES = -I../include -I$(LIBM3LPATH) -I$(LSIPDXPATH) -I../src_bridges


TARGET = flowpsi

FFILES= $(notdir $(wildcard *.loci))
OBJS=$(FFILES:%.loci=%.o)


LOCAL_LIBS = $(LIBXML2_LIBS) -L$(LOCI_BASE)/lib ../src_bridges/libbridges.a

JUNK = *~  core debug .schedule .schedule-*

default:
	$(MAKE) $(TARGET)

all:
	$(MAKE) $(TARGET)


$(TARGET): $(OBJS) 
	$(LD) -o  $(TARGET) $(RPATH)$(FLOWPSI_RPATH) $(OBJS) $(LOCAL_LIBS) $(LIBS) $(LIBSLM3LLSPDX) $(LIBOPENMPI) $(LDFLAGS)
	mv flowpsi ../bin/flowpsi

FRC : 

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

clean:
	rm -fr $(OBJS) $(TARGET) $(JUNK)

distclean: 
	rm -f $(DEPEND_FILES)
	rm -fr $(OBJS) $(TARGET) $(JUNK) $(LOCI_LPP_FILES)


DEPEND_FILES=$(subst .o,.d,$(OBJS))


#include automatically generated dependencies
-include $(DEPEND_FILES)
