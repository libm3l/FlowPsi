# Set FLOWPSI_BASE to point to the FLOWPSI installation directory to compile
# this module

# To use a module, compile the module and then put the resulting shared object
# file either in the installed lib directory with flowpsi or some other place 
# that is contained in the LD_LIBRARY_PATH.  Once this is done, the module 
# can be loaded into flowpsi using the -load_module option.  For example
# chem -load_module extensions case_name
# where extensions above is the MODULE_NAME given below
# multiple modules can be loaded with multiple -load_module options
include $(FLOWPSI_BASE)/flowpsi.conf
include $(LOCI_BASE)/Loci.conf

#COPT= -g -O0
#COPTLESS= -g -O0 

INCLUDES = -DUSE_LOCI_TYPES -I$(FLOWPSI_BASE)/include -I./include

# Put objects in the module here
OBJS =	solverGust.o  

# List your object files here
VOGDEFOBJS =  solverGust.o 

MODULE_NAME = Gust

LOCAL_LIBS = 

JUNK = *~  core ti_files ii_files rii_files

LIB_OBJS=$(OBJS:.o=_lo.o)

$(MODULE_NAME)FPSI_m.so: $(LIB_OBJS)
	$(SHARED_LD) $(SHARED_LD_FLAGS) $(MODULE_NAME)FPSI_m.so $(LIB_FLAGS) $(LIB_OBJS) $(LIBS)
	mv $(MODULE_NAME)FPSI_m.so ../../lib/$(MODULE_NAME)FPSI_m.so


all: $(MODULE_NAME)FPSI_m.so 

FRC : 

clean:
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)FPSI_m.so  $(JUNK)

LOCI_FILES = $(wildcard *.loci)
LOCI_LPP_FILES = $(LOCI_FILES:.loci=.cc)

distclean: 
	rm -fr $(DEPEND_FILES)
	rm -fr $(OBJS) $(LIB_OBJS) $(MODULE_NAME)FPSI_m.so $(VOGDEFOBJS) $(JUNK) $(LOCI_LPP_FILES)

DEPEND_FILES=$(subst .o,.d,$(OBJS)) 


#include automatically generated dependencies
-include $(DEPEND_FILES)
