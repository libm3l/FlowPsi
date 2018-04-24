PATHL  = $(shell pwd)

include ../flowpsi.conf
INCLUDES = -I$(LIBM3LPATH) -I$(LSIPDXPATH)

all:
	gcc -O -o Client1_FakedSimulink.out Client1_FakedSimulink.c -lpthread -lm $(LIBSLM3LLSPDX) $(INCLUDES)
	gcc -O -o Client1_FakedSimulink_Quaternion.out Client1_FakedSimulink_Quaternion.c -lpthread -lm $(LIBSLM3LLSPDX) $(INCLUDES)
clean:
	rm *.o *.out
