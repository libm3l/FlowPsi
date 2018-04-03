PATHL  = $(shell pwd)

include ../flowpsi.conf

all:

	ln -sf  $(LIBM3LPATH)/libm3l.* .
	ln -sf  $(LSIPDXPATH)/liblsipdx.so* .
	ln -sf  $(LSIPDXPATH)/lsipdx.h .

	gcc -O -o Client1_FakedSimulink.out Client1_FakedSimulink.c -L$(LSIPDXPATH)  -llsipdx -Wl,-rpath=$(LSIPDXPATH) -L$(LIBM3LPATH)  -lm3l -Wl,-rpath=$(LIBM3LPATH) -lpthread -lm

	gcc -O -o Client1_FakedSimulinkIC1.out Client1_FakedSimulinkIC1.c -L$(LSIPDXPATH)  -llsipdx -Wl,-rpath=$(LSIPDXPATH) -L$(LIBM3LPATH)  -lm3l -Wl,-rpath=$(LIBM3LPATH) -lpthread -lm

	gcc -O -o Client1_FakedSimulinkIC2.out Client1_FakedSimulinkIC2.c -L$(LSIPDXPATH)  -llsipdx -Wl,-rpath=$(LSIPDXPATH) -L$(LIBM3LPATH)  -lm3l -Wl,-rpath=$(LIBM3LPATH) -lpthread -lm

	gcc -O -o Client1_FakedSimulinkIC3.out Client1_FakedSimulinkIC3.c -L$(LSIPDXPATH)  -llsipdx -Wl,-rpath=$(LSIPDXPATH) -L$(LIBM3LPATH)  -lm3l -Wl,-rpath=$(LIBM3LPATH) -lpthread -lm
clean:
	rm *.o *.out
