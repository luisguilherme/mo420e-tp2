# makefile para compilar programa "cutting-stock.c"
# que usa geração de colunas no Xpress para resolver
# a relaxação linear do cutting-stock.
# para compilar outro programa, faca 
#    make TARGETS=outro.o
# plataformas: linux solaris hpux aix
#
# Preparado por Cid de Souza (IC-UNICAMP) e
#               Edna Hoshino (DTC-UFMS)
OSTYPE=linux



CFLAGS-linux=-D_REENTRANT
CFLAGS-linux-gnu=$(CFLAGS-linux)
CFLAGS-solaris=-D_REENTRANT
CFLAGS-hpux=+Z +DAportable -D_POSIX_C_SOURCE=199506L
CFLAGS-aix=-q32 -D_THREAD_SAFE
CFLAGS=-g $(CFLAGS-$(OSTYPE)) -I$(XPRESSDIR)/include 

LDFLAGS-hpux=-Wl,+s
LDFLAGS-aix=-brtl -bnolibpath
LDFLAGS=$(LDFLAGS-$(OSTYPE)) -L$(XPRESSDIR)/lib


LOADLIBS-aix=-lC
LOADLIBS=$(LOADLIBS-$(OSTYPE)) -lm -lxprs

.c:
	$(CC) $(CFLAGS) $(LDFLAGS)  $< $(LOADLIBS) -o $*


TARGETS=cutting

All: $(TARGETS)
