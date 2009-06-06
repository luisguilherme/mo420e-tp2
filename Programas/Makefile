CXX = c++
CFLAGS += -Wall -g -D_REENTRANT -I$(XPRESSDIR)/include
LDFLAGS = -L$(XPRESSDIR)/lib -lxprs -lm  

headers = milsp_instance.H instance.H column_generation.H uls_instance.H
objects = milsp_instance.o milsp.o column_generation.o uls_instance.o
sources = $(objects:.o=.cc)

all: $(sources) $(headers) $(objects)

.cc.o:
	$(CXX) $(CFLAGS) $< -c -o $@

milsp: $(objects) $(headers)
	$(CXX) $(CFLAGS) $(LDFLAGS) $(objects) -o $@

.PHONY: clean
clean: 
	rm -Rf $(objects) milsp


