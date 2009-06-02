CXX = c++
CFLAGS = -Wall -g -I$(XPRESSDIR)/include 

headers = milsp_instance.H instance.H column_generation.H
objects = milsp_instance.o milsp.o column_generation.o 
sources = $(objects:.o=.cc)

all: $(objects)

.cc.o:
	$(CXX) $(CFLAGS) $< -c -o $@

milsp: $(objects) $(headers)
	$(CXX) $(CFLAGS) $(objects) -o $@

.PHONY: clean
clean: 
	rm -Rf $(objects) milsp


