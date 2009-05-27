CXX = c++
CFLAGS = -Wall -g

headers = milsp_instance.H
objects = milsp_instance.o milsp.o
sources = $(objects:.o=.cc)

all: $(objects)

.cc.o:
	$(CXX) $(CFLAGS) $< -c -o $@

milsp: $(objects) $(headers)
	$(CXX) $(CFLAGS) $(objects) -o $@

.PHONY: clean
clean: 
	rm -Rf $(objects) milsp


