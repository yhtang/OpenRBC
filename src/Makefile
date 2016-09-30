VECINFO ?= 0
AALLOC  ?= 0
CXX     ?= g++
ARCH    ?= x86

ifeq ($(CXX),g++)
	ifeq ($(ARCH),x86)
		CFLAGS  = -std=c++11 -g -O4 -march=native -mtune=native -ffast-math -flax-vector-conversions -fopenmp
	else ifeq ($(ARCH),power8)
		CFLAGS	= -std=c++11 -g -O4 -mcpu=power8 -ffast-math -fopenmp
	endif
endif

ifeq ($(CXX),icpc)
CFLAGS	= -std=c++11 -g -Ofast -march=native -mtune=native -fopenmp
endif

ifeq ($(CXX),xlC_r)
	CFLAGS	= -std=c++11 -g -O2 -qhot=level=1 -qarch=pwr8 -qtune=pwr8:smt8 -fopenmp
endif

ifeq ($(DEBUG),1) 
CFLAGS	+= -Og
endif

ifeq ($(VECINFO),1)
CFLAGS += -fopt-info-vec -fopt-info-vec-missed
endif

ifeq ($(AALLOC),1)
CFLAGS += -D_USE_CUSTOM_ALIGNED_ALLOC
endif

#CPP	= $(wildcard *.cpp) 
CPP	= openrbc.cpp
HPP	= $(wildcard *.h) 
BIN	= $(CPP:.cpp=)
 
%: %.cpp $(HPP) 
	$(CXX) $(CFLAGS) $< -o $@ 

default: $(BIN)
 
clean: 
	rm -f $(BIN) 
 
