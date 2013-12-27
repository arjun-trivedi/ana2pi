ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
F77 = gfortran
F_FLAGS =
C++ = g++
CXX = g++
C_FLAGS =
CXXFLAGS =      -O2 -fPIC -w -fmessage-length=0 $(shell root-config --cflags) -Wno-deprecated
INCS =          -I$(CLAS6INC) -I$(HOME)/include -I. -I$(shell root-config --incdir)
OBJS =          data_h10.o 
LIBS =          $(shell root-config --glibs)
LIBOUT =        libAna2pi.so
TARGET = ana2pi

%.o: %.cpp
	$(CXX) $(INCS) $(CXXFLAGS) -c $< -o $@

$(TARGET):	lib ana2pi.o 
	$(CXX) -o $(TARGET) $(CXXFLAGS) ana2pi.o $(INCS) $(LIBS) -L. libAna2pi.so

all:	dict lib $(TARGET)

clean:
	-rm -f $(OBJS) $(LIBOUT) ana2pi ana2pi.o

lib:	$(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) -o $(LIBOUT)

dict: data_h10.h particle_constants.h 
	-rm ep_dict.*
	@echo ${LD_LIBRARY_PATH}
	rootcint ep_dict.cpp -c data_h10.h particle_constants.h LinkDef.h
	$(CXX) $(INCS) $(CXXFLAGS) -c ep_dict.cpp -o ep_dict.o
