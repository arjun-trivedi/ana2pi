#SVNREV = -D'SVN_REV="$(shell svnversion -n .)"'
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
OBJS =          epconfig.o data_h10.o data_ana.o data_eid.o data_efid.o data_skim_q.o data_mom.o data_pid.o data_ekin.o data_top.o eid.o ep_processor.o sel_h10.o ep_dict.o
LIBS =          $(shell root-config --glibs) -lProofPlayer
LIBOUT =        libAna2pi.so
#TARGET =	      test
TARGET = ana2pi

%.o: %.cpp
	$(CXX) $(INCS) $(CXXFLAGS) -c $< -o $@

#test: lib test.o
#	$(CXX) -o $(TARGET) $(CXXFLAGS) test.o -L. -lAna2pi $(INCS) $(LIBS)

$(TARGET):	lib ana2pi.o 
	$(CXX) -o $(TARGET) $(CXXFLAGS) ana2pi.o $(INCS) $(LIBS) -L. libAna2pi.so

all:	dict lib $(TARGET)

clean:
	-rm -f $(OBJS) $(LIBOUT) ana2pi ana2pi.o

cleansel:
	rm -f sel_h10.o $(LIBOUT) ana2pi ana2pi.o

lib:	$(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) -o $(LIBOUT)

dict: data_h10.h data_ana.h data_eid.h data_efid.h data_skim_q.h data_mom.h data_pid.h data_ekin.h data_top.h sel_h10.h particle_constants.h 
	-rm ep_dict.*
	@echo ${LD_LIBRARY_PATH}
	rootcint ep_dict.cpp -c data_h10.h data_ana.h data_eid.h data_efid.h data_skim_q.h data_mom.h data_pid.h data_ekin.h data_top.h sel_h10.h particle_constants.h LinkDef.h
	$(CXX) $(INCS) $(CXXFLAGS) -c ep_dict.cpp -o ep_dict.o
