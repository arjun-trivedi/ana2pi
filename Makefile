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
SRC = epconfig.cpp h10looper.cpp ep_processor.cpp data_h10.cpp data_ana.cpp eid.cpp data_eid.cpp data_ekin.cpp data_efid.cpp data_eeff.cpp data_skim_q.cpp data_skim_q_elast.cpp data_mom.cpp data_pid.cpp data_peff.cpp data_pid_elast.cpp data_2pi.cpp data_elastic.cpp cut_eff.cpp
OBJS = $(patsubst %.cpp,obj/%.o,$(SRC)) 
LIBS =          $(shell root-config --glibs)
LIBOUT =        $(WORKSPACE)/ana2pi/sobj/lib_proc_h10.so
TARGET = proc_h10

obj/%.o: %.cpp
	$(CXX) $(INCS) $(CXXFLAGS) -c $< -o $@

$(TARGET):	lib obj/proc_h10.o 
	$(CXX) -o $(TARGET) $(CXXFLAGS) obj/proc_h10.o $(INCS) $(LIBS) -L. $(LIBOUT)

all: dict lib $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) proc_h10 obj/proc_h10.o

$(OBJS): | obj

obj:
	@mkdir -p $@

sobj:
	@mkdir -p $@

lib:	sobj $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) ep_dict.o $(LIBS) -o $(LIBOUT)

dict: data_h10.h data_ana.h data_eid.h data_efid.h data_eeff.h data_skim_q.h data_skim_q_elast.h data_mom.h data_pid.h data_peff.h data_pid_elast.h data_ekin.h data_2pi.h h10looper.h particle_constants.h data_elastic.h
	-rm ep_dict.*
	@echo ${LD_LIBRARY_PATH}
	rootcint ep_dict.cpp -c data_h10.h data_ana.h data_eid.h data_efid.h data_eeff.h data_skim_q.h data_skim_q_elast.h data_mom.h data_pid.h data_peff.h data_pid_elast.h data_ekin.h data_2pi.h h10looper.h particle_constants.h data_elastic.h LinkDef.h
	$(CXX) $(INCS) $(CXXFLAGS) -c ep_dict.cpp -o ep_dict.o
