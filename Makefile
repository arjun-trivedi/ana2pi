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
#OBJDIR = obj
#OBJS =          $(OBJDIR)/epconfig.o $(OBJDIR)/h10looper.o $(OBJDIR)/ep_processor.o $(OBJDIR)/data_h10.o $(OBJDIR)/data_ana.o $(OBJDIR)/eid.o $(OBJDIR)/data_eid.o $(OBJDIR)/data_ekin.o $(OBJDIR)/data_efid.o $(OBJDIR)/data_skim_q.o $(OBJDIR)/data_mom.o $(OBJDIR)/data_pid.o $(OBJDIR)/data_top.o 
SRC = epconfig.cpp h10looper.cpp ep_processor.cpp data_h10.cpp data_ana.cpp eid.cpp data_eid.cpp data_ekin.cpp data_efid.cpp data_skim_q.cpp data_mom.cpp data_pid.cpp data_top.cpp
OBJS = $(patsubst %.cpp,obj/%.o,$(SRC)) 
LIBS =          $(shell root-config --glibs)
LIBOUT =        $(WORKSPACE)/ana2pi/sobj/libana.so
TARGET = ana

obj/%.o: %.cpp
	$(CXX) $(INCS) $(CXXFLAGS) -c $< -o $@

$(TARGET):	lib obj/ana.o 
	$(CXX) -o $(TARGET) $(CXXFLAGS) obj/ana.o $(INCS) $(LIBS) -L. $(LIBOUT)

all: lib $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) ana obj/ana.o

$(OBJS): | obj

obj:
	@mkdir -p $@

sobj:
	@mkdir -p $@

lib:	sobj $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) -o $(LIBOUT)

# dict: particle_constants.h data_h10.h data_ana.h data_eid.h data_ekin.h 
# 	-rm ep_dict.*
# 	@echo ${LD_LIBRARY_PATH}
# 	rootcint particle_constants.h ep_dict.cpp -c data_h10.h LinkDef.h
# 	$(CXX) $(INCS) $(CXXFLAGS) -c ep_dict.cpp -o ep_dict.o
