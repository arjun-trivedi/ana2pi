ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
C++ = g++
CXX = g++
C_FLAGS =
CXXFLAGS =      -O2 -fPIC -w -fmessage-length=0 $(shell root-config --cflags) -Wno-deprecated
INCS =          -I$(CLAS6INC) -I$(HOME)/include -I. -I$(shell root-config --incdir)
SRC = pid.cpp h10looper_e1f.cpp h10looper_2pi.cpp 
OBJS = $(patsubst %.cpp,obj/%.o,$(SRC)) 
LIBS =          $(shell root-config --glibs) -lgfortran
LIBOUT =        $(ELAST_LITE)/sobj/lib_proc_h10_lite.so
TARGET = proc_h10_lite

# The following lines are needed for Fortran compiler (used to compile cut_fid_e16.f)
COMPILER = /usr/bin/gfortran
F77OPT   = -O2 -fno-automatic -w
F77OPT  += -ffixed-line-length-none -fno-second-underscore -DLinux

obj/%.o: %.cpp
	$(CXX) $(INCS) $(CXXFLAGS)  -c $< -o $@

$(TARGET): lib_cut_fid_e16  lib_cut_theta_vs_p_e16 lib_cut_sc_pd_e16 lib_cut_sc_pd_e16 lib_e_corr_sub lib_pi_corr_sub lib_pr lib_vertex_e16 obj/proc_h10_lite.o 
	$(CXX) -o $(TARGET) $(CXXFLAGS) obj/cut_fid_e16.o obj/cut_theta_vs_p_e16.o obj/cut_sc_pd_e16.o obj/e_corr_sub.o obj/pi_corr_sub.o obj/pr.o obj/vertex_e16.o obj/proc_h10_lite.o $(INCS) $(LIBS) -L. $(LIBOUT)

all: $(OBJS) lib $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) obj/cut_fid_e16.o obj/cut_theta_vs_p_e16.o obj/cut_sc_pd_e16.o obj/e_corr_sub.o obj/pi_corr_sub.o obj/pr.o obj/vertex_e16.o obj/proc_h10_lite.o $(TARGET)

$(OBJS): | obj

obj:
	@mkdir -p $@

sobj:
	@mkdir -p $@

lib_cut_fid_e16: cut_fid_e16.f
	$(COMPILER) $(F77OPT) -c $< -o obj/cut_fid_e16.o

lib_cut_theta_vs_p_e16: cut_theta_vs_p_e16.f
	$(COMPILER) $(F77OPT) -c $< -o obj/cut_theta_vs_p_e16.o

lib_cut_sc_pd_e16: cut_sc_pd_e16.f
	$(COMPILER) $(F77OPT) -c $< -o obj/cut_sc_pd_e16.o

lib_e_corr_sub: e_corr_sub.f
	$(COMPILER) $(F77OPT) -c $< -o obj/e_corr_sub.o

lib_pi_corr_sub: pi_corr_sub.f
	$(COMPILER) $(F77OPT) -c $< -o obj/pi_corr_sub.o

lib_pr: pr.f
	$(COMPILER) $(F77OPT) -c $< -o obj/pr.o

lib_vertex_e16: vertex_e16.f
	$(COMPILER) $(F77OPT) -c $< -o obj/vertex_e16.o

lib:	sobj $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) -o $(LIBOUT)
