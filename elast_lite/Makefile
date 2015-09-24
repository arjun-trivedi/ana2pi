ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
C++ = g++
CXX = g++
C_FLAGS =
CXXFLAGS =      -O2 -fPIC -w -fmessage-length=0 $(shell root-config --cflags) -Wno-deprecated
INCS =          -I$(CLAS6INC) -I$(HOME)/include -I. -I$(shell root-config --incdir)
SRC = pid.cpp h10looper_e1f.cpp  
OBJS = $(patsubst %.cpp,obj/%.o,$(SRC)) 
LIBS =          $(shell root-config --glibs)
LIBOUT =        $(ELAST_LITE)/sobj/lib_proc_h10_lite.so
TARGET = proc_h10_lite

obj/%.o: %.cpp
	$(CXX) $(INCS) $(CXXFLAGS)  -c $< -o $@

$(TARGET): obj/proc_h10_lite.o 
	$(CXX) -o $(TARGET) $(CXXFLAGS) obj/proc_h10_lite.o $(INCS) $(LIBS) -L. $(LIBOUT)

all: $(OBJS) lib $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) obj/proc_h10_lite.o $(TARGET)

$(OBJS): | obj

obj:
	@mkdir -p $@

sobj:
	@mkdir -p $@

lib:	sobj $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) -o $(LIBOUT)
