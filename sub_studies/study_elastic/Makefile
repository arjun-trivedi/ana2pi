ATEID   = ${HOME}/CLAS/workspace/at-eid
EPCUTS  = ${HOME}/CLAS/workspace/epana/inc
MYUTILS = ${HOME}/CLAS/workspace/utils

DIR_OBJ  =	./obj
DIR_SOBJ = 
DIR_BIN  =	.

TARGET      =   test#anaCCEC_sim anaCCEC_exp verifyeCuts_sim verifyeCuts_exp elasYield_sim elasYield_exp
OBJS        =   
SOBJS       =  
OBJ_PATHS   =   $(addprefix $(DIR_OBJ)/,$(OBJS))
SOBJ_PATHS  =   $(addprefix $(DIR_SOBJ)/,$(SOBJS))
TARGET_PATH =	$(addprefix $(DIR_BIN)/,$(TARGET))

LIBS     =       $(shell root-config --glibs) -L${MYUTILS}/sobj -lutils -L${ATEID}/sobj -lateid -lgfortran
INCS     =      -I$(shell root-config --incdir) -I./include -I${ATEID} -I${EPCUTS} -I${MYUTILS}/include
ROOTSYS	 =	 $(shell root-config --exec-prefix)
F77      = f95
F_FLAGS  =  
C++      = g++
C_FLAGS  = 
CXXFLAGS =      -O2 -fmessage-length=0 $(shell root-config --cflags) -Wno-deprecated

COMPILER = /usr/bin/gfortran
F77OPT   = -O2 -fno-automatic -w
F77OPT  += -ffixed-line-length-none -fno-second-underscore -DLinux

all: clean ${TARGET}

%.o: %.cpp
	$(CXX) -c $(INCS) $(CXXFLAGS) $< -o $(DIR_OBJ)/$@

elaslib.o: elaslib.f        
	$(COMPILER) $(F77OPT) -c $< -o $(DIR_OBJ)/$@

test: dirs main.o elaslib.o epTree_Looper.o
	$(CXX) -o $(DIR_BIN)/test $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/elaslib.o $(DIR_OBJ)/main.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)

anaCCEC_sim: dirs setup_anaCCEC_sim.o epTree_Looper.o anaCCEC.o 
	$(CXX) -o $(DIR_BIN)/anaCCEC_sim $(DIR_OBJ)/setup_anaCCEC_sim.o $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/anaCCEC.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)

anaCCEC_exp: dirs setup_anaCCEC_exp.o epTree_Looper.o anaCCEC.o 
	$(CXX) -o $(DIR_BIN)/anaCCEC_exp $(DIR_OBJ)/setup_anaCCEC_exp.o $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/anaCCEC.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)
	
verifyeCuts_sim: dirs setup_verifyeCuts_sim.o epTree_Looper.o anaCCEC.o 
	$(CXX) -o $(DIR_BIN)/verifyeCuts_sim $(DIR_OBJ)/setup_verifyeCuts_sim.o $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/anaCCEC.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)

verifyeCuts_exp: dirs setup_verifyeCuts_exp.o epTree_Looper.o anaCCEC.o
	$(CXX) -o $(DIR_BIN)/verifyeCuts_exp $(DIR_OBJ)/setup_verifyeCuts_exp.o $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/anaCCEC.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)

elasYield_th: dirs setup_elasYield_th.o epTree_Looper.o anaCCEC.o 
#	$(CXX) -o $(DIR_BIN)/elasYield_th $(DIR_OBJ)/setup_elasYield_th.o $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/anaCCEC.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)

elasYield_sim: dirs setup_elasYield_sim.o epTree_Looper.o anaCCEC.o 
	$(CXX) -o $(DIR_BIN)/elasYield_sim $(DIR_OBJ)/setup_elasYield_sim.o $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/anaCCEC.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)

elasYield_exp: dirs setup_elasYield_exp.o epTree_Looper.o anaCCEC.o
	$(CXX) -o $(DIR_BIN)/elasYield_exp $(DIR_OBJ)/setup_elasYield_exp.o $(DIR_OBJ)/epTree_Looper.o $(DIR_OBJ)/anaCCEC.o $(SOBJ_PATHS) $(LIBS) $(CXXFLAGS)

dirs:
	-mkdir -p $(DIR_OBJ)
	-mkdir -p $(DIR_BIN)

clean:
	-rm -rf $(DIR_OBJ)/*.o  $(DIR_BIN)/$(TARGET)
