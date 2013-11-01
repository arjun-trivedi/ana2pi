#include <TString.h>

const Int_t nTops    = 5;
const Int_t nVarsets = 3;
TString at_topNames[nTops] = {"1", "2",   "3",   "4",   "1:2:3:4"};
TString vm_topNames[nTops] = {"1", "2:1", "3:1", "4:1", "1:2:3:4"};

Int_t nQ2Wbins_e1fs1    = 8; 
const Int_t nSims_e1fs1 = 6;
TString q2wbng_e1fs1    = "1-1.400-1.500__8-1.600-1.800";
TString simDir_e1fs1[nSims_e1fs1] = 
                          {"/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.4-1.5__1.6-1.8/Q2W__1.4-1.5__1.6-1.8__100M",
                           "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.4-1.5__1.6-1.8/Q2W__1.4-1.5__1.6-1.8__100M__200M",
                           "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.4-1.5__1.6-1.8/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii",
                           "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.4-1.5__1.6-1.8/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii__400M",
                           "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.4-1.5__1.6-1.8/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii__400M__400Mii",
                           "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.4-1.5__1.6-1.8/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii__400M__400Mii__400Miii"
                           };
TString expDir_e1fs1 = ("/data/trivedia/e1f/ana2pi/experiment/e1f.goldenruns.qskim/Q2W__1.4-1.5__1.6-1.8");

Int_t nQ2Wbins_e1fs2    = 24;
const Int_t nSims_e1fs2 = 4;
TString q2wbng_e1fs2    = "1-2.000-2.400__24-1.300-1.900";
TString simDir_e1fs2[nSims_e1fs2] = 
                           {"/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.9-2.5__1.3-1.9/Q2W__1.9-2.5__1.3-1.9__400M",
                            "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.9-2.5__1.3-1.9/Q2W__1.9-2.5__1.3-1.9__400M__200M",
                            "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.9-2.5__1.3-1.9/Q2W__1.9-2.5__1.3-1.9__400M__200M__200Mii",
                            "/data/trivedia/e1f/ana2pi/simulation_2pi/qskim/Q2W__1.9-2.5__1.3-1.9/Q2W__1.9-2.5__1.3-1.9__400M__200M__200Mii__200Miii" 
                           };
TString expDir_e1fs2 = "/data/trivedia/e1f/ana2pi/experiment/e1f.goldenruns.qskim/Q2W__1.9-2.5__1.3-1.9";


