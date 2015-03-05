//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 22 12:23:30 2012 by ROOT version 5.31/01
// from TTree epTree/Tree containing variables for Electron Identification
// found on file: ep_skim.root
//////////////////////////////////////////////////////////

#ifndef epTree_Looper_h
#define epTree_Looper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <vector>
#include <string>
using namespace std;

class epTree_Looper {
public :
   epTree_Looper();
   virtual ~epTree_Looper();
   void fill_th_elasYield(double nbins, double xlow, double xhigh);
};

#endif

#ifdef epTree_Looper_cxx
epTree_Looper::epTree_Looper()
{
}

epTree_Looper::~epTree_Looper()
{
}

#endif // #ifdef epTree_Looper_cxx
