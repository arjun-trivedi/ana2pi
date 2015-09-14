#define epTree_Looper_cxx
#include "epTree_Looper.h"

#include "elaslib.h"

#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THnSparse.h>

#include <fstream>

epTree_Looper::epTree_Looper()
{
}

void epTree_Looper::fill_hTTnorm(float nbins, float xlow, float xhigh){
   gStyle->SetOptStat("neMR");

   float BE=5.499;
   float E1F_TARGET_RAD_LENGTH=0.00562;
   float WCUT=1.028;

   TFile* f_out = new TFile("TT.root", "RECREATE");
   TH1D* hTTnorm=new TH1D("hTTnorm","hTTnorm", nbins, xlow, xhigh);
   TH1D* hTTnorm_norad=new TH1D("hTTnorm_norad","hTTnorm_norad", nbins, xlow, xhigh);
   
   for (int ibin = 0; ibin < nbins; ibin ++){
	   double theta=hTTnorm->GetBinLowEdge(ibin+1);
	   double binc=elasrad(BE,theta,E1F_TARGET_RAD_LENGTH,WCUT);
	   double binc_norad=elas(BE,theta);
	   hTTnorm->SetBinContent(ibin+1,binc);
           hTTnorm_norad->SetBinContent(ibin+1,binc_norad);
           hTTnorm->SetBinError(ibin+1,0);
           hTTnorm_norad->SetBinError(ibin+1,0);
   }
   f_out->Write();
}//end fill_th_elasYield
