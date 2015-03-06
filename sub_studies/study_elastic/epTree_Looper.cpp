#define epTree_Looper_cxx
#include "epTree_Looper.h"

#include "elaslib.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THnSparse.h>

#include <fstream>

void epTree_Looper::fill_th_elasYield(double nbins, double xlow, double xhigh){
   gStyle->SetOptStat("neMR");

   TFile* f_out = new TFile("th_elasYield.root", "RECREATE");

   TH1F* h_nVtheta;
   TH1F* h_norad_nVtheta;
   TH1F* h_nVthetax;
   TH1F* h_norad_nVthetax;

   double binwidth = (xhigh-xlow)/nbins;
   h_nVtheta         = new TH1F("nVtheta",         "Number of Elastic Reactions(#theta)", nbins, xlow, xhigh);
   h_norad_nVtheta   = new TH1F("norad_nVtheta",   "Number of Elastic Reactions(#theta); no Rad",  nbins, xlow, xhigh);
   h_nVthetax        = new TH1F("nVthetax",        "Number of Elastic Reactions(#theta)", 100, 10, 60);
   h_norad_nVthetax  = new TH1F("norad_nVthetax",  "Number of Elastic Reactions(#theta); no Rad.",  100, 10, 60);

   //Fill th_theta
   for (int iBin = 0; iBin < nbins; iBin ++){
	   double theta = xlow + (iBin*binwidth);
	   double xsec_elas    = diffxsec_elas   (5.499, theta);                  //d[sigma]/d[omega]
	   double xsec_elasrad = diffxsec_elasrad(5.499, theta, 2.5*0.00577, 1.1);//d[sigma]/d[omega]
	   xsec_elas    = xsec_elas    * TMath::Sin((TMath::Pi()/180)*theta) * (2*TMath::Pi());//d[sigma]/d[theta]
	   xsec_elasrad = xsec_elasrad * TMath::Sin((TMath::Pi()/180)*theta) * (2*TMath::Pi());//d[sigma]/d[theta]
	   unsigned int norad_yield = int((xsec_elas    * binwidth * 6250000) + 0.5);//n(theta)
	   unsigned int wrad_yield  = int((xsec_elasrad * binwidth * 6250000) + 0.5);//n(theta)
	   h_nVtheta      ->SetBinContent(iBin+1,wrad_yield);
	   h_norad_nVtheta->SetBinContent(iBin+1,norad_yield);
	   //cout << "theta: norad_Yield = " << theta << " : " << norad_yield << endl;
   	   //cout << "theta: wrad_Yield = " << theta << " : " << wrad_yield << endl;
   }
   //Fill th_thetax
   double nbinsx = h_nVthetax->GetNbinsX();
   double xlowx = h_nVthetax->GetXaxis()->GetXmin();
   double xmaxx = h_nVthetax->GetXaxis()->GetXmax();
   double binwidthx= (xmaxx-xlowx)/nbinsx;
   for (int iBin = 0; iBin < nbinsx; iBin ++){
	   double theta = xlowx + (iBin*binwidthx);
	   double xsec_elas    = diffxsec_elas   (5.479, theta);                  //d[sigma]/d[omega]
	   double xsec_elasrad = diffxsec_elasrad(5.479, theta, 2.5*0.00577, 1.1);//d[sigma]/d[omega]
	   xsec_elas    = xsec_elas    * TMath::Sin((TMath::Pi()/180)*theta) * (2*TMath::Pi());//d[sigma]/d[theta]
	   xsec_elasrad = xsec_elasrad * TMath::Sin((TMath::Pi()/180)*theta) * (2*TMath::Pi());//d[sigma]/d[theta]
	   unsigned int norad_yield = int((xsec_elas    * binwidthx * 6250000) + 0.5);//n(theta)
	   unsigned int wrad_yield  = int((xsec_elasrad * binwidthx * 6250000) + 0.5);//n(theta)
	   h_nVthetax      ->SetBinContent(iBin+1,wrad_yield);
	   h_norad_nVthetax->SetBinContent(iBin+1,norad_yield);
   }
   f_out->Write();
}//end fill_th_elasYield*/
