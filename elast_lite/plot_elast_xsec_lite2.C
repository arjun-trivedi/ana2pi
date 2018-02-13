#include <TF1.h>
#include <TH1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>

void plot_elast_xsec_lite2(TString obsdir, TString sector_bng="full", TString binw="binw-1"){
  printf("In plot_elast_xsec_lite2(). obsdir=%s,sector_bng=%s\n",obsdir.Data(),sector_bng.Data());
  /*TString obsname(getenv("obsname"));
  printf("obsname=%s\n",obsname.Data());*/

  //! Check to see if input args are OK
  if (sector_bng!="full" && sector_bng!="central"){
    printf("sector_bng=full or central only. Exiting.\n");
    return;
  }
  if (binw!="binw-1" && binw!="binw-05"){
    printf("binw=binw-1 or binw-05 only. Exiting.\n");
    return;
  }
  //! Get input data
  TFile* fER=new TFile(TString::Format("%s/delastER.root",obsdir.Data()));
  TFile* fSR=new TFile(TString::Format("%s/delastSR.root",obsdir.Data()));
  TFile* fST=new TFile(TString::Format("%s/delastST.root",obsdir.Data()));
  TFile* fTT=new TFile(TString::Format("%s/elast_lite/thrtcl_xsec/TT_theta-binw-1.root",getenv("ANA2PI")));
  if (!fER->IsOpen() || !fSR->IsOpen() || !fST->IsOpen() || !fTT->IsOpen()){
    printf("One or more of the input files do not exists. Exiting.\n");
    return;
  }

  //! Prepare normalization histogram 
  double LUM=28.18;// #fb^-1
  double LUM_INVFB_TO_INVMICROB=1000000000.0;
  double DPHI=0;
  if (sector_bng=="full"){
    DPHI=60*TMath::DegToRad();
  }else if (sector_bng=="central"){
    DPHI=2*TMath::DegToRad();
  }
  int nbins;
  if      (binw=="binw-1")nbins=32;
  else if (binw=="binw-05")nbins=64;
  double theta_min=14;
  double theta_max=46;
  TH1D* hlumnorm=new TH1D("lumnorm","lumnorm",nbins,theta_min,theta_max);
  for (int ibin=0;ibin<nbins;ibin++){
    double theta_min=hlumnorm->GetBinLowEdge(ibin+1);
    double theta_max=theta_min+hlumnorm->GetBinWidth(ibin+1);
    double Dcostheta=fabs(TMath::Cos(theta_max*TMath::DegToRad())-TMath::Cos(theta_min*TMath::DegToRad()));
    double Domega=Dcostheta*DPHI;
    hlumnorm->SetBinContent(ibin+1,LUM*LUM_INVFB_TO_INVMICROB*Domega);
    hlumnorm->SetBinError(ibin+1,0);
  }

  //! Output
  TString fout_name="";
  if (sector_bng=="full"){
    fout_name="yield_full_sector.root";
  }else if (sector_bng=="central"){
    fout_name="yield_central_sector.root";
  }
  TFile* fout=new TFile(TString::Format("%s/%s",obsdir.Data(),fout_name.Data()),"RECREATE");

  //! Get histograms
  //! hTTnorm
  TH1D* hTTnorm=(TH1D*)fTT->Get("hTTnorm");
  //! hSR,hST,hSA,hSC,hER,hEC and normalize them as needed
  const int NSCTR=6;
  TH1D* hSR[NSCTR];
  TH1D* hST[NSCTR];
  TH1D* hSA[NSCTR];
  TH1D* hSC[NSCTR];
  TH1D* hER[NSCTR];
  TH1D* hEC[NSCTR];
  TH1D* hEClumnorm[NSCTR];
  for (int isctr=0;isctr<NSCTR;isctr++){
    TString hname="";
    if (sector_bng=="full"){
      hname=TString::Format("delast/hf_s%d",isctr+1);
    }else if (sector_bng=="central"){
      hname=TString::Format("delast/hc_s%d",isctr+1);
    }
    //! 1. Get SR
    //hSR[isctr]=(TH1D*)fSR[LTE]->Get(hname.Data());
    //hSR[isctr]->SetName(TString::Format("hSR_s%d",isctr+1));
    cout <<hname.Data()<<endl;
    TH1D* h=(TH1D*)fSR->Get(hname.Data());
    hSR[isctr]=(TH1D*)h->Clone(TString::Format("hSR_s%d",isctr+1));
    //! 2. Get ST
    //hST[isctr]=(TH1D*)fST[LTE]->Get(hname.Data());
    //hST[isctr]->SetName(TString::Format("hST_s%d",isctr+1));	
    h=(TH1D*)fST->Get(hname.Data());
    hST[isctr]=(TH1D*)h->Clone(TString::Format("hST_s%d",isctr+1));
    //! 3. Calculate SA
    //! [02-13-18] Since in EI's sim's SR data, sectors 4,5,6 have no data (ST is OK),
    //!            for their acceptance, I am getting data from sectors 1,2,3 respectively
    if ((isctr+1)<4){ //!s=1,2,3 => calculate SA
      hSR[isctr]->Sumw2();
      hST[isctr]->Sumw2();
      hSA[isctr]=(TH1D*)hSR[isctr]->Clone(TString::Format("hSA_s%d",isctr+1));
      hSA[isctr]->Divide(hST[isctr]);
    }else{//!s=4,5,6 => Use SA from 1,2,3 respectively
      if (isctr==3){
        hSA[isctr]=(TH1D*)hSA[0]->Clone(TString::Format("hSA_s%d",isctr+1));
      }else if (isctr==4){
        hSA[isctr]=(TH1D*)hSA[1]->Clone(TString::Format("hSA_s%d",isctr+1));
      }else if (isctr==5){
        hSA[isctr]=(TH1D*)hSA[2]->Clone(TString::Format("hSA_s%d",isctr+1));
      }
    }
    if (sector_bng=="full"){
      hSA[isctr]->SetMinimum(0);
      hSA[isctr]->SetMaximum(0.8);
    }else if (sector_bng=="central"){
      hSA[isctr]->SetMinimum(0);
      hSA[isctr]->SetMaximum(1.2);
    }
    //! Calculate SC
    hSC[isctr]=(TH1D*)hSR[isctr]->Clone(TString::Format("hSC_s%d",isctr+1));
    hSC[isctr]->Divide(hSA[isctr]);
    //! 4. Get ER
    //hER[isctr]=(TH1D*)fER[LTE]->Get(hname.Data());
    //hER[isctr]->SetName(TString::Format("hER_s%d",isctr));
    h=(TH1D*)fER->Get(hname.Data());
    hER[isctr]=(TH1D*)h->Clone(TString::Format("hER_s%d",isctr+1));
    //! 5. Calculate EC
    hEC[isctr]=(TH1D*)hER[isctr]->Clone(TString::Format("hEC_s%d",isctr+1));
    hEC[isctr]->Divide(hSA[isctr]); 
    //! 6. Calculate EClumnorm
    hEClumnorm[isctr]=(TH1D*)hEC[isctr]->Clone(TString::Format("hEClumnorm_s%d",isctr+1));
    hEClumnorm[isctr]->Divide(hlumnorm);
  }

  //! Now plot histograms
  gStyle->SetOptStat(0);

  //! Use the following TLine in ratio plots
  double xmin=hlumnorm->GetXaxis()->GetXmin();
  double xmax=hlumnorm->GetXaxis()->GetXmax();
  TLine* ln= new TLine(xmin,1,xmax,1);

  TH1D* hrto_EClumnorm_TTnorm[NSCTR];
  TCanvas *crto_EClumnorm_TTnorm=new TCanvas("crto_EClumnorm_TTnorm","crto_EClumnorm_TTnorm",1000,800);
  crto_EClumnorm_TTnorm->Divide(3,2);
  for (int isctr=0;isctr<NSCTR;isctr++){
    crto_EClumnorm_TTnorm->cd(isctr+1);
    hrto_EClumnorm_TTnorm[isctr]=(TH1D*)hEClumnorm[isctr]->Clone(TString::Format("hrto_EClumnorm_TTnorm_s%d",isctr+1));
    hrto_EClumnorm_TTnorm[isctr]->Divide(hTTnorm);
    hrto_EClumnorm_TTnorm[isctr]->SetMinimum(0);
    hrto_EClumnorm_TTnorm[isctr]->SetMaximum(2);
    //! Axis labels and aesthetics
    //! + Optmized for .pdf (NOT for .png)
    //! x-axis
    hrto_EClumnorm_TTnorm[isctr]->GetXaxis()->SetTitleOffset(0.7); //0.7
    hrto_EClumnorm_TTnorm[isctr]->GetXaxis()->SetLabelSize(0.03); //0.05
    hrto_EClumnorm_TTnorm[isctr]->GetXaxis()->SetTitleSize(0.07); //0.10
    hrto_EClumnorm_TTnorm[isctr]->SetXTitle("#theta [#circ]");
    //! y-axis
    hrto_EClumnorm_TTnorm[isctr]->GetYaxis()->SetTitleOffset(0.78); //0.7
    hrto_EClumnorm_TTnorm[isctr]->GetYaxis()->SetLabelSize(0.03); //0.05
    hrto_EClumnorm_TTnorm[isctr]->GetYaxis()->SetTitleSize(0.06); //0.10
    hrto_EClumnorm_TTnorm[isctr]->SetYTitle("#Delta#sigma/#Delta#Omega^{exp}/#Delta#sigma/#Delta#Omega^{Bosted}");
    hrto_EClumnorm_TTnorm[isctr]->Draw(); 
    /*//! if sector=1, then draw TPavetext which is useful for viewing .png of TCanvas
    if (isctr+1==1){
      TPaveText *pt = new TPaveText(.30,.80,.90,.90,"NDC");
      pt->AddText(TString::Format("%s_%s",obsname.Data(),sector_bng.Data()));
      pt->Draw("same");
    }*/
    ln->Draw("same");
  }
  crto_EClumnorm_TTnorm->Write();	
  crto_EClumnorm_TTnorm->SaveAs(TString::Format("%s/crto_EClumnorm_TTnorm_%s.png",obsdir.Data(),sector_bng.Data()));
  crto_EClumnorm_TTnorm->SaveAs(TString::Format("%s/crto_EClumnorm_TTnorm_%s.pdf",obsdir.Data(),sector_bng.Data()));

  TCanvas *cSA=new TCanvas("cSA","cSA",1000,800);
  cSA->Divide(3,2);
  for (int isctr=0;isctr<NSCTR;isctr++){
    cSA->cd(isctr+1);
    hSA[isctr]->Draw();
    /*if (isctr+1==1){
      TPaveText *pt = new TPaveText(.30,.80,.90,.90,"NDC");
      pt->AddText(TString::Format("%s_%s",obsname.Data(),sector_bng.Data()));
      pt->Draw("same");
    }*/
  }
  cSA->Write();
  cSA->SaveAs(TString::Format("%s/cSA_%s.png",obsdir.Data(),sector_bng.Data()));
  cSA->SaveAs(TString::Format("%s/cSA_%s.pdf",obsdir.Data(),sector_bng.Data()));

  TCanvas *cEC=new TCanvas("cEC","cEC");
  cEC->Divide(3,2);
  for (int isctr=0;isctr<NSCTR;isctr++){
    cEC->cd(isctr+1);
    hEC[isctr]->Draw();
  }
  cEC->Write();

  TCanvas *cER=new TCanvas("cER","cER");
  cER->Divide(3,2);
  for (int isctr=0;isctr<NSCTR;isctr++){
    cER->cd(isctr+1);
    hER[isctr]->Draw();
  }
  cER->Write();

  TCanvas *crto_ER_SR=new TCanvas("crto_ER_SR","crto_ER_SR");
  crto_ER_SR->Divide(3,2);
  TH1D* hrto_ER_SR[NSCTR];
  for (int isctr=0;isctr<NSCTR;isctr++){
    crto_ER_SR->cd(isctr+1);
    TH1D* hERn=(TH1D*)hER[isctr]->DrawNormalized("",1000);
    TH1D* hSRn=(TH1D*)hSR[isctr]->DrawNormalized("sames",1000);
    hrto_ER_SR[isctr]=(TH1D*)hERn->Clone(TString::Format("rto_ER_SR_s%d",isctr+1));
    hrto_ER_SR[isctr]->Divide(hSRn);
    hrto_ER_SR[isctr]->Draw();
    ln->Draw("same");
  }
  crto_ER_SR->Write();

  TCanvas *crto_EC_SC=new TCanvas("crto_EC_SC","crto_EC_SC");
  crto_EC_SC->Divide(3,2);
  TH1D* hrto_EC_SC[NSCTR];
  for (int isctr=0;isctr<NSCTR;isctr++){
    crto_EC_SC->cd(isctr+1);
    TH1D* hECn=(TH1D*)hEC[isctr]->DrawNormalized("",1000);
    TH1D* hSCn=(TH1D*)hSC[isctr]->DrawNormalized("sames",1000);
    hrto_EC_SC[isctr]=(TH1D*)hECn->Clone(TString::Format("rto_EC_SC_s%d",isctr+1));
    hrto_EC_SC[isctr]->Divide(hSCn);
    hrto_EC_SC[isctr]->Draw();
    ln->Draw("same");
  }
  crto_EC_SC->Write();

  fout->Write();
}
