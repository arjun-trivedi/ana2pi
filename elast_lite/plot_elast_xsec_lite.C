#include <TF1.h>
#include <TH1.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>

void plot_elast_xsec_lite(TString sector_bng="full"){
  if (sector_bng!="full" && sector_bng!="central"){
    printf("sector_bng=full or central only. Exiting.\n");
    return;
  }
  //! Get input data
  TFile* fER=new TFile(TString::Format("%s/lite_vs_old/delastR_lite.root",getenv("DELASTDIR_EXP")));
  TFile* fSR=new TFile(TString::Format("%s/lite_vs_old/delastR_lite.root",getenv("DELASTDIR_SIM")));
  TFile* fST=new TFile(TString::Format("%s/lite_vs_old/delastT_lite.root",getenv("DELASTDIR_SIM")));

  TFile* fTT=new TFile(TString::Format("%s/elast_lite/thrtcl_xsec/TT.root",getenv("ANA2PI")));
  TH1D* hTTnorm=(TH1D*)fTT->Get("hTTnorm");

  //! Prepare theoretical xsec
  double LUM=19.844;// #fb^-1
  double LUM_INVFB_TO_INVMICROB=1000000000.0;
  double DPHI=0;
  if (sector_bng=="full"){
    DPHI=60*TMath::DegToRad();
  }else if (sector_bng=="central"){
    DPHI=2*TMath::DegToRad();
  }
  TH1D* hlumnorm=new TH1D("lumnorm","lumnorm",32,14,46);
  for (int ibin=0;ibin<32;ibin++){
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
  TFile* fout=new TFile(TString::Format("%s/lite_vs_old/lite/%s",getenv("OBSDIR_ELASTIC"),fout_name.Data()),"RECREATE");
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
      hname=TString::Format("hf_s%d",isctr+1);
    }else if (sector_bng=="central"){
      hname=TString::Format("hc_s%d",isctr+1);
    }
    //! 1. Get SR
    //hSR[isctr]=(TH1D*)fSR[LTE]->Get(hname.Data());
    //hSR[isctr]->SetName(TString::Format("hSR_s%d",isctr+1));
    TH1D* h=(TH1D*)fSR->Get(hname.Data());
    hSR[isctr]=(TH1D*)h->Clone(TString::Format("hSR_s%d",isctr+1));
    //! 2. Get ST
    //hST[isctr]=(TH1D*)fST[LTE]->Get(hname.Data());
    //hST[isctr]->SetName(TString::Format("hST_s%d",isctr+1));	
    h=(TH1D*)fST->Get(hname.Data());
    hST[isctr]=(TH1D*)h->Clone(TString::Format("hST_s%d",isctr+1));
    //! 3. Calculate SA
    hSR[isctr]->Sumw2();
    hST[isctr]->Sumw2();
    hSA[isctr]=(TH1D*)hSR[isctr]->Clone(TString::Format("hSA_s%d",isctr+1));
    hSA[isctr]->Divide(hST[isctr]);
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
  //! Use the following TLine in ratio plots
  double xmin=hlumnorm->GetXaxis()->GetXmin();
  double xmax=hlumnorm->GetXaxis()->GetXmax();
  TLine* ln= new TLine(xmin,1,xmax,1);

  TH1D* hrto_EClumnorm_TTnorm[NSCTR];
  TCanvas *crto_EClumnorm_TTnorm=new TCanvas("crto_EClumnorm_TTnorm","crto_EClumnorm_TTnorm");
  crto_EClumnorm_TTnorm->Divide(3,2);
  for (int isctr=0;isctr<NSCTR;isctr++){
    crto_EClumnorm_TTnorm->cd(isctr+1);
    hrto_EClumnorm_TTnorm[isctr]=(TH1D*)hEClumnorm[isctr]->Clone(TString::Format("hrto_EClumnorm_TTnorm_s%d",isctr+1));
    hrto_EClumnorm_TTnorm[isctr]->Divide(hTTnorm);
    hrto_EClumnorm_TTnorm[isctr]->SetMinimum(0);
    hrto_EClumnorm_TTnorm[isctr]->SetMaximum(2);
    hrto_EClumnorm_TTnorm[isctr]->Draw(); 
    ln->Draw("same");
  }
  crto_EClumnorm_TTnorm->Write();	

  TCanvas *cSA=new TCanvas("cSA","cSA");
  cSA->Divide(3,2);
  for (int isctr=0;isctr<NSCTR;isctr++){
    cSA->cd(isctr+1);
    hSA[isctr]->Draw();
  }
  cSA->Write();

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
