/*#include <TCut.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTree.h>
#include <TF1.h>
#include <TString.h>

#include <iostream>
using namespace std;*/

TString prtcl[]={"e","p","pip","pim"};

Int_t theta_min[]={0, 0, 0,  0};
Int_t theta_max[]={60,60,120,120};
Int_t theta_bins[]={60,60,120,120};
Int_t p_min[]={1,0,0,0};
Int_t p_max[]={5,4,3,3};

Int_t sctr_phi_min[]={0, 30, 90,  150, 210, 270};
Int_t sctr_phi_max[]={0, 90 ,150, 210, 270, 330};
Int_t phi_bins[]={60,60,60,60,60,60};

//TString prtcl[]={"e","p","pip","pim"};
void plot_phi_V_theta(int top=1,bool exp=kTRUE,bool draw_cut=kFALSE){
  cout<<"top="<<top<<endl;
  TFile* fin;
  TFile* fout;
  TString datadir;
  if (exp==kTRUE){
    cout<<"exp"<<endl;
    datadir=gSystem->ExpandPathName("$D2PIDIR_EXP");
  }else{
    cout<<"sim"<<endl;
    datadir=gSystem->ExpandPathName("$D2PIDIR_SIM");
  }
  fin=new TFile(TString::Format("%s/mon_det_ineff/d2piRmon.root",datadir.Data()));
  TTree* t=(TTree*)fin->Get("d2piR/tR");
  //cout<<"Tree name:"<<t->GetTitle()<<endl;
  fout=new TFile(TString::Format("%s/mon_det_ineff/phiVtheta_top%d.root",datadir.Data(),top),"RECREATE");

  TCanvas *c[4];
  TCanvas *clog[4];
  for (int i=0;i<4;i++){
    c[i]=new TCanvas(TString::Format("c%s",prtcl[i].Data()),TString::Format("%s",prtcl[i].Data()));
    c[i]->Divide(3,2);
    clog[i]=new TCanvas(TString::Format("c%s_log",prtcl[i].Data()),TString::Format("%s",prtcl[i].Data()));
    clog[i]->Divide(3,2);
  }

  gStyle->SetOptStat(0);
  TCut cut_top(TString::Format("top==%d",top));
  for (int i=0;i<4;i++){
    for (int isctr=0;isctr<6;isctr++){
      TCut cut_sctr(TString::Format("sector_%s==%d",prtcl[i].Data(),isctr+1));
      c[i]->cd(isctr+1);
      TString hname=TString::Format("hphiVtheta_%s_sctr%d",prtcl[i].Data(),isctr+1);
      TString cmd=TString::Format("phi_%s:theta_%s>>%s(%d,%d,%d,%d,%d,%d)",
                                   prtcl[i].Data(),prtcl[i].Data(),hname.Data(),
                                   theta_bins[i],theta_min[i],theta_max[i],
                                   phi_bins[isctr],sctr_phi_min[isctr],sctr_phi_max[isctr]);
      cout <<"cmd="<<cmd.Data()<<endl;
      t->Draw(cmd,cut_sctr&&cut_top,"colz");
      TH2F* h=(TH2F*)gDirectory->Get(hname);
      cout<<"h = "<<h->GetName()<<endl;     
      //! Draw histogram in log scale
      TPad* pad=(TPad*)clog[i]->cd(isctr+1);
      pad->SetLogz(); 
      h->Draw("colz");
      //! Write hist to file
      fout->Write(h->GetName(),TObject::kWriteDelete);
    }
  }
  fout->cd();
  for (int i=0;i<4;i++){
    c[i]->Write(); 
    clog[i]->Write();
  }
}
