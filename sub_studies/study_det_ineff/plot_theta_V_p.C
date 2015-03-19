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
TF1* cut_lw[4][6];
TF1* cut_hg[4][6];
void setup_cuts(){
  for (int i=0;i<4;i++){
    for (int isctr=0;isctr<6;isctr++){
      TString name=TString::Format("%s_%d",prtcl[i].Data(),isctr+1).Data();
      if (prtcl[i]=="e"){
        if (isctr+1==1||isctr+1==2||isctr+1==5||isctr+1==6){//! no cuts
          cut_lw[i][isctr]=new TF1(TString::Format("%s_lw",name.Data()),"0", 0,5);
          cut_hg[i][isctr]=new TF1(TString::Format("%s_hg",name.Data()),"60",0,5);
        }else if (isctr+1==3){
          cut_lw[i][isctr]=new TF1(TString::Format("%s_lw",name.Data()),"30-1*x",2.0,3.5);
          cut_hg[i][isctr]=new TF1(TString::Format("%s_hg",name.Data()),"32-1*x",2.0,3.5);        
        }else if (isctr+1==4){
          cut_lw[i][isctr]=new TF1(TString::Format("%s_lw",name.Data()),"21.0-1*x",2.0,4.5);
          cut_hg[i][isctr]=new TF1(TString::Format("%s_hg",name.Data()),"24.5-1*x",2.0,4.5);
        }
      }
    }
  }
}
//TString prtcl[]={"e","p","pip","pim"};
void plot_theta_V_p(int top=1,bool exp=kTRUE,bool draw_cut=kFALSE){
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
  fout=new TFile(TString::Format("%s/mon_det_ineff/thetaVp_top%d.root",datadir.Data(),top),"RECREATE");

  if (draw_cut){
    setup_cuts();
  }

  //TString prtcl[]={"e","p","pip","pim"};
  TCanvas *c[4];
  TCanvas *clog[4];
  for (int i=0;i<4;i++){
    c[i]=new TCanvas(TString::Format("c%s",prtcl[i].Data()),TString::Format("%s",prtcl[i].Data()));
    c[i]->Divide(3,2);
    clog[i]=new TCanvas(TString::Format("c%s_log",prtcl[i].Data()),TString::Format("%s",prtcl[i].Data()));
    clog[i]->Divide(3,2);
  }

  Int_t sctr_phi_min[]={0, 30, 90,  150, 210, 270};
  Int_t sctr_phi_max[]={0, 90 ,150, 210, 270, 330};

  Int_t theta_min[]={0, 0, 0,  0};
  Int_t theta_max[]={60,60,120,120};
  Int_t p_min[]={1,0,0,0};
  Int_t p_max[]={5,4,3,3};

  gStyle->SetOptStat(0);
  TCut cut_top(TString::Format("top==%d",top));
  for (int i=0;i<4;i++){
    for (int isctr=0;isctr<6;isctr++){
      TCut cut_sctr(TString::Format("sector_%s==%d",prtcl[i].Data(),isctr+1));
      c[i]->cd(isctr+1);
      TString hname=TString::Format("hthetaVp_%s_sctr%d",prtcl[i].Data(),isctr+1);
      TString cmd=TString::Format("theta_%s:p_%s>>%s(100,%d,%d,150,%d,%d)",
                                   prtcl[i].Data(),prtcl[i].Data(),hname.Data(),
                                   p_min[i],p_max[i],theta_min[i],theta_max[i]);
      /*TString cmd=TString::Format("theta_%s:p_%s>>hthetaVp_%s_sctr%d(100,%d,%d,100,%d,%d)",
                                   prtcl[i].Data(),prtcl[i].Data(),prtcl[i].Data(),isctr+1,
                                   p_min[i],p_max[i],theta_min[i],theta_max[i]);*/
      cout <<"cmd="<<cmd.Data()<<endl;
      t->Draw(cmd,cut_sctr&&cut_top,"colz");
      if (draw_cut){
        if (prtcl[i]=="e"){
          cut_lw[i][isctr]->Draw("same");
          cut_hg[i][isctr]->Draw("same");
        }
      }
      TH2F* h=(TH2F*)gDirectory->Get(hname);
      cout<<"h = "<<h->GetName()<<endl;      
      //! Draw histogram in log scale
      TPad* pad=(TPad*)clog[i]->cd(isctr+1);
      pad->SetLogz(); 
      h->Draw("colz");
      if (draw_cut){
        if (prtcl[i]=="e"){
          cut_lw[i][isctr]->Draw("same");
          cut_hg[i][isctr]->Draw("same");
        }
      }
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
