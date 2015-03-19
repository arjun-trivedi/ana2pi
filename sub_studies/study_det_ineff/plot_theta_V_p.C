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
Int_t p_min[]={1,0,0,0};
Int_t p_max[]={5,4,3,3};

TF1** cut_lw[4][6];
TF1** cut_hg[4][6];
void setup_cuts(){
  for (int i=0;i<4;i++){
    Double_t pmin=p_min[i];
    Double_t pmax=p_max[i];
    Int_t thetamin=theta_min[i];
    Int_t thetamax=theta_max[i];
    printf("pmin,pax,thetamin,thetamax=%.1f,%.1f,%d,%d\n",pmin,pmax,thetamin,thetamax);
    for (int isctr=0;isctr<6;isctr++){
      TString name=TString::Format("%s_%d",prtcl[i].Data(),isctr+1).Data();
      if (prtcl[i]=="e"){
        if (isctr+1==1||isctr+1==2||isctr+1==5||isctr+1==6){//! no cuts
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),TString::Format("%d",thetamin),pmin,pmax);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),TString::Format("%d",thetamax),pmin,pmax);
        }else if (isctr+1==3){
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];          
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"30-1*x",2.0,3.5);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"32-1*x",2.0,3.5);        
        }else if (isctr+1==4){
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"21.0-1*x",2.0,4.5);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"24.5-1*x",2.0,4.5);
        }
      }
      if (prtcl[i]=="p"){
        if (isctr+1==1||isctr+1==4||isctr+1==6){//! no cuts
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),TString::Format("%d",thetamin),pmin,pmax);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),TString::Format("%d",thetamax),pmin,pmax);
        }else if (isctr+1==2){
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"16+8*x-1*x*x",0.5,3.0);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"19+8*x-1*x*x",0.5,3.0); 
        }else if (isctr+1==3){
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"5+8*x-1*x*x",0.5,3.5);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"8+8*x-1*x*x",0.5,3.5);
        }else if (isctr+1==5){
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"19+5*x-1*x*x",0.5,3.0);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"22+5*x-1*x*x",0.5,3.0);
        }
      }
      if (prtcl[i]=="pip"){
        if (isctr+1==1||isctr+1==4||isctr+1==5||isctr+1==6){//! no cuts
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),TString::Format("%d",thetamin),pmin,pmax);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),TString::Format("%d",thetamax),pmin,pmax);
        }else if (isctr+1==2){
          cut_lw[i][isctr]=new TF1*[1];
          cut_hg[i][isctr]=new TF1*[1];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"50+105*x-90*x*x",0.03,0.7);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),TString::Format("%d",thetamax),0.03,0.7);
        }else if (isctr+1==3){
          cut_lw[i][isctr]=new TF1*[4];
          cut_hg[i][isctr]=new TF1*[4];
          cut_lw[i][isctr][0]=new TF1(TString::Format("%s_%d_lw",name.Data(),1),"80+37*x-20*x*x",0.0,0.7);
          cut_hg[i][isctr][0]=new TF1(TString::Format("%s_%d_hg",name.Data(),1),"85+40*x-15*x*x",0.0,0.7);
          cut_lw[i][isctr][1]=new TF1(TString::Format("%s_%d_lw",name.Data(),2),"64.5+32*x-20*x*x",0.0,0.7);
          cut_hg[i][isctr][1]=new TF1(TString::Format("%s_%d_hg",name.Data(),2),"68+35*x-18*x*x",0.0,0.7);
          cut_lw[i][isctr][2]=new TF1(TString::Format("%s_%d_lw",name.Data(),3),"15+30*x-7*x*x",0.1,1.5);
          cut_hg[i][isctr][2]=new TF1(TString::Format("%s_%d_hg",name.Data(),3),"20.5+30*x-7*x*x",0.1,1.5);
          cut_lw[i][isctr][3]=new TF1(TString::Format("%s_%d_lw",name.Data(),4),"-13+25*x-5*x*x",0.6,2.5);
          cut_hg[i][isctr][3]=new TF1(TString::Format("%s_%d_hg",name.Data(),4),"-8.0+25*x-5*x*x",0.6,2.5);
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

  /*Int_t theta_min[]={0, 0, 0,  0};
  Int_t theta_max[]={60,60,120,120};
  Int_t p_min[]={1,0,0,0};
  Int_t p_max[]={5,4,3,3};*/

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
        if (prtcl[i]=="e"||prtcl[i]=="p"){
          cut_lw[i][isctr][0]->Draw("same");
          cut_hg[i][isctr][0]->Draw("same");
        }else if (prtcl[i]=="pip"){
          Int_t ncuts=1;
          if (isctr+1==3) ncuts=4;
          for (int j=0;j<ncuts;j++){
            cut_lw[i][isctr][j]->Draw("same");
            cut_hg[i][isctr][j]->Draw("same");
          }
        }
      }
      TH2F* h=(TH2F*)gDirectory->Get(hname);
      cout<<"h = "<<h->GetName()<<endl;     
      if (draw_cut){
        if (prtcl[i]=="e"||prtcl[i]=="p"){
          h->GetListOfFunctions()->Add(cut_lw[i][isctr][0]);
          h->GetListOfFunctions()->Add(cut_hg[i][isctr][0]); 
        }else if (prtcl[i]=="pip"){
          Int_t ncuts=1;
          if (isctr+1==3) ncuts=4;
          for (int j=0;j<ncuts;j++){
            h->GetListOfFunctions()->Add(cut_lw[i][isctr][j]);
            h->GetListOfFunctions()->Add(cut_hg[i][isctr][j]);
          }
        }
      }
      //! Draw histogram in log scale
      TPad* pad=(TPad*)clog[i]->cd(isctr+1);
      pad->SetLogz(); 
      h->Draw("colz");
      if (draw_cut){
        if (prtcl[i]=="e"||prtcl[i]=="p"){
          cut_lw[i][isctr][0]->Draw("same");
          cut_hg[i][isctr][0]->Draw("same");
        }else if (prtcl[i]=="pip"){
          Int_t ncuts=1;
          if (isctr+1==3) ncuts=4;
          for (int j=0;j<ncuts;j++){
            cut_lw[i][isctr][j]->Draw("same");
            cut_hg[i][isctr][j]->Draw("same");
          }
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
