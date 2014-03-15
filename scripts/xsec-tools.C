#include "xsec-tools.h"

#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TKey.h>
#include <THStack.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TPaveStats.h>
#include <TText.h>

#include "myTHnTool.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

Bool_t setup(TString xsectype/*= "vm"*/, TString q2w_bng/*=""*/, bool pol/*=false*/){
  Bool_t ret = kTRUE;

  //gStyle->SetOptStat(0);

  TString pwd = gSystem->pwd();
  cout <<"pwd = "<<pwd << endl;

  TString path;
  if(q2w_bng=="")
  {
    if      (pwd.Contains("1.4-1.5__1.6-1.8"))  _q2w_bng = Q2W_BNG1;    
    else if (pwd.Contains("1.9-2.5__1.3-1.9"))  _q2w_bng = Q2W_BNG2;

    /*if      (pwd.Contains("experiment"))      _dtype_exp = "exp";
    else if (pwd.Contains("simulation"))      _dtype_sim = "sim";*/
    _dtype_exp = "exp";
    _dtype_sim = "sim";

    path=".";
  }else{
    _q2w_bng=q2w_bng;
    printf("q2wbng = %s\n", _q2w_bng.Data());
    if(_q2w_bng.Contains(Q2W_BNG1)){
      path="/data/trivedia/e1f/ana2pi/experiment/e1f.goldenruns.qskim_071513/Q2W__1.4-1.5__1.6-1.8";
    }else if (_q2w_bng.Contains(Q2W_BNG2)){
      path="/data/trivedia/e1f/ana2pi/experiment/e1f.goldenruns.qskim_071513/Q2W__1.9-2.5__1.3-1.9";
    }else{
      printf("Could not determine path\n");
      return kFALSE; 
    }
    _dtype_exp="exp";
    _dtype_sim="sim";
  }

  if (pol) {
    _dtype_exp = "pol__exp";
    _dtype_sim = "pol__sim";
  }
  

  cout << "xsectype = " << xsectype << endl;
  if (xsectype=="at"){
    copy(AT_TOPNAMES,AT_TOPNAMES+5, _topNames);
    _fyexp[0] = TFile::Open(TString::Format("%s/1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[1] = TFile::Open(TString::Format("%s/2__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[2] = TFile::Open(TString::Format("%s/3__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[3] = TFile::Open(TString::Format("%s/4__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[4] = TFile::Open(TString::Format("%s/1:2:3:4__%s__%s.root", path.Data(), _q2w_bng.Data(), _dtype_exp.Data()));
    _fysim[0] = TFile::Open(TString::Format("%s/simdir/1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[1] = TFile::Open(TString::Format("%s/simdir/2__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[2] = TFile::Open(TString::Format("%s/simdir/3__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[3] = TFile::Open(TString::Format("%s/simdir/4__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[4] = TFile::Open(TString::Format("%s/simdir/1:2:3:4__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
  }else if (xsectype=="vm"){
    copy(VM_TOPNAMES,VM_TOPNAMES+5, _topNames);
    _fyexp[0] = TFile::Open(TString::Format("%s/1__%s__%s.root",   path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[1] = TFile::Open(TString::Format("%s/2:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[2] = TFile::Open(TString::Format("%s/3:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[3] = TFile::Open(TString::Format("%s/4:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fyexp[4] = TFile::Open(TString::Format("%s/1:2:3:4__%s__%s.root",  path.Data(),_q2w_bng.Data(), _dtype_exp.Data()));
    _fysim[0] = TFile::Open(TString::Format("%s/simdir/1__%s__%s.root",   path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[1] = TFile::Open(TString::Format("%s/simdir/2:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[2] = TFile::Open(TString::Format("%s/simdir/3:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[3] = TFile::Open(TString::Format("%s/simdir/4:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
    _fysim[4] = TFile::Open(TString::Format("%s/simdir/1:2:3:4__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype_sim.Data()));
  }else {
    printf("Could not determine xsectype\n");
    ret=kFALSE;
  }

  TString q2bng = _q2w_bng.Tokenize("__")->At(0)->GetName();
  _q2min = ((TString)q2bng.Tokenize("-")->At(1)->GetName()).Atof();
  _q2max = ((TString)q2bng.Tokenize("-")->At(2)->GetName()).Atof();
  _dq2 = _q2max - _q2min;
  
  return ret;
}

TH1F* normalizeYield(TH1F* hYW)
{
  TH1F* hYWnorm = (TH1F*)hYW->Clone("hYWnorm");
  hYWnorm->SetTitle( TString::Format("Q2W = %s",_q2w_bng.Data()) );
  hYWnorm->SetYTitle(XSECTITLE.Data());  
  for (int ibin = 0; ibin < hYW->GetNbinsX(); ibin++)
    {
      float wmin = hYW->GetBinLowEdge(ibin+1);
      float dw   = hYW->GetBinWidth(ibin+1);
      float vgflux = getvgflux(wmin,_q2min);
      float factor = 1000000000;
      float norm = LUM*vgflux*_dq2*dw*factor;
      /*printf("[wmin,q2min] = %f:%f\n",wmin,_q2min);
      printf("lum:vgflux:dw:dq2:factor = %f:%f:%f:%f:%.0E\n",
              LUM,vgflux,dw,_dq2,factor);
      printf("norm = %f\n",norm);*/

      float yield     = hYW->GetBinContent(ibin+1);
      float normYield = yield/norm;
      
      hYWnorm->SetBinContent(ibin+1,normYield);
      hYWnorm->SetBinError(ibin+1,0); //tmp till errors are correctly propagated
    }
  return hYWnorm;
}

void plotxsec(TString xsectype/*= "vm"*/, bool ploty/*=kFALSE*/,bool sim/*=kFALSE*/,bool comp/*=kFALSE*/){
  if (setup(xsectype)==kFALSE) return;

  gStyle->SetOptStat(0);

  TString methodB = TString::Format("Yield - Method B {%s}",_topNames[4].Data());
  TString methodA = TString::Format("Yield - Method A");

  TH1F* hYW_Dir[5];
  TH1F* hYW_5D[5];
  TH1F* hYW_TH[5];
  TH1F* hYWnorm_Dir[5];
  TH1F* hYWnorm_5D[5];
  TLegend* l_Dir = new TLegend(0.1,0.7,0.5,0.9);
  TLegend* l_5D  = new TLegend(0.1,0.7,0.5,0.9);
  TString hname_hYW_Dir = TString::Format("hYW_Dir/Varset1/hYW_ACC_CORR");
  TString hname_hYW_5D  = TString::Format("hYW/Varset1/hYW");
  TString hname_hYW_TH  = TString::Format("hYW_Dir/Varset1/hYW_TH");

  TCanvas* c_Dir = new TCanvas(TString::Format("Dir-%s:%s", xsectype.Data(), _q2w_bng.Data()), 
                               TString::Format("Dir-%s:%s", xsectype.Data(), _q2w_bng.Data()) );
  c_Dir->SetGrid(); 
  TCanvas* c_5D = new TCanvas(TString::Format("5D-%s:%s", xsectype.Data(), _q2w_bng.Data()), 
                              TString::Format("5D-%s:%s", xsectype.Data(), _q2w_bng.Data()) );
  c_5D->SetGrid();
  //! For each Top: hYW_Dir,hYW_5D --> hYWnorm_Dir,hYWnorm_5D & then draw
  TFile* fy[5];
  if (sim) memcpy(fy,_fysim,sizeof(_fysim));
  else     memcpy(fy,_fyexp,sizeof(_fyexp));
  for (int iTop=0;iTop<5;iTop++){
    if (fy[iTop]==NULL) continue;
    cout << "Top = " << iTop << endl;
    hYW_Dir[iTop] = (TH1F*)fy[iTop]->Get(hname_hYW_Dir);
    hYW_5D[iTop]  = (TH1F*)fy[iTop]->Get(hname_hYW_5D);
    if(sim) hYW_TH[iTop]  = (TH1F*)fy[iTop]->Get(hname_hYW_TH);
    if (hYW_Dir[iTop]==NULL || hYW_5D[iTop]==NULL) continue;
    hYW_Dir[iTop]->SetMarkerStyle(20+iTop);
    hYW_Dir[iTop]->SetMarkerColor(2+iTop);
    //hYW_Dir[iTop]->SetMaximum(150000);
    hYW_5D[iTop]->SetMarkerStyle(20+iTop); //+5
    hYW_5D[iTop]->SetMarkerColor(2+iTop);  //+5
    //hYW_5D[iTop]->SetMaximum(500000);
    /*if (sim) hYW_TH[iTop]->SetMarkerColor(kFullStar);  
    if (sim) hYW_TH[iTop]->SetMarkerColor(kBlack);*/
    if (sim) hYW_TH[iTop]->SetFillColor(kRed);  
    if (sim) hYW_TH[iTop]->SetFillStyle(3001);
    //Normalize yields to LUM*vgflux*binw_w*binw_q2
    hYWnorm_Dir[iTop] = normalizeYield(hYW_Dir[iTop]);
    hYWnorm_Dir[iTop]->SetMaximum(20);
    hYWnorm_5D[iTop]  = normalizeYield(hYW_5D[iTop]);
    hYWnorm_5D[iTop]->SetMaximum(20);

   if (iTop==0){
     c_5D->cd();
     l_5D->SetHeader( TString::Format("%s - 5D Method", XSECTITLE.Data()) );
     l_5D->AddEntry( hYWnorm_5D[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     hYWnorm_5D[iTop]->SetMinimum(0);
     ploty?hYW_5D[iTop]->Draw("ep"):hYWnorm_5D[iTop]->Draw("ep");
     if (sim)
     {
      hYW_TH[iTop]->Draw("sames");
      l_5D->AddEntry( hYW_TH[iTop], TString::Format("Thrown") );
     } 

     c_Dir->cd();
     l_Dir->SetHeader( TString::Format("%s - Direct Method", XSECTITLE.Data()) );
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     hYWnorm_Dir[iTop]->SetMinimum(0);
     ploty?hYW_Dir[iTop]->Draw("ep"):hYWnorm_Dir[iTop]->Draw("ep");
     if (sim)
     {
      hYW_TH[iTop]->Draw("bar sames");
      l_Dir->AddEntry( hYW_TH[iTop], TString::Format("Thrown") );
     }    
   }else{
     c_5D->cd();
     l_5D->AddEntry( hYWnorm_5D[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     ploty?hYW_5D[iTop]->Draw("ep sames"):hYWnorm_5D[iTop]->Draw("ep sames");
     if (sim)
     {
      hYW_TH[iTop]->Draw("bar sames");
     }

     c_Dir->cd();
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     ploty?hYW_Dir[iTop]->Draw("ep sames"):hYWnorm_Dir[iTop]->Draw("ep sames");
     if (sim)
     {
      hYW_TH[iTop]->Draw("bar sames");
     }
   }
   c_5D->cd();
   l_5D->Draw("same");
   c_Dir->cd();
   l_Dir->Draw("same"); 
  }
  //new - comp 5D & Dir for Top2 & 5
  if (comp){
    TCanvas* c_comp = new TCanvas(TString::Format("5D_Dir_Cpmp-%s", xsectype.Data()), TString::Format("5D & Dir Comp %s", xsectype.Data()));
    TLegend* l_comp = new TLegend(0.1,0.7,0.5,0.9);
    for (int iTop=0;iTop<5;iTop++){
      if (iTop!=1 && iTop!=4) continue;
      hYW_5D[iTop]->SetMarkerStyle(1);
      hYW_5D[iTop]->SetFillStyle(3001+iTop);
      hYW_5D[iTop]->SetFillColor(hYW_5D[iTop]->GetMarkerColor());
      hYW_5D[iTop]->SetLineColor(hYW_5D[iTop]->GetMarkerColor());
      l_comp->AddEntry(hYW_5D[iTop], TString::Format("5D Method-%s", _topNames[iTop].Data()));
      l_comp->AddEntry(hYW_Dir[iTop], TString::Format("Direct Method-%s", _topNames[iTop].Data()));
      if (iTop==1){
        hYW_5D[iTop]->Draw();
        hYW_Dir[iTop]->Draw("sames");
      }else{
        hYW_5D[iTop]->Draw("sames");
        hYW_Dir[iTop]->Draw("sames");
      }
    }
    l_comp->Draw("same");     
  }
}

void plotxsec_dnp(TString xsectype/*= "vm"*/, bool ploty/*=kFALSE*/,bool sim/*=kFALSE*/,bool comp/*=kFALSE*/){
  if (setup(xsectype)==kFALSE) return;

  gStyle->SetOptStat(0);

  TString methodB = TString::Format("Yield - Method B {%s}",_topNames[4].Data());
  TString methodA = TString::Format("Yield - Method A");

  TH1F* hYW[5];
  TH1F* hYWnorm[5];
  TLegend* l = new TLegend(0.1,0.7,0.5,0.9);
  TString hname_hYW  = TString::Format("hYW/Varset1/hYW");
  
  TCanvas* c = new TCanvas("itgyield","itgyield");
  //! For each Top: hYW_Dir,hYW_5D --> hYWnorm_Dir,hYWnorm_5D & then draw
  TFile* fy[5];
  if (sim) memcpy(fy,_fysim,sizeof(_fysim));
  else     memcpy(fy,_fyexp,sizeof(_fyexp));
  for (int iTop=0;iTop<5;iTop++){
    if (iTop!=4) continue;
    if (fy[iTop]==NULL) continue;
    cout << "Top = " << iTop << endl;
    hYW[iTop] = (TH1F*)fy[iTop]->Get(hname_hYW);
    if (hYW[iTop]==NULL) continue;
    /*hYW[iTop]->SetMarkerStyle(20+iTop);
    hYW[iTop]->SetMarkerColor(2+iTop);*/
    //Normalize yields to LUM*vgflux*binw_w*binw_q2
    hYWnorm[iTop] = normalizeYield(hYW[iTop]);
    hYWnorm[iTop]->SetMarkerStyle(kFullCircle);
    //hYWnorm[iTop]->SetMaximum(20);
    
   if (iTop==4){
     c->cd();
     //l->SetHeader("[2.0-2.4]_[1.3-1.9]");
     TPaveText* pt = new TPaveText(0.1,0.75,0.5,0.9,"NDC");
     TText* q2wt = pt->AddText("[2.0-2.4]_[1.3-1.9]");
     q2wt->SetTextColor(kBlue);
     pt->AddText("Normalized Yield");
     //l->AddEntry( hYWnorm[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     hYWnorm[iTop]->SetMinimum(0);
     hYWnorm[iTop]->SetTitle("");
     hYWnorm[iTop]->SetYTitle("Yield [a.u.]");
     hYWnorm[iTop]->Draw("ep");
     pt->Draw();
     //hYWnorm[iTop]->SetMinimum(0);
     //ploty?hYW[iTop]->Draw("ep"):hYWnorm[iTop]->Draw("ep");*/
   }
   TString outdir = "itgyield";
   gSystem->mkdir(outdir,1);
   c->SaveAs(TString::Format("%s/%s.eps", outdir.Data(),c->GetName())); 
  }
    //new - comp 5D & Dir for Top2 & 5
}

void plot1Dxsec(seq_t seq/*=FULL*/){
  if (setup("vm")==kFALSE) return;

  //! COSMETICS for this method
  gStyle->SetOptStat(0);

  //! INPUT DATA for this method
  
  //! OBJECTS for this method
  //* Use any of the input files to get the number of Q2Wbins
  const int kNumQ2Wbins = _fyexp[0]->GetNkeys() - 2; //subtract 2 for hYW_Dir and hYW dirs
  TCanvas* cexp[kNumQ2Wbins][nVARSET];
  TCanvas* csim[kNumQ2Wbins][nVARSET];
  for(int iq2wbin=0;iq2wbin<kNumQ2Wbins;iq2wbin++){
    for(int iVarset=0;iVarset<nVARSET;iVarset++){
      TString cname_exp = TString::Format("chY1Dexp_%02d_%d",iq2wbin+1,iVarset+1);
      cexp[iq2wbin][iVarset] = new TCanvas(cname_exp,cname_exp,800,600);
      cexp[iq2wbin][iVarset]->Divide(2,2);
      TString cname_sim = TString::Format("chY1Dsim_%02d_%d",iq2wbin+1,iVarset+1);
      csim[iq2wbin][iVarset] = new TCanvas(cname_sim,cname_sim,800,600);
      csim[iq2wbin][iVarset]->Divide(2,2);
    }
  }

  //! GET all the OBJECTS from INPUT DATA 
  for (int iTop=0;iTop<5;iTop++){//nTOP loop
    if (_fyexp[iTop]==NULL) continue;
    TIter nextkey(_fyexp[iTop]->GetListOfKeys());
    TKey *key;
    int iq2wbin = 0;
    while (key = (TKey*)nextkey()) {//Q2W loop
      TString Q2Wdirname = key->GetName();
      if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
      printf("iq2wbin,Q2Wdirname = %d,%s\n", iq2wbin, Q2Wdirname.Data());
      for (int iVarset=0; iVarset<nVARSET; iVarset++){//nVARSET loop
        for (int iVar=0; iVar<nVAR;iVar++){//nVAR loop
          if (iVar==ALPHA) continue;
          TString hname = TString::Format("%s/hY1D/Varset%d/hY1D_%s_%s",
                                Q2Wdirname.Data(),iVarset+1,seqTitle[seq].Data(),varName[iVar].Data());
          TString hname_th = TString::Format("%s/hY1D/Varset%d/hY1D_%s_%s",
                                Q2Wdirname.Data(),iVarset+1,"TH",varName[iVar].Data());
          cout << hname << endl;
          TH1F* h1Dexp=(TH1F*)_fyexp[iTop]->Get(hname);
          TH1F* h1Dsim=(TH1F*)_fysim[iTop]->Get(hname);
          TH1F* h1Dsim_th=(TH1F*)_fysim[iTop]->Get(hname_th);
          h1Dexp->SetMarkerStyle(20+iTop);
          h1Dexp->SetMarkerColor(2+iTop);
          h1Dexp->SetLineColor(2+iTop);
          h1Dsim->SetMarkerStyle(20+iTop);
          h1Dsim->SetMarkerColor(2+iTop);
          h1Dsim->SetLineColor(2+iTop);
          h1Dsim_th->SetFillColor(2+iTop);
          h1Dsim_th->SetFillStyle(3001+iTop);
          h1Dsim_th->SetLineColor(kBlack);
          /*TCanvas* ct = new TCanvas("t","t");
          h1Dsim_th->Draw("hist");
          printf("histogram options = %s,%s\n", h1Dexp->GetOption(), h1Dsim->GetOption());*/

          cexp[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==0){
            if (seq==RECO) h1Dexp->SetMaximum(10*h1Dexp->GetMaximum());
            h1Dexp->Draw();
          }else h1Dexp->Draw("same");

          csim[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==0){
            if (seq==RECO) h1Dsim->SetMaximum(10*h1Dsim->GetMaximum());
            h1Dsim->Draw();
            if (seq==FULL) h1Dsim_th->Draw("hist same");
          }else {
            h1Dsim->Draw("same");
            if (seq==FULL) h1Dsim_th->Draw("hist same");
          }
        }//end nVAR loop
      }//end nVARSET loop
      iq2wbin+=1;
    }//end Q2W loop
  }//end nTOP loop

  //! DRAW OBJECTS
  TString outdir_exp = "1Dxsec/exp";
  TString outdir_sim = "1Dxsec/sim";
  gSystem->mkdir(outdir_exp,1);
  gSystem->mkdir(outdir_sim,1);
  for(int iq2wbin=0;iq2wbin<kNumQ2Wbins;iq2wbin++){
    for(int iVarset=0;iVarset<nVARSET;iVarset++){
      TString csavename_exp = TString::Format("%s/h1D_%02d_%d_%s.png",outdir_exp.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      cexp[iq2wbin][iVarset]->SaveAs(csavename_exp);
      cexp[iq2wbin][iVarset]->Close();
      TString csavename_sim = TString::Format("%s/h1D_%02d_%d_%s.png",outdir_sim.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      csim[iq2wbin][iVarset]->SaveAs(csavename_sim);
      csim[iq2wbin][iVarset]->Close();
    }
  }
}

void plot1Dxsec_dnp(seq_t seq/*=FULL*/){
  if (setup("vm")==kFALSE) return;

  //! COSMETICS for this method
  gStyle->SetOptStat(0);

  //! INPUT DATA for this method
  
  //! OBJECTS for this method
  //* Use any of the input files to get the number of Q2Wbins
  const int kNumQ2Wbins = _fyexp[0]->GetNkeys() - 2; //subtract 2 for hYW_Dir and hYW dirs
  TCanvas* c[kNumQ2Wbins][nVARSET];
  /*TCanvas* cexp[kNumQ2Wbins][nVARSET];
  TCanvas* csim[kNumQ2Wbins][nVARSET];*/
  for(int iq2wbin=0;iq2wbin<kNumQ2Wbins;iq2wbin++){
    for(int iVarset=0;iVarset<nVARSET;iVarset++){
      TString cname = TString::Format("chY1D_%02d_%d",iq2wbin+1,iVarset+1);
      c[iq2wbin][iVarset] = new TCanvas(cname,cname,800,600);
      c[iq2wbin][iVarset]->Divide(2,2);
      /*TString cname_exp = TString::Format("chY1Dexp_%02d_%d",iq2wbin+1,iVarset+1);
      cexp[iq2wbin][iVarset] = new TCanvas(cname_exp,cname_exp,800,600);
      cexp[iq2wbin][iVarset]->Divide(2,2);
      TString cname_sim = TString::Format("chY1Dsim_%02d_%d",iq2wbin+1,iVarset+1);
      csim[iq2wbin][iVarset] = new TCanvas(cname_sim,cname_sim,800,600);
      csim[iq2wbin][iVarset]->Divide(2,2);*/
    }
  }
  //TCanvas *cs= TCanvas();

  //! GET all the OBJECTS from INPUT DATA 
  for (int iTop=0;iTop<5;iTop++){//nTOP loop
    if (iTop!=4) continue;
    if (_fyexp[iTop]==NULL) continue;
    TIter nextkey(_fyexp[iTop]->GetListOfKeys());
    TKey *key;
    int iq2wbin = 0;
    while (key = (TKey*)nextkey()) {//Q2W loop
      TString Q2Wdirname = key->GetName();
      if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
      printf("iq2wbin,Q2Wdirname = %d,%s\n", iq2wbin, Q2Wdirname.Data());
      for (int iVarset=0; iVarset<nVARSET; iVarset++){//nVARSET loop
        for (int iVar=0; iVar<nVAR;iVar++){//nVAR loop
          if (iVar==ALPHA) continue;
          TString hname = TString::Format("%s/hY1D/Varset%d/hY1D_%s_%s",
                                Q2Wdirname.Data(),iVarset+1,seqTitle[seq].Data(),varName[iVar].Data());
          TString hname_th = TString::Format("%s/hY1D/Varset%d/hY1D_%s_%s",
                                Q2Wdirname.Data(),iVarset+1,"TH",varName[iVar].Data());
          cout << hname << endl;
          TH1F* h1Dexp=(TH1F*)_fyexp[iTop]->Get(hname);
          TH1F* h1Dsim=(TH1F*)_fysim[iTop]->Get(hname);
          TH1F* h1Dsim_th=(TH1F*)_fysim[iTop]->Get(hname_th);
          h1Dexp->SetMarkerStyle(kFullCircle);
          h1Dexp->SetMarkerColor(kGreen);
          h1Dexp->SetLineColor(kGreen);
          h1Dsim->SetMarkerStyle(kFullCircle);
          h1Dsim->SetMarkerColor(kRed);
          h1Dsim->SetLineColor(kRed);

          h1Dsim_th->SetFillColor(2+iTop);
          h1Dsim_th->SetFillStyle(3001+iTop);
          h1Dsim_th->SetLineColor(kBlack);
          /*TCanvas* ct = new TCanvas("t","t");
          h1Dsim_th->Draw("hist");
          printf("histogram options = %s,%s\n", h1Dexp->GetOption(), h1Dsim->GetOption());*/

          //! Make Titles nice
          TPaveText* pt = new TPaveText(0.5, 0.6, 1.0, 0.75, "NDC");
          TText* q2wt = pt->AddText(TString::Format("[Q^{2}][W] = %s",Q2Wdirname.Data()));
          q2wt->SetTextColor(kBlue);
          TText* vart = pt->AddText(TString::Format("Varset 1: %s,%s,%s,%s",
                                       varTitle[iVarset][0].Data(),varTitle[iVarset][1].Data(),
                                       varTitle[iVarset][2].Data(),varTitle[iVarset][3].Data()));
          vart->SetTextSize(0.05);
          int kNorm = 900;
          c[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==4){
            if (seq==RECO) h1Dexp->SetMaximum(10*h1Dexp->GetMaximum());
            h1Dexp->SetTitle("");
            h1Dsim->SetTitle("");
            h1Dexp->SetTitleSize(0.05);
            h1Dsim->SetTitleSize(0.05);
            h1Dexp->DrawNormalized("X0 E0",kNorm);
            h1Dsim->DrawNormalized("X0 E0 same",kNorm);
            if (iVar+1==1){
              TLegend* l = new TLegend(0.7,0.75,1,0.9);
              l->AddEntry(h1Dexp,"exp","p");
              l->AddEntry(h1Dsim,"sim","p");
              l->Draw("same");
            }
          }else {
            h1Dexp->DrawNormalized("same",kNorm);
            h1Dsim->DrawNormalized("same",kNorm);
          }
          if (iVar+1==1) pt->Draw();
          
          
          /*cexp[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==0){
            if (seq==RECO) h1Dexp->SetMaximum(10*h1Dexp->GetMaximum());
            h1Dexp->Draw();
          }else h1Dexp->Draw("same");

          csim[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==0){
            if (seq==RECO) h1Dsim->SetMaximum(10*h1Dsim->GetMaximum());
            h1Dsim->Draw();
            if (seq==FULL) h1Dsim_th->Draw("hist same");
          }else {
            h1Dsim->Draw("same");
            if (seq==FULL) h1Dsim_th->Draw("hist same");
          }*/
        }//end nVAR loop
      }//end nVARSET loop
      iq2wbin+=1;
    }//end Q2W loop
  }//end nTOP loop

  //! DRAW OBJECTS
  TString outdir = "1Dxsec_dnp";
  gSystem->mkdir(outdir,1);
  /*TString outdir_exp = "1Dxsec_dnp/exp";
  TString outdir_sim = "1Dxsec_dnp/sim";*/
  /*gSystem->mkdir(outdir_exp,1);
  gSystem->mkdir(outdir_sim,1);*/
  for(int iq2wbin=0;iq2wbin<kNumQ2Wbins;iq2wbin++){
    for(int iVarset=0;iVarset<nVARSET;iVarset++){
      TString csavename = TString::Format("%s/h1D_%02d_%d_%s.png",outdir.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      c[iq2wbin][iVarset]->SaveAs(csavename);
      csavename = TString::Format("%s/h1D_%02d_%d_%s.eps",outdir.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      c[iq2wbin][iVarset]->SaveAs(csavename);
      c[iq2wbin][iVarset]->Close();
      /*TString csavename_exp = TString::Format("%s/h1D_%02d_%d_%s.png",outdir_exp.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      cexp[iq2wbin][iVarset]->SaveAs(csavename_exp);
      cexp[iq2wbin][iVarset]->Close();
      TString csavename_sim = TString::Format("%s/h1D_%02d_%d_%s.png",outdir_sim.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      csim[iq2wbin][iVarset]->SaveAs(csavename_sim);
      csim[iq2wbin][iVarset]->Close();*/
    }
  }
}

void plot1Dxsec_prop14(seq_t seq/*=FULL*/){
  if (setup("vm")==kFALSE) return;

  //! COSMETICS for this method
  gStyle->SetOptStat(0);

  //! INPUT DATA for this method
  
  //! OBJECTS for this method
  //* Use any of the input files to get the number of Q2Wbins
  const int kNumQ2Wbins = _fyexp[0]->GetNkeys() - 2; //subtract 2 for hYW_Dir and hYW dirs
  TCanvas* c[kNumQ2Wbins][nVARSET];
  /*TCanvas* cexp[kNumQ2Wbins][nVARSET];
  TCanvas* csim[kNumQ2Wbins][nVARSET];*/
  for(int iq2wbin=0;iq2wbin<kNumQ2Wbins;iq2wbin++){
    for(int iVarset=0;iVarset<nVARSET;iVarset++){
      TString cname = TString::Format("chY1D_%02d_%d",iq2wbin+1,iVarset+1);
      c[iq2wbin][iVarset] = new TCanvas(cname,cname,800,600);
      c[iq2wbin][iVarset]->Divide(2,2);
      /*TString cname_exp = TString::Format("chY1Dexp_%02d_%d",iq2wbin+1,iVarset+1);
      cexp[iq2wbin][iVarset] = new TCanvas(cname_exp,cname_exp,800,600);
      cexp[iq2wbin][iVarset]->Divide(2,2);
      TString cname_sim = TString::Format("chY1Dsim_%02d_%d",iq2wbin+1,iVarset+1);
      csim[iq2wbin][iVarset] = new TCanvas(cname_sim,cname_sim,800,600);
      csim[iq2wbin][iVarset]->Divide(2,2);*/
    }
  }
  //TCanvas *cs= TCanvas();

  //! GET all the OBJECTS from INPUT DATA 
  for (int iTop=0;iTop<5;iTop++){//nTOP loop
    if (iTop!=4) continue;
    if (_fyexp[iTop]==NULL) continue;
    TIter nextkey(_fyexp[iTop]->GetListOfKeys());
    TKey *key;
    int iq2wbin = 0;
    while (key = (TKey*)nextkey()) {//Q2W loop
      TString Q2Wdirname = key->GetName();
      if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
      printf("iq2wbin,Q2Wdirname = %d,%s\n", iq2wbin, Q2Wdirname.Data());
      for (int iVarset=0; iVarset<nVARSET; iVarset++){//nVARSET loop
        for (int iVar=0; iVar<nVAR;iVar++){//nVAR loop
          if (iVar==ALPHA) continue;
          TString hname = TString::Format("%s/hY1D/Varset%d/hY1D_%s_%s",
                                Q2Wdirname.Data(),iVarset+1,seqTitle[seq].Data(),varName[iVar].Data());
          TString hname_th = TString::Format("%s/hY1D/Varset%d/hY1D_%s_%s",
                                Q2Wdirname.Data(),iVarset+1,"TH",varName[iVar].Data());
          cout << hname << endl;
          TH1F* h1Dexp=(TH1F*)_fyexp[iTop]->Get(hname);
          TH1F* h1Dsim=(TH1F*)_fysim[iTop]->Get(hname);
          TH1F* h1Dsim_th=(TH1F*)_fysim[iTop]->Get(hname_th);
          h1Dexp->SetMarkerStyle(kFullCircle);
          h1Dexp->SetMarkerColor(kGreen);
          h1Dexp->SetLineColor(kGreen);
          h1Dsim->SetMarkerStyle(kFullCircle);
          h1Dsim->SetMarkerColor(kRed);
          h1Dsim->SetLineColor(kRed);

          h1Dsim_th->SetFillColor(2+iTop);
          h1Dsim_th->SetFillStyle(3001+iTop);
          h1Dsim_th->SetLineColor(kBlack);
          /*TCanvas* ct = new TCanvas("t","t");
          h1Dsim_th->Draw("hist");
          printf("histogram options = %s,%s\n", h1Dexp->GetOption(), h1Dsim->GetOption());*/

          //! Make Titles nice
          TPaveText* pt = new TPaveText(0.5, 0.6, 1.0, 0.75, "NDC");
          TText* q2wt = pt->AddText(TString::Format("[Q^{2}][W] = %s",Q2Wdirname.Data()));
          q2wt->SetTextColor(kBlue);
          TText* vart = pt->AddText(TString::Format("Varset 1: %s,%s,%s,%s",
                                       varTitle[iVarset][0].Data(),varTitle[iVarset][1].Data(),
                                       varTitle[iVarset][2].Data(),varTitle[iVarset][3].Data()));
          vart->SetTextSize(0.05);
          int kNorm = 900;
          c[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==4){
            if (seq==RECO) h1Dexp->SetMaximum(10*h1Dexp->GetMaximum());
            h1Dexp->SetTitle("");
            h1Dsim->SetTitle("");
            h1Dexp->SetTitleSize(0.05);
            h1Dsim->SetTitleSize(0.05);
            h1Dexp->DrawNormalized("X0 E0",kNorm);
            h1Dsim->DrawNormalized("X0 E0 same",kNorm);
            if (iVar+1==1){
              TLegend* l = new TLegend(0.7,0.75,1,0.9);
              l->AddEntry(h1Dexp,"exp","p");
              l->AddEntry(h1Dsim,"sim","p");
              l->Draw("same");
            }
          }else {
            h1Dexp->DrawNormalized("same",kNorm);
            h1Dsim->DrawNormalized("same",kNorm);
          }
          if (iVar+1==1) pt->Draw();
          
          
          /*cexp[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==0){
            if (seq==RECO) h1Dexp->SetMaximum(10*h1Dexp->GetMaximum());
            h1Dexp->Draw();
          }else h1Dexp->Draw("same");

          csim[iq2wbin][iVarset]->cd(iVar+1);
          if(iTop==0){
            if (seq==RECO) h1Dsim->SetMaximum(10*h1Dsim->GetMaximum());
            h1Dsim->Draw();
            if (seq==FULL) h1Dsim_th->Draw("hist same");
          }else {
            h1Dsim->Draw("same");
            if (seq==FULL) h1Dsim_th->Draw("hist same");
          }*/
        }//end nVAR loop
      }//end nVARSET loop
      iq2wbin+=1;
    }//end Q2W loop
  }//end nTOP loop

  //! DRAW OBJECTS
  TString outdir = "Prop14";
  gSystem->mkdir(outdir,1);
  /*TString outdir_exp = "1Dxsec_dnp/exp";
  TString outdir_sim = "1Dxsec_dnp/sim";*/
  /*gSystem->mkdir(outdir_exp,1);
  gSystem->mkdir(outdir_sim,1);*/
  for(int iq2wbin=0;iq2wbin<kNumQ2Wbins;iq2wbin++){
    for(int iVarset=0;iVarset<nVARSET;iVarset++){
      TString csavename = TString::Format("%s/h1D_%02d_%d_%s.png",outdir.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      c[iq2wbin][iVarset]->SaveAs(csavename);
      csavename = TString::Format("%s/h1D_%02d_%d_%s.eps",outdir.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      c[iq2wbin][iVarset]->SaveAs(csavename);
      c[iq2wbin][iVarset]->Close();
      /*TString csavename_exp = TString::Format("%s/h1D_%02d_%d_%s.png",outdir_exp.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      cexp[iq2wbin][iVarset]->SaveAs(csavename_exp);
      cexp[iq2wbin][iVarset]->Close();
      TString csavename_sim = TString::Format("%s/h1D_%02d_%d_%s.png",outdir_sim.Data(),iq2wbin+1,iVarset+1,seqTitle[seq].Data());
      csim[iq2wbin][iVarset]->SaveAs(csavename_sim);
      csim[iq2wbin][iVarset]->Close();*/
    }
  }
}

void plotcompxsec(TString xsectype/*= "vm"*/, bool ploty/*=kFALSE*/,bool comp/*=kFALSE*/){
  TString methodB = TString::Format("Yield - Method B {%s}",_topNames[4].Data());
  TString methodA = TString::Format("Yield - Method A");

  TH1F* hYW_Dir[5];
  TH1F* hYW_5D[5];
  TH1F* hYWnorm_Dir[5];
  TH1F* hYWnorm_5D[5];
  TLegend* l_Dir = new TLegend(0.1,0.7,0.5,0.9);
  TLegend* l_5D  = new TLegend(0.1,0.7,0.5,0.9);
  TString hname_hYW_Dir  = TString::Format("hYW_Dir/Varset1/hYW_ACC_CORR");
  TString hname_hYW_5D   = TString::Format("hYW/Varset1/hYW");

  TCanvas* c_Dir = new TCanvas(TString::Format("Dir-%s:%s", xsectype.Data(), _q2w_bng.Data()), 
                               TString::Format("Dir-%s:%s", xsectype.Data(), _q2w_bng.Data()) );
  c_Dir->SetGrid(); 
  TCanvas* c_5D = new TCanvas(TString::Format("5D-%s:%s", xsectype.Data(), _q2w_bng.Data()), 
                              TString::Format("5D-%s:%s", xsectype.Data(), _q2w_bng.Data()) );
  c_5D->SetGrid();

  if (setup(xsectype,Q2W_BNG1)==kFALSE) return;

  //! For each Top: hYW_Dir,hYW_5D --> hYWnorm_Dir,hYWnorm_5D & then draw
  for (int iTop=0;iTop<5;iTop++){
    if (iTop!=1 && iTop!=4) continue;
    if (_fyexp[iTop]==NULL) continue;
    cout << "Top = " << iTop << endl;
    hYW_Dir[iTop] = (TH1F*)_fyexp[iTop]->Get(hname_hYW_Dir);
    hYW_5D[iTop]  = (TH1F*)_fyexp[iTop]->Get(hname_hYW_5D);
    if (hYW_Dir[iTop]==NULL || hYW_5D[iTop]==NULL) continue;
    hYW_Dir[iTop]->SetMarkerStyle(20+iTop);
    hYW_Dir[iTop]->SetMarkerColor(2+iTop);
    hYW_Dir[iTop]->SetMaximum(150000);
    hYW_5D[iTop]->SetMarkerStyle(20+iTop); //+5
    hYW_5D[iTop]->SetMarkerColor(2+iTop);  //+5
    hYW_5D[iTop]->SetMaximum(500000);
    //Normalize yields to LUM*vgflux*binw_w*binw_q2
    hYWnorm_Dir[iTop] = normalizeYield(hYW_Dir[iTop]);
    hYWnorm_Dir[iTop]->SetMaximum(10);
    hYWnorm_5D[iTop]  = normalizeYield(hYW_5D[iTop]);
    hYWnorm_5D[iTop]->SetMaximum(10);

   if (iTop==1){
     c_5D->cd();
     l_5D->SetHeader( TString::Format("%s - 5D", XSECTITLE.Data()) );
     l_5D->AddEntry( hYWnorm_5D[iTop], TString::Format("5D-%s", _topNames[iTop].Data()) );
     hYWnorm_5D[iTop]->SetMinimum(0);
     ploty?hYW_5D[iTop]->Draw("p"):hYWnorm_5D[iTop]->Draw("p");

     c_Dir->cd();
     l_Dir->SetHeader( TString::Format("%s - 5D", XSECTITLE.Data()) );
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("Direct-%s", _topNames[iTop].Data()) );
     hYWnorm_Dir[iTop]->SetMinimum(0);
     ploty?hYW_Dir[iTop]->Draw("p"):hYWnorm_Dir[iTop]->Draw("p");
   }else{
     c_5D->cd();
     l_5D->AddEntry( hYWnorm_5D[iTop], TString::Format("5D-%s", _topNames[iTop].Data()) );
     ploty?hYW_5D[iTop]->Draw("p sames"):hYWnorm_5D[iTop]->Draw("p sames");

     c_Dir->cd();
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("Direct-%s", _topNames[iTop].Data()) );
     ploty?hYW_Dir[iTop]->Draw("p sames"):hYWnorm_Dir[iTop]->Draw("p sames");
   }
   c_5D->cd();
   l_5D->Draw("same");
   c_Dir->cd();
   l_Dir->Draw("same"); 
  }

  if (setup(xsectype,Q2W_BNG2)==kFALSE) return;

  //! For each Top: hYW_Dir,hYW_5D --> hYWnorm_Dir,hYWnorm_5D & then draw
  for (int iTop=0;iTop<5;iTop++){
    if (iTop!=1 && iTop!=4) continue;
    if (_fyexp[iTop]==NULL) continue;
    cout << "Top = " << iTop << endl;
    hYW_Dir[iTop] = (TH1F*)_fyexp[iTop]->Get(hname_hYW_Dir);
    hYW_5D[iTop]  = (TH1F*)_fyexp[iTop]->Get(hname_hYW_5D);
    if (hYW_Dir[iTop]==NULL || hYW_5D[iTop]==NULL) continue;
    hYW_Dir[iTop]->SetMarkerStyle(20+iTop);
    hYW_Dir[iTop]->SetMarkerColor(2+iTop);
    hYW_Dir[iTop]->SetMaximum(150000);
    hYW_5D[iTop]->SetMarkerStyle(20+iTop); //+5
    hYW_5D[iTop]->SetMarkerColor(2+iTop);  //+5
    hYW_5D[iTop]->SetMaximum(500000);
    //Normalize yields to LUM*vgflux*binw_w*binw_q2
    hYWnorm_Dir[iTop] = normalizeYield(hYW_Dir[iTop]);
    hYWnorm_Dir[iTop]->SetMaximum(10);
    hYWnorm_5D[iTop]  = normalizeYield(hYW_5D[iTop]);
    hYWnorm_5D[iTop]->SetMaximum(10);

   /*if (iTop==0){
     c_5D->cd();
     l_5D->SetHeader( TString::Format("%s - 5D", XSECTITLE.Data()) );
     l_5D->AddEntry( hYWnorm_5D[iTop], TString::Format("5D-%s", _topNames[iTop].Data()) );
     hYWnorm_5D[iTop]->SetMinimum(0);
     ploty?hYW_5D[iTop]->Draw("p"):hYWnorm_5D[iTop]->Draw("p");

     c_Dir->cd();
     l_Dir->SetHeader( TString::Format("%s - 5D", XSECTITLE.Data()) );
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("Direct-%s", _topNames[iTop].Data()) );
     hYWnorm_Dir[iTop]->SetMinimum(0);
     ploty?hYW_Dir[iTop]->Draw("p"):hYWnorm_Dir[iTop]->Draw("p");
   }else{*/
     c_5D->cd();
     l_5D->AddEntry( hYWnorm_5D[iTop], TString::Format("5D-%s", _topNames[iTop].Data()) );
     ploty?hYW_5D[iTop]->Draw("p sames"):hYWnorm_5D[iTop]->Draw("p sames");

     c_Dir->cd();
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("Direct-%s", _topNames[iTop].Data()) );
     ploty?hYW_Dir[iTop]->Draw("p sames"):hYWnorm_Dir[iTop]->Draw("p sames");
   /*}*/
   /*c_5D->cd();
   l_5D->Draw("same");
   c_Dir->cd();
   l_Dir->Draw("same"); */
  }


  //new - comp 5D & Dir for Top2 & 5
  if (comp){
    TCanvas* c_comp = new TCanvas(TString::Format("5D_Dir_Cpmp-%s", xsectype.Data()), TString::Format("5D & Dir Comp %s", xsectype.Data()));
    TLegend* l_comp = new TLegend(0.1,0.7,0.5,0.9);
    for (int iTop=0;iTop<5;iTop++){
      if (iTop!=1 && iTop!=4) continue;
      hYW_5D[iTop]->SetMarkerStyle(1);
      hYW_5D[iTop]->SetFillStyle(3001+iTop);
      hYW_5D[iTop]->SetFillColor(hYW_5D[iTop]->GetMarkerColor());
      hYW_5D[iTop]->SetLineColor(hYW_5D[iTop]->GetMarkerColor());
      l_comp->AddEntry(hYW_5D[iTop], TString::Format("5D-%s", _topNames[iTop].Data()));
      l_comp->AddEntry(hYW_Dir[iTop], TString::Format("Dir-%s", _topNames[iTop].Data()));
      if (iTop==1){
        hYW_5D[iTop]->Draw();
        hYW_Dir[iTop]->Draw("sames");
      }else{
        hYW_5D[iTop]->Draw("sames");
        hYW_Dir[iTop]->Draw("sames");
      }
    }
    l_comp->Draw("same");     
  }
}

/****************************************************************************************/
/* plotxsec_CommonBins(), integrates the yield in hY5D common to all Tops. ("common-bins")
/* for every Q2Wbin and plots this integrated yield vs Q2Wbin.
/*
/* - seq    = seq_t, Enumeration type defined in constants.h 
/* - sim    = if sim-yields are to be used
/* - Q2Wbin = If this option is specified, only the specified Q2Wbin will be processed AND 
/*            additionally, a linearized version of the hY5D for every Top will be drawn.
/*
/******************************************************************************************/

void plotxsec_CommonBins(seq_t seq/*=ACC_CORR*/, bool sim/*=kFALSE*/, int Q2Wbin/*=0*/){
  if (setup("vm")==kFALSE) return;

  //! Stats. Box option for this function
  if (Q2Wbin>0) gStyle->SetOptStat("nemMrRi");
  else          gStyle->SetOptStat(0);
  
  //!sim or data
  TFile* fy[5];
  if (sim) memcpy(fy,_fysim,sizeof(_fysim));
  else     memcpy(fy,_fyexp,sizeof(_fyexp));
 
  //!data objects for 5D yields for each Top
  THnSparse* hY5D[nTOP];

  //![Q2,W] binning information
  TString wbng   = _q2w_bng.Tokenize("__")->At(1)->GetName();
  TString nwbins = wbng.Tokenize("-")->At(0)->GetName();
  TString wmin   = wbng.Tokenize("-")->At(1)->GetName();
  TString wmax   = wbng.Tokenize("-")->At(2)->GetName();

  //! Some constants needed to define good variable names
  //! hY5D[seq]_all-bins = Yield[seq]_com-bins + Yield[seq]_noncom-bins
  //!   * seq = seq_t defined on constants.h; user specified
  //!   * com-bins    = bins common to all Tops
  //!   * noncom-bins = bins unique to each Top
  //!   * all-bins    = com-bins + noncom-bins
  const int kNumBinGrps = 3;
  enum BinGrps {kAll=0,kCom=1,kNonCom=2};
  TString bingrp_name[kNumBinGrps] ={"all_bins","com_bins","noncom_bins"};
  TString bingrp_title[kNumBinGrps]={"All Bins", "Common Bins", "NonCommon Bins"};

  //! For every Q2w bin in seq: Data Objects per Top and each BinGrp:
  //!    Yield Vs. W 
  //!    Nbins Vs. W ; where entries in each of the Nbins contribute to the Yield
  TH1F* hYvW[nTOP][kNumBinGrps];
  TH1I* hNbinsvW[nTOP][kNumBinGrps];
  
  TLegend* leg[kNumBinGrps];
  for (int ibingrp=0;ibingrp<kNumBinGrps;ibingrp++){
      leg[ibingrp] = new TLegend(0.1,0.7,0.5,0.9);
      leg[ibingrp]->SetHeader(TString::Format("%s Yield(%s)", seqTitle[seq].Data(), bingrp_title[ibingrp].Data()));
  }

  for(int iTop=0;iTop<nTOP;iTop++){
    for (int ibingrp=0;ibingrp<kNumBinGrps;ibingrp++){
      //! First create hYvsW
      TString hname = TString::Format("hYvW_%s", bingrp_name[ibingrp].Data());
      TString htitle = hname;
      hYvW[iTop][ibingrp] = new TH1F(hname,htitle, nwbins.Atoi(), wmin.Atof(), wmax.Atof());
      hYvW[iTop][ibingrp]->SetXTitle("W[GeV]");
      hYvW[iTop][ibingrp]->SetMarkerStyle(20+iTop); //+5
      hYvW[iTop][ibingrp]->SetMarkerColor(2+iTop);
      hYvW[iTop][ibingrp]->SetLineColor(2+iTop);
      leg[ibingrp]->AddEntry( hYvW[iTop][ibingrp], TString::Format("%s", _topNames[iTop].Data()) );

      //! Now create hNbinsVW (not TLegend entries refer ony hYvsW, but should be representative for hNumBinsvsW too)
      hname = TString::Format("hNbinsvsW_%s", bingrp_name[ibingrp].Data());
      htitle = hname;
      hNbinsvW[iTop][ibingrp] = new TH1I(hname,htitle, nwbins.Atoi(), wmin.Atof(), wmax.Atof());
      hNbinsvW[iTop][ibingrp]->SetXTitle("W[GeV]");
      hNbinsvW[iTop][ibingrp]->SetMarkerStyle(24+iTop); //+5
      hNbinsvW[iTop][ibingrp]->SetMarkerColor(2+iTop);
      hNbinsvW[iTop][ibingrp]->SetLineColor(2+iTop);
    }
  }
  
 
  //!Loop over Q2W dirs, get h5Ds and their yields only common bins!
  TIter nextkey(fy[0]->GetListOfKeys());
  TKey *key;
  int counterQ2Wbin=0;
  bool draw=kFALSE;
  if (Q2Wbin>0) draw=kTRUE;

  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    printf("##### Processing %s #####\n", Q2Wdirname.Data());
    counterQ2Wbin+=1;
    
    //!if Q2Wbin option > 0, then process only the specified Q2Wbin
    if (Q2Wbin>0 && (counterQ2Wbin != Q2Wbin)) continue; 

    TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
    TString wlow = wrange.Tokenize(",")->At(0)->GetName();
    wlow.Remove(0,1); //remove "["
    double w = wlow.Atof();
      
    char hname[200];                 
    sprintf(hname, "%s/hY5D/Varset1/hY5D_%s", Q2Wdirname.Data(), seqTitle[seq].Data());
    for(int iTop=0;iTop<nTOP;iTop++){
      hY5D[iTop] = (THnSparse*)fy[iTop]->Get(hname);
    }

    //! Create Data structures 
    float    intg_hY5D[nTOP][kNumBinGrps] = {0.0};
    float intgErr_hY5D[nTOP][kNumBinGrps] = {0.0};
    int     nbins_hY5D[nTOP][kNumBinGrps] = {0};
    
    myTHnTool hntool(kFALSE);

    //! 1. Get integral over all bins
    printf("-----Integral over All Bins-----\n");
    for(int iTop=0;iTop<nTOP;iTop++){
      if (iTop==2 || iTop==3) continue; //not doing for now
      int nbins = hntool.GetIntegral(hY5D[iTop], intg_hY5D[iTop][kAll], intgErr_hY5D[iTop][kAll]);
      
      /*** Reason for why nbins = hntool.GetNbinsNotEq0() rather than nbins = hntool.GetIntegral(): ***/
      /* This is due to a Technical Caveat that affects on the hY5D_HOLE:
              Since HOLE = TH - ACC_CORR, the bins that are set to 0, are also included in THnSparse:GetNbins()
              method. However, the Integral is not even meant to be calculated over these 0-bins. I should 
              workaround this Technical Caveat, but till I implement this workaround, this is the solution to
              get the correct number of bins over which the integral is calculated for hY5D_HOLE.*/
      //nbins_hY5D[iTop][kAll]=nbins;
      nbins_hY5D[iTop][kAll]=hntool.GetNbinsNotEq0(hY5D[iTop]); 

      printf("Top=%d: Yield=%f; nbins=%d\n",iTop+1,intg_hY5D[iTop][kAll], nbins_hY5D[iTop][kAll] );
    }

    //! 2. Get integral over only common-bins
    //NOTE, in hntool.GetIntegralCommonBins(), currently, hY5D[2] & hY5D[3] are not considered
    printf("-----Integral over Common Bins-----\n");
    int num_combins = hntool.GetIntegralCommonBins(
                  hY5D[0],              hY5D[1],              hY5D[2],              hY5D[3],
             intg_hY5D[0][kCom],   intg_hY5D[1][kCom],   intg_hY5D[2][kCom],   intg_hY5D[3][kCom],
          intgErr_hY5D[0][kCom],intgErr_hY5D[1][kCom],intgErr_hY5D[2][kCom],intgErr_hY5D[3][kCom],
          draw);
    nbins_hY5D[0][kCom]=nbins_hY5D[1][kCom]=num_combins;
    printf("Number of common bins = %d\n",num_combins);
    printf("Top=%d: Yield=%f; nbins=%d\n",1,intg_hY5D[0][kCom], nbins_hY5D[0][kCom]);
    printf("Top=%d: Yield=%f; nbins=%d\n",2,intg_hY5D[1][kCom], nbins_hY5D[1][kCom]);
    
    //!2. Get integral over non-common-bins
    printf("-----Integral over NonCommon Bins-----\n");
    for(int iTop=0;iTop<nTOP;iTop++){
      if (iTop==2 || iTop==3) continue; //not doing for now
      intg_hY5D[iTop][kNonCom] = intg_hY5D[iTop][kAll] - intg_hY5D[iTop][kCom];
      intgErr_hY5D[iTop][kNonCom] = TMath::Sqrt(TMath::Power(intgErr_hY5D[iTop][kAll],2)+
                                                TMath::Power(intgErr_hY5D[iTop][kCom],2));
      nbins_hY5D[iTop][kNonCom] = nbins_hY5D[iTop][kAll] - nbins_hY5D[iTop][kCom];
    }
    printf("Top=%d: Yield=%f; nbins=%d\n",1,intg_hY5D[0][kNonCom], nbins_hY5D[0][kNonCom]);
    printf("Top=%d: Yield=%f; nbins=%d\n",2,intg_hY5D[1][kNonCom], nbins_hY5D[1][kNonCom]);

    //! Now fill Histogram Objects
    for(int iTop=0;iTop<nTOP;iTop++){
      if (iTop==2 || iTop==3) continue; //not doing for now
      for (int ibingrp=0; ibingrp<kNumBinGrps;ibingrp++){
        hYvW[iTop][ibingrp]->SetBinContent(hYvW[iTop][ibingrp]->FindBin(w+.005), intg_hY5D[iTop][ibingrp]); //+5Mev = _instrinsic.Wbinw
        hYvW[iTop][ibingrp]->  SetBinError(hYvW[iTop][ibingrp]->FindBin(w+.005), intgErr_hY5D[iTop][ibingrp]); //+5Mev = _instrinsic.Wbinw
        hNbinsvW[iTop][ibingrp]->SetBinContent(hNbinsvW[iTop][ibingrp]->FindBin(w+.005), nbins_hY5D[iTop][ibingrp]);
      }
    }
    printf("##### Done #####\n");
    printf("\n");
  }

  //!Now Draw the hYvW Histogram Objects
  TString cname = TString::Format("%s Yield Analysis", seqTitle[seq].Data());
  TCanvas* cy = new TCanvas(cname, cname);
  cy->Divide(kNumBinGrps,1);
  for (int ibingrp=0; ibingrp<kNumBinGrps;ibingrp++){
    cy->cd(ibingrp+1);
    for(int iTop=0;iTop<nTOP;iTop++){
      if (iTop==2 || iTop==3) continue; //not doing for now
      if (iTop==0) {
        if (seq==HOLE) {
          gPad->SetLogy();
          hYvW[iTop][ibingrp]->SetMinimum(.000001);
        }else{
          hYvW[iTop][ibingrp]->SetMinimum(0);
          //set max of histograms as per [Top=2][kAll]
          float max_t1 = hYvW[0][kAll]->GetBinContent(hYvW[0][kAll]->GetMaximumBin()) + 100;
          float max_t2 = hYvW[1][kAll]->GetBinContent(hYvW[1][kAll]->GetMaximumBin()) + 100;
          if (max_t2>max_t1) hYvW[iTop][ibingrp]->SetMaximum(max_t2);
          //printf("max_t2 = %f\n",max_t2);
        }
        hYvW[iTop][ibingrp]->Draw("p");
      }else {
        hYvW[iTop][ibingrp]->Draw("p sames");
      }
    }
  }

  //}

  //!Now Draw the hNbinsvW Histogram Objects
  /*for(int iTop=0;iTop<nTOP;iTop++){
    if (iTop==2 || iTop==3) continue; //not doing for now*/

    for (int ibingrp=0; ibingrp<kNumBinGrps;ibingrp++){
      cy->cd(ibingrp+1);
      gPad->Update();
      float pad_ymin, pad_ymax;
      TString axopt;
      if (gPad->GetLogy()){
        pad_ymin = TMath::Power(10, gPad->GetUymin());
        pad_ymax = TMath::Power(10, gPad->GetUymax());
        axopt = "+LG"; 
      } else {
        pad_ymin = gPad->GetUymin();
        pad_ymax = gPad->GetUymax();
        axopt = "+L";
      }
      //set rightaxis and scale as per Top whose hNbinsvW maximum is the largest in this bingrp
      vector<int> v_nbins;
      v_nbins.push_back(hNbinsvW[0][ibingrp]->GetMaximum());
      v_nbins.push_back(hNbinsvW[1][ibingrp]->GetMaximum());
      /*v_nbins.push_back(hNbinsvW[2][ibingrp]->GetMaximum());
      v_nbins.push_back(hNbinsvW[3][ibingrp]->GetMaximum());*/
      sort(v_nbins.begin(),v_nbins.end());
      int max_nbins = v_nbins.at(v_nbins.size()-1);

      float rightmax = max_nbins;
      float scale = pad_ymax/rightmax;
      hNbinsvW[0][ibingrp]->Scale(scale);
      hNbinsvW[1][ibingrp]->Scale(scale);
      /*hNbinsvW[2][ibingrp]->Scale(scale);
      hNbinsvW[3][ibingrp]->Scale(scale);*/

      hNbinsvW[0][ibingrp]->Draw("p same");
      hNbinsvW[1][ibingrp]->Draw("p same");
      /*hNbinsvW[2][ibingrp]->Draw("p same");
      hNbinsvW[3][ibingrp]->Draw("p same");*/

      TGaxis *axis = new TGaxis(gPad->GetUxmax(),pad_ymin, gPad->GetUxmax(), pad_ymax,
                                pad_ymin,rightmax,510,axopt);
      axis->SetLineColor(kMagenta);
      axis->SetLabelColor(kMagenta);
      axis->SetTitle("nbins");
      axis->SetTitleColor(kMagenta);
      axis->Draw();
      
    }

  //}

  //!Finally draw the Legend
  for (int ibingrp=0; ibingrp<kNumBinGrps;ibingrp++){
    cy->cd(ibingrp+1);
    leg[ibingrp]->Draw("same");
  }
}

void plotasym(int top){
  if (setup("vm","","pol__")==kFALSE) return;

  THStack* hs = new THStack(TString::Format("asym_top%d",top), TString::Format("asym_top%d",top));

  int itop = top-1;
  TFile* fy = _fyexp[itop];
  TIter nextkey(fy->GetListOfKeys());
  TKey *key;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;
    TString hname = TString::Format("%s/hAsym/Varset1/hAsym_ACC_CORR_phi",Q2Wdirname.Data());
    cout << "hname = " << hname << endl;
    TH1D* h = (TH1D*)fy->Get(hname);
    if (h==NULL) cout << "histogram not found" << endl;
    //h->Draw();
    hs->Add(h,"e1");
  }
  TCanvas *c = new TCanvas(hs->GetName(),hs->GetTitle());
  hs->Draw("pads");

}

void tnail(){
  TPad *sel = (TPad*)gPad->GetSelectedPad();
  int px = gPad->GetEventX();
  int py = gPad->GetEventY();
  if (sel && sel != c  && sel != ct) {
    ct->cd();
    TPad *newpad = (TPad*)sel->Clone();
    ct->GetListOfPrimitives()->Add(newpad);
    newpad->SetPad(0,0,1,1);
    selold_tnail = newpad;
    ct->Update();
    ct->cd();
  }
}

/*void plotr(int top, TString q2w_bng=""){
  if (setup("vm",q2w_bng,"pol__")==kFALSE) return;

  c = new TCanvas("R", "R");
  ct = new TCanvas("ct","ct");
  c->AddExec("tnail","tnail()");

  int itop = top-1;
  TFile* fy = _fyexp[itop];

  int nq2wbins = fy->GetNkeys();
  int nrows = 4;
  int ncols = nq2wbins/nrows;
  c->Divide(nrows, ncols);
  TIter nextkey(fy->GetListOfKeys());
  TKey *key;
  int iq2wbin = 0;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;
    TString hname_pos = TString::Format("%s/hPhi_POS/Varset1/theta/hRvVar",Q2Wdirname.Data());
    TString hname_neg = TString::Format("%s/hPhi_NEG/Varset1/theta/hRvVar",Q2Wdirname.Data());
    cout << "hname_pos = " << hname_pos << endl;
    cout << "hname_neg = " << hname_neg << endl;
    TH1D* hphi_pos = (TH1D*)fy->Get(hname_pos);
    TH1D* hphi_neg = (TH1D*)fy->Get(hname_neg);
    if (hphi_pos==NULL || hphi_neg==NULL) cout << "histogram not found" << endl;

    //cosmetics and display
    hphi_pos->SetLineColor(kRed);
    hphi_neg->SetLineColor(kBlue);
    c->cd(iq2wbin+1);
    hphi_pos->Draw();
    hphi_neg->Draw("sames");
    iq2wbin+=1;
  }
  c->Update();
}*/

//void plotphi(int top, bool phi=kFALSE){
void plotphi(int top,bool sim=kFALSE){
  //! SETUP   
  if (setup("vm","","pol__")==kFALSE) return;
  bool exp = !sim;
  //! COSMETICS
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111); //pcev = 1111
  //gStyle->SetTitleW(1.5);

  //! INPUT data
  int itop = top-1;
  TFile* fy;
  if (sim) fy = _fysim[itop];
  else fy = _fyexp[itop];

  //!OUTPUT data
  TString outdir_root;
  if (sim) outdir_root = "./polobs/sim";
  else outdir_root = "./polobs/exp";

  TF1 *fphi = new TF1("fphi", "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0, 360);
  fphi->SetParameter(0,1);
  fphi->SetParameter(1,10);
  fphi->SetParameter(2,20);
  fphi->SetParameter(3,100);
  fphi->SetParName(0, "A");
  fphi->SetParName(1, "B");
  fphi->SetParName(2, "C");
  fphi->SetParName(3, "hPR");

  TF1 *fsinphi = new TF1("fsinphi", "[0] + [1]*sin(x*TMath::DegToRad())",0, 360);
  fsinphi->SetParameter(0,1);
  fsinphi->SetParameter(1,10);
  fsinphi->SetParName(0, "Offset");
  fsinphi->SetParName(1, "2hPR");
  
  int nq2wBins = fy->GetNkeys();
  
  TIter nextkey(fy->GetListOfKeys());
  TKey *key;
  int iq2wbin = 0;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;

    //! Prepare outdir for Q2W
    TString outdir = TString::Format("%s/top%d/%s/Varset1/theta", outdir_root.Data(), top, Q2Wdirname.Data());
    gSystem->mkdir(outdir,1);

    //! Get W information for this Q2Wbin. NOTE Q2 information already obtained in setup
    TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
    TString wmin_str = wrange.Tokenize(",")->At(0)->GetName();
    TString wmax_str = wrange.Tokenize(",")->At(1)->GetName();
    wmin_str.Remove(TString::kLeading,'[');
    wmax_str.Remove(TString::kTrailing,'[');  
    float wmin = wmin_str.Atof();
    float wmax = wmax_str.Atof();
    float dw = wmax-wmin;
        
    //! Get Normalization for this Q2Wbin
    double vgflux = getvgflux(W,_q2min);
    double factor = 1000000000;
    double norm = LUM*vgflux*_dq2*dw*factor;
    //printf("wmin:wmax=%f:%f",wmin, wmax);


    //! Prepare:
    //!   - hsinphi
    //!   - hnorm
    //!   - hRvVar   

    THnSparse* hY5D = (THnSparse*)fy->Get(TString::Format("%s/hY5D/Varset1/hY5D_FULL",Q2Wdirname.Data()) );

    //if (hsinphi==NULL){
      int nphibins = hY5D->GetAxis(PHI)->GetNbins();
      double phimin = hY5D->GetAxis(PHI)->GetXmin();
      double phimax = hY5D->GetAxis(PHI)->GetXmax();
      TH1D hsinphi("hsinphi", "hsinphi", nphibins, phimin, phimax);
      for (int ibin = 0; ibin < hsinphi.GetNbinsX(); ibin++){
        double phi = hsinphi.GetBinLowEdge(ibin+1) * TMath::DegToRad();
        hsinphi.SetBinContent(ibin+1, TMath::Sin(phi));
        hsinphi.SetBinError(ibin+1, 0);
      }
    //}

    TH1D hnorm (TString::Format("hnorm_%s",Q2Wdirname.Data()),TString::Format("hnorm_%s",Q2Wdirname.Data()), 
                nphibins, phimin, phimax);
    for (int ibin = 0; ibin < hsinphi.GetNbinsX(); ibin++){
        hnorm.SetBinContent(ibin+1, norm);
        hnorm.SetBinError(ibin+1, 0);
     }

    //! Create hRvVar
    TH1F hRvVar[nHEL];
    int nvarbins = hY5D->GetAxis(THETA)->GetNbins();
    float varmin = hY5D->GetAxis(THETA)->GetXmin();
    float varmax = hY5D->GetAxis(THETA)->GetXmax();
    //printf("numphibins:numvarbins = %d, %d\n", nphibins,nvarbins);
    TString xtitle = TString::Format("%s", varTitle[0][THETA].Data());
    TString ytitle = TString::Format("hPR");
    TString title_unp = TString::Format("%s vs %s [h=%s] %s", 
                                        ytitle.Data(), xtitle.Data(), helTitle[UNPOL].Data(), Q2Wdirname.Data());
    TString title_pos = TString::Format("%s vs %s [h=%s] %s", 
                                        ytitle.Data(), xtitle.Data(), helTitle[POS].Data(), Q2Wdirname.Data());
    TString title_neg = TString::Format("%s vs %s [h=%s] %s", 
                                        ytitle.Data(), xtitle.Data(), helTitle[NEG].Data(), Q2Wdirname.Data());
    hRvVar[UNPOL] = TH1F("hRvVar",title_unp, nvarbins, varmin, varmax);
    hRvVar[UNPOL].SetXTitle(xtitle);
    hRvVar[UNPOL].SetYTitle(ytitle);
    hRvVar[UNPOL].SetLineColor(kBlue);
    hRvVar[POS] = TH1F("hRvVar",title_pos, nvarbins, varmin, varmax);
    hRvVar[POS].SetXTitle(xtitle);
    hRvVar[POS].SetYTitle(ytitle);
    hRvVar[POS].SetLineColor(kRed);
    hRvVar[NEG] = TH1F("hRvVar",title_neg, nvarbins, varmin, varmax);
    hRvVar[NEG].SetXTitle(xtitle);
    hRvVar[NEG].SetYTitle(ytitle);
    hRvVar[NEG].SetLineColor(kBlack);

    //! Create Canvas for displaying hphi_diff
    TCanvas cphi_diff = TCanvas("hphi_diff", "hphi_diff", 900, 600);
    cphi_diff.Divide(5,2);

    for (int i=0; i<nvarbins; i++){
      TString hname_unp = TString::Format("%s/hPhi/Varset1/theta/hphi_proj_%02d",Q2Wdirname.Data(),i+1);
      TString hname_pos = TString::Format("%s/hPhi_POS/Varset1/theta/hphi_proj_%02d",Q2Wdirname.Data(),i+1);
      TString hname_neg = TString::Format("%s/hPhi_NEG/Varset1/theta/hphi_proj_%02d",Q2Wdirname.Data(),i+1);
      cout << "hname_unp = " << hname_unp << endl;
      cout << "hname_pos = " << hname_pos << endl;
      cout << "hname_neg = " << hname_neg << endl;
      TH1D* hphi[nHEL];//_unp, hphi_pos, hphi_neg;
      TH1D* hphi_norm[nHEL];//, hphi_pos_norm, hphi_neg_norm;
      TH1D* hphi_norm_sphi[nHEL];//, hphi_pos_norm_sphi, hphi_neg_norm_sphi;
      hphi[UNPOL] = (TH1D*)fy->Get(hname_unp);
      if (exp){
        hphi[POS] = (TH1D*)fy->Get(hname_pos);
        hphi[NEG] = (TH1D*)fy->Get(hname_neg);
      }
      if ( (sim && hphi[UNPOL]==NULL) || (exp && (hphi[UNPOL]==NULL || hphi[POS]==NULL || hphi[NEG]==NULL)) ) cout << "histogram not found" << endl;
      hphi[UNPOL]->SetLineColor(kBlue);
      if (exp) {
        hphi[POS]->SetLineColor(kRed);
        hphi[NEG]->SetLineColor(kBlack);
      }

      hphi_norm[UNPOL] = (TH1D*)hphi[UNPOL]->Clone((TString)hphi[UNPOL]->GetName()+"_norm");
      hphi_norm[UNPOL]->Divide(&hnorm);
      if (exp){
        hphi_norm[POS] = (TH1D*)hphi[POS]->Clone((TString)hphi[POS]->GetName()+"_norm");
        hphi_norm[POS]->Divide(&hnorm);
        hphi_norm[NEG] = (TH1D*)hphi[NEG]->Clone((TString)hphi[NEG]->GetName()+"_norm");
        hphi_norm[NEG]->Divide(&hnorm);
      }

      hphi_norm_sphi[UNPOL] = (TH1D*)hphi_norm[UNPOL]->Clone((TString)hphi_norm[UNPOL]->GetName()+"_sphi");
      hphi_norm_sphi[UNPOL]->Multiply(&hsinphi);
      if (exp){
        hphi_norm_sphi[POS] = (TH1D*)hphi_norm[POS]->Clone((TString)hphi_norm[POS]->GetName()+"_sphi");
        hphi_norm_sphi[POS]->Multiply(&hsinphi);
        hphi_norm_sphi[NEG] = (TH1D*)hphi_norm[NEG]->Clone((TString)hphi_norm[NEG]->GetName()+"_sphi");
        hphi_norm_sphi[NEG]->Multiply(&hsinphi);
      }

      //! Apply Method 2
      double integerr, integ = 0.0;
      //! Note how Underflow and Overflow bins are explicity avoided when calling TH1::IntegralAndError()
      integ = hphi_norm_sphi[UNPOL]->IntegralAndError(1, hphi_norm_sphi[UNPOL]->GetNbinsX(), integerr);
      hRvVar[UNPOL].SetBinContent(i+1, integ/TMath::Pi());
      hRvVar[UNPOL].SetBinError(i+1, integerr/TMath::Pi());
      if (exp){
        integ = hphi_norm_sphi[POS]->IntegralAndError(1, hphi_norm_sphi[POS]->GetNbinsX(), integerr);
        hRvVar[POS].SetBinContent(i+1, integ/TMath::Pi());
        hRvVar[POS].SetBinError(i+1, integerr/TMath::Pi());
        integ = hphi_norm_sphi[NEG]->IntegralAndError(1, hphi_norm_sphi[NEG]->GetNbinsX(), integerr);
        hRvVar[NEG].SetBinContent(i+1, integ/TMath::Pi());
        hRvVar[NEG].SetBinError(i+1, integerr/TMath::Pi());
      }

      //! Apply Method 3
      TH1D* hphi_diff=NULL;
      fsinphi->SetParameters(1,10);
      if (exp){
        hphi_diff = (TH1D*)hphi_norm[POS]->Clone((TString)hphi_norm[UNPOL]->GetName()+"_hdiff");
        hphi_diff->Add(hphi_norm[NEG],-1);
        hphi_diff->SetLineColor(kGreen);
        cphi_diff.cd(i+1);
        hphi_diff->Fit(fsinphi,"");
        hphi_diff->Draw();
        //fsinphi->Draw();
      }

      //! Apply Method 1
      fphi->SetParameters(1,10,20,100);
      TCanvas cfit = TCanvas("cfit","cfit");
      hphi_norm[UNPOL]->Fit(fphi,"");
      if (exp){
        fphi->SetParameters(1,10,20,100);
        hphi_norm[POS]->Fit(fphi,"");
        fphi->SetParameters(1,10,20,100);
        hphi_norm[NEG]->Fit(fphi,"");
      }

      //! COSMETIC and DISPLAY
      //! Modify titles
      TObjArray *tarr;
      char t[100];
      TPaveText *pt[nHEL];

      pt[UNPOL] = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
      TString tunp = hphi[UNPOL]->GetTitle();
      tarr = tunp.Tokenize("|");
      sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
      pt[UNPOL]->AddText(tarr->At(0)->GetName());
      pt[UNPOL]->AddText(t);
      hphi[UNPOL]->SetTitle("");
      hphi_norm[UNPOL]->SetTitle("");
      //cout << "here2" << endl;
      if (exp){
        pt[POS] = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
        TString tpos = hphi[POS]->GetTitle();
        tarr = tpos.Tokenize("|");
        sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
        pt[POS]->AddText(tarr->At(0)->GetName());
        pt[POS]->AddText(t);
        hphi[POS]->SetTitle("");
        hphi_norm[POS]->SetTitle("");

        pt[NEG] = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
        TString tneg = hphi[NEG]->GetTitle();
        tarr = tneg.Tokenize("|");
        sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
        pt[NEG]->AddText(tarr->At(0)->GetName());
        pt[NEG]->AddText(t);
        hphi[NEG]->SetTitle("");
        hphi_norm[NEG]->SetTitle("");
      }

      //! Display Normalized Phi distributions for bh=0/+/-
      TCanvas cphi_norm = TCanvas(hphi_norm[UNPOL]->GetName(),hphi_norm[UNPOL]->GetName(),900, 600);
      cphi_norm.Divide(3,1);
      cphi_norm.cd(1);
      hphi_norm[UNPOL]->Draw();
      gPad->Update();
      TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
      s->SetTextSize(0.05);  
      s->SetX1NDC(0.2);
      s->SetX2NDC(0.8);
      s->SetY1NDC(0.7);
      s->SetY2NDC(0.9);
      s->SetFitFormat(".1e");
      s->Draw();
      gPad->Update();
      pt[UNPOL]->Draw();
      //cout << "here3" << endl;
      if (exp){
        cphi_norm.cd(2);
        hphi_norm[POS]->Draw();
        gPad->Update();
        s = (TPaveStats*) gPad->GetPrimitive("stats");
        s->SetTextSize(0.05);  
        s->SetX1NDC(0.2);
        s->SetX2NDC(0.8);
        s->SetY1NDC(0.7);
        s->SetY2NDC(0.9);
        s->SetFitFormat(".1e");
        s->Draw();
        gPad->Update();
        pt[POS]->Draw();
        cphi_norm.cd(3);
        hphi_norm[NEG]->Draw();
        gPad->Update();
        s = (TPaveStats*) gPad->GetPrimitive("stats");
        s->SetTextSize(0.05);  
        s->SetX1NDC(0.2);
        s->SetX2NDC(0.8);
        s->SetY1NDC(0.7);
        s->SetY2NDC(0.9);
        s->SetFitFormat(".1e");
        s->Draw();
        gPad->Update();
        pt[NEG]->Draw();
      }
      cphi_norm.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cphi_norm.GetName()));
      cphi_norm.Close();

      //Diplay phi distribution for bh+ = bh-
      /*if (exp){
        cphi_diff.cd(i+1);
        hphi_diff->Draw();
      }*/
    }//! end varbin loop
    //! Draw and save hRvVar & hsinphi
    TCanvas csinphi("sinphi","sinphi");
    hsinphi.Draw();
    csinphi.SaveAs(TString::Format("%s/%s.png", outdir.Data(), csinphi.GetName()));

    TCanvas cnorm("norm","norm");
    hnorm.Draw();
    cnorm.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cnorm.GetName()));

    TCanvas cRvVar("RvVar", "RvVar");
    cRvVar.SetGridy();
    TLine l(0,0,180,0);
    if (exp && top==5){
      hRvVar[UNPOL].SetMinimum(-0.01);
      hRvVar[UNPOL].SetMaximum(0.01);
    }
    hRvVar[UNPOL].Draw();
    l.Draw("same");
    if (exp){
      hRvVar[POS].Draw("same");
      hRvVar[NEG].Draw("same");
    }
    cRvVar.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cRvVar.GetName()));

    if (exp){
      cphi_diff.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cphi_diff.GetName()));
      cphi_diff.Close();
    }

    iq2wbin+=1;   
  }//end loop over Q2W bins
}//end plotphi()

void plotphi4(int top,bool sim=kFALSE){
  //! SETUP   
  if (setup("vm","","pol__")==kFALSE) return;
  bool exp = !sim;
  //! COSMETICS
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111); //pcev = 1111
  //gStyle->SetTitleW(1.5);

  //! INPUT data
  int itop = top-1;
  TFile* fy;
  if (sim) fy = _fysim[itop];
  else fy = _fyexp[itop];

  //!OUTPUT data
  TString outdir_root;
  if (sim) outdir_root = "./polobs4/sim";
  else outdir_root = "./polobs4/exp";

  TF1 *fphi = new TF1("fphi", "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0, 360);
  fphi->SetParameter(0,1);
  fphi->SetParameter(1,10);
  fphi->SetParameter(2,20);
  fphi->SetParameter(3,100);
  fphi->SetParName(0, "A");
  fphi->SetParName(1, "B");
  fphi->SetParName(2, "C");
  fphi->SetParName(3, "hPR");

  TF1 *fsinphi = new TF1("fsinphi", "[0] + [1]*sin(x*TMath::DegToRad())",0, 360);
  fsinphi->SetParameter(0,1);
  fsinphi->SetParameter(1,10);
  fsinphi->SetParName(0, "Offset");
  fsinphi->SetParName(1, "2hPR");
  
  int nq2wBins = fy->GetNkeys();
  
  TIter nextkey(fy->GetListOfKeys());
  TKey *key;
  int iq2wbin = 0;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;

    //! Prepare outdir for Q2W
    TString outdir = TString::Format("%s/top%d/%s/Varset1/theta", outdir_root.Data(), top, Q2Wdirname.Data());
    gSystem->mkdir(outdir,1);

    //! Get W information for this Q2Wbin. NOTE Q2 information already obtained in setup
    TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
    TString wmin_str = wrange.Tokenize(",")->At(0)->GetName();
    TString wmax_str = wrange.Tokenize(",")->At(1)->GetName();
    wmin_str.Remove(TString::kLeading,'[');
    wmax_str.Remove(TString::kTrailing,'[');  
    float wmin = wmin_str.Atof();
    float wmax = wmax_str.Atof();
    float dw = wmax-wmin;
        
    //! Get Normalization for this Q2Wbin
    double vgflux = getvgflux(W,_q2min);
    double factor = 1000000000;
    double norm = LUM*vgflux*_dq2*dw*factor;
    //printf("wmin:wmax=%f:%f",wmin, wmax);


    //! Prepare:
    //!   - hsinphi
    //!   - hnorm
    //!   - hRvVar   

    THnSparse* hY5D = (THnSparse*)fy->Get(TString::Format("%s/hY5D/Varset1/hY5D_FULL",Q2Wdirname.Data()) );

    //if (hsinphi==NULL){
      int nphibins = hY5D->GetAxis(PHI)->GetNbins();
      double phimin = hY5D->GetAxis(PHI)->GetXmin();
      double phimax = hY5D->GetAxis(PHI)->GetXmax();
      TH1D hsinphi("hsinphi", "hsinphi", nphibins, phimin, phimax);
      TH1D hmsinphi("hmsinphi", "hmsinphi", nphibins, phimin, phimax);
      for (int ibin = 0; ibin < hsinphi.GetNbinsX(); ibin++){
        double phi = hsinphi.GetBinLowEdge(ibin+1) * TMath::DegToRad();
        hsinphi.SetBinContent(ibin+1, TMath::Sin(phi));
        hsinphi.SetBinError(ibin+1, 0);
        hmsinphi.SetBinContent(ibin+1, -TMath::Sin(phi));
        hmsinphi.SetBinError(ibin+1, 0);
      }
    //}

    TH1D hnorm (TString::Format("hnorm_%s",Q2Wdirname.Data()),TString::Format("hnorm_%s",Q2Wdirname.Data()), 
                nphibins, phimin, phimax);
    for (int ibin = 0; ibin < hsinphi.GetNbinsX(); ibin++){
        hnorm.SetBinContent(ibin+1, norm);
        hnorm.SetBinError(ibin+1, 0);
     }

    //! Create hRvVar
    TH1F hRvVar[nHEL];
    int nvarbins = hY5D->GetAxis(THETA)->GetNbins();
    float varmin = hY5D->GetAxis(THETA)->GetXmin();
    float varmax = hY5D->GetAxis(THETA)->GetXmax();
    //printf("numphibins:numvarbins = %d, %d\n", nphibins,nvarbins);
    TString xtitle = TString::Format("%s", varTitle[0][THETA].Data());
    TString ytitle = TString::Format("hPR");
    TString title_unp = TString::Format("%s vs %s [h=%s] %s", 
                                        ytitle.Data(), xtitle.Data(), helTitle[UNPOL].Data(), Q2Wdirname.Data());
    TString title_pos = TString::Format("%s vs %s [h=%s] %s", 
                                        ytitle.Data(), xtitle.Data(), helTitle[POS].Data(), Q2Wdirname.Data());
    TString title_neg = TString::Format("%s vs %s [h=%s] %s", 
                                        ytitle.Data(), xtitle.Data(), helTitle[NEG].Data(), Q2Wdirname.Data());
    hRvVar[UNPOL] = TH1F("hRvVar",title_unp, nvarbins, varmin, varmax);
    hRvVar[UNPOL].SetXTitle(xtitle);
    hRvVar[UNPOL].SetYTitle(ytitle);
    hRvVar[UNPOL].SetLineColor(kBlue);
    hRvVar[POS] = TH1F("hRvVar",title_pos, nvarbins, varmin, varmax);
    hRvVar[POS].SetXTitle(xtitle);
    hRvVar[POS].SetYTitle(ytitle);
    hRvVar[POS].SetLineColor(kRed);
    hRvVar[NEG] = TH1F("hRvVar",title_neg, nvarbins, varmin, varmax);
    hRvVar[NEG].SetXTitle(xtitle);
    hRvVar[NEG].SetYTitle(ytitle);
    hRvVar[NEG].SetLineColor(kBlack);

    //! Create Canvas for displaying hphi_diff
    TCanvas cphi_diff = TCanvas("hphi_diff", "hphi_diff", 900, 600);
    cphi_diff.Divide(5,2);

    for (int i=0; i<nvarbins; i++){
      TString hname_unp = TString::Format("%s/hPhi/Varset1/theta/hphi_proj_%02d",Q2Wdirname.Data(),i+1);
      TString hname_pos = TString::Format("%s/hPhi_POS/Varset1/theta/hphi_proj_%02d",Q2Wdirname.Data(),i+1);
      TString hname_neg = TString::Format("%s/hPhi_NEG/Varset1/theta/hphi_proj_%02d",Q2Wdirname.Data(),i+1);
      cout << "hname_unp = " << hname_unp << endl;
      cout << "hname_pos = " << hname_pos << endl;
      cout << "hname_neg = " << hname_neg << endl;
      TH1D* hphi[nHEL];//_unp, hphi_pos, hphi_neg;
      TH1D* hphi_norm[nHEL];//, hphi_pos_norm, hphi_neg_norm;
      TH1D* hphi_norm_sphi[nHEL];//, hphi_pos_norm_sphi, hphi_neg_norm_sphi;
      hphi[UNPOL] = (TH1D*)fy->Get(hname_unp);
      if (exp){
        hphi[POS] = (TH1D*)fy->Get(hname_pos);
        hphi[NEG] = (TH1D*)fy->Get(hname_neg);
      }
      if ( (sim && hphi[UNPOL]==NULL) || (exp && (hphi[UNPOL]==NULL || hphi[POS]==NULL || hphi[NEG]==NULL)) ) cout << "histogram not found" << endl;
      hphi[UNPOL]->SetLineColor(kBlue);
      if (exp) {
        hphi[POS]->SetLineColor(kRed);
        hphi[NEG]->SetLineColor(kBlack);
      }

      hphi_norm[UNPOL] = (TH1D*)hphi[UNPOL]->Clone((TString)hphi[UNPOL]->GetName()+"_norm");
      hphi_norm[UNPOL]->Divide(&hnorm);
      if (exp){
        hphi_norm[POS] = (TH1D*)hphi[POS]->Clone((TString)hphi[POS]->GetName()+"_norm");
        hphi_norm[POS]->Divide(&hnorm);
        hphi_norm[NEG] = (TH1D*)hphi[NEG]->Clone((TString)hphi[NEG]->GetName()+"_norm");
        hphi_norm[NEG]->Divide(&hnorm);
      }

      hphi_norm_sphi[UNPOL] = (TH1D*)hphi_norm[UNPOL]->Clone((TString)hphi_norm[UNPOL]->GetName()+"_sphi");
      hphi_norm_sphi[UNPOL]->Multiply(&hmsinphi);
      if (exp){
        hphi_norm_sphi[POS] = (TH1D*)hphi_norm[POS]->Clone((TString)hphi_norm[POS]->GetName()+"_sphi");
        hphi_norm_sphi[POS]->Multiply(&hsinphi);
        hphi_norm_sphi[NEG] = (TH1D*)hphi_norm[NEG]->Clone((TString)hphi_norm[NEG]->GetName()+"_sphi");
        hphi_norm_sphi[NEG]->Multiply(&hmsinphi);
      }

      //! Apply Method 2
      double integerr, integ = 0.0;
      //! Note how Underflow and Overflow bins are explicity avoided when calling TH1::IntegralAndError()
      integ = hphi_norm_sphi[UNPOL]->IntegralAndError(1, hphi_norm_sphi[UNPOL]->GetNbinsX(), integerr);
      hRvVar[UNPOL].SetBinContent(i+1, integ/TMath::Pi());
      hRvVar[UNPOL].SetBinError(i+1, integerr/TMath::Pi());
      if (exp){
        integ = hphi_norm_sphi[POS]->IntegralAndError(1, hphi_norm_sphi[POS]->GetNbinsX(), integerr);
        hRvVar[POS].SetBinContent(i+1, integ/TMath::Pi());
        hRvVar[POS].SetBinError(i+1, integerr/TMath::Pi());
        integ = hphi_norm_sphi[NEG]->IntegralAndError(1, hphi_norm_sphi[NEG]->GetNbinsX(), integerr);
        hRvVar[NEG].SetBinContent(i+1, integ/TMath::Pi());
        hRvVar[NEG].SetBinError(i+1, integerr/TMath::Pi());
      }

      //! Apply Method 3
      TH1D* hphi_diff=NULL;
      fsinphi->SetParameters(1,10);
      if (exp){
        hphi_diff = (TH1D*)hphi_norm[POS]->Clone((TString)hphi_norm[UNPOL]->GetName()+"_hdiff");
        hphi_diff->Add(hphi_norm[NEG],-1);
        hphi_diff->SetLineColor(kGreen);
        cphi_diff.cd(i+1);
        hphi_diff->Fit(fsinphi,"");
        hphi_diff->Draw();
        //fsinphi->Draw();
      }

      //! Apply Method 1
      fphi->SetParameters(1,10,20,100);
      TCanvas cfit = TCanvas("cfit","cfit");
      hphi_norm[UNPOL]->Fit(fphi,"");
      if (exp){
        fphi->SetParameters(1,10,20,100);
        hphi_norm[POS]->Fit(fphi,"");
        fphi->SetParameters(1,10,20,100);
        hphi_norm[NEG]->Fit(fphi,"");
      }

      //! COSMETIC and DISPLAY
      //! Modify titles
      TObjArray *tarr;
      char t[100];
      TPaveText *pt[nHEL];

      pt[UNPOL] = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
      TString tunp = hphi[UNPOL]->GetTitle();
      tarr = tunp.Tokenize("|");
      sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
      pt[UNPOL]->AddText(tarr->At(0)->GetName());
      pt[UNPOL]->AddText(t);
      hphi[UNPOL]->SetTitle("");
      hphi_norm[UNPOL]->SetTitle("");
      //cout << "here2" << endl;
      if (exp){
        pt[POS] = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
        TString tpos = hphi[POS]->GetTitle();
        tarr = tpos.Tokenize("|");
        sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
        pt[POS]->AddText(tarr->At(0)->GetName());
        pt[POS]->AddText(t);
        hphi[POS]->SetTitle("");
        hphi_norm[POS]->SetTitle("");

        pt[NEG] = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
        TString tneg = hphi[NEG]->GetTitle();
        tarr = tneg.Tokenize("|");
        sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
        pt[NEG]->AddText(tarr->At(0)->GetName());
        pt[NEG]->AddText(t);
        hphi[NEG]->SetTitle("");
        hphi_norm[NEG]->SetTitle("");
      }

      //! Display Normalized Phi distributions for bh=0/+/-
      TCanvas cphi_norm = TCanvas(hphi_norm[UNPOL]->GetName(),hphi_norm[UNPOL]->GetName(),900, 600);
      cphi_norm.Divide(3,1);
      cphi_norm.cd(1);
      hphi_norm[UNPOL]->Draw();
      gPad->Update();
      TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
      if (s!=NULL){
        s->SetTextSize(0.05);  
        s->SetX1NDC(0.2);
        s->SetX2NDC(0.8);
        s->SetY1NDC(0.7);
        s->SetY2NDC(0.9);
        s->SetFitFormat(".1e");
        s->Draw();
      }
      gPad->Update();
      pt[UNPOL]->Draw();
      //cout << "here3" << endl;
      if (exp){
        cphi_norm.cd(2);
        hphi_norm[POS]->Draw();
        gPad->Update();
        s = (TPaveStats*) gPad->GetPrimitive("stats");
        if (s!=NULL){
          s->SetTextSize(0.05);  
          s->SetX1NDC(0.2);
          s->SetX2NDC(0.8);
          s->SetY1NDC(0.7);
          s->SetY2NDC(0.9);
          s->SetFitFormat(".1e");
          s->Draw();
        }
        gPad->Update();
        pt[POS]->Draw();
        cphi_norm.cd(3);
        hphi_norm[NEG]->Draw();
        gPad->Update();
        s = (TPaveStats*) gPad->GetPrimitive("stats");
        if (s!=NULL){
          s->SetTextSize(0.05);  
          s->SetX1NDC(0.2);
          s->SetX2NDC(0.8);
          s->SetY1NDC(0.7);
          s->SetY2NDC(0.9);
          s->SetFitFormat(".1e");
          s->Draw();
        }
        gPad->Update();
        pt[NEG]->Draw();
      }
      cphi_norm.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cphi_norm.GetName()));
      cphi_norm.Close();

      //Diplay phi distribution for bh+ = bh-
      /*if (exp){
        cphi_diff.cd(i+1);
        hphi_diff->Draw();
      }*/
    }//! end varbin loop
    //! Draw and save hRvVar & hsinphi
    TCanvas csinphi("sinphi","sinphi");
    hsinphi.Draw();
    csinphi.SaveAs(TString::Format("%s/%s.png", outdir.Data(), csinphi.GetName()));

    TCanvas cmsinphi("msinphi","msinphi");
    hmsinphi.Draw();
    cmsinphi.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cmsinphi.GetName()));

    TCanvas cnorm("norm","norm");
    hnorm.Draw();
    cnorm.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cnorm.GetName()));

    TCanvas cRvVar("RvVar", "RvVar");
    cRvVar.SetGridy();
    TLine l(0,0,180,0);
    if (exp && top==5){
      hRvVar[UNPOL].SetMinimum(-0.01);
      hRvVar[UNPOL].SetMaximum(0.01);
    }
    hRvVar[UNPOL].Draw();
    l.Draw("same");
    if (exp){
      hRvVar[POS].Draw("same");
      hRvVar[NEG].Draw("same");
    }
    cRvVar.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cRvVar.GetName()));

    if (exp){
      cphi_diff.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cphi_diff.GetName()));
      cphi_diff.Close();
    }

    iq2wbin+=1;   
  }//end loop over Q2W bins
}//end plotphi4()

void plotr(int top,bool sim=kFALSE){
  //! SETUP   
  if (setup("vm","","pol__")==kFALSE) return;
  bool exp = !sim;
  //! COSMETICS
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111); //pcev = 1111
  //gStyle->SetTitleW(1.5);

  //! INPUT data
  int itop = top-1;
  TFile* fy;
  if (sim) fy = _fysim[itop];
  else fy = _fyexp[itop];

  //!OUTPUT data
  TString outdir_root;
  if (sim) outdir_root = "./polobs2/sim";
  else outdir_root = "./polobs2/exp";

  int nq2wBins = fy->GetNkeys();
  
  TIter nextkey(fy->GetListOfKeys());
  TKey *key;
  int iq2wbin = 0;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;

    //! PREPARE outdir for Q2W
    TString outdir = TString::Format("%s/top%d/%s/Varset1/theta", outdir_root.Data(), top, Q2Wdirname.Data());
    gSystem->mkdir(outdir,1);

    //! GET W information for this Q2Wbin. NOTE Q2 information already obtained in setup
    TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
    TString wmin_str = wrange.Tokenize(",")->At(0)->GetName();
    TString wmax_str = wrange.Tokenize(",")->At(1)->GetName();
    wmin_str.Remove(TString::kLeading,'[');
    wmax_str.Remove(TString::kTrailing,'[');  
    float wmin = wmin_str.Atof();
    float wmax = wmax_str.Atof();
    float dw = wmax-wmin;
        
    //! GET Normalization for this Q2Wbin
    double vgflux = getvgflux(W,_q2min);
    double factor = 1000000000;
    double norm = LUM*vgflux*_dq2*dw*factor;
    //printf("wmin:wmax=%f:%f",wmin, wmax);


    //! GET needed hY5Ds for this Q2W bin   
    THnSparse* hY5D[nSEQ][nHEL];
    //! 1. for extracting h=POS R (hack method)
    hY5D[FULL][POS]     = (THnSparse*)fy->Get(TString::Format("%s/hY5D_POS/Varset1/hY5D_FULL",Q2Wdirname.Data()) );
    hY5D[ACC_CORR][NEG] = (THnSparse*)fy->Get(TString::Format("%s/hY5D_NEG/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    //! 2. for extracting h=UNPOL R
    hY5D[FULL][UNPOL]     = (THnSparse*)fy->Get(TString::Format("%s/hY5D/Varset1/hY5D_FULL",Q2Wdirname.Data()) );
    
    //! for extracting h=POS R (hack method!)
    //! 1. MULTIPLY hY5D[FULL][POS] by sin(phi)
    //! 2. MULTIPLY hY5D[ACC_CORR][NEG] by -sin(phi)
    //! 3. ADD hY5D[FULL][POS]*sin(phi) & hY5D[ACC_CORR][NEG]*(-sin(phi))
    //! 4. PROJECT on to THETA
    int nbins = hY5D[FULL][POS]->GetNbins();
    TAxis* axphi = hY5D[FULL][POS]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc = hY5D[FULL][POS]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[FULL][POS]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[FULL][POS]->SetBinContent(ibin, TMath::Sin(phi)*binc);
      hY5D[FULL][POS]->SetBinError(ibin, TMath::Sin(phi)*binerr);
    }
    nbins = hY5D[ACC_CORR][NEG]->GetNbins();
    axphi = hY5D[ACC_CORR][NEG]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc = hY5D[ACC_CORR][NEG]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[ACC_CORR][NEG]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[ACC_CORR][NEG]->SetBinContent(ibin, -TMath::Sin(phi)*binc);
      hY5D[ACC_CORR][NEG]->SetBinError(ibin, -TMath::Sin(phi)*binerr);
    }
    THnSparse* hY5D_pos_sphi = (THnSparse*)hY5D[FULL][POS]->Clone("hY5D_pos_sphi");
    hY5D_pos_sphi->Add(hY5D[ACC_CORR][NEG],1);
    TH1D* hRvVar_pos = (TH1D*)hY5D_pos_sphi->Projection(THETA, "E");
    hRvVar_pos->Scale(1/TMath::Pi());
    hRvVar_pos->Scale(1/norm);
    hRvVar_pos->SetLineColor(kRed);
    
    //! for extracing h=UNPOL R
    //! 1. MULTIPLY hY5D[FULL][UNPOL] by sin(phi)
    //! 2. PROJECT on to THETA
    nbins = hY5D[FULL][UNPOL]->GetNbins();
    axphi = hY5D[FULL][UNPOL]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc = hY5D[FULL][UNPOL]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[FULL][UNPOL]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[FULL][UNPOL]->SetBinContent(ibin, TMath::Sin(phi)*binc);
      hY5D[FULL][UNPOL]->SetBinError(ibin, TMath::Sin(phi)*binerr);
    }
    TH1D* hRvVar_unpol = (TH1D*)hY5D[FULL][UNPOL]->Projection(THETA, "E");
    hRvVar_unpol->Scale(1/TMath::Pi());
    hRvVar_unpol->Scale(1/norm);
    hRvVar_unpol->SetLineColor(kBlue);
    
    //! DRAW & SAVE hRvVar_pos, hRvVar_unpol
    TCanvas cRvVar("RvVar", "RvVar");
    cRvVar.SetGridy();
    TLine l(0,0,180,0);
    if (exp && top==5){
      hRvVar_unpol->SetMinimum(-0.01);
      hRvVar_unpol->SetMaximum(0.01);
    }
    hRvVar_unpol->Draw();
    l.Draw("same");
    if (exp){
      hRvVar_pos->Draw("same");
    }
    cRvVar.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cRvVar.GetName()));

    iq2wbin+=1;   
  }//end loop over Q2W bins
}//end plotr()

void plotr3(int top,bool sim=kFALSE){
  //! SETUP   
  if (setup("vm","","pol__")==kFALSE) return;
  bool exp = !sim;
  myTHnTool hntool(kFALSE);
  //! COSMETICS
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111); //pcev = 1111
  //gStyle->SetTitleW(1.5);

  //! INPUT data
  int itop = top-1;
  TFile* fyexp;
  TFile* fysim;
  fyexp = _fyexp[itop];
  fysim = _fysim[itop];

  //!OUTPUT data
  TString outdir_exp_root = "./polobs3/exp";
  TString outdir_sim_root = "./polobs3/sim";
  
  int nq2wBins = fyexp->GetNkeys();
  
  TIter nextkey(fyexp->GetListOfKeys());
  TKey *key;
  int iq2wbin = 0;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;

    //! PREPARE outdir for Q2W
    TString outdir_exp = TString::Format("%s/top%d/%s/Varset1/theta", outdir_exp_root.Data(), top, Q2Wdirname.Data());
    TString outdir_sim = TString::Format("%s/top%d/%s/Varset1/theta", outdir_sim_root.Data(), top, Q2Wdirname.Data());
    gSystem->mkdir(outdir_exp,1);
    gSystem->mkdir(outdir_sim,1);

    //! GET W information for this Q2Wbin. NOTE Q2 information already obtained in setup
    TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
    TString wmin_str = wrange.Tokenize(",")->At(0)->GetName();
    TString wmax_str = wrange.Tokenize(",")->At(1)->GetName();
    wmin_str.Remove(TString::kLeading,'[');
    wmax_str.Remove(TString::kTrailing,'[');  
    float wmin = wmin_str.Atof();
    float wmax = wmax_str.Atof();
    float dw = wmax-wmin;
        
    //! GET Normalization for this Q2Wbin
    double vgflux = getvgflux(W,_q2min);
    double factor = 1000000000;
    double norm = LUM*vgflux*_dq2*dw*factor;
    //printf("wmin:wmax=%f:%f",wmin, wmax);


    //! GET needed hY5Ds for this Q2W bin   
    const int kDTypes = 2;
    enum {kExp,kSim};
    THnSparse* hY5D[kDTypes][nSEQ][nHEL];
    //! 1. for extracting R_lt from h=POS/NEG/UNP  Dtype= exp/sim data
    hY5D[kExp][ACC_CORR][POS]  = (THnSparse*)fyexp->Get(TString::Format("%s/hY5D_POS/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    hY5D[kExp][ACC_CORR][NEG]  = (THnSparse*)fyexp->Get(TString::Format("%s/hY5D_NEG/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    hY5D[kExp][ACC_CORR][UNPOL] = (THnSparse*)fyexp->Get(TString::Format("%s/hY5D/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    /*hY5D[kExp][ACC_CORR][UNPOL] = (THnSparse*)hY5D[kExp][ACC_CORR][POS]->Clone("unpol");
    hY5D[kExp][ACC_CORR][UNPOL]->Add(hY5D[kExp][ACC_CORR][NEG],1);*/

    hY5D[kSim][ACC_CORR][UNPOL] = (THnSparse*)fysim->Get(TString::Format("%s/hY5D/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    //! OUTPUT objects for this Q2W bin
    TH1D* hRvVar[kDTypes][nHEL];
    TH1D* hRvVar_avg;// avg of hRvVar[kExp][POS] & hRvVar[kExp][NEG]
    
    //! Extract R_lt from [kExp][ACC_CORR][POS] data and compare with [kSim][ACC_CORR][UNPOL](ER PS only)
    //! 1. First make a clone of hY5D[kSim][ACC_CORR][UNPOL] 
    //! 2. Set contents in non-ER PS of hY5D[kSim][ACC_CORR][UNPOL] to zero
    //! 3. MULTIPLY hY5D[kExp][ACC_CORR][POS] & hY5D[kSim][ACC_CORR][UNPOL] by sin(phi) 
    //! 4. PROJECT  hY5D[kExp][ACC_CORR][POS] & hY5D[kSim][ACC_CORR][UNPOL] on to THETA
    hY5D[kSim][ACC_CORR][POS] = (THnSparse*)hY5D[kSim][ACC_CORR][UNPOL]->Clone("pos");
    hntool.SetBinContentsToZero(hY5D[kSim][ACC_CORR][POS], hY5D[kExp][ACC_CORR][POS]);

    int nbins    = hY5D[kExp][ACC_CORR][POS]->GetNbins();
    TAxis* axphi = hY5D[kExp][ACC_CORR][POS]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc   = hY5D[kExp][ACC_CORR][POS]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[kExp][ACC_CORR][POS]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[kExp][ACC_CORR][POS]->SetBinContent(ibin, TMath::Sin(phi)*binc);
      hY5D[kExp][ACC_CORR][POS]->SetBinError(ibin, TMath::Sin(phi)*binerr);
    }
    nbins = hY5D[kSim][ACC_CORR][POS]->GetNbins();
    axphi = hY5D[kSim][ACC_CORR][POS]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc   = hY5D[kSim][ACC_CORR][POS]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[kSim][ACC_CORR][POS]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[kSim][ACC_CORR][POS]->SetBinContent(ibin, TMath::Sin(phi)*binc);
      hY5D[kSim][ACC_CORR][POS]->SetBinError(ibin, TMath::Sin(phi)*binerr);
    }
    hRvVar[kExp][POS] = (TH1D*)hY5D[kExp][ACC_CORR][POS]->Projection(THETA, "E");
    hRvVar[kExp][POS]->Scale(1/TMath::Pi());
    hRvVar[kExp][POS]->Scale(1/norm);
    hRvVar[kExp][POS]->SetLineColor(kRed);
    hRvVar[kSim][POS] = (TH1D*)hY5D[kSim][ACC_CORR][POS]->Projection(THETA, "E");
    hRvVar[kSim][POS]->Scale(1/TMath::Pi());
    hRvVar[kSim][POS]->Scale(1/norm);
    hRvVar[kSim][POS]->SetLineColor(kRed);

    //! Extract R_lt from [kExp][ACC_CORR][NEG] data and compare with [kSim][ACC_CORR][UNPOL](ER PS only)
    //! 1. First make a clone of hY5D[kSim][ACC_CORR][UNPOL] 
    //! 2. Set contents in non-ER PS of hY5D[kSim][ACC_CORR][UNPOL] to zero
    //! 3. MULTIPLY hY5D[kExp][ACC_CORR][NEG] & hY5D[kSim][ACC_CORR][UNPOL] by sin(phi) 
    //! 4. PROJECT  hY5D[kExp][ACC_CORR][NEG] & hY5D[kSim][ACC_CORR][UNPOL] on to THETA
    hY5D[kSim][ACC_CORR][NEG] = (THnSparse*)hY5D[kSim][ACC_CORR][UNPOL]->Clone("neg");
    hntool.SetBinContentsToZero(hY5D[kSim][ACC_CORR][NEG], hY5D[kExp][ACC_CORR][NEG]);

    nbins = hY5D[kExp][ACC_CORR][NEG]->GetNbins();
    axphi = hY5D[kExp][ACC_CORR][NEG]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc   = hY5D[kExp][ACC_CORR][NEG]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[kExp][ACC_CORR][NEG]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[kExp][ACC_CORR][NEG]->SetBinContent(ibin, -TMath::Sin(phi)*binc);
      hY5D[kExp][ACC_CORR][NEG]->SetBinError(ibin, -TMath::Sin(phi)*binerr);
    }
    nbins = hY5D[kSim][ACC_CORR][NEG]->GetNbins();
    axphi = hY5D[kSim][ACC_CORR][NEG]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc   = hY5D[kSim][ACC_CORR][NEG]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[kSim][ACC_CORR][NEG]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[kSim][ACC_CORR][NEG]->SetBinContent(ibin, -TMath::Sin(phi)*binc);
      hY5D[kSim][ACC_CORR][NEG]->SetBinError(ibin, -TMath::Sin(phi)*binerr);
    }
    hRvVar[kExp][NEG] = (TH1D*)hY5D[kExp][ACC_CORR][NEG]->Projection(THETA, "E");
    hRvVar[kExp][NEG]->Scale(1/TMath::Pi());
    hRvVar[kExp][NEG]->Scale(1/norm);
    hRvVar[kExp][NEG]->SetLineColor(kBlack);
    hRvVar[kSim][NEG] = (TH1D*)hY5D[kSim][ACC_CORR][NEG]->Projection(THETA, "E");
    hRvVar[kSim][NEG]->Scale(1/TMath::Pi());
    hRvVar[kSim][NEG]->Scale(1/norm);
    hRvVar[kSim][NEG]->SetLineColor(kBlack);

    //!Compute hRvVar_avg;
    hRvVar_avg = (TH1D*)hRvVar[kExp][POS]->Clone("avg");
    hRvVar_avg->Add(hRvVar[kExp][NEG],1);
    hRvVar_avg->Scale(0.5);
    hRvVar_avg->SetLineColor(kGreen);
    hRvVar_avg->SetMarkerStyle(kFullCircle);

    //! Extract R_lt from [kExp][ACC_CORR][UNPOL] data 
    //! 1. MULTIPLY hY5D[kExp][ACC_CORR][UNPOL] by sin(phi) 
    //! 4. PROJECT  hY5D[kExp][ACC_CORR][UNPOL] on to THETA
    nbins = hY5D[kExp][ACC_CORR][UNPOL]->GetNbins();
    axphi = hY5D[kExp][ACC_CORR][UNPOL]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc   = hY5D[kExp][ACC_CORR][UNPOL]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[kExp][ACC_CORR][UNPOL]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[kExp][ACC_CORR][UNPOL]->SetBinContent(ibin, -TMath::Sin(phi)*binc);
      hY5D[kExp][ACC_CORR][UNPOL]->SetBinError(ibin, -TMath::Sin(phi)*binerr);
    }
    hRvVar[kExp][UNPOL] = (TH1D*)hY5D[kExp][ACC_CORR][UNPOL]->Projection(THETA, "E");
    hRvVar[kExp][UNPOL]->Scale(1/TMath::Pi());
    hRvVar[kExp][UNPOL]->Scale(1/norm);
    hRvVar[kExp][UNPOL]->SetLineColor(kBlue);
        
    //! DRAW & SAVE hRvVar_exp, hRvVar_sim
    TLine l(0,0,180,0);
    TCanvas cRvVar_exp("RvVar", "RvVar");
    cRvVar_exp.SetGridy();
    /*hRvVar[kExp][POS]->SetMinimum(-0.008);
    hRvVar[kExp][POS]->SetMaximum(0.008);
    hRvVar[kExp][POS]->Draw();
    hRvVar[kExp][NEG]->Draw("same");*/
    hRvVar_avg->SetMinimum(-0.003);
    hRvVar_avg->SetMaximum(0.003);
    hRvVar_avg->Draw("ep");
    //hRvVar[kExp][UNPOL]->Draw("same");
    l.Draw("same");
    cRvVar_exp.SaveAs(TString::Format("%s/%s.png", outdir_exp.Data(), cRvVar_exp.GetName()));

    TCanvas cRvVar_sim("RvVar", "RvVar");
    cRvVar_sim.SetGridy();
    hRvVar[kSim][POS]->SetMinimum(-1);
    hRvVar[kSim][POS]->SetMaximum(1);
    hRvVar[kSim][POS]->Draw();
    hRvVar[kSim][NEG]->Draw("same");
    l.Draw("same");
    cRvVar_sim.SaveAs(TString::Format("%s/%s.png", outdir_sim.Data(), cRvVar_sim.GetName()));

    iq2wbin+=1;   
  }//end loop over Q2W bins
}//end plotr3()

void plotr3_dnp(int top,bool sim=kFALSE){
  //! SETUP   
  if (setup("vm","","pol__")==kFALSE) return;
  bool exp = !sim;
  myTHnTool hntool(kFALSE);
  //! COSMETICS
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111); //pcev = 1111
  //gStyle->SetTitleW(1.5);

  //! INPUT data
  int itop = top-1;
  TFile* fyexp;
  TFile* fysim;
  fyexp = _fyexp[itop];
  fysim = _fysim[itop];

  //!OUTPUT data
  TString outdir_root = "./polobs3_dnp";
    
  int nq2wBins = fyexp->GetNkeys();
  
  TIter nextkey(fyexp->GetListOfKeys());
  TKey *key;
  int iq2wbin = 0;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;

    //! PREPARE outdir for Q2W
    TString outdir = TString::Format("%s/top%d/%s/Varset1/theta", outdir_root.Data(), top, Q2Wdirname.Data());
    TString outdir_eps = TString::Format("%s/top%d/%d/Varset1/theta", outdir_root.Data(), top, iq2wbin+1);
    //TString outdir_sim = TString::Format("%s/top%d/%s/Varset1/theta", outdir_sim_root.Data(), top, Q2Wdirname.Data());
    gSystem->mkdir(outdir,1);
    gSystem->mkdir(outdir_eps,1);
    //gSystem->mkdir(outdir_sim,1);

    //! GET W information for this Q2Wbin. NOTE Q2 information already obtained in setup
    TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
    TString wmin_str = wrange.Tokenize(",")->At(0)->GetName();
    TString wmax_str = wrange.Tokenize(",")->At(1)->GetName();
    wmin_str.Remove(TString::kLeading,'[');
    wmax_str.Remove(TString::kTrailing,'[');  
    float wmin = wmin_str.Atof();
    float wmax = wmax_str.Atof();
    float dw = wmax-wmin;
        
    //! GET Normalization for this Q2Wbin
    /*double vgflux = getvgflux(W,_q2min);
    double factor = 1000000000;
    double norm = LUM*vgflux*_dq2*dw*factor;
    cout << "norm =" << norm << endl;*/
    double norm = 50000; //at 11-12-13: arbitrary norm
    //printf("wmin:wmax=%f:%f",wmin, wmax);


    //! GET needed hY5Ds for this Q2W bin   
    const int kDTypes = 2;
    enum {kExp,kSim};
    THnSparse* hY5D[kDTypes][nSEQ][nHEL];
    //! 1. for extracting R_lt from h=POS/NEG/UNP  Dtype= exp/sim data
    hY5D[kExp][ACC_CORR][POS]  = (THnSparse*)fyexp->Get(TString::Format("%s/hY5D_POS/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    hY5D[kExp][ACC_CORR][NEG]  = (THnSparse*)fyexp->Get(TString::Format("%s/hY5D_NEG/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    //hY5D[kExp][ACC_CORR][UNPOL] = (THnSparse*)fyexp->Get(TString::Format("%s/hY5D/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    /*hY5D[kExp][ACC_CORR][UNPOL] = (THnSparse*)hY5D[kExp][ACC_CORR][POS]->Clone("unpol");
    hY5D[kExp][ACC_CORR][UNPOL]->Add(hY5D[kExp][ACC_CORR][NEG],1);*/

    //hY5D[kSim][ACC_CORR][UNPOL] = (THnSparse*)fysim->Get(TString::Format("%s/hY5D/Varset1/hY5D_ACC_CORR",Q2Wdirname.Data()) );
    //! OUTPUT objects for this Q2W bin
    TH1D* hRvVar[kDTypes][nHEL];
    TH1D* hRvVar_avg;// avg of hRvVar[kExp][POS] & hRvVar[kExp][NEG]
    
    //! Extract R_lt from [kExp][ACC_CORR][POS] 
    //! 1. MULTIPLY hY5D[kExp][ACC_CORR][POS] by sin(phi) 
    //! 2. PROJECT  hY5D[kExp][ACC_CORR][POS] on to THETA
    int nbins    = hY5D[kExp][ACC_CORR][POS]->GetNbins();
    TAxis* axphi = hY5D[kExp][ACC_CORR][POS]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc   = hY5D[kExp][ACC_CORR][POS]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[kExp][ACC_CORR][POS]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[kExp][ACC_CORR][POS]->SetBinContent(ibin, TMath::Sin(phi)*binc);
      hY5D[kExp][ACC_CORR][POS]->SetBinError(ibin, TMath::Sin(phi)*binerr);
    }
    hRvVar[kExp][POS] = (TH1D*)hY5D[kExp][ACC_CORR][POS]->Projection(THETA, "E");
    hRvVar[kExp][POS]->Scale(1/TMath::Pi());
    hRvVar[kExp][POS]->Scale(1/norm);
    hRvVar[kExp][POS]->SetLineColor(kRed);
    
    //! Extract R_lt from [kExp][ACC_CORR][NEG] data
    //! 1. MULTIPLY hY5D[kExp][ACC_CORR][NEG] by sin(phi) 
    //! 2. PROJECT  hY5D[kExp][ACC_CORR][NEG] on to THETA
    nbins = hY5D[kExp][ACC_CORR][NEG]->GetNbins();
    axphi = hY5D[kExp][ACC_CORR][NEG]->GetAxis(PHI);
    for (int ibin=0;ibin<nbins;ibin++){
      int bincoord[5] = {0};
      float binc   = hY5D[kExp][ACC_CORR][NEG]->GetBinContent(ibin, bincoord);
      float binerr = hY5D[kExp][ACC_CORR][NEG]->GetBinError(ibin);
      double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
      hY5D[kExp][ACC_CORR][NEG]->SetBinContent(ibin, -TMath::Sin(phi)*binc);
      hY5D[kExp][ACC_CORR][NEG]->SetBinError(ibin, -TMath::Sin(phi)*binerr);
    }
    hRvVar[kExp][NEG] = (TH1D*)hY5D[kExp][ACC_CORR][NEG]->Projection(THETA, "E");
    hRvVar[kExp][NEG]->Scale(1/TMath::Pi());
    hRvVar[kExp][NEG]->Scale(1/norm);
    hRvVar[kExp][NEG]->SetLineColor(kBlack);
    
    //!Compute hRvVar_avg;
    hRvVar_avg = (TH1D*)hRvVar[kExp][POS]->Clone("avg");
    hRvVar_avg->Add(hRvVar[kExp][NEG],1);
    hRvVar_avg->Scale(0.5);
    hRvVar_avg->SetLineColor(kMagenta);
    hRvVar_avg->SetMarkerStyle(kFullCircle);
    hRvVar_avg->SetTitle(TString::Format("D^{%s} vs. %s",varTitle[0][THETA].Data(),varTitle[0][THETA].Data()));
    hRvVar_avg->SetXTitle(varTitle[0][THETA]);
    hRvVar_avg->SetTitleSize(0.05);
    //! Make Titles nice
    hRvVar_avg->SetTitle("");
    TPaveText* pt = new TPaveText(0.3, 0.85, 0.7, 1.0, "NDC");
    TText* q2wt = pt->AddText(TString::Format("[Q^{2}][W] = %s",Q2Wdirname.Data()));
    q2wt->SetTextColor(kBlue);
    TText* vart = pt->AddText(TString::Format("D^{%s} vs. %s",varTitle[0][THETA].Data(),varTitle[0][THETA].Data()));
    vart->SetTextSize(0.05);

    //! DRAW & SAVE hRvVar_exp, hRvVar_sim
    TLine l(0,0,180,0);
    TCanvas cRvVar("RvVar", "RvVar");
    //cRvVar.SetGridy();
    /*hRvVar[kExp][POS]->SetMinimum(-0.008);
    hRvVar[kExp][POS]->SetMaximum(0.008);
    hRvVar[kExp][POS]->Draw();
    hRvVar[kExp][NEG]->Draw("same");*/
    hRvVar_avg->SetMinimum(-0.003);
    hRvVar_avg->SetMaximum(0.003);
    hRvVar_avg->Draw("ep");
    pt->Draw();
    //hRvVar[kExp][UNPOL]->Draw("same");
    l.Draw("same");
    cRvVar.SaveAs(TString::Format("%s/%s.png", outdir.Data(), cRvVar.GetName()));
    cRvVar.SaveAs(TString::Format("%s/%s_%d.eps", outdir_root.Data(), cRvVar.GetName(),iq2wbin+1));

    iq2wbin+=1;   
  }//end loop over Q2W bins
}//end plotr3_dnp()

