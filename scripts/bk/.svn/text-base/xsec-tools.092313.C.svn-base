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

#include "myTHnTool.h"

#include <string>
#include <iostream>
#include <algorithm>
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

    if      (pwd.Contains("experiment"))      _dtype = "exp";
    else if (pwd.Contains("simulation"))      _dtype = "sim";
    path=".";
  }else{
    _q2w_bng=q2w_bng;
    printf("q2wbng = %s\n", _q2w_bng.Data());
    if(_q2w_bng.Contains(Q2W_BNG1))
    {
      path="/data/trivedia/e1f/ana2pi/experiment/e1f.goldenruns.qskim_071513/Q2W__1.4-1.5__1.6-1.8";
    }else if (_q2w_bng.Contains(Q2W_BNG2)){
      path="/data/trivedia/e1f/ana2pi/experiment/e1f.goldenruns.qskim_071513/Q2W__1.9-2.5__1.3-1.9";
    }else 
    {
      printf("Could not determine path\n");
      return kFALSE; 
    }
    _dtype="exp";
  }

  if (pol) _dtype = "pol__exp";
  

  cout << "xsectype = " << xsectype << endl;
  if (xsectype=="at"){
    copy(AT_TOPNAMES,AT_TOPNAMES+5, _topNames);
    _fyexp[0] = TFile::Open(TString::Format("%s/1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[1] = TFile::Open(TString::Format("%s/2__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[2] = TFile::Open(TString::Format("%s/3__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[3] = TFile::Open(TString::Format("%s/4__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[4] = TFile::Open(TString::Format("%s/1:2:3:4__%s__%s.root", path.Data(), _q2w_bng.Data(), _dtype.Data()));
    _fysim[0] = TFile::Open(TString::Format("%s/simdir/1__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[1] = TFile::Open(TString::Format("%s/simdir/2__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[2] = TFile::Open(TString::Format("%s/simdir/3__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[3] = TFile::Open(TString::Format("%s/simdir/4__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[4] = TFile::Open(TString::Format("%s/simdir/1:2:3:4__%s__%s.root", _q2w_bng.Data(), "sim"));
  }else if (xsectype=="vm"){
    copy(VM_TOPNAMES,VM_TOPNAMES+5, _topNames);
    _fyexp[0] = TFile::Open(TString::Format("%s/1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[1] = TFile::Open(TString::Format("%s/2:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[2] = TFile::Open(TString::Format("%s/3:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[3] = TFile::Open(TString::Format("%s/4:1__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fyexp[4] = TFile::Open(TString::Format("%s/1:2:3:4__%s__%s.root", path.Data(),_q2w_bng.Data(), _dtype.Data()));
    _fysim[0] = TFile::Open(TString::Format("%s/simdir/1__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[1] = TFile::Open(TString::Format("%s/simdir/2:1__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[2] = TFile::Open(TString::Format("%s/simdir/3:1__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[3] = TFile::Open(TString::Format("%s/simdir/4:1__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
    _fysim[4] = TFile::Open(TString::Format("%s/simdir/1:2:3:4__%s__%s.root", path.Data(),_q2w_bng.Data(), "sim"));
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
     ploty?hYW_5D[iTop]->Draw("p"):hYWnorm_5D[iTop]->Draw("p");
     if (sim)
     {
      hYW_TH[iTop]->Draw("sames");
      l_5D->AddEntry( hYW_TH[iTop], TString::Format("Thrown") );
     } 

     c_Dir->cd();
     l_Dir->SetHeader( TString::Format("%s - Direct Method", XSECTITLE.Data()) );
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     hYWnorm_Dir[iTop]->SetMinimum(0);
     ploty?hYW_Dir[iTop]->Draw("p"):hYWnorm_Dir[iTop]->Draw("p");
     if (sim)
     {
      hYW_TH[iTop]->Draw("bar sames");
      l_Dir->AddEntry( hYW_TH[iTop], TString::Format("Thrown") );
     }    
   }else{
     c_5D->cd();
     l_5D->AddEntry( hYWnorm_5D[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     ploty?hYW_5D[iTop]->Draw("p sames"):hYWnorm_5D[iTop]->Draw("p sames");
     if (sim)
     {
      hYW_TH[iTop]->Draw("bar sames");
     }

     c_Dir->cd();
     l_Dir->AddEntry( hYWnorm_Dir[iTop], TString::Format("%s", _topNames[iTop].Data()) );
     ploty?hYW_Dir[iTop]->Draw("p sames"):hYWnorm_Dir[iTop]->Draw("p sames");
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

  gStyle->SetOptStat("nemMrRi");

  //sim or data
  TFile* fy[5];
  if (sim) memcpy(fy,_fysim,sizeof(_fysim));
  else     memcpy(fy,_fyexp,sizeof(_fyexp));
 
  //data objects for 5D yields for each Top
  int nTops=4;
  THnSparse* hY5D[nTops];

  //data object yield vs W for each top (only in common bins!)
  TH1F* hYW[nTops];
  TLegend* l = new TLegend(0.1,0.7,0.5,0.9);
  l->SetHeader("Yield - 5D Method(Common AccCorr Bins)");
  TString wbng   = _q2w_bng.Tokenize("__")->At(1)->GetName();
  TString nwbins = wbng.Tokenize("-")->At(0)->GetName();
  TString wmin   = wbng.Tokenize("-")->At(1)->GetName();
  TString wmax   = wbng.Tokenize("-")->At(2)->GetName();
  for(int iTop=0;iTop<nTops;iTop++){
    hYW[iTop] = new TH1F(TString::Format("hYW_%d",iTop+1),TString::Format("hYW_%d",iTop+1), nwbins.Atoi(), wmin.Atof(), wmax.Atof());
    hYW[iTop]->SetXTitle("W[GeV]");
    hYW[iTop]->SetMarkerStyle(20+iTop); //+5
    hYW[iTop]->SetMarkerColor(2+iTop);
    l->AddEntry( hYW[iTop], TString::Format("%s", _topNames[iTop].Data()) );
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
    printf("Processing %s\n", Q2Wdirname.Data());
    counterQ2Wbin+=1;
    
    //!if Q2Wbin option > 0, then process only the specified Q2Wbin
    if (Q2Wbin>0 && (counterQ2Wbin != Q2Wbin)) continue; 

    TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
    TString wlow = wrange.Tokenize(",")->At(0)->GetName();
    wlow.Remove(0,1); //remove "["
    double w = wlow.Atof();
      
    char hname[200];                 
    sprintf(hname, "%s/hY5D/Varset1/hY5D_%s", Q2Wdirname.Data(), seqTitle[seq].Data());
    for(int iTop=0;iTop<nTops;iTop++){
      hY5D[iTop] = (THnSparse*)fy[iTop]->Get(hname);
    }
    float intg_hY5D[nTops];// = {0,0,0,0};
    float intgErr_hY5D[nTops];

    myTHnTool hntool(kFALSE);
    int nComBins = hntool.GetIntegralCommonBins(hY5D[0],hY5D[1],hY5D[2],hY5D[3],
                                      intg_hY5D[0],intg_hY5D[1],intg_hY5D[2],intg_hY5D[3],
                                      intgErr_hY5D[0],intgErr_hY5D[1],intgErr_hY5D[2],intgErr_hY5D[3],
                                      draw);
    printf("Number of common bins for q2w bin %s = %d\n", Q2Wdirname.Data(), nComBins);
    //printf("intgs =  %.2f,%.2f,%.2f,%.2f\n",intg_hY5D[0],intg_hY5D[1],intg_hY5D[2],intg_hY5D[3]);
    //printf("intgs Errors =  %.2f,%.2f,%.2f,%.2f\n",intgErr_hY5D[0],intgErr_hY5D[1],intgErr_hY5D[2],intgErr_hY5D[3]);

    for(int iTop=0;iTop<nTops;iTop++){
      hYW[iTop]->SetBinContent(hYW[iTop]->FindBin(w+.005), intg_hY5D[iTop]); //+5Mev = _instrinsic.Wbinw
      hYW[iTop]->  SetBinError(hYW[iTop]->FindBin(w+.005), intgErr_hY5D[iTop]); //+5Mev = _instrinsic.Wbinw
    }
  }
  TCanvas* cy_combins = new TCanvas();
  if (seq==HOLE) cy_combins->SetLogy();
  for(int iTop=0;iTop<nTops;iTop++){
    if (iTop==0) {
      if (seq==HOLE) hYW[iTop]->SetMinimum(1);
      hYW[iTop]->Draw("p");
    }else hYW[iTop]->Draw("p sames");
  }
  l->Draw("same");
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

void plotr(int top){
  if (setup("vm","","pol__")==kFALSE) return;

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
    TH1D* hpos = (TH1D*)fy->Get(hname_pos);
    TH1D* hneg = (TH1D*)fy->Get(hname_neg);
    if (hpos==NULL || hneg==NULL) cout << "histogram not found" << endl;

    //cosmetics and display
    hpos->SetLineColor(kRed);
    hneg->SetLineColor(kBlue);
    c->cd(iq2wbin+1);
    hpos->Draw();
    hneg->Draw("sames");
    iq2wbin+=1;
  }
  c->Update();
}

void plotphi(int top, bool phi=kFALSE){
  if (setup("vm","","pol__")==kFALSE) return;
  gStyle->SetOptFit(1111); //pcev = 1111
  //gStyle->SetTitleW(1.5);

  int itop = top-1;
  TFile* fy = _fyexp[itop];

  int nq2wbins = fy->GetNkeys();
  TIter nextkey(fy->GetListOfKeys());
  TKey *key;
  int iq2wbin = 0;
  while (key = (TKey*)nextkey()) {
    TString Q2Wdirname = key->GetName();
    if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
    cout << "Q2Wdirname = " << Q2Wdirname << endl;

    THStack* hs[10];
    TF1 *fphi = new TF1("fphi", "([0] + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()) + [3]*sin(x*TMath::DegToRad()))",0, 360);
    fphi->SetParameter(0,1);
    fphi->SetParameter(1,10);
    fphi->SetParameter(2,20);
    fphi->SetParameter(3,100);
    fphi->SetParName(0, "A");
    fphi->SetParName(1, "B");
    fphi->SetParName(2, "C");
    fphi->SetParName(3, "hPR");
    for (int i = 0; i < 10; i++)
    {
      TString hname_unp = TString::Format("%s/hPhi/Varset1/theta/hphi_%02d",Q2Wdirname.Data(),i+1);
      TString hname_pos = TString::Format("%s/hPhi_POS/Varset1/theta/hphi_%02d",Q2Wdirname.Data(),i+1);
      TString hname_neg = TString::Format("%s/hPhi_NEG/Varset1/theta/hphi_%02d",Q2Wdirname.Data(),i+1);
      cout << "hname_unp = " << hname_unp << endl;
      cout << "hname_pos = " << hname_pos << endl;
      cout << "hname_neg = " << hname_neg << endl;
      TH1D* hunp = (TH1D*)fy->Get(hname_unp);
      TH1D* hpos = (TH1D*)fy->Get(hname_pos);
      TH1D* hneg = (TH1D*)fy->Get(hname_neg);
      if (hunp==NULL || hpos==NULL || hneg==NULL) cout << "histogram not found" << endl;
      hunp->SetLineColor(kBlack);
      hpos->SetLineColor(kRed);
      hneg->SetLineColor(kBlue);
      hunp->Fit(fphi);
      hpos->Fit(fphi);
      hneg->Fit(fphi);

      //! Modify titles
      TObjArray *tarr;
      char t[100];

      TPaveText *ptunp = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
      TString tunp = hunp->GetTitle();
      tarr = tunp.Tokenize("|");
      sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
      ptunp->AddText(tarr->At(0)->GetName());
      ptunp->AddText(t);
      hunp->SetTitle("");

      TPaveText *ptpos = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
      TString tpos = hpos->GetTitle();
      tarr = tpos.Tokenize("|");
      sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
      ptpos->AddText(tarr->At(0)->GetName());
      ptpos->AddText(t);
      hpos->SetTitle("");

      TPaveText *ptneg = new TPaveText(0.0, 0.9, 1.0, 1, "NDC");
      TString tneg = hneg->GetTitle();
      tarr = tneg.Tokenize("|");
      sprintf(t, "%s:%s:%s", tarr->At(1)->GetName(), tarr->At(2)->GetName(), tarr->At(3)->GetName());
      ptneg->AddText(tarr->At(0)->GetName());
      ptneg->AddText(t);
      hneg->SetTitle("");

      //! Draw histograms
      TCanvas *c = new TCanvas(hpos->GetName(),hpos->GetName(),900, 600);
      c->Divide(3,1);
      c->cd(1);
      hunp->Draw();
      ptunp->Draw();
      c->cd(2);
      hpos->Draw();
      ptpos->Draw();
      c->cd(3);
      hneg->Draw();
      ptneg->Draw();
      TString dir = TString::Format("./polobs/top%d/%s/Varset1/theta", top, Q2Wdirname.Data());
      gSystem->mkdir(dir,1);
      c->SaveAs(TString::Format("%s/%s.png", dir.Data(), c->GetName()));
      c->Close();
      c->Delete();
    }
  }
}



