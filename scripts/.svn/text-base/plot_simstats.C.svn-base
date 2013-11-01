#include <TFile.h>
#include <THnSparse.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>

#include <vector>
#include <iostream>
using namespace std;

void saveGraph(TGraph* gr[2], TLegend* l);
void seth4AxesRange(THnSparse* h4, int Top, int Q2Wbin, int Varset);

void plot_simstats(TString xsectype="vm", Int_t Varset=1){
  TFile* f = new TFile(TString::Format("simstats_%s.root",xsectype.Data()));  

  THnSparse* h4_n0sr;
  THnSparse* h4_n0acc;
  THnSparse* h4_nFsr;
  THnSparse* h4_nFth;
  
  h4_n0sr  = (THnSparse*)f->Get("h4_n0sr");
  h4_n0acc = (THnSparse*)f->Get("h4_n0acc");
  h4_nFsr  = (THnSparse*)f->Get("h4_nFsr");
  h4_nFth  = (THnSparse*)f->Get("h4_nFth");
   

  //Does not matter which of the h4s used to get the following data
  Int_t nTops    = h4_n0sr->GetAxis(0)->GetNbins();
  Int_t nSims    = h4_n0sr->GetAxis(1)->GetNbins();
  Int_t nQ2Wbins = h4_n0sr->GetAxis(2)->GetNbins();
  Int_t nVarsets = h4_n0sr->GetAxis(3)->GetNbins();
  printf("[nTops][nSims][nQ2Wbins][nVarSets] = [%d][%d][%d][%d]\n",nTops,nSims,nQ2Wbins,nVarsets);
 
  TH1D*** h_n0sr = new TH1D**[nQ2Wbins];
  TH1D*** h_Rn0sr = new TH1D**[nQ2Wbins]; //ratio: n0sr[iTop]/n0sr[Top=2]
  TH1D*** h_Fn0sr = new TH1D**[nQ2Wbins]; //fract: n0sr[iTop]/nFth
  TH1D*** h_n0acc = new TH1D**[nQ2Wbins];
  //TH2D*** h_n0srVn0acc = new TH2D**[nQ2Wbins];
  TH1D*** h_nFsr = new TH1D**[nQ2Wbins];
  TH1D*** h_nFth = new TH1D**[nQ2Wbins];

  for (int iQ2Wbin = 0; iQ2Wbin < nQ2Wbins; ++iQ2Wbin)
  {
    h_n0sr[iQ2Wbin] = new TH1D*[nTops];
    h_Rn0sr[iQ2Wbin] = new TH1D*[nTops];
    h_Fn0sr[iQ2Wbin] = new TH1D*[nTops];
    h_n0acc[iQ2Wbin] = new TH1D*[nTops];
    //h_n0srVn0acc[iQ2Wbin] = new TH2D*[nTops];
    h_nFsr[iQ2Wbin] = new TH1D*[nTops];
    h_nFth[iQ2Wbin] = new TH1D*[nTops]; 
  }

  for (int iQ2Wbin = 0; iQ2Wbin < nQ2Wbins; ++iQ2Wbin)
  {
    for (int iTop = 0; iTop < nTops; ++iTop)
    {
      seth4AxesRange(h4_n0sr,iQ2Wbin+1, iTop+1, Varset);
      seth4AxesRange(h4_n0acc,iQ2Wbin+1, iTop+1, Varset);
      seth4AxesRange(h4_nFsr,iQ2Wbin+1, iTop+1, Varset);
      seth4AxesRange(h4_nFth,iQ2Wbin+1, iTop+1, Varset);

      h_n0sr[iQ2Wbin][iTop] = (TH1D*)h4_n0sr->Projection(1);
      h_n0acc[iQ2Wbin][iTop] = (TH1D*)h4_n0acc->Projection(1);
      h_nFsr[iQ2Wbin][iTop] = (TH1D*)h4_nFsr->Projection(1);
      h_nFth[iQ2Wbin][iTop] = (TH1D*)h4_nFth->Projection(1);

      if (iTop+1==2 && iQ2Wbin+1==24)
      {
        h_n0sr[iQ2Wbin][iTop]->Draw();
      }
    }
  }

  //make R/F histograms for holes
  for (int iQ2Wbin = 0; iQ2Wbin < nQ2Wbins; ++iQ2Wbin)
  {
    gStyle->SetOptStat(0);
    TLegend* lR = new TLegend(0.7,0.3,0.9,0.5);
    TCanvas* cR = new TCanvas(TString::Format("cR_%02d", iQ2Wbin+1));
    TLegend* lF = new TLegend(0.7,0.3,0.9,0.5);
    TCanvas* cF = new TCanvas(TString::Format("cF_%02d", iQ2Wbin+1));
    TLegend* lnFth = new TLegend(0.7,0.3,0.9,0.5);
    TCanvas* cnFth = new TCanvas(TString::Format("cnFth_%02d", iQ2Wbin+1));
    TLegend* lnFsr = new TLegend(0.7,0.3,0.9,0.5);
    TCanvas* cnFsr = new TCanvas(TString::Format("cnFsr_%02d", iQ2Wbin+1));
    for (int iTop = 0; iTop < nTops; ++iTop)
    {
      h_Rn0sr[iQ2Wbin][iTop] = (TH1D*)h_n0sr[iQ2Wbin][iTop]->Clone();
      h_Rn0sr[iQ2Wbin][iTop]->Divide(h_n0sr[iQ2Wbin][1]);
      h_Rn0sr[iQ2Wbin][iTop]->SetMinimum(0);
      h_Rn0sr[iQ2Wbin][iTop]->SetTitle(TString::Format("h[Top]/h[2]:%s",h4_n0sr->GetAxis(2)->GetBinLabel(iQ2Wbin+1)));
      h_Rn0sr[iQ2Wbin][iTop]->SetXTitle("nSim");
      h_Rn0sr[iQ2Wbin][iTop]->SetMarkerStyle(iTop+20);
      h_Rn0sr[iQ2Wbin][iTop]->SetMarkerColor(iTop+2);
      lR->AddEntry(h_Rn0sr[iQ2Wbin][iTop],TString::Format("top%d",iTop+1));
      cR->cd();
      if (iTop==0) h_Rn0sr[iQ2Wbin][iTop]->Draw("p");
      else h_Rn0sr[iQ2Wbin][iTop]->Draw("p sames");
      
      h_Fn0sr[iQ2Wbin][iTop] = (TH1D*)h_n0sr[iQ2Wbin][iTop]->Clone();
      h_Fn0sr[iQ2Wbin][iTop]->Divide(h_nFth[iQ2Wbin][iTop]);
      h_Fn0sr[iQ2Wbin][iTop]->SetMinimum(0);
      h_Fn0sr[iQ2Wbin][iTop]->SetTitle(TString::Format("h[Top]/nThBins:%s",h4_n0sr->GetAxis(2)->GetBinLabel(iQ2Wbin+1)));
      h_Fn0sr[iQ2Wbin][iTop]->SetXTitle("nSim");
      h_Fn0sr[iQ2Wbin][iTop]->SetMarkerStyle(iTop+20);
      h_Fn0sr[iQ2Wbin][iTop]->SetMarkerColor(iTop+2);
      lF->AddEntry(h_Rn0sr[iQ2Wbin][iTop],TString::Format("top%d",iTop+1));
      cF->cd();
      if (iTop==0) h_Fn0sr[iQ2Wbin][iTop]->Draw("p");
      else h_Fn0sr[iQ2Wbin][iTop]->Draw("p sames");

      h_nFth[iQ2Wbin][iTop]->SetTitle(TString::Format("nFth:%s",h4_nFth->GetAxis(2)->GetBinLabel(iQ2Wbin+1)));
      h_nFth[iQ2Wbin][iTop]->SetXTitle("nSim");
      h_nFth[iQ2Wbin][iTop]->SetYTitle("nBins");
      h_nFth[iQ2Wbin][iTop]->SetMarkerStyle(iTop+20);
      h_nFth[iQ2Wbin][iTop]->SetMarkerColor(iTop+2);
      lnFth->AddEntry(h_nFth[iQ2Wbin][iTop],TString::Format("top%d",iTop+1));
      cnFth->cd();
      if (iTop==0) h_nFth[iQ2Wbin][iTop]->Draw("p");
      else h_nFth[iQ2Wbin][iTop]->Draw("p sames");

      h_nFsr[iQ2Wbin][iTop]->SetTitle(TString::Format("nFsr:%s",h4_nFsr->GetAxis(2)->GetBinLabel(iQ2Wbin+1)));
      h_nFsr[iQ2Wbin][iTop]->SetXTitle("nSim");
      h_nFsr[iQ2Wbin][iTop]->SetYTitle("nBins");
      h_nFsr[iQ2Wbin][iTop]->SetMaximum(20000);
      h_nFsr[iQ2Wbin][iTop]->SetMarkerStyle(iTop+20);
      h_nFsr[iQ2Wbin][iTop]->SetMarkerColor(iTop+2);
      lnFsr->AddEntry(h_nFsr[iQ2Wbin][iTop],TString::Format("top%d",iTop+1));
      cnFsr->cd();
      if (iTop==0) h_nFsr[iQ2Wbin][iTop]->Draw("p");
      else h_nFsr[iQ2Wbin][iTop]->Draw("p sames");

    } 

    cR->cd();
    lR->Draw();
    cR->SaveAs(TString::Format("%s/%s.jpg", "simstats",cR->GetName()));
    cF->cd();
    lF->Draw();
    cF->SaveAs(TString::Format("%s/%s.jpg", "simstats",cF->GetName()));
    cnFth->cd();
    lnFth->Draw();
    cnFth->SaveAs(TString::Format("%s/%s.jpg", "simstats",cnFth->GetName()));
    cnFsr->cd();
    lnFsr->Draw();
    cnFsr->SaveAs(TString::Format("%s/%s.jpg", "simstats",cnFsr->GetName()));
  }
  delete h_n0sr;
  delete h_n0acc;
  delete h_nFsr;
  delete h_nFth;
}

void seth4AxesRange(THnSparse* h4, int Q2Wbin, int Top, int Varset)
{
  h4->GetAxis(0)->SetRange(Top,Top);
  h4->GetAxis(2)->SetRange(Q2Wbin,Q2Wbin);
  h4->GetAxis(3)->SetRange(Varset, Varset);
}

void saveGraph(TGraph* gr[2], TLegend* l)
{
  TString cname  = gr[0]->GetName();
  TString ctitle = gr[0]->GetTitle();
  TCanvas* c = new TCanvas(cname, ctitle, 1000, 500);
  c->Divide(2,1);
  c->cd(1);
  gr[0]->Draw("AP");
  l->Draw();
  c->cd(2);
  gr[1]->Draw("AP");
  l->Draw();
  c->SaveAs(TString::Format("%s/%s.jpg","simstats",cname.Data()));
  c->Close();
}

/*void tnailZoom() {
   TPad *sel = (TPad*)gPad->GetSelectedPad();
   int px = gPad->GetEventX();
   int py = gPad->GetEventY();
   if (sel && 
       sel!=c_n0acc[0] && sel!=c_n0acc[1] && sel!=c_n0acc[2] && sel!=c_n0acc[3] && sel!=c_n0acc[4] &&
       sel!=c_n0sr[0]  && sel!=c_n0sr[1]  && sel!=c_n0sr[2]  && sel!=c_n0sr[3]  && sel!=c_n0sr[4] &&
       sel!=c_n0srVn0acc[0][0] && sel!=c_n0srVn0acc[0][1] && sel!=c_n0srVn0acc[0][2] && sel!=c_n0srVn0acc[0][3] && sel!=c_n0srVn0acc[0][4] && 
       sel!=c_n0srVn0acc[1][0] && sel!=c_n0srVn0acc[1][1] && sel!=c_n0srVn0acc[1][2] && sel!=c_n0srVn0acc[1][3] && sel!=c_n0srVn0acc[1][4] && 
       sel!=c_nFsr[0] && sel!=c_nFsr[1] && sel!=c_nFsr[2] && sel!=c_nFsr[3] && sel!=c_nFsr[4] && 
       sel!=c_nFth[0] && sel!=c_nFth[1] && sel!=c_nFth[2] && sel!=c_nFth[3] && sel!=c_nFth[4]) {
      if (selold) delete selold;
      c_zoom->cd();
      TPad *newpad = (TPad*)sel->Clone();
      c_zoom->GetListOfPrimitives()->Add(newpad);
      newpad->SetPad(0,0,1,1);
      selold = newpad;
      c_zoom->Update();
   }
}*/
