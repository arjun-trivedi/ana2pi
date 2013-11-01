#include "TFile.h"
#include "TKey.h"
#include "TH2.h"
#include "THStack.h"
#include "TCanvas.h"

#include "CanvasManager.h"

TH2* normXslices(TH2 *h) {
  TH2D *h2 = (TH2D*)h->Clone();
  Int_t nby = h2->GetYaxis()->GetNbins();
  for (Int_t bx = 1; bx <= h2->GetXaxis()->GetNbins(); bx++) {
    h2->GetXaxis()->SetRange(bx,bx);
    TH1 *h1 = h2->ProjectionY();
    Double_t max = h1->GetMaximum();
    for (Int_t by = 1; by <= nby; by++) {
      Double_t bc = h2->GetBinContent(bx,by);
      if (bc > 0) {
Int_t global_bin = h2->GetBin(bx,by);
h2->SetBinContent(global_bin,bc/max);
      }
    }
  }
  h2->GetXaxis()->SetRange(0,h2->GetXaxis()->GetNbins()+1);
  h2->SetMaximum(1);
  return h2;
}

TH2* normYslices(TH2 *h) {
  TH2D *h2 = (TH2D*)h->Clone();
  Int_t nbx = h2->GetXaxis()->GetNbins();
  for (Int_t by = 1; by <= h2->GetYaxis()->GetNbins(); by++) {
    h2->GetYaxis()->SetRange(by,by);
    TH1 *h1 = h2->ProjectionX();
    Double_t max = h1->GetMaximum();
    for (Int_t bx = 1; bx <= nbx; bx++) {
      Double_t bc = h2->GetBinContent(bx,by);
      if (bc > 0) {
Int_t global_bin = h2->GetBin(bx,by);
h2->SetBinContent(global_bin,bc/max);
      }
    }
  }
  h2->GetYaxis()->SetRange(0,h2->GetYaxis()->GetNbins()+1);
  h2->SetMaximum(1);
  return h2;
}

//coption=TCanvas Draw option; hoption= TH2 draw optopn
void plot(TString dirname, char* hoption="", char* coption=""){
  TFile* fyield;
  if (gDirectory->GetFile())  fyield = gDirectory->GetFile();
  else fyield = TFile::Open("yield.root", "r");

  TDirectory* dir = fyield->GetDirectory(dirname);
  if (dir==NULL){
    printf("%s not found in %s\n", dirname.Data(), fyield->GetName());
    return;
  }
 
  TString dirpath = dir->GetPath();
  printf("dirpath = %s\n", dirpath.Data());
  TObjArray* dirpathArray = dirpath.Tokenize("/");
  TString proc, procmode, sector;
  if (dirpathArray->At(1) != NULL )proc = dirpathArray->At(1)->GetName();
  if (dirpathArray->At(2) != NULL )procmode = dirpathArray->At(2)->GetName();
  if (dirpathArray->At(3) != NULL) sector = dirpathArray->At(3)->GetName();
  printf("proc/procmode/sector = %s/%s/%s\n", proc.Data(), procmode.Data(), sector.Data()); 

  Int_t nHistsEid = 5;
  TString hists_eid[5] = {"heoutVein", "hsfoutVsfin", "hsftotVp", "hnphe", "hCCthetaVCCseg"};
  
  Int_t nHistsEFid = 3;
  TString hists_efid[3] = {"hephiVetheta", "hescyVescx", "heecyVeecx"};
  
  Int_t nHistsPid = 21;
  TString hists_pid[21] = {"hP_betaVp", "hP_dtVp", "hP_eoutVein", "hP_sfoutVsfin", "hP_sftotVp", "hP_nphe", "hP_CCthetaVCCseg",
                          "hPip_betaVp", "hPip_dtVp", "hPip_eoutVein", "hPip_sfoutVsfin", "hPip_sftotVp", "hPip_nphe", "hPip_CCthetaVCCseg",
                          "hPim_betaVp", "hPim_dtVp", "hPim_eoutVein", "hPim_sfoutVsfin", "hPim_sftotVp", "hPim_nphe", "hPim_CCthetaVCCseg"};

  Int_t nHists;
  TString* hists = NULL;

  if(strcmp(proc,"eid")==0 || strcmp(proc,"eid2")==0){
    nHists = nHistsEid;
    hists = new TString[nHists];
    hists = hists_eid;
  }else if (strcmp(proc,"fid")==0 || strcmp(proc,"fid2")==0){
    nHists = nHistsEFid;
    hists = new TString[nHists];
    hists = hists_efid;
  }else if (strcmp(proc,"pid")==0 || strcmp(proc,"pid2")==0){
    nHists = nHistsPid;
    hists = new TString[nHists];
    hists = hists_pid;
  }else{
    printf("anaMonTools not implemented for proc %s\n", proc.Data());
    return;
  }

  Int_t _cww = 350;
  Int_t _cwh = 250;

  if( strcmp(sector,"")!=0) {
    printf("making plots for %s\n", sector.Data());
    CanvasManager* cm = new CanvasManager(_cww, _cwh);
    THStack *hs = new THStack(TString::Format("hs%s", sector.Data()), TString::Format("hs%s", sector.Data()));
    for (Int_t iHist=0; iHist<nHists; iHist++){
      TObject* h = dir->FindKey(hists[iHist].Data())->ReadObj();
      if ( (strcmp(proc,"fid")==0 || strcmp(proc,"fid2")==0) && strcmp(h->GetName(), "hephiVetheta")==0 && strcmp(hoption,"norm")==0 ){
        TObject* hnorm = (TObject*)normXslices((TH2*)h);
        hs->Add((TH2*)hnorm, "colz");
      }else{
      	if (strcmp(h->ClassName(),"TH2F")==0) hs->Add((TH2*)h, "colz");
      	else hs->Add((TH1*)h);
      }
    }
    cm->add(hs);
    cm->draw(coption);
  }else{
    printf("making plots for all sectors\n");
    CanvasManager* cm = new CanvasManager(_cww, _cwh);
    THStack *hs[7];
    for (Int_t iSector=0;iSector<7;iSector++){
      if (iSector==0) continue; //not for sector0
      dir->cd(TString::Format("sector%d", iSector));
      hs[iSector] = new THStack(TString::Format("hs%d", iSector), TString::Format("hs%d", iSector));  
      for (Int_t iHist=0; iHist<nHists; iHist++){
        TObject* h = gDirectory->FindKey(hists[iHist].Data())->ReadObj();
        if ( (strcmp(proc,"fid")==0 || strcmp(proc,"fid2")==0) && strcmp(h->GetName(), "hephiVetheta")==0 && strcmp(hoption,"norm")==0 ){
          TObject* hnorm = (TObject*)normXslices((TH2*)h);
          hs[iSector]->Add((TH2*)hnorm, "colz");
        }else{
          if (strcmp(h->ClassName(),"TH2F")==0) hs[iSector]->Add((TH2*)h, "colz");
          else hs[iSector]->Add((TH1*)h);
        }
      }
      cm->add(hs[iSector]);
    }
    cm->draw(coption);
  }
   
}
