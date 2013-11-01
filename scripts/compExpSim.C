void compExpSim(char* fexpname, char* fsimname, char* type = "recocorr", bool domc=false){
  gStyle->SetOptStat("MRei");
 
  TFile *fexp =   TFile::Open(fexpname);
  TFile *fsim = TFile::Open(fsimname);

  TCanvas* cc[100];
  TCanvas* ccmc[100];

  //determine 2piEvt
  TString evtSel2PiName = TString("evt2Pi") + TString(fexpname[0]);
  //create direction in pwd for pngs
  Bool_t isRecursive=kTRUE;
  TString outdirname = TString::Format("pngs/%s/%s/", evtSel2PiName.Data(), type);
  gSystem->mkdir(outdirname, isRecursive);

  fexp->GetListOfKeys()->Print(); //all Q2W dirs
  TIter nextkey(fexp->GetListOfKeys());
  TKey *key;

  
  Int_t i = 0;
  TH1* hexpnorm;
  TH1* hsimnorm;
  while (key = (TKey*)nextkey()) {
     TString Q2Wdirname = key->GetName();
     printf("Q2Wdir = %s\n", Q2Wdirname.Data());

     TString cnamemc       =  Q2Wdirname + TString("/") + TString("LowStats(I)/h1Da_MC");
     TString cnamereco     =  Q2Wdirname + TString("/") + TString("LowStats(I)/h1Da_RECO");
     TString cnamerecocorr =  Q2Wdirname + TString("/") + TString("LowStats(I)/h1Da_RECOCORR");
     //printf("cnamereco = %s\n", cnamereco.Data());

	 TCanvas* cexp;
	 if (type=="reco") cexp = (TCanvas*)fexp->Get(cnamereco.Data());
     else cexp = (TCanvas*)fexp->Get(cnamerecocorr.Data());
     cexp->SetName("cexp");
     cexp->SetTitle("cexp");
     cexp->Draw();
     cexp->Modified();
     cexp->Update();
     cout << " cexp size = " << ((TList*)cexp->GetListOfPrimitives())->GetSize() << endl; 

     TCanvas* csim;
     if (type=="reco") csim = (TCanvas*)fsim->Get(cnamereco.Data());
     else csim = (TCanvas*)fsim->Get(cnamerecocorr.Data());
     csim->SetName("csim");
     csim->SetTitle("csim");
     csim->Draw();
     csim->Modified();
     csim->Update();
     cout << " csim size = " << ((TList*)csim->GetListOfPrimitives())->GetSize() << endl;

     //if(domc){
       TCanvas* cmc = (TCanvas*)fsim->Get(cnamemc.Data());
       cmc->SetName("cmc");
       cmc->SetTitle("cmc");
       cmc->Draw();
       cmc->Modified();
       cmc->Update();
       cout << " cmc size = " << ((TList*)cmc->GetListOfPrimitives())->GetSize() << endl;
     
     if(domc){
       ccmc[i] = new TCanvas((TString("mc")+  TString("_") + Q2Wdirname).Data(), (TString("mc")+  TString("_") + Q2Wdirname).Data(), 1400, 900);
       ccmc[i]->Divide(4,3);
     }


     cc[i]   = new TCanvas((evtSel2PiName +  TString("_") + Q2Wdirname).Data(), (evtSel2PiName + TString("_") + Q2Wdirname).Data(), 1400, 900);
     cc[i]->Divide(4,3);
     
     for (int iPad=1;iPad<13;iPad++){
       TH1D* hexp = (TH1D*)cexp->GetPad(iPad)->GetListOfPrimitives()->At(1);
       hexp->SetTitleSize(0.05,"x");
       hexp->SetLineWidth(2);
       hexp->SetLineColor(kBlue);
       TH1D* hsim = (TH1D*)csim->GetPad(iPad)->GetListOfPrimitives()->At(1);
       hsim->SetTitleSize(0.05, "x");
       hsim->SetLineWidth(2);
       hsim->SetLineColor(kRed);
       //if(domc){
         TH1D* hmc = (TH1D*)cmc->GetPad(iPad)->GetListOfPrimitives()->At(1);
         hmc->SetTitleSize(0.05, "x");
         hmc->SetLineWidth(2);
         hmc->SetLineColor(kGreen+2);
       //}

       cc[i]->cd(iPad);
       gPad->SetGridx();
       hexpnorm = hexp->DrawNormalized("", 1000);
       TPaveStats *psexp = (TPaveStats*)gPad->GetPrimitive("stats");
       psexp->SetName("expstats");
       psexp->SetFillStyle(0);
       psexp->SetLineColor(kBlue-7);
       psexp->SetTextColor(kBlue-7);
       psexp->SetX1NDC(0.7);
       psexp->SetX2NDC(1);
       psexp->SetY1NDC(0.7);
       psexp->SetY2NDC(0.9);
       psexp->Draw();

       hsimnorm = hsim->DrawNormalized("sames", 1000);
       gPad->Update();
       TPaveStats *pssim= (TPaveStats*)gPad->GetPrimitive("stats");
       pssim->SetName("simstats");
       pssim->SetFillStyle(0);
       pssim->SetLineColor(kRed-7);
       pssim->SetTextColor(kRed-7);
       pssim->SetX1NDC(0.7);
       pssim->SetX2NDC(1);
       pssim->SetY1NDC(0.5);
       pssim->SetY2NDC(0.7);
       pssim->Draw();
       
       if (type=="mc"){
         hmc->DrawNormalized("sames", 1000);
         gPad->Update();
         TPaveStats *psmc= (TPaveStats*)gPad->GetPrimitive("stats");
         psmc->SetName("mcstats");
         psmc->SetFillStyle(0);
         psmc->SetLineColor(kGreen+2);
         psmc->SetTextColor(kGreen+2);
         psmc->SetX1NDC(0.7);
         psmc->SetX2NDC(1);
         psmc->SetY1NDC(0.3);
         psmc->SetY2NDC(0.5);
         psmc->Draw();
  
         //gPad->Modified();
         //gPad->Update();
       }
       
       Double_t max = 0;
       hexpnorm->GetMaximum()>hsimnorm->GetMaximum()?max=hexpnorm->GetMaximum():max=hsimnorm->GetMaximum();
       hexpnorm->SetMaximum(max+20);
       gPad->Update();

       //gPad->Modified();
       //gPad->Update();

       //TPaveText *pt = (TPaveText*)(gPad->GetPrimitive("title")); 
       //pt->SetTextSize(0.06); 
       //gPad->Modified();

       if (domc){
         ccmc[i]->cd(iPad);
         gPad->SetGridx();
         hmc->DrawNormalized("", 1000);
         gPad->Update();
         TPaveStats *psmc= (TPaveStats*)gPad->GetPrimitive("stats");
         psmc->SetName("mcstats");
         psmc->SetFillStyle(0);
         psmc->SetLineColor(kGreen+2);
         psmc->SetTextColor(kGreen+2);
         psmc->SetX1NDC(0.7);
         psmc->SetX2NDC(1);
         psmc->SetY1NDC(0.3);
         psmc->SetY2NDC(0.5);
         psmc->Draw();
  
         //gPad->Modified();
         //gPad->Update();
       }
     }
     //TString csavename = TString("pngs/") + evtSel2PiName + TString("/") + TString(cc[i]->GetName()) + TString(".png");
     TString csavename = outdirname + TString(cc[i]->GetName()) + TString(".png");
     cout << "csavename = " << csavename.Data() << endl;
     cc[i]->SaveAs(csavename.Data());

     if(domc){
       //TString csimsavename = TString("pngs/") +  TString("mc") + TString("/") + TString(ccmc[i]->GetName()) + TString(".png");
       TString csimsavename = outdirname + TString(cc[i]->GetName()) + TString(".png");
       cout << "csimsavename = " << csimsavename.Data() << endl;
       ccmc[i]->SaveAs(csimsavename.Data());
       
       cmc->Close();
       ccmc[i]->Close();
     }
 
     cexp->Close();
     csim->Close();
     cc[i]->Close();
     i+=1;
     //if (i == 1) break;
  }
  
}
