{
TFile *fn=TFile::Open("delastR_npcorr_e16.root");
TFile *fw=TFile::Open("delastR_wpcorr_e16.root");
TH1F* hn=fn->Get("delast/helast");
TH1F* hw=fw->Get("delast/helast");
hn->SetName("npcorr");
hw->SetName("wpcorr");

//! Draw and save plots
TString tag="cW-elast_e16_check";
TFile* fout=new TFile(TString::Format("%s.root",tag.Data()),"RECREATE");

gStyle->SetOptFit(1111);
TCanvas *c=new TCanvas("c","c",500,650);
c->Divide(1,2);

c->cd(1);
hn->Draw();
hn->Fit("gaus","","",0.9,1);

c->cd(2);
hw->Draw();
hw->Fit("gaus","","",0.9,0.97);

c->Update();
c->Write();
c->SaveAs(TString::Format("%s.png",tag.Data()));
}
