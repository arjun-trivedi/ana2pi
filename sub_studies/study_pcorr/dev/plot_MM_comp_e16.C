{
TFile *fn=TFile::Open("d2piR_cutset1_npcorr_e16.root");
TFile *fw=TFile::Open("d2piR_cutset1_wpcorr_e16.root");
TH1F* hn=(TH1F*)fn->Get("d2pi/MM/pre_cut/hmmppip_prec_fW");
TH1F* hw=(TH1F*)fw->Get("d2pi/MM/pre_cut/hmmppip_prec_fW");
hn->SetLineColor(kBlack);
hw->SetLineColor(kBlue);
hn->SetName("hmmppip_npcorr");
hw->SetName("hmmppip_wpcorr");

//! Draw and save plots
TString tag="cMM-top2_e16_check";
TFile* fout=new TFile(TString::Format("%s.root",tag.Data()),"RECREATE");

TCanvas *c=new TCanvas("c","c",500,500);
TLegend* l=new TLegend(0.6,0.7,0.9,0.9);
hn->Draw();
hw->Draw("same");
l->AddEntry(hn,"npcorr");
l->AddEntry(hw,"wpcorr");
l->Draw("same");

c->Update();
c->Write();
c->SaveAs(TString::Format("%s.png",tag.Data()));
}
