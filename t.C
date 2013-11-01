{
TH1F* hh = new TH1F("th","th",10,0,10);
//hh->Fill(0,100);
hh->SetBinContent(0+1,100);
hh->SetBinError(0+1,100);
TCanvas *c = new TCanvas();
c->cd();
hh->SetMaximum(10);
hh->Draw("p");
}
