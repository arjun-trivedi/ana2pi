void ts(){
  gStyle->SetOptFit(1111);
  TH1F* hh = new TH1F("t", "t", 10, 0, 10);
  //hh->SetBinContent(1,0.0000001);
  //hh->SetBinContent(2,0.0001);
  hh->FillRandom("gaus",1000);
  hh->Fit("gaus");
  TCanvas* c = new TCanvas();
  hh->Draw();
  //c->Update();
  gPad->Update();
  TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
  //s->Draw();
  s->SetTextSize(0.05);
  s->SetFitFormat(".2e");
  s->SetX1NDC(0.5);
  s->SetX2NDC(0.9);
  s->SetY1NDC(0.5);
  s->SetY2NDC(0.7);
  s->Draw();
  gPad->Update(); 
}
