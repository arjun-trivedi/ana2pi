{
  TH1F* h[6];
  for (int i=0;i<6;i++){
	h[i]=(TH1F*)_file0->Get(TString::Format("/ER/sector%d/phibinnum3/hTHETA",i+1));
  }

  TH1F* hn[6];
  //! Take ratios of ER w.r.t. sector1
  for (int i=0;i<6;i++){
	hn[i]=(TH1F*)h[i]->Clone();
	hn[i]->Divide(h[0]);
  }

  TCanvas *c =new TCanvas();
  c.Divide(2,3);
  for (int i=0;i<6;i++){
    c.cd(i+1);
    hn[i]->Draw();
  }

}
