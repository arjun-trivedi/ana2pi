TH2F* norm(TH2F* h2){
  TH2F* h2c=h2->Clone(TString::Format("%s_norm",h2->GetName()));
  int nxbins=h2c->GetNbinsX();
  int nybins=h2c->GetNbinsY();

  for (int ixbin=0;ixbin<nxbins;ixbin++){
    TH1F* hy=(TH1F*)h2c->ProjectionY("_py",ixbin+1,ixbin+1);
    float max=hy->GetMaximum();
    if (max==0) continue;
    //! Now normalize h2 for this projection
    for (int iybin=0;iybin<nybins;iybin++){
      float binc=h2c->GetBinContent(ixbin+1,iybin+1);
      binc=binc/max;
      h2c->SetBinContent(ixbin+1,iybin+1,binc);
    }
  }
  return h2c;
}
