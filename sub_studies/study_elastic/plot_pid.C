{
  //gStyle->SetOptStat(0);
  TFile* _file0= new TFile(gSystem->ExpandPathName("$DELASTDIR_EXP/mon/delastRmon.root"));
  _file0->cd("delast/cut");
  TCanvas c_betavp;
  c_betavp.Divide(3,2);
  for (int i=0;i<6;i++){
    TCut cut_sctr(TString::Format("sector_e==%d",i+1));
    c_betavp.cd(i+1);
    t->Draw(TString::Format("beta_p:p_p>>hbetavp_%d(100,0,3,100,0,1)",i+1),cut_sctr,"colz");
  }
}
