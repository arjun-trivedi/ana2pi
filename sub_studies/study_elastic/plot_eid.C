{
  //gStyle->SetOptStat(0);
  TFile* _file0= new TFile(gSystem->ExpandPathName("$DELASTDIR_EXP/mon/delastRmon.root"));
  _file0->cd("delast/cut");

  //! Over all EID
  TCanvas* c_nphe=new TCanvas("nphe","nphe");
  c_nphe->Divide(3,2);
  TCanvas* c_sf=new TCanvas("sf","sf");
  c_sf->Divide(3,2);
  for (int i=0;i<6;i++){
    TCut cut_sctr(TString::Format("sector_e==%d",i+1));
    c_nphe->cd(i+1);
    t->Draw(TString::Format("nphe>>hnphe_%d(100,0,300)",i+1),cut_sctr,"");
    c_sf->cd(i+1);
    t->Draw(TString::Format("etot/p_e_eid:p_e_eid>>hsfvsp_%d(100,0,5,100,0,0.5)",i+1),cut_sctr,"colz");
  }

  //! See dependence on theta_e
  TCut *cut_theta[2];
  TCanvas* c_nphe_cut_theta[2];
  TCanvas* c_sf_cut_theta[2];
  for (int i=0;i<2;i++){
    if (i==0) cut_theta[i]=new TCut("theta_e<25");
    else      cut_theta[i]=new TCut("theta_e>25");

    c_nphe_cut_theta[i]=new TCanvas(TString::Format("nphe%d",i+1),TString::Format("nphe%d",i+1));	
    c_nphe_cut_theta[i]->Divide(3,2);
    c_sf_cut_theta[i]=new TCanvas(TString::Format("sf%d",i+1),TString::Format("sf%d",i+1));          
    c_sf_cut_theta[i]->Divide(3,2);
  }

  for (int i=0;i<2;i++){
    for (int j=0;j<6;j++){
      TCut cut_sctr(TString::Format("sector_e==%d",j+1));
      c_nphe_cut_theta[i]->cd(j+1);
      t->Draw(TString::Format("nphe>>hnphe_cut_theta%d_sector%d(100,0,300)",i+1,j+1),cut_sctr&&(*cut_theta[i]),"");
      c_sf_cut_theta[i]->cd(j+1);
      t->Draw(TString::Format("etot/p_e_eid:p_e_eid>>hsfvsp_cut_theta%d_sector%d(100,0,5,100,0,0.5)",i+1,j+1),cut_sctr&&(*cut_theta[i]),"colz");
    }	
  }
}
