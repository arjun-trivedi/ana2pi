void plot_d2pi_Q2W(int top=1){
  const int NDTYP=2;//3;
  TFile* f[NDTYP];//0=ER,1=SR,2=ST
  enum{ER,SR};//,ST};
  TString dtyp_name[]={"ER","SR"};//,"ST"};
  int clr[]={kBlue,kRed};//,kGreen};
  f[0]=TFile::Open("$STUDY_BG_DATADIR/d2pi_exp/d2piR_hvy.root");
  f[1]=TFile::Open("$STUDY_BG_DATADIR/d2pi_sim/d2piR_hvy.root");
  TTree* t[NDTYP];
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    t[idtyp]=(TTree*)f[idtyp]->Get("d2piR/tR");
  }

  TCut cut_top(TString::Format("top==%d",top));
  TCut cut_mm2;
  if (top==1){
    cut_mm2="mm2ppippim>-0.0005 && mm2ppippim<0.0005";
  }else if (top==2){
    cut_mm2="mm2ppip>0.00 && mm2ppip<0.04";
  }else{
      printf("Not implemented for top=%d",top);
      return;
  }
  TCut cut=TCut(cut_top);//&&cut_mm2;
  cut+=cut_mm2;


  TH2F* h_q2w[NDTYP];
  TCanvas*c = new TCanvas("c","c");//default canvas
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    cout<<f[idtyp]->GetName()<<" " <<t[idtyp]->GetName()<<endl;
    c->cd();
    t[idtyp]->Draw("Q2:W>>h(100,0,3.5,150,0,6)",cut);

    gStyle->SetOptStat("ne");//mri"); 
    h_q2w[idtyp]=(TH2F*)gDirectory->Get("h");
    h_q2w[idtyp]->SetName(TString::Format("q2w_%s_top%d",dtyp_name[idtyp].Data(),top));
    h_q2w[idtyp]->SetTitle(TString::Format("q2w_%s_top%d",dtyp_name[idtyp].Data(),top));
  }
  c->Close();

  char* outdir="/data/trivedia/e1f/study_d2pi_kinematics";
  TCanvas* cER=new TCanvas("q2w_ER","q2w_ER");
  h_q2w[ER]->Draw("colz");
  cER->SaveAs(TString::Format("%s/%s.eps",outdir,h_q2w[ER]->GetName()));

  TCanvas* cSR=new TCanvas("q2w_SR","q2w_SR");
  h_q2w[SR]->Draw("colz");
  cSR->SaveAs(TString::Format("%s/%s.eps",outdir,h_q2w[SR]->GetName()));
}
