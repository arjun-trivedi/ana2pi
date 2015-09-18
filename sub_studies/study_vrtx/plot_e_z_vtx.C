void plot_e_z_vtx(Long64_t nentries=1000000000){
  TFile *fexp=TFile::Open(TString::Format("%s/study_elast_xsec_081715/delastR_tree.root",getenv("DELASTDIR_EXP")));
  TFile *fsim=TFile::Open(TString::Format("%s/study_elast_xsec_081715/recon/delastR_tree.root",getenv("DELASTDIR_SIM")));
  TTree* texp=(TTree*)fexp->Get("delast/monitor/t");
  TTree* tsim=(TTree*)fsim->Get("delast/monitor/t");

  TH1F* hexp[6];
  TH1F* hsim[6];
  for (int isctr=0;isctr<6;isctr++){
    printf("Processign sector %d\n",isctr+1);
    TCut cut_sctr(TString::Format("sector==%d",isctr+1));
    //! exp
    printf("Processign exp\n");
    texp->Draw("vz_e>>hexp(100,-40,-10)",cut_sctr,"",nentries);
    TH1F* htmp_exp=gDirectory->Get("hexp");
    hexp[isctr]=(TH1F*)htmp_exp->Clone(TString::Format("hexp%d",isctr+1));
    hexp[isctr]->SetLineColor(isctr+1);
    //! sim
    printf("Processign sim\n");
    tsim->Draw("vz_e>>hsim(100,-40,-10)",cut_sctr,"",nentries);
    TH1F* htmp_sim=gDirectory->Get("hsim");
    hsim[isctr]=(TH1F*)htmp_sim->Clone(TString::Format("hsim%d",isctr+1));
    hsim[isctr]->SetLineColor(isctr+1);
  }
  TCanvas *cexp=new TCanvas("exp","exp");
  TLegend* leg_exp=new TLegend(0.8,0.8,1,1);
  for (int isctr=0;isctr<6;isctr++){
    if (isctr==0) hexp[isctr]->Draw();
    else hexp[isctr]->Draw("same");
    leg_exp->AddEntry(hexp[isctr],hexp[isctr]->GetName(),"l");
  }
  leg_exp->Draw();
  
  TCanvas *csim=new TCanvas("sim","sim");
  TLegend* leg_sim=new TLegend(0.8,0.8,1,1);
  for (int isctr=0;isctr<6;isctr++){
    if (isctr==0) hsim[isctr]->Draw();
    else hsim[isctr]->Draw("same");
    leg_sim->AddEntry(hsim[isctr],hsim[isctr]->GetName(),"l");
  }
  leg_sim->Draw();

  //sector:exp vs sim
  gStyle->SetOptStat("nmMrR");
  Float_t z_vtx_cut_min_exp[6]={-28.25,-27.50,-27.50,-27.50,-28.25,-28.75};
  Float_t z_vtx_cut_max_exp[6]={-23.00,-22.50,-22.25,-22.50,-23.00,-23.50};
  Float_t z_vtx_cut_min_sim[6]={-27.50,-27.50,-27.50,-27.50,-27.50,-27.50};
  Float_t z_vtx_cut_max_sim[6]={-22.50,-22.50,-22.50,-22.50,-22.50,-22.50};
  TLine* lvtxcutmin_exp[6];
  TLine* lvtxcutmax_exp[6];
  TLine* lvtxcutmin_sim[6];
  TLine* lvtxcutmax_sim[6]; 
  TCanvas* csctr=new TCanvas("sctr","sctr");
  TLegend* leg_sctr=new TLegend(0.8,0.8,1,1);
  csctr->Divide(3,2);
  for (int isctr=0;isctr<6;isctr++){
    csctr->cd(isctr+1);
    hexp[isctr]->SetLineColor(kBlue);
    hsim[isctr]->SetLineColor(kRed);
    hexp[isctr]->Sumw2();
    hsim[isctr]->Sumw2();
    TH1F* hexpn=hexp[isctr]->DrawNormalized("hist E",1000);
    TH1F* hsimn=hsim[isctr]->DrawNormalized("hist E sames",1000);
    float ymax=hexpn->GetMaximum();
    ymax=ymax>hsimn->GetMaximum()?ymax:hsimn->GetMaximum();
    ymax+=10;
    hexpn->SetMaximum(ymax);
    hsimn->SetMaximum(ymax);
    gPad->Update();
    //! Draw z vtx cut
    lvtxcutmin_exp[isctr]=new TLine(z_vtx_cut_min_exp[isctr],0,z_vtx_cut_min_exp[isctr],ymax);
    lvtxcutmax_exp[isctr]=new TLine(z_vtx_cut_max_exp[isctr],0,z_vtx_cut_max_exp[isctr],ymax);
    lvtxcutmin_sim[isctr]=new TLine(z_vtx_cut_min_sim[isctr],0,z_vtx_cut_min_sim[isctr],ymax);
    lvtxcutmax_sim[isctr]=new TLine(z_vtx_cut_max_sim[isctr],0,z_vtx_cut_max_sim[isctr],ymax);
    lvtxcutmin_exp[isctr]->SetLineColor(kBlue);
    lvtxcutmax_exp[isctr]->SetLineColor(kBlue);
    lvtxcutmin_sim[isctr]->SetLineColor(kRed);
    lvtxcutmax_sim[isctr]->SetLineColor(kRed);
    lvtxcutmin_exp[isctr]->SetLineStyle(5);
    lvtxcutmax_exp[isctr]->SetLineStyle(5);
    lvtxcutmin_sim[isctr]->SetLineStyle(5);
    lvtxcutmax_sim[isctr]->SetLineStyle(5);
    lvtxcutmin_exp[isctr]->Draw("same");
    lvtxcutmax_exp[isctr]->Draw("same");
    lvtxcutmin_sim[isctr]->Draw("same");
    lvtxcutmax_sim[isctr]->Draw("same");
    if (isctr==0){
      leg_sctr->AddEntry(hexp[isctr],"exp","l");
      leg_sctr->AddEntry(hsim[isctr],"sim","l");
      leg_sctr->Draw();
    }
  }
}
