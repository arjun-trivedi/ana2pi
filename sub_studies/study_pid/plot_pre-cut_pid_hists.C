{
  TFile* f[2];// 0=Exp,1=Sim
  TString name[]={"Exp","Sim"};
  f[0]=TFile::Open("$D2PIDIR_EXP/mon_pid_new/dpid.root");
  f[1]=TFile::Open("$D2PIDIR_SIM/mon_pid_new/dpid.root");
 
  TCanvas* c_bVp[2]; 
  TCanvas* c_dtVp[2];
  TH2F* h2_bVp_pos[2];
  TH2F* h2_bVp_neg[2];
  TH2F* h2_dtVp_pos_p[2];
  TH2F* h2_dtVp_pos_pip[2];
  TH2F* h2_dtVp_neg_pim[2];

  TCanvas *c=new TCanvas("c","c");//default canvas
  for (int i=0;i<2;i++){
    TTree* t=(TTree*)f[i]->Get("pid/tree");
    c->cd();
    //bVp
    t->Draw("b:p>>h_bVp_pos(100,0,5,100,0,1.2)","q==+1&&dc>0&&sc>0","colz");
    t->Draw("b:p>>h_bVp_neg(100,0,5,100,0,1.2)","q==-1&&dc>0&&sc>0","colz");
    //dtVp under particle assumption that all q==1 = proton or pip
    t->Draw("dt_p:p>>h_dtVp_pos_p(100,0,5,100,-5,5)","q==+1&&dc>0&&sc>0","colz");
    t->Draw("dt_pip:p>>h_dtVp_pos_pip(100,0,5,100,-5,5)","q==+1&&dc>0&&sc>0","colz");
    //dtVp under particle assumption that all q==-1 = pim
    t->Draw("dt_pim:p>>h_dtVp_neg_pim(100,0,5,100,-5,5)","q==-1&&dc>0&&sc>0","colz");
   
    gStyle->SetOptStat("ne");//mri"); 
    h2_bVp_pos[i]=(TH2F*)gDirectory->Get("h_bVp_pos");
    h2_bVp_pos[i]->SetName(TString::Format("h_bVp_pos_%s",name[i].Data()));
    h2_bVp_neg[i]=(TH2F*)gDirectory->Get("h_bVp_neg");
    h2_bVp_neg[i]->SetName(TString::Format("h_bVp_neg_%s",name[i].Data()));
    c_bVp[i]=new TCanvas(TString::Format("bVp_%s",name[i].Data()),TString::Format("bVp_%s",name[i].Data()));
    c_bVp[i].Divide(1,2);
    c_bVp[i].cd(1);
    h2_bVp_pos[i]->Draw("colz");
    c_bVp[i].cd(2);
    h2_bVp_neg[i]->Draw("colz");

    h2_dtVp_pos_p[i]=(TH2F*)gDirectory->Get("h_dtVp_pos_p");
    h2_dtVp_pos_p[i]->SetName(TString::Format("h_dtVp_pos_p_%s",name[i].Data()));
    h2_dtVp_pos_pip[i]=(TH2F*)gDirectory->Get("h_dtVp_pos_pip");
    h2_dtVp_pos_pip[i]->SetName(TString::Format("h_dtVp_pos_pip_%s",name[i].Data()));
    h2_dtVp_neg_pim[i]=(TH2F*)gDirectory->Get("h_dtVp_neg_pim");
    h2_dtVp_neg_pim[i]->SetName(TString::Format("h_dtVp_neg_pim_%s",name[i].Data()));
    c_dtVp[i]=new TCanvas(TString::Format("dtVp_%s",name[i].Data()),TString::Format("dtVp_%s",name[i].Data()));
    c_dtVp[i].Divide(2,2);
    c_dtVp[i].cd(1);
    h2_dtVp_pos_p[i]->Draw("colz");
    c_dtVp[i].cd(3);
    h2_dtVp_pos_pip[i]->Draw("colz");
    c_dtVp[i].cd(2);
    h2_dtVp_neg_pim[i]->Draw("colz");
  }
}
