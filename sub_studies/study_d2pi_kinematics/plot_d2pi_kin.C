void plot_d2pi_kin(int top=1){
  const int NDTYP=3;
  TFile* f[NDTYP];//0=ER,1=SR,2=ST
  TString dtyp_name[]={"ER","SR","ST"};
  int clr[]={kBlue,kRed,kGreen};
  f[0]=TFile::Open("$D2PIDIR_EXP/mon_d2pi/d2piR.root");
  f[1]=TFile::Open("$D2PIDIR_SIM/mon_d2pi/d2piR.root");
  f[2]=TFile::Open("$D2PIDIR_SIM/mon_d2pi/d2piT.root");
  TTree* t[NDTYP];
  for (int i=0;i<NDTYP;i++){
   if ( i==0 || i==1 )t[i]=(TTree*)f[i]->Get("d2piR/tR");
   if (i==2)          t[i]=(TTree*)f[i]->Get("d2piT/tT");
  }

  const int NPART=3;//0=p,1=pip,2=pim
  TString part_name[]={"p","pip","pim"};
  TH1F* h_p[NDTYP][NPART]; 

  TCanvas*c = new TCanvas("c","c");//default canvas
  for (int i=0;i<NDTYP;i++){
    c->cd();
    //plot momentum hists
    TString cut=TString::Format("top==%d",top);
    if (i==2) cut="";//no cut for ST
    t[i]->Draw("p_p>>h_p(100,0,5)",cut);
    t[i]->Draw("p_pip>>h_pip(100,0,5)",cut);
    t[i]->Draw("p_pim>>h_pim(100,0,5)",cut);

    gStyle->SetOptStat("ne");//mri"); 
    h_p[i][0]=(TH1F*)gDirectory->Get("h_p");
    h_p[i][1]=(TH1F*)gDirectory->Get("h_pip");
    h_p[i][2]=(TH1F*)gDirectory->Get("h_pim");
    h_p[i][0]->SetName(TString::Format("h_p_%s_%s",dtyp_name[i].Data(),part_name[0].Data()));
    h_p[i][0]->SetLineColor(clr[i]);
    h_p[i][1]->SetName(TString::Format("h_p_%s_%s",dtyp_name[i].Data(),part_name[1].Data()));
    h_p[i][1]->SetLineColor(clr[i]);
    h_p[i][2]->SetName(TString::Format("h_p_%s_%s",dtyp_name[i].Data(),part_name[2].Data()));
    h_p[i][2]->SetLineColor(clr[i]);
  }
  c->Close();

  TCanvas* c_p=new TCanvas(TString::Format("c_p_%d",top),TString::Format("c_p_%d",top));
  //h_p[2][2]->Draw();
  c_p->Divide(2,2);
  for (int j=0;j<NPART;j++){
    c_p->cd(j+1);
    int nhist=0;
    for (int i=0;i<NDTYP;i++){
      //if (i==0||i==1) continue;
      if (nhist==0) h_p[i][j]->DrawNormalized("",1000);
      else          h_p[i][j]->DrawNormalized("sames",1000);
      nhist+=1;
    }
  }
  
}
