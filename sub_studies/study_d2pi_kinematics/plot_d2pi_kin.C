void plot_d2pi_kin(int top=1){
  const int NDTYP=3;
  TFile* f[NDTYP];//0=ER,1=SR,2=ST
  enum{ER,SR,ST};
  TString dtyp_name[]={"ER","SR","ST"};
  int clr[]={kBlue,kRed,kGreen};
  f[0]=TFile::Open("$D2PIDIR_EXP/mon_d2pi/d2piR.root");
  f[1]=TFile::Open("$D2PIDIR_SIM/mon_d2pi/d2piR.root");
  f[2]=TFile::Open("$D2PIDIR_SIM/mon_d2pi/d2piT.root");
  TTree* t[NDTYP];
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
   if (idtyp==ER || idtyp==SR)t[idtyp]=(TTree*)f[idtyp]->Get("d2piR/tR");
   if (idtyp==ST)             t[idtyp]=(TTree*)f[idtyp]->Get("d2piT/tT");
  }

  const int NPART=4;//0=e,1=p,2=pip,3=pim
  TString part_name[]={"e","p","pip","pim"};
  TH1F* h_p[NDTYP][NPART]; 

  TCanvas*c = new TCanvas("c","c");//default canvas
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    c->cd();
    //plot momentum hists
    TString cut=TString::Format("top==%d",top);
    if (idtyp==ST) cut="";//no cut for ST
    t[idtyp]->Draw("p_e>>h_e(100,0,5)",cut);
    t[idtyp]->Draw("p_p>>h_p(100,0,5)",cut);
    t[idtyp]->Draw("p_pip>>h_pip(100,0,5)",cut);
    t[idtyp]->Draw("p_pim>>h_pim(100,0,5)",cut);

    gStyle->SetOptStat("ne");//mri"); 
    h_p[idtyp][0]=(TH1F*)gDirectory->Get("h_e");
    h_p[idtyp][1]=(TH1F*)gDirectory->Get("h_p");
    h_p[idtyp][2]=(TH1F*)gDirectory->Get("h_pip");
    h_p[idtyp][3]=(TH1F*)gDirectory->Get("h_pim");
    h_p[idtyp][0]->SetName(TString::Format("h_p_%s_%s",dtyp_name[idtyp].Data(),part_name[0].Data()));
    h_p[idtyp][0]->SetLineColor(clr[idtyp]);
    h_p[idtyp][1]->SetName(TString::Format("h_p_%s_%s",dtyp_name[idtyp].Data(),part_name[1].Data()));
    h_p[idtyp][1]->SetLineColor(clr[idtyp]);
    h_p[idtyp][2]->SetName(TString::Format("h_p_%s_%s",dtyp_name[idtyp].Data(),part_name[2].Data()));
    h_p[idtyp][2]->SetLineColor(clr[idtyp]);
    h_p[idtyp][3]->SetName(TString::Format("h_p_%s_%s",dtyp_name[idtyp].Data(),part_name[3].Data()));
    h_p[idtyp][3]->SetLineColor(clr[idtyp]);
  }
  c->Close();

  TCanvas* c_p=new TCanvas(TString::Format("c_p_%d",top),TString::Format("c_p_%d",top));
  //h_p[2][2]->Draw();
  c_p->Divide(2,2);
  for (int iprt=0;iprt<NPART;iprt++){
    c_p->cd(iprt+1);
    int nhist=0;
    for (int idtyp=0;idtyp<NDTYP;idtyp++){
      //if (i==0||i==1) continue;
      if (nhist==0) h_p[idtyp][iprt]->DrawNormalized("",1000);
      else          h_p[idtyp][iprt]->DrawNormalized("sames",1000);
      nhist+=1;
    }
  }
  
}
