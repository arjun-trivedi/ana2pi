void plot_elast_kinematics(){
  const int NDTYP=3;
  TFile* f[NDTYP];//0=ER,1=SR,2=ST
  enum{ER,SR,ST};
  TString dtyp_name[]={"ER","SR","ST"};
  int clr[]={kBlue,kRed,kGreen};
  f[0]=TFile::Open("$DELASTDIR_EXP/study_elast_xsec_081715/delastR_tree.root");
  f[1]=TFile::Open("$DELASTDIR_SIM/study_elast_xsec_081715/recon/delastR_tree.root");
  f[2]=TFile::Open("$DELASTDIR_SIM/study_elast_xsec_081715/evtgen/delastT_tree.root");
  TTree* t[NDTYP];
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
   if (idtyp==ER || idtyp==SR)t[idtyp]=(TTree*)f[idtyp]->Get("delast/cut/t");
   if (idtyp==ST)             t[idtyp]=(TTree*)f[idtyp]->Get("delast/tT");//! This is creater after cut for ST!
  }

  const int NPART=2;//0=e,1=p
  TString part_name[]={"e","p"};
  TH1F* h_p[NDTYP][NPART]; 

  TCanvas*c = new TCanvas("c","c");//default canvas
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    c->cd();
    //plot momentum hists
    t[idtyp]->Draw("p_e>>h_e(100,0,5)");
    t[idtyp]->Draw("p_p>>h_p(100,0,5)");

    gStyle->SetOptStat("ne");//mri"); 
    h_p[idtyp][0]=(TH1F*)gDirectory->Get("h_e");
    h_p[idtyp][1]=(TH1F*)gDirectory->Get("h_p");
    h_p[idtyp][0]->SetName(TString::Format("h_p_%s_%s",dtyp_name[idtyp].Data(),part_name[0].Data()));
    h_p[idtyp][0]->SetLineColor(clr[idtyp]);
    h_p[idtyp][1]->SetName(TString::Format("h_p_%s_%s",dtyp_name[idtyp].Data(),part_name[1].Data()));
    h_p[idtyp][1]->SetLineColor(clr[idtyp]);
  }
  c->Close();

  TCanvas* c_p=new TCanvas(TString::Format("c_p"),TString::Format("c_p"));
  //h_p[2][2]->Draw();
  c_p->Divide(2,1);
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
