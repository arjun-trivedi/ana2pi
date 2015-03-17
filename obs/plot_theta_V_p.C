void plot_theta_V_p(bool logz=kTRUE){
  TFile* _file0= new TFile(gSystem->ExpandPathName("$D2PIDIR_EXP/mon_det_ineff/d2piRmon.root"));
  _file0->cd("d2piR");

  gStyle->SetOptStat(0);
  TString prtcl[]={"e","p","pip","pim"};
  TCanvas *c[4];
  for (int i=0;i<4;i++){
    c[i]=new TCanvas(TString::Format("c%s",prtcl[i].Data()),TString::Format("%s",prtcl[i].Data()));
    c[i]->Divide(3,2);
  }

  Int_t sctr_phi_min[]={0, 30, 90,  150, 210, 270};
  Int_t sctr_phi_max[]={0, 90 ,150, 210, 270, 330};

  Int_t theta_min[]={0, 0, 0,  0};
  Int_t theta_max[]={60,60,120,120};

  TCut cut_top("top==1");
  for (int i=0;i<4;i++){
    for (int isctr=0;isctr<6;isctr++){
      //TCut* cut_sctr;
      /*if (isctr==0){
        cut_sctr= new TCut(TString::Format("(phi_%s>0&&phi_%s<30)||(phi_%s>330&&phi_%s<360)",
                                           prtcl[i].Data(),prtcl[i].Data(),prtcl[i].Data(),prtcl[i].Data()));
        cut_sctr= new TCut(TString::Format("(phi_%s>-180&&phi_%s<-150)||(phi_%s>150&&phi_%s<180)",
                                           prtcl[i].Data(),prtcl[i].Data(),prtcl[i].Data(),prtcl[i].Data()));
      }else{
        cut_sctr= new TCut(TString::Format("phi_%s>%d&&phi_%s<%d",
                                           prtcl[i].Data(),sctr_phi_min[isctr],prtcl[i].Data(),sctr_phi_max[isctr]));
      }*/
      TCut cut_sctr(TString::Format("sector_%s==%d",prtcl[i].Data(),isctr+1));
      TPad* pad=c[i]->cd(isctr+1);
      if (logz)
        pad->SetLogz();      
      TString cmd=TString::Format("theta_%s:p_%s>>hthetaVp_%s_sctr%d(100,0,5,100,%d,%d)",
                                   prtcl[i].Data(),prtcl[i].Data(),prtcl[i].Data(),isctr+1,theta_min[i],theta_max[i]);
      cout <<"cmd="<<cmd.Data()<<endl;
      tR->Draw(cmd,cut_sctr&&cut_top,"colz");
    }
  }
}
