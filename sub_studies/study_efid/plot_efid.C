void plot_efid(TString slctn,TString cutlvl){
  const int NDTYP=2;
  enum {EXP,SIM};
  TFile* f[NDTYP];// 0=Exp,1=Sim
  TString dtyp_name[]={"Exp","Sim"};

  if(slctn=="2pi"){
    //f[EXP]=TFile::Open("$D2PIDIR_EXP/data_efid_071415/defid_after_d2piR.root");
    f[EXP]=TFile::Open("$D2PIDIR_EXP/data_efid_071415/defid.root");
    f[SIM]=TFile::Open("$D2PIDIR_SIM/data_efid_071415/defid.root");
  }else if (slctn=="elastic"){
    //f[EXP]=TFile::Open("$D2PIDIR_EXP/data_eid_070915/elastic/deid.root");
    f[EXP]=TFile::Open("$DELASTDIR_EXP/test_efid_elast_070715/defid_after_delastR.root");
    f[SIM]=TFile::Open("$D2PIDIR_SIM/data_eid_070915/elastic/deid.root");
  }
  else{
    cout <<"slctn not recognized. Going to Exit"<<endl;
  }
  printf("%s\n",slctn.Data());
 
  if (cutlvl!="monitor" && cutlvl!="cut"){
    printf("cutlvl=%s not recognized. Only mon and cut allowed\n",cutlvl.Data());
    return;
  }
  
  //! Prepare OUTDIR & fout.root
  TString wspace=getenv("WORKSPACE");
  TString OUTDIR=TString::Format("%s/ana2pi/sub_studies/study_efid/hists/%s/%s",wspace.Data(),slctn.Data(),cutlvl.Data());
  gSystem->mkdir(OUTDIR,1);
  TFile* fout[NDTYP];
  fout[EXP]=new TFile(TString::Format("%s/fexp.root",OUTDIR.Data()),"RECREATE");
  fout[SIM]=new TFile(TString::Format("%s/fsim.root",OUTDIR.Data()),"RECREATE");
  

  const int NHST=3;
  TString hst_name[]= {"hephiVetheta","hescyVescx","heecyVeecx"};
  TString hst_title[]={"hephiVetheta","hescyVescx","heecyVeecx"};
  
  const int NSCTR=6;
 
  const int NVW=2;
  enum {LIN,LOG};
  TH2F* h[NDTYP][NHST][NSCTR][NVW];
  TCanvas* c[NDTYP][NHST][NVW];

  gROOT->ProcessLine(".L misc-CINT-funcs.C");

  //! Start getting hists indexed by dtyp,hst,sctr
  gStyle->SetOptStat("ne");
  for (int idtyp=0;idtyp<NDTYP;idtyp++){ //begin dtyp loop
    //if (idtyp>0)continue;
    cout<<dtyp_name[idtyp]<<endl;
    for (int ihst=0;ihst<NHST;ihst++){// begin nhst loop
      //if (ihst>1)continue;
      for (int isctr=0;isctr<NSCTR;isctr++){ //begin sctr loop
        TString hist=TString::Format("efid/%s/sector%d/%s",cutlvl.Data(),isctr+1,hst_name[ihst].Data());
        cout<<hist<<endl;
        h[idtyp][ihst][isctr][LIN]=(TH2F*)f[idtyp]->Get(hist);
        h[idtyp][ihst][isctr][LIN]->SetName(TString::Format("%s_%s_s%d",dtyp_name[idtyp].Data(),hst_name[ihst].Data(),isctr+1));
      }//end sctr loop
    }// end nhst loop
  }//end dtyp loop

  //! Now plot hists  
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    //if (idtyp>0)continue;
    fout[idtyp]->cd();
    for (int ihst=0;ihst<NHST;ihst++){
      TString outdir=TString::Format("%s/%s/%s",OUTDIR.Data(),dtyp_name[idtyp].Data(),hst_name[ihst].Data());
      gSystem->mkdir(outdir,1);
      TString cname=TString::Format("c_%s_%s",dtyp_name[idtyp].Data(),hst_name[ihst].Data());
      c[idtyp][ihst][LIN]=new TCanvas(TString::Format("%s_lin",cname.Data()));
      c[idtyp][ihst][LOG]=new TCanvas(TString::Format("%s_log",cname.Data()));
      c[idtyp][ihst][LIN]->Divide(3,2);
      c[idtyp][ihst][LOG]->Divide(3,2);
      for (int isctr=0;isctr<NSCTR;isctr++){
        c[idtyp][ihst][LIN]->cd(isctr+1);
        h[idtyp][ihst][isctr][LIN]->Draw("colz");
        c[idtyp][ihst][LOG]->cd(isctr+1);
        h[idtyp][ihst][isctr][LOG]=norm2D(h[idtyp][ihst][isctr][LIN]);
        h[idtyp][ihst][isctr][LOG]->Draw("colz");
        //TH2F* htmp=norm2D(h[idtyp][ihst][isctr]);
        //htmp->Draw("colz");
      }
      c[idtyp][ihst][LIN]->SaveAs(TString::Format("%s/%s.jpg",outdir.Data(),c[idtyp][ihst][LIN]->GetName()));
      c[idtyp][ihst][LIN]->Write();
      c[idtyp][ihst][LOG]->SaveAs(TString::Format("%s/%s.jpg",outdir.Data(),c[idtyp][ihst][LOG]->GetName()));
      c[idtyp][ihst][LOG]->Write();
    }
  }
  //! Close fouts
  fout[EXP]->Close();
  fout[SIM]->Close();
}
