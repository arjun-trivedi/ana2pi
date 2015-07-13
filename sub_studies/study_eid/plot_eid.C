TH2F* norm(TH2F* h2){
  TH2F* h2c=h2->Clone(TString::Format("%s_norm",h2->GetName()));
  int nxbins=h2c->GetNbinsX();
  int nybins=h2c->GetNbinsY();

  for (int ixbin=0;ixbin<nxbins;ixbin++){
    TH1F* hy=(TH1F*)h2c->ProjectionY("_py",ixbin+1,ixbin+1);
    float max=hy->GetMaximum();
    if (max==0) continue;
    //! Now normalize h2 for this projection
    for (int iybin=0;iybin<nybins;iybin++){
      float binc=h2c->GetBinContent(ixbin+1,iybin+1);
      binc=binc/max;
      h2c->SetBinContent(ixbin+1,iybin+1,binc);
    }
  }
  return h2c;
}

void plot_eid(TString slctn,TString seq,int nentries=1000000000){
  const int NDTYP=2;
  enum {EXP,SIM};
  TFile* f[NDTYP];// 0=Exp,1=Sim
  TString dtyp_name[]={"Exp","Sim"};

  TCut slctn_cut;
  if(slctn=="all"){
    slctn_cut="";
    f[EXP]=TFile::Open("$D2PIDIR_EXP/data_eid_070915/deid.root");
    f[SIM]=TFile::Open("$D2PIDIR_SIM/data_eid_070915/deid.root");
  }else if (slctn=="ana2pi"){
    slctn_cut="W<2.125";
    f[EXP]=TFile::Open("$D2PIDIR_EXP/data_eid_070915/deid.root");
    f[SIM]=TFile::Open("$D2PIDIR_SIM/data_eid_070915/deid.root");
  }
  else if (slctn=="elastic"){
    slctn_cut="W<1";
    f[EXP]=TFile::Open("$D2PIDIR_EXP/data_eid_070915/elastic/deid.root");
    f[SIM]=TFile::Open("$D2PIDIR_SIM/data_eid_070915/elastic/deid.root");
  }
  else{
    cout <<"slctn not recognized. Going to Exit"<<endl;
  }
  printf("%s=%s\n",slctn.Data(),slctn_cut.GetTitle());
  
  TTree* t[NDTYP];
  if (seq=="mon"){
    t[EXP]=(TTree*)f[EXP]->Get("eid/monitor/t");
    t[SIM]=(TTree*)f[SIM]->Get("eid/monitor/t");
  }else if (seq=="cut"){
    t[EXP]=(TTree*)f[EXP]->Get("eid/cut/t");
    t[SIM]=(TTree*)f[SIM]->Get("eid/cut/t");
  }else if (seq=="evtsel"){
    t[EXP]=(TTree*)f[EXP]->Get("eid2/cut/t");
    t[SIM]=(TTree*)f[SIM]->Get("eid2/cut/t");
  }else{
    cout <<"seq not recognized. Going to Exit"<<endl;
  }
  printf("%s=%s\n",seq.Data(),t[EXP]->GetTitle());

  TString wspace=getenv("WORKSPACE");
  gROOT->ProcessLine(".L ${WORKSPACE}/ana2pi/epconfig.cpp++");
  gROOT->ProcessLine(".L ${WORKSPACE}/ana2pi/eid.cpp++");
  Eid* eid_tool[NDTYP];
  eid_tool[EXP]=new Eid(TString::Format("%s/ana2pi/eid/eid.exp.out",wspace.Data()).Data());
  eid_tool[SIM]=new Eid(TString::Format("%s/ana2pi/eid/eid.mc.out",wspace.Data()).Data());

  const int NHST=4;//removed "hsfoutVsfin","hetotVp","ecU","ecV","ecW"
  TString hst_name[]={"heoutVein","hsfoutVsfin","hsftotVp",             "hnphe","hCCthetaVCCseg"};
  TString hst_title[]={"#E_{out} vs. #E_{in}","hsfoutVsfin","hsftotVp", "hnphe","hCCthetaVCCseg"};
  TString hst_drw_cmd[]={
                         TString::Format("ec_eo:ec_ei>>%s(150,0,1,150,0,1)",        hst_name[0].Data()),
                         TString::Format("ec_eo/p:ec_ei/p>>%s(150,0,0.5,150,0,0.5)",hst_name[1].Data()),
                         TString::Format("etot/p:p>>%s(160,0,5,150,0,0.5)",         hst_name[2].Data()),
                         //TString::Format("etot:p>>%s(160,0,5,150,0,2)",             hst_name[3].Data()),
                         //TString::Format("ecU>>%s(100,0,450)",                        hst_name[2].Data()), 
                         //TString::Format("ecV>>%s(100,0,450)",                        hst_name[3].Data()),
                         //TString::Format("ecW>>%s(100,0,450)",                        hst_name[4].Data()),
                         TString::Format("nphe>>%s(100,0,300)",                     hst_name[3].Data()),
                         TString::Format("cc_theta:cc_segm>>%s(20,0,20,100,0,50)",   hst_name[4].Data())
                         };
  
  const int NSCTR=6;
 
  TH2F* h[NDTYP][NHST][NSCTR];
  TCanvas* c[NDTYP][NHST];

  //! Start making plots indexed by dtyp,hst,sctr
  TCut pre_cut("etot>0.001 && ec_ei>0.001 && ec_eo>0.001");
  TCanvas *c1=new TCanvas("c1","c1");//default canvas
  for (int idtyp=0;idtyp<NDTYP;idtyp++){ //begin dtyp loop
    //if (idtyp>0)continue;
    cout<<dtyp_name[idtyp]<<endl;
    for (int ihst=0;ihst<NHST;ihst++){// begin nhst loop
      //if (ihst>1)continue;
      for (int isctr=0;isctr<NSCTR;isctr++){ //begin sctr loop
        TString outdir=TString::Format("%s/ana2pi/sub_studies/study_eid/hists/%s/%s/%s/%s",wspace.Data(),slctn.Data(),seq.Data(),dtyp_name[idtyp].Data(),hst_name[ihst].Data());
        gSystem->mkdir(outdir,1);
        c1->cd();
        TCut sctr_cut=TString::Format("sector==%d",isctr+1);
        TCut cut=pre_cut&&sctr_cut&&slctn_cut;
        cout<<hst_drw_cmd[ihst]<<endl;
        t[idtyp]->Draw(hst_drw_cmd[ihst],cut,"colz",nentries);
   
        gStyle->SetOptStat("ne");//mri"); 
        if (ihst<=2||ihst==4){//! i.e. 2D hists
          TH2F* h2tmp=(TH2F*)gDirectory->Get(hst_name[ihst]);
          h[idtyp][ihst][isctr]=(TH2F*)h2tmp->Clone();
        }else{//! i.e. 1D hists
          TH1F* h1tmp=(TH1F*)gDirectory->Get(hst_name[ihst]);
          h[idtyp][ihst][isctr]=(TH2F*)h1tmp->Clone();
        }
        h[idtyp][ihst][isctr]->SetName(TString::Format("%s_%s_s%d",       hst_name[ihst].Data(),dtyp_name[idtyp].Data(),isctr+1));
        h[idtyp][ihst][isctr]->SetTitle(TString::Format("%s_%s", pre_cut.GetTitle(),slctn.Data()));
        if (ihst<=2||ihst==4){//! i.e. 2D hists
          h[idtyp][ihst][isctr]->Draw("colz");
        }else{//! i.e. 1D hists
          h[idtyp][ihst][isctr]->Draw();
        }
        if (ihst==2)eid_tool[idtyp]->DrawSFcuts(isctr+1);
        c1->SaveAs(TString::Format("%s/%s.jpg",outdir.Data(),h[idtyp][ihst][isctr]->GetName()));
      }//end sctr loop
    }// end nhst loop
  }//end dtyp loop

  //! Now plot hists for intereactive viewing 
  //! Seperate canvases for each dtyp for 2D
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    for (int ihst=0;ihst<NHST;ihst++){
      //if (ihst>1) continue;
      if (ihst>2&&ihst!=4) continue;
      cout<<hst_name[ihst]<<endl;
      TString cname=TString::Format("c_%s_%s",dtyp_name[idtyp].Data(),hst_name[ihst].Data());
      c[idtyp][ihst]=new TCanvas(cname,cname);
      c[idtyp][ihst]->Divide(3,2);
      for (int isctr=0;isctr<NSCTR;isctr++){
        c[idtyp][ihst]->cd(isctr+1);
        //h[idtyp][ihst][isctr]=norm(h[idtyp][ihst][isctr]);
        h[idtyp][ihst][isctr]->Draw("colz");
        if (ihst==2)eid_tool[idtyp]->DrawSFcuts(isctr+1);
      }
    }
  }

  //! Same canvas for each for 1D
  for (int ihst=0;ihst<NHST;ihst++){
    if (ihst<=2||ihst==4) continue;
    TString cname=TString::Format("c_%s",hst_name[ihst].Data());
    c[EXP][ihst]=new TCanvas(cname,cname);
    c[EXP][ihst]->Divide(3,2);
    for (int isctr=0;isctr<NSCTR;isctr++){
      c[EXP][ihst]->cd(isctr+1);
      h[EXP][ihst][isctr]->Draw();
      /*h[EXP][ihst][isctr]->SetLineColor(kBlue);
      h[SIM][ihst][isctr]->SetLineColor(kRed);
      h[EXP][ihst][isctr]->DrawNormalized("",1000);
      h[SIM][ihst][isctr]->DrawNormalized("sames",1000);*/
    }
  }
}
