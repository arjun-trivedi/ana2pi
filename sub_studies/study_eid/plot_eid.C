void plot_eid(TString slctn,int nentries=1000000){
  TCut slctn_cut;
  if      (slctn=="ep_enrchmnt") slctn_cut="nphe>30 && ec_ei>0.06";//removed etot>0.16
  else if (slctn=="ana2pi")      slctn_cut="W<2.125";
  else if (slctn=="elastic")     slctn_cut="W<1";
  else{
    cout <<"slctn not recognized. Going to Exit"<<endl;
  }
  printf("%s=%s\n",slctn.Data(),slctn_cut.GetTitle());

  const int NDTYP=2;
  enum {EXP,SIM};
  TFile* f[NDTYP];// 0=Exp,1=Sim
  TString dtyp_name[]={"Exp","Sim"};
  //f[EXP]=TFile::Open("$D2PIDIR_EXP/data_eid_070915/deid.root");
  //f[SIM]=TFile::Open("$D2PIDIR_SIM/data_eid_070915/deid.root");
  f[EXP]=TFile::Open("$D2PIDIR_EXP/data_eid_070915/elastic/deid.root");
  f[SIM]=TFile::Open("$D2PIDIR_SIM/data_eid_070915/elastic/deid.root");

  const int NHST=7;//removed "hsfoutVsfin","hetotVp"
  TString hst_name[]={"heoutVein","hsftotVp",             "hU","hV","hW","hnphe","hCCthetaVCCseg"};
  TString hst_title[]={"#E_{out} vs. #E_{in}","hsftotVp", "hU","hV","hW","hnphe","hCCthetaVCCseg"};
  TString hst_drw_cmd[]={
                         TString::Format("ec_eo:ec_ei>>%s(150,0,1,150,0,1)",        hst_name[0].Data()),
                         //TString::Format("ec_eo/p:ec_ei/p>>%s(150,0,0.5,150,0,0.5)",hst_name[1].Data()),
                         TString::Format("etot/p:p>>%s(160,0,5,150,0,0.5)",         hst_name[1].Data()),
                         //TString::Format("etot:p>>%s(160,0,5,150,0,2)",             hst_name[3].Data()),
                         TString::Format("ecU>>%s(100,0,450)",                        hst_name[2].Data()), 
                         TString::Format("ecV>>%s(100,0,450)",                        hst_name[3].Data()),
                         TString::Format("ecW>>%s(100,0,450)",                        hst_name[4].Data()),
                         TString::Format("nphe>>%s(100,0,450)",                     hst_name[5].Data()),
                         TString::Format("cc_theta:cc_segm>>%s(20,0,20,100,0,50)",   hst_name[6].Data())
                         };
  
  const int NSCTR=6;
 
  TH2F* h[NDTYP][NHST][NSCTR];
  TCanvas* c[NDTYP][NHST];

  TString wspace=getenv("WORKSPACE");
  gROOT->ProcessLine(".L ${WORKSPACE}/ana2pi/epconfig.cpp++");
  gROOT->ProcessLine(".L ${WORKSPACE}/ana2pi/eid.cpp++");
  Eid* eid_tool[NDTYP];

  //! Start making plots indexed by dtyp,hst,sctr
  TCut pre_cut("etot>0.001 && ec_ei>0.001 && ec_eo>0.001");
  TCanvas *c1=new TCanvas("c1","c1");//default canvas
  for (int idtyp=0;idtyp<NDTYP;idtyp++){ //begin dtyp loop
    //if (idtyp>0)continue;
    cout<<dtyp_name[idtyp]<<endl;
    //! Prepare eid_tool
    TString cut_file_name;
    if      (strcmp(dtyp_name[idtyp].Data(),"Exp")==0) cut_file_name="eid.exp.out";
    else if (strcmp(dtyp_name[idtyp].Data(),"Sim")==0) cut_file_name="eid.mc.out";
    eid_tool[idtyp]=new Eid(TString::Format("%s/ana2pi/eid/%s",wspace.Data(),cut_file_name.Data()).Data());
    //! Get hists from Tree
    //TTree* t=(TTree*)f[idtyp]->Get("eid/monitor/t");
    TTree* t=(TTree*)f[idtyp]->Get("eid2/monitor/t");
    for (int ihst=0;ihst<NHST;ihst++){// begin nhst loop
      //if (ihst>2)continue;
      for (int isctr=0;isctr<NSCTR;isctr++){ //begin sctr loop
        c1->cd();
        //TString sctr_cut=TString::Format("sector==%d");// && W<2.125 && 60<ecU && ecU<400 && ecV<360 && ecW <395",isctr+1);
        TCut sctr_cut=TString::Format("sector==%d",isctr+1);
        TCut cut=pre_cut&&sctr_cut&&slctn_cut;
        cout<<hst_drw_cmd[ihst]<<endl;
        t->Draw(hst_drw_cmd[ihst],cut,"colz",nentries);
   
        gStyle->SetOptStat("ne");//mri"); 
        if (ihst<=1||ihst==6){//! i.e. 2D hists
          TH2F* h2tmp=(TH2F*)gDirectory->Get(hst_name[ihst]);
          h[idtyp][ihst][isctr]=(TH2F*)h2tmp->Clone();
        }else{//! i.e. 1D hists
          TH1F* h1tmp=(TH1F*)gDirectory->Get(hst_name[ihst]);
          h[idtyp][ihst][isctr]=(TH2F*)h1tmp->Clone();
        }
        h[idtyp][ihst][isctr]->SetName(TString::Format("%s_%s_s%d",       hst_name[ihst].Data(),dtyp_name[idtyp].Data(),isctr+1));
        //h[idtyp][ihst][isctr]->SetTitle(TString::Format("%s(%s,s%d)[%s]", hst_name[ihst].Data(),dtyp_name[idtyp].Data(),isctr+1,cut.GetTitle()));
        h[idtyp][ihst][isctr]->SetTitle(TString::Format("%s_%s", pre_cut.GetTitle(),slctn.Data()));
      }//end sctr loop
      //! Now plot hists
      /*TString cname=TString::Format("c_%s_%s",hst_name[ihst].Data(),dtyp_name[idtyp].Data());
      c[idtyp][ihst]=new TCanvas(cname,cname);
      c[idtyp][ihst]->Divide(3,2);
      for (int isctr=0;isctr<NSCTR;isctr++){
        c[idtyp][ihst]->cd(isctr+1);
        h[idtyp][ihst][isctr]->Draw("colz");
        if (ihst==2)eid_tool[idtyp]->DrawSFcuts(isctr+1);
      }*/
    }// end nhst loop
  }//end dtyp loop
  cout<<"here"<<endl;
  //! Now plot hists
  //! Seperate canvases for each dtyp for 2D
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    for (int ihst=0;ihst<NHST;ihst++){
      if (ihst>1&&ihst!=6) continue;
      cout<<hst_name[ihst]<<endl;
      TString cname=TString::Format("c_%s_%s",dtyp_name[idtyp].Data(),hst_name[ihst].Data());
      cout<<"cname="<<cname<<endl;
      c[idtyp][ihst]=new TCanvas(cname,cname);
      cout<<"here2"<<endl;
      c[idtyp][ihst]->Divide(3,2);
      for (int isctr=0;isctr<NSCTR;isctr++){
        c[idtyp][ihst]->cd(isctr+1);
        h[idtyp][ihst][isctr]->Draw("colz");
        if (ihst==1)eid_tool[idtyp]->DrawSFcuts(isctr+1);
      }
    }
  }

  //! Same canvas for each for 1D
  for (int ihst=0;ihst<NHST;ihst++){
    if (ihst<=1||ihst==6) continue;
    TString cname=TString::Format("c_%s",hst_name[ihst].Data());
    c[EXP][ihst]=new TCanvas(cname,cname);
    c[EXP][ihst]->Divide(3,2);
    for (int isctr=0;isctr<NSCTR;isctr++){
      c[EXP][ihst]->cd(isctr+1);
      h[EXP][ihst][isctr]->SetLineColor(kBlue);
      h[SIM][ihst][isctr]->SetLineColor(kRed);
      h[EXP][ihst][isctr]->DrawNormalized("",1000);
      h[SIM][ihst][isctr]->DrawNormalized("sames",1000);
    }
  }
}
