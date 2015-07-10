{
  int NDTYP=2;
  enum {EXP,SIM};
  TFile* f[NDTYP];// 0=Exp,1=Sim
  TString dtyp_name[]={"Exp","Sim"};
  f[EXP]=TFile::Open("$D2PIDIR_EXP/data_eid_070915/deid.root");
  f[SIM]=TFile::Open("$D2PIDIR_SIM/data_eid_070915/deid.root");

  int NSCTR=6;

  /*TObjArray** hists;
  for (int isctr=0;isctr<NSCTR;isctr++){
    hists[isctr] = new TObjArray(9);
    hists[isctr]->Add(new TH2F("heoutVein",      TString::Format("#E_{out} vs. #E_{in}(sector%d)", isctr),     150,0,1,150,0,1));
    hists[isctr]->Add(new TH2F("hsfoutVsfin",    TString::Format("SF_{out} vs. SF_{in}(sector%d)", isctr),     100,0,0.5,100,0,0.5));
    hists[isctr]->Add(new TH2F("hsftotVp",       TString::Format("SF_{tot} vs. p(sector%d)", isctr),           160,0,5,100,0,0.5));
    hists[isctr]->Add(new TH2F("hetotVp",        TString::Format("E_{tot} vs. p(sector%d)", isctr),            160,0,5,100,0,2));
    hists[isctr]->Add(new TH1F("hU",             TString::Format("U(sector%d)", isctr),                        100,0,450));
    hists[isctr]->Add(new TH1F("hV",             TString::Format("V(sector%d)", isctr),                        100,0,450));
    hists[isctr]->Add(new TH1F("hW",             TString::Format("V(sector%d)", isctr),                        100,0,450));
    hists[isctr]->Add(new TH1F("hnphe",          TString::Format("nphe(sector%d)", isctr),                     100,0,300));
    hists[isctr]->Add(new TH2F("hCCthetaVCCseg" ,TString::Format("CC_{#theta} vs. CC_{seg}(sector%d)", isctr), 20,0,20,100,0,50)); 
  }*/

  int NHST=9;
  TString hst_name[]={"heoutVein","hsfoutVsfin","hsftotVp","hetotVp","hU","hV","hW","hnphe","hCCthetaVCCseg"};
  TString hst_title[]={"#E_{out} vs. #E_{in}","hsfoutVsfin","hsftotVp","hetotVp","hU","hV","hW","hnphe","hCCthetaVCCseg"};
  TString hst_drw_cmd[]={TString::Format("ec_eo:ec_ei>>%s(150,0,1,150,0,1)",    hst_name[0].Data()),
                         TString::Format("ec_eo/p:ec_ei/p>>%s(150,0,0.5,150,0,0.5)",hst_name[1].Data())};
   
  TH2F* h[NHST][NDTYP][NSCTR];
  TCanvas* c[NHST][NDTYP];

  TCanvas *c1=new TCanvas("c1","c1");//default canvas
  for (int idtyp=0;idtyp<NDTYP;idtyp++){
    //! Get hists from Tree
    TTree* t=(TTree*)f[idtyp]->Get("eid/monitor/t");
    for (int ihst=0;ihst<NHST;ihst++){
      if (ihst>1)continue;
      for (int isctr=0;isctr<NSCTR;isctr++){
         c1->cd();
        //TString drw_cmd=TString::Format("ec_eo:ec_ei>>heoutVein(150,0,1,150,0,1)");//TString::Format("ec_eo:ec_ei>>heoutVein_s%d(150,0,1,150,0,1)",isctr+1);
        TString sctr_cut=TString::Format("sector==%d",isctr+1);
        //t->Draw(drw_cmd,sctr_cut,"colz");
        t->Draw(hst_drw_cmd[ihst],sctr_cut,"colz");
   
        gStyle->SetOptStat("ne");//mri"); 
        //h[idtyp][isctr]=(TH2F*)gDirectory->Get(TString::Format("heoutVein_s%d",isctr+1));
        TH2F* htmp=(TH2F*)gDirectory->Get(hst_name[ihst]);
        h[ihst][idtyp][isctr]=(TH2F*)htmp->Clone();
        h[ihst][idtyp][isctr]->SetName(TString::Format("%s_%s_s%d",       hst_name[ihst].Data(),dtyp_name[idtyp].Data(),isctr+1));
        h[ihst][idtyp][isctr]->SetTitle(TString::Format("%s(%s,sector%d)",hst_name[ihst].Data(),dtyp_name[idtyp].Data(),isctr+1));
      }//sctr loop
      //! Now plot hists
      TString cname=TString::Format("c_%s_%s",hst_name[ihst].Data(),dtyp_name[idtyp].Data());
      c[ihst][idtyp]=new TCanvas(cname,cname);
      c[ihst][idtyp]->Divide(3,2);
      for (int isctr=0;isctr<NSCTR;isctr++){
        c[ihst][idtyp]->cd(isctr+1);
        h[ihst][idtyp][isctr]->Draw("colz");
      }
    }//nhst loop
    //! Now plot hists
  }//dtyp loop
}
