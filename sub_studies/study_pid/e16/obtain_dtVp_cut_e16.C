{
  TFile* f[2];// 0=exp,1=sim
  TString dtyp_name[]={"exp","sim"};
  TString part_name[3]={"p","pip","pim"};
  TString part_name_latex[3]={"p","#pi^{+}","#pi^{-}"};
  f[0]=TFile::Open("$D2PIDIR_EXP_E16/data_pid_dtvp_020117/dpid_dtvp_2pi.root");
  f[1]=TFile::Open("$D2PIDIR_SIM_E16/data_pid_dtvp_020117/dpid_dtvp_2pi.root");
 
  TH2F* h2_dtVp[2][3];

  //! Create directory where output (in .root file and as .pngs and .pdf) is stored
  //! First create date in mmddyy format
  TDatime dt=TDatime();
  TString yyyy=TString::Format("%d",dt.GetYear());
  TString yy(yyyy(2,2));
  TString dd=TString::Format("%02d",dt.GetDay());
  TString mm=TString::Format("%02d",dt.GetMonth());
  TString mmddyy=TString::Format("%s%s%s",mm.Data(),dd.Data(),yy.Data());
  TString OUTDIR=TString::Format("%s/dtVp_cuts_%s",gSystem->ExpandPathName("$STUDY_PID_DTVP_E16_DATADIR"),mmddyy.Data());
  gSystem->mkdir(OUTDIR); //! Will create dir only if it does not exist
  //! 'touch' OUTDIR to reflect that it has been worked upon
  TString cmd=TString::Format(".!touch %s",OUTDIR.Data());
  gROOT->ProcessLine(cmd);
  //! .root output
  TFile* fout=new TFile(TString::Format("%s/fout.root",OUTDIR.Data()),"RECREATE");
  
  TCanvas *c=new TCanvas("c","c");//default canvas
  for (int i=0;i<2;i++){
    TTree* t=(TTree*)f[i]->Get("pid/tree");;
    c->cd();
    //! Note upper momenta limits for p,pip,pim determined from $ANANOTE/figures/Acceptance/kin_comp_ER_SR_ST/hQ2/c_p_top1.pdf
    //dtVp under particle assumption that all q==1 = proton
    t->Draw("dt_p:p>>h_dtVp_pos_p(80,0,4.5,250,-5,5)","q==+1&&dc>0&&sc>0","colz");
    //dtVp under particle assumption that all q==1 = pip
    t->Draw("dt_pip:p>>h_dtVp_pos_pip(50,0,2.5,250,-5,5)","q==+1&&dc>0&&sc>0","colz");
    //dtVp under particle assumption that all q==-1 = pim
    t->Draw("dt_pim:p>>h_dtVp_neg_pim(50,0,2.5,250,-5,5)","q==-1&&dc>0&&sc>0","colz");
   
    gStyle->SetOptStat("ne");//mri"); 
    gStyle->SetOptFit(1111);
    //! Make dtVp histograms
    h2_dtVp[i][0]=(TH2F*)gDirectory->Get("h_dtVp_pos_p");
    h2_dtVp[i][0]->SetName(TString::Format("h_dtVp_pos_p_%s",dtyp_name[i].Data()));
    h2_dtVp[i][0]->GetYaxis()->SetRangeUser(-5,5);
    h2_dtVp[i][1]=(TH2F*)gDirectory->Get("h_dtVp_pos_pip");
    h2_dtVp[i][1]->SetName(TString::Format("h_dtVp_pos_pip_%s",dtyp_name[i].Data()));
    h2_dtVp[i][1]->GetYaxis()->SetRangeUser(-5,5);
    h2_dtVp[i][2]=(TH2F*)gDirectory->Get("h_dtVp_neg_pim");
    h2_dtVp[i][2]->SetName(TString::Format("h_dtVp_neg_pim_%s",dtyp_name[i].Data()));
    h2_dtVp[i][2]->GetYaxis()->SetRangeUser(-5,5);

    //! Now fit dt's momentum projections for each particle and get cut parameters
    //! + dt fits are performed in bins where there are high statistics
    //! + As stats get low at high momentum limits, "Method to obtain fit cure to dt vs. p for high momenta" is used:

    //! Method to obtain fit cure to dt vs. p for high momenta:
    //! =======================================================
    //! + Since data beyond 'pmax' (see below) is either not available or sparse, fits to dt cannot be made and
    //! therefore, the fit curve to dt vs. p points cannot be determined in the high momentum bins.
    //!
    //! + In order to get fit curves in p region where no data exists, I am going to use dt-mu/sg from the last bin
    //! in p where there was data to obtain dt-mu/sg
    //!
    //! + max-momenta(p,pip,pim)    ={4.50,2.50,2.50} ($ANANOTE/figures/Acceptance/kin_comp_ER_SR_ST/hQ2/c_p_top1.pdf)
    //! + max-fit-momenta(p,pip,pim)={3.25,2.00,2.75} (obtained by visual inspection)
    //! (Ignore any values quoted for pim because after a while, since its not used in the analysis, I stopped fine tuning its parameters)
    

    //! First define momentum bins where dt fits can be made
    //! + Bins defined are optimized visually by making sure that there are enough momentum bins with good stats
    //! + Note, because of dt fits cannot be made upto PMAX[3]={4.50,2.50,2.50} and these limits need to be truncated to {3.25,2.00,2.75}
    //! + This truncation is OK because bins ignored have insignificantly low statistics (see $ANANOTE/figures/Acceptance/kin_comp_ER_SR_ST/hQ2/c_p_top1.pdf)
    //! + For these high momentum bins "Method to obtain fit cure to dt vs. p for high momenta" is used
    int nbins[3]={10,4,10}; //!0=p,1=pip,2=pim
    float pmin[3][10]=
    {
      {0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,3.00},//0.00},
      {0.00,0.50,1.00,1.50,0.00,0.00,0.00,0.00,0.00,0.00},//2.50,1.50,1.75,2.00,2.25},//2.50},//{0.00,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50},
      {0.00,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.50}//0.00}
    };
    float pmax[3][10]=
    {
      {0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.25},//,0.00},
      {0.50,1.00,1.50,2.00,0.50,0.00,0.00,0.00,0.00,0.00},//3.00,1.75,2.00,2.25,2.50},//,3.00},//{0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,3.00},
      {0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.75}//,0.00}
    };
    //! Now define momentum bins where "Method to obtain fit cure to dt vs. p for high momenta" will be used
    //! + Currently, this is used only for the:
    //!   + proton: There are 5 high momentum bins with low stats in which the dt vs p fit needs to constrained
    //!   + pip: There is 1 high momentum bin in which dt vs p fit needs to constrained
    //!     (Note, for pip, this is only done for experiment because in the simulation the fit directly extended i.e. without being
    //!      constrained by a dt-par in the high momentum, actually works better)     
    //!   + pim: Ignored becaue it is not being used in the analysis 
    int nbins_high_p[3]={5,1,2};// chosen so as to use minimum number of bins whose binning can be expressed within 2 decimal places. 
    float pmin_high_p[3][5]=//! low bin for each particle should be pmax for each particle!
    {
      {3.25,3.50,3.75,4.00,4.25},//{3.25,3.60,3.95,4.30,4.65},
      {2.00,0.00,0.00,0.00,0.00},
      {2.50,3.75.0.00,0.00,0.00} 
    }; 
    float pmax_high_p[3][5]=
    {
      {3.50,3.75,4.00,4.25,4.50},//{3.60,3.95,4.30,4.65,5.00},
      {2.50,0.00,0.00,0.00,0.00},
      {3.75,5.00,0.00,0.00,0.00}
    };
    //! For each particle, get dt-mu/sg for each pbin and 
    //! put values into histograms(hcut_h,hcut_l) that will later be fit
    for (int j=0;j<3;j++){
      //! create outdir pngs
      TString outdir=TString::Format("%s/%s/%s",OUTDIR.Data(),dtyp_name[i].Data(),part_name[j].Data());
      bool rcrsv=kTRUE;
      gSystem->mkdir(outdir,rcrsv);
      //! create outdir in root file
      if (fout->GetDirectory(dtyp_name[i].Data())!=NULL){//! i.e. dtyp dir already exists; cd() there and create prtcl dir
         fout->cd(dtyp_name[i].Data());
         gDirectory->mkdir(part_name[j].Data())->cd();
      }else{//! directly create dtyp/prtcl dir
         TString outdir_root=TString::Format("%s/%s",dtyp_name[i].Data(),part_name[j].Data());
         fout->mkdir(outdir_root);
         fout->cd(outdir_root);
      }
      /*TString outdir_root=TString::Format("%s/%s",dtyp_name[i].Data(),part_name[j].Data());
      if (fout->mkdir(outdir_root)==NULL){//! i.e. dtyp_name dir already exists
        fout->cd(dtyp_name[i].Data())->mkdir(part_name[j].Data());
        TString outdir_root=TString::Format("%s/%s",dtyp_name[i].Data(),part_name[j].Data());;
      }
      fout->cd(outdir_root);     
      return;*/

      //clear current cut data for jth particle
      TString cmd=TString::Format(".!rm %s/*",outdir.Data());
      gROOT->ProcessLine(cmd);

      TH1F* hmean=new TH1F("hmean","hmean",100,0,5);
      TH1F* hcut_h=new TH1F("hcut_h","hcut_h",100,0,5);
      TH1F* hcut_l=new TH1F("hcut_l","hcut_l",100,0,5);
      hmean->SetMarkerStyle(kStar);
      hmean->SetMarkerColor(kRed);
      hcut_h->SetMarkerStyle(kStar);
      hcut_h->SetMarkerColor(kRed);
      hcut_l->SetMarkerStyle(kStar);
      hcut_l->SetMarkerColor(kRed);
      //! fill hcut_l/h with dt obtained from fits in momentum bins
      for (int ibin=0;ibin<nbins[j];ibin++){
        int bmin=h2_dtVp[i][j]->GetXaxis()->FindBin(pmin[j][ibin]);
        int bmax=h2_dtVp[i][j]->GetXaxis()->FindBin(pmax[j][ibin]);
        TH1F* h=(TH1F*)h2_dtVp[i][j]->ProjectionY(TString::Format("_py_%d",ibin+1),bmin,bmax);
        //h->SetTitle(TString::Format("%s_%s_p_%.2f-%.2f",dtyp_name[i].Data(),part_name[j].Data(),pmin[j][ibin],pmax[j][ibin]));
        h->SetTitle(TString::Format("p=[%.2f GeV, %.2f GeV)",pmin[j][ibin],pmax[j][ibin]));
        //TCanvas* c=new TCanvas("c","c");
        TString cname=TString::Format("c_%02d_%s_%s",ibin+1,dtyp_name[i].Data(),part_name[j].Data());
        TCanvas* c=new TCanvas(cname,cname);
        h->Draw();
        //! Axes title
        h->SetXTitle("#Deltat (ns)");
        c->SetLeftMargin(0.20);
        h->GetYaxis().SetTitleOffset(1.5);
        h->SetYTitle("N_{entries}");
        //! Fit
        if (part_name[j]=="pip"){
          //! 100715
          //h->Fit("gaus","","",-0.20,0.30);
          //! 100715.II
          /*if (dtyp_name[i]=="exp"){
            h->Fit("gaus","","",-0.60,0.40);
          }else if (dtyp_name[i]=="sim"){
            h->Fit("gaus","","",-0.40,0.40);
          }*/
          //! [02-13-16] E16
          if (dtyp_name[i]=="exp"){
            if      ((ibin+1)<4)  {h->Fit("gaus","","",-0.50,0.40);}
            else if ((ibin+1)==4) {h->Fit("gaus","","",-0.40,0.40);}
            else if ((ibin+1)==5) {h->Fit("gaus","","",-0.20,0.40);}
          }else if (dtyp_name[i]=="sim"){
            if      ((ibin+1)<5)  h->Fit("gaus","","",-0.50,0.50);
            else if ((ibin+1)==5) h->Fit("gaus","","",-0.20,0.50);
          }
        }else if(part_name[j]=="pim"){
          //! 100715
          //h->Fit("gaus","","",-0.20,0.20); 
          //! 100715.II
          //h->Fit("gaus","","",-0.40,0.40);
          //! [02-13-16] E16
          if (dtyp_name[i]=="exp"){
            h->Fit("gaus","","",-0.50,0.40);
          }else if (dtyp_name[i]=="sim"){
            h->Fit("gaus","","",-0.50,0.50);
          }
        } else if(part_name[j]=="p"){
          //! [02-13-16] E16
          if (dtyp_name[i]=="exp"){
            h->Fit("gaus","","",-0.50,0.40);
          }else if (dtyp_name[i]=="sim"){
            h->Fit("gaus","","",-0.50,0.50);
          }
        }/*else{
          //! 100715
          //h->Fit("gaus","","",-0.20,0.20);
          //! 100715.II
          h->Fit("gaus","","",-0.40,0.40);
        }*/
        //c->SaveAs(TString::Format("%s/c%02d.jpg",outdir.Data(),ibin+1));
        //c->Write(TString::Format("c_%s_%s_pbin%02d",dtyp_name[i].Data(),part_name[j].Data(),ibin+1));
        c->SaveAs(TString::Format("%s/%s.jpg",outdir.Data(),cname.Data()));
	      c->SaveAs(TString::Format("%s/%s.pdf",outdir.Data(),cname.Data()));
        c->Write(cname.Data());
        c->Close();
        //! Get fit parameters
        TF1 *fitf=h->GetFunction("gaus");
        Float_t constant=fitf->GetParameter(0);
        float mean=fitf->GetParameter(1);
        float sigma=fitf->GetParameter(2);
        float mean_plus_3sigma=mean+3*sigma;
        float mean_minus_3sigma=mean-3*sigma;
        //Fill hmean,hsigma3
        float pavg=(pmin[j][ibin]+pmax[j][ibin])/2;
        int bavg=hmean->GetXaxis()->FindBin(pavg);
        hmean->SetBinContent(bavg,mean);
        hcut_h->SetBinContent(bavg,mean_plus_3sigma);
        hcut_l->SetBinContent(bavg,mean_minus_3sigma);
      }//end p-bin loop
      //! Now fill hcut_l/h for high momentum bins where there is no data available to fit 
      //! using the method noted earlier
      //! Get bin content of the the last p bin where dt-mu/sg was obtained using a fit
      if(part_name[j]=="p" || (part_name[j]=="pip" && dtyp_name[i]=="exp")){
        float pavg_last_fitted_bin=(pmin[j][nbins[j]-1]+pmax[j][nbins[j]-1])/2;
        float bavg_last_fitted_bin=hmean->GetXaxis()->FindBin(pavg_last_fitted_bin);
        dt_mean_last_fitted_bin=hmean->GetBinContent(bavg_last_fitted_bin);
        dt_h_last_fitted_bin=hcut_h->GetBinContent(bavg_last_fitted_bin);
        dt_l_last_fitted_bin=hcut_l->GetBinContent(bavg_last_fitted_bin);
        //! Now fill hcut_h/l
        for (int ibin=0;ibin<nbins_high_p[j];ibin++){
          float pavg=(pmin_high_p[j][ibin]+pmax_high_p[j][ibin])/2;
          int bavg=hmean->GetXaxis()->FindBin(pavg);
          hmean->SetBinContent(bavg,dt_mean_last_fitted_bin);
          hcut_h->SetBinContent(bavg,dt_h_last_fitted_bin);
          hcut_l->SetBinContent(bavg,dt_l_last_fitted_bin);
        }// end high_p-bin loop
      }
      
      //Fit hcut_l and hcut_h
      TF1 *fl = new TF1("fl","pol3",-3,3);
      TF1 *fh = new TF1("fh","pol3",-3,3);
      hcut_l->Fit("fl","","",pmin[j][0],pmax_high_p[j][nbins_high_p[j]-1]);
      hcut_h->Fit("fh","","",pmin[j][0],pmax_high_p[j][nbins_high_p[j]-1]);
      //c->Close();
      FILE* fcut;
      fcut=fopen(TString::Format("%s/cutpars.txt",outdir.Data()).Data(), "w");
      //fprintf(fcut,"cut_l(pol3:p0:p1:p2:p3) %.5f %.5f %.5f %.5f\n",fl->GetParameter(0),fl->GetParameter(1),fl->GetParameter(2),fl->GetParameter(3)); 
      //fprintf(fcut,"cut_h(pol3:p0:p1:p2:p3) %.5f %.5f %.5f %.5f\n",fh->GetParameter(0),fh->GetParameter(1),fh->GetParameter(2),fh->GetParameter(3));
      fprintf(fcut,"cut_l(pol3:p0:p1:p2:p3) %f %f %f %f\n",fl->GetParameter(0),fl->GetParameter(1),fl->GetParameter(2),fl->GetParameter(3));
      fprintf(fcut,"cut_h(pol3:p0:p1:p2:p3) %f %f %f %f\n",fh->GetParameter(0),fh->GetParameter(1),fh->GetParameter(2),fh->GetParameter(3));
      fclose(fcut);
      //! Draw fit
      //TCanvas* ccut=new TCanvas("cut","cut");
      cname=TString::Format("c_cut_%s_%s",dtyp_name[i].Data(),part_name[j].Data());
      TCanvas* ccut=new TCanvas(cname,cname);
      h2_dtVp[i][j]->Draw("colz");
      //! Axes title
      h2_dtVp[i][j]->SetXTitle("p [GeV]");
      ccut->SetLeftMargin(0.20);
      h2_dtVp[i][j]->GetYaxis().SetTitleOffset(1.5);
      h2_dtVp[i][j]->SetYTitle("#Deltat [ns]");
      //! Adjust title from TTree::Draw()
      h2_dtVp[i][j]->SetTitle(TString::Format("#Deltat versus momentum under %s mass assumption",part_name_latex[j].Data()));
      //! cut functions
      hmean->Draw("P same");
      hcut_h->Draw("P same");
      hcut_l->Draw("P same");
      //ccut->SaveAs(TString::Format("%s/ccut.jpg",outdir.Data()));
      //ccut->Write(TString::Format("c_cut_%s_%s",dtyp_name[i].Data(),part_name[j].Data()));
      ccut->SaveAs(TString::Format("%s/%s.jpg",outdir.Data(),cname.Data()));
      ccut->SaveAs(TString::Format("%s/%s.pdf",outdir.Data(),cname.Data()));
      ccut->Write(cname.Data());
      ccut->Close();
      //gPad->Update();
    }//end j loop
  }//end i loop
}
