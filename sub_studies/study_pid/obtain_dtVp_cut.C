{
  TFile* f[2];// 0=exp,1=sim
  TString dtyp_name[]={"exp","sim"};
  TString part_name[3]={"p","pip","pim"};
  f[0]=TFile::Open("$D2PIDIR_EXP/mon_pid_new/dpid.root");
  f[1]=TFile::Open("$D2PIDIR_SIM/mon_pid_new/dpid.root");
 
  TCanvas* c_cut[2];
  TH2F* h2_dtVp[2][3];

  TCanvas *c=new TCanvas("c","c");//default canvas
  for (int i=0;i<2;i++){
    TTree* t=(TTree*)f[i]->Get("pid/tree");;
    c->cd();
    //dtVp under particle assumption that all q==1 = proton
    t->Draw("dt_p:p>>h_dtVp_pos_p(100,0,5,250,-5,5)","q==+1&&dc>0&&sc>0","colz");
    //dtVp under particle assumption that all q==1 = pip
    t->Draw("dt_pip:p>>h_dtVp_pos_pip(100,0,5,250,-5,5)","q==+1&&dc>0&&sc>0","colz");
    //dtVp under particle assumption that all q==-1 = pim
    t->Draw("dt_pim:p>>h_dtVp_neg_pim(100,0,5,250,-5,5)","q==-1&&dc>0&&sc>0","colz");
   
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

    //! Now fit momentum projections for each particle
    //! and get cut parameters
    c_cut[i]=new TCanvas(TString::Format("cut_dtVp_%s",dtyp_name[i].Data()),TString::Format("cut_dtVp_%s",dtyp_name[i].Data()));
    c_cut[i]->Divide(2,2);
    int nbins[3]={10,5,10}; //!0=p,1=pip,2=pim
    float pmin[3][11]=
    {
      {0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,3.00,0.00},
      {0.00,0.50,1.00,1.50,2.00,2.50,1.50,1.75,2.00,2.25,2.50},//{0.00,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50},
      {0.00,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.50,0.00}
    };
    float pmax[3][11]=
    {
      {0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.25,0.00},
      {0.50,1.00,1.50,2.00,2.50,3.00,1.75,2.00,2.25,2.50,3.00},//{0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,3.00},
      {0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.75,0.00}
    };
    for (int j=0;j<3;j++){
      //clear current cut data for jth particle
      TString cmd=TString::Format(".!rm dtVp_cuts/%s/%s/*",dtyp_name[i].Data(),part_name[j].Data());
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
      for (int ibin=0;ibin<nbins[j];ibin++){
        int bmin=h2_dtVp[i][j]->GetXaxis()->FindBin(pmin[j][ibin]);
        int bmax=h2_dtVp[i][j]->GetXaxis()->FindBin(pmax[j][ibin]);
        TH1F* h=(TH1F*)h2_dtVp[i][j]->ProjectionY(TString::Format("_py_%d",ibin+1),bmin,bmax);
        h->SetTitle(TString::Format("p_%.2f-%.2f",pmin[j][ibin],pmax[j][ibin]));
        TCanvas* c=new TCanvas("c","c");
        h->Draw();
        if (part_name[j]=="pip"){
          h->Fit("gaus","","",-0.20,0.30);//-0.40,0.20
        }else if(part_name[j]=="pim"){
          h->Fit("gaus","","",-0.20,0.20); 
        }else{
          h->Fit("gaus","","",-0.20,0.20);
        }
        c->SaveAs(TString::Format("dtVp_cuts/%s/%s/c%02d.jpg",dtyp_name[i].Data(),part_name[j].Data(),ibin+1));
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
      //Fit hcut_l and hcut_h
      TF1 *fl = new TF1("fl","pol3",-3,3);
      TF1 *fh = new TF1("fh","pol3",-3,3);
      hcut_l->Fit("fl","","",pmin[j][0],pmax[j][nbins[j]-1]); 
      hcut_h->Fit("fh","","",pmin[j][0],pmax[j][nbins[j]-1]);
      //c->Close();
      FILE* fcut;
      fcut=fopen(TString::Format("dtVp_cuts/%s/%s/cutpars.txt",dtyp_name[i].Data(),part_name[j].Data()).Data(), "w");
      fprintf(fcut,"cut_l(pol3:p0:p1:p2:p3) %.2f %.2f %.2f %.2f\n",fl->GetParameter(0),fl->GetParameter(1),fl->GetParameter(2),fl->GetParameter(3)); 
      fprintf(fcut,"cut_h(pol3:p0:p1:p2:p3) %.2f %.2f %.2f %.2f\n",fh->GetParameter(0),fh->GetParameter(1),fh->GetParameter(2),fh->GetParameter(3));
      fclose(fcut);
      //! Draw fit
      TCanvas* ccut=new TCanvas("cut","cut");
      h2_dtVp[i][j]->Draw("colz");
      hmean->Draw("P same");
      hcut_h->Draw("P same");
      hcut_l->Draw("P same");
      ccut->SaveAs(TString::Format("dtVp_cuts/%s/%s/ccut.jpg",dtyp_name[i].Data(),part_name[j].Data()));	
      //ccut->Close();
      //! For interactive display
      c_cut[i]->cd(j+1);
      h2_dtVp[i][j]->Draw("colz");
      hmean->Draw("P same");
      hcut_h->Draw("P same");
      hcut_l->Draw("P same");
    }//end j loop
  }//end i loop
}
