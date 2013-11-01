{
  Int_t nSim = 6;
  Int_t nQ2Wbin = 8;
  FILE* f[nSim];
  FILE* f[0] = fopen(gSystem->ExpandPathName("$(E1F_SIM2PI_ANADIR)/Q2W__1.4-1.5__1.6-1.8__100M/simstats.txt"), "r");
  FILE* f[1] = fopen(gSystem->ExpandPathName("$(E1F_SIM2PI_ANADIR)/Q2W__1.4-1.5__1.6-1.8__100M__200M/simstats.txt"), "r");
  FILE* f[2] = fopen(gSystem->ExpandPathName("$(E1F_SIM2PI_ANADIR)/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii/simstats.txt"), "r");
  FILE* f[3] = fopen(gSystem->ExpandPathName("$(E1F_SIM2PI_ANADIR)/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii_400M/simstats.txt"), "r");
  FILE* f[4] = fopen(gSystem->ExpandPathName("$(E1F_SIM2PI_ANADIR)/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii_400M_400Mii/simstats.txt"), "r");
  FILE* f[5] = fopen(gSystem->ExpandPathName("$(E1F_SIM2PI_ANADIR)/Q2W__1.4-1.5__1.6-1.8__100M__200M__200Mii_400M_400Mii_400Miii/simstats.txt"), "r"); 

  TString Q2Wdirname[nQ2Wbin]; 
  Int_t n0simRECObins[nSim][nQ2Wbin];
  Int_t nFsimRECObins[nSim][nQ2Wbin];
  Int_t nRECO[nSim][nQ2Wbin];
  Int_t nMC[nSim][nQ2Wbin];

 
  for (Int_t iSim = 0; iSim < nSim; iSim++){ 
    cout << "Sim = " << endl;
    if (f[iSim]){
      char str[500];
      char tmpQ2Wdirname[200];
      Int_t iQ2Wbin = 0;
      while (fgets(str, sizeof(str), f[iSim])){
        sscanf(str, "%s %d %d %d %d", 
               tmpQ2Wdirname, &n0simRECObins[iSim][iQ2Wbin], &nFsimRECObins[iSim][iQ2Wbin], &nRECO[iSim][iQ2Wbin], &nMC[iSim][iQ2Wbin]);
       Q2Wdirname[iQ2Wbin] = tmpQ2Wdirname;
       printf("%s	%d	%d	%d	%d\n",
             Q2Wdirname[iQ2Wbin].Data(), n0simRECObins[iSim][iQ2Wbin], nFsimRECObins[iSim][iQ2Wbin], nRECO[iSim][iQ2Wbin], nMC[iSim][iQ2Wbin]);
       //cout << Q2Wdirname[iQ2Wbin] << "	" << n0simRECObins[iSim][iQ2Wbin] << "	" << nFsimRECObins[iSim][iQ2Wbin] << "	" << nRECO[iSim][iQ2Wbin] << "	" << nMC[iSim][iQ2Wbin] << endl;
        iQ2Wbin += 1;
      }
    }
  }

  TCanvas* c_0binStatsStudy = new TCanvas("0binStatsStudy", "0binStatsStudy", 1500, 400);
  c_0binStatsStudy->Divide(4,2);
  TGraph* g_n0simRECObinsVnMC[nQ2Wbin];

  TCanvas* c_FbinStatsStudy = new TCanvas("FbinStatsStudy", "FbinStatsStudy", 1500, 400);
  c_FbinStatsStudy->Divide(4,2);
  TGraph* g_nFsimRECObinsVnMC[nQ2Wbin];

  TCanvas* c_nRECOVnMC = new TCanvas("nRECOVnMC", "nRECO_V_nMC", 1500, 400);
  c_nRECOVnMC->Divide(4,2);
  TGraph* g_nRECObinsVnMC[nQ2Wbin]; 

  for (Int_t iQ2Wbin =0; iQ2Wbin < nQ2Wbin; iQ2Wbin++){
    Int_t v_nMC[nSim];
    Int_t v_n0simRECObins[nSim];
    Int_t v_nFsimRECObins[nSim];
    Int_t v_nRECO[nSim];
    for (Int_t iSim = 0; iSim < nSim; iSim++){
      v_nMC[iSim] = nMC[iSim][iQ2Wbin];
      v_n0simRECObins[iSim] = n0simRECObins[iSim][iQ2Wbin];
      v_nFsimRECObins[iSim] = nFsimRECObins[iSim][iQ2Wbin];
      v_nRECO[iSim] = nRECO[iSim][iQ2Wbin];
    }
    g_n0simRECObinsVnMC[iQ2Wbin] = new TGraph(nSim, v_nMC, v_n0simRECObins);
    g_n0simRECObinsVnMC[iQ2Wbin]->SetLineWidth(2);
    g_n0simRECObinsVnMC[iQ2Wbin]->SetLineColor(2);
    g_n0simRECObinsVnMC[iQ2Wbin]->SetNameTitle(Q2Wdirname[iQ2Wbin].Data(), Q2Wdirname[iQ2Wbin].Data());
    g_nFsimRECObinsVnMC[iQ2Wbin] = new TGraph(nSim, v_nMC, v_nFsimRECObins);
    g_nFsimRECObinsVnMC[iQ2Wbin]->SetLineWidth(2);
    g_nFsimRECObinsVnMC[iQ2Wbin]->SetLineColor(4);
    g_nFsimRECObinsVnMC[iQ2Wbin]->SetNameTitle(Q2Wdirname[iQ2Wbin].Data(), Q2Wdirname[iQ2Wbin].Data());
    g_nRECObinsVnMC[iQ2Wbin] = new TGraph(nSim, v_nMC, v_nRECO);
    g_nRECObinsVnMC[iQ2Wbin]->SetNameTitle(Q2Wdirname[iQ2Wbin].Data(), Q2Wdirname[iQ2Wbin].Data());
   
    c_0binStatsStudy->cd(iQ2Wbin+1);
    g_n0simRECObinsVnMC[iQ2Wbin]->Draw("AL*");
    g_n0simRECObinsVnMC[iQ2Wbin]->GetXaxis()->SetTitle("nThrown");
    g_n0simRECObinsVnMC[iQ2Wbin]->GetYaxis()->SetTitle("n0simRECObins");
   
    c_FbinStatsStudy->cd(iQ2Wbin+1); 
    g_nFsimRECObinsVnMC[iQ2Wbin]->Draw("AL*");
    g_nFsimRECObinsVnMC[iQ2Wbin]->GetXaxis()->SetTitle("nThrown");
    g_nFsimRECObinsVnMC[iQ2Wbin]->GetYaxis()->SetTitle("nFsimRECObins");

    c_nRECOVnMC->cd(iQ2Wbin+1);
    g_nRECObinsVnMC[iQ2Wbin]->Draw("AL*");
    g_nRECObinsVnMC[iQ2Wbin]->GetXaxis()->SetTitle("nThrown");
    g_nRECObinsVnMC[iQ2Wbin]->GetYaxis()->SetTitle("nReco");
    
  }

}
