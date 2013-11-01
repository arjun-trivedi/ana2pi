/*********************************************************
10-15-13:
This program verfies the following:
1. That the 2 techniques to obtain R_lt(and error) using Method 2 are the same.
  Technique 1. 
    i.   hY5D-->hphiproj_thetaslice<i>
    ii.  hphiproj_theta_slice<i> = hphiproj_theta_slice<i>*hsinphi
    iii. R_theta<i> = Int(hphiproj_theta_slice<i>*hsinphi)
  Technique 2.  
    i.   hY5D = hY5D*sin(phi)
    ii.  hY5D-->htheta_proj
    iii. R_theta<i> = htheta_proj.GetBinContent(i)
**********************************************************/
{
 //Following only for UNPOL data
 //1. Get hY5D_FULL
 THnSparse* hY5Dtmp = (THnSparse*)gDirectory->Get("hY5D/Varset1/hY5D_FULL");
 THnSparse* hY5D = (THnSparse*)hY5Dtmp->Clone(hY5Dtmp->GetName());
 //2. Get hPhi projections for THETA bins made from hY5D_FULL
 TH1D* hphiproj_theta[10];
 for (int i=0;i<10;i++){
   //cout << TString::Format("hPhi/Varset1/theta/hphi_proj_%02d",i+1) << endl;
   TH1D* htmp = (TH1D*)gDirectory->Get(TString::Format("hPhi/Varset1/theta/hphi_proj_%02d",i+1));
   hphiproj_theta[i] = (TH1D*)htmp->Clone(htmp->GetName());
 }
 /*TCanvas *c1 = new TCanvas();
 hphiproj_theta[0]->Draw();*/

 //!3. Verify if (hY5D_FULL-->)h1D_THETA[i] =  int(hphiproj_theta[i]) 
 //!3.1 First make h1D_THETA;
 TH1D* h1D_THETA = (TH1D*)hY5D->Projection(THETA,"E");
 //!3.2 Now Verify
 for (int i=0;i<10;i++){
   double itgerr;
   double itg = hphiproj_theta[i]->IntegralAndError(1,hphiproj_theta[i]->GetNbinsX(),itgerr);
   
   printf("[%.2f:%.2f][[%.2f:%.2f]\n",h1D_THETA->GetBinContent(i+1),h1D_THETA->GetBinError(i+1),itg,itgerr);
 }

 //!4. Verify if (hY5D_FULL*sin(phi)-->)h1D_THETA[i]/norm = int(hphiproj_theta[i]*hsinphi)/norm
 Double_t norm = 1000;
 //!4.1 First CREATE (hY5D_FULL*sin(phi)-->)h1D_THETA[i]
 int nbins = hY5D->GetNbins();
 TAxis* axphi = hY5D->GetAxis(PHI);
 for (int ibin=0;ibin<nbins;ibin++){
   int bincoord[5] = {0};
   float binc = hY5D->GetBinContent(ibin, bincoord);
   float binerr = hY5D->GetBinError(ibin);
   double phi = axphi->GetBinLowEdge(bincoord[PHI]) * TMath::DegToRad();
   hY5D->SetBinContent(ibin, TMath::Sin(phi)*binc);
   hY5D->SetBinError(ibin, TMath::Sin(phi)*binerr);
 }
 h1D_THETA = (TH1D*)hY5D->Projection(THETA,"E");
 h1D_THETA->Scale(1/TMath::Pi());
 h1D_THETA->Scale(1/norm); 
 //!4.2 Now CREATE hphiproj_theta[i]*hsinphi 
 //!4.2.1 CREATE hsinphi
 int nphibins  = hY5D->GetAxis(PHI)->GetNbins();
 double phimin = hY5D->GetAxis(PHI)->GetXmin();
 double phimax = hY5D->GetAxis(PHI)->GetXmax();
 TH1D hsinphi("hsinphi", "hsinphi", nphibins, phimin, phimax);
 for (int ibin = 0; ibin < hsinphi.GetNbinsX(); ibin++){
   double phi = hsinphi.GetBinLowEdge(ibin+1) * TMath::DegToRad();
   hsinphi.SetBinContent(ibin+1, TMath::Sin(phi));
   hsinphi.SetBinError(ibin+1, 0);
 }
 //!4.2.2. CREATE hphiproj_theta[i]*sin(phi)
 for (int i=0;i<10;i++){
   hphiproj_theta[i]->Multiply(&hsinphi);
 }
 //!4.3 VERIFY
 for (int i=0;i<10;i++){
   double itgerr;
   double itg = hphiproj_theta[i]->IntegralAndError(1,hphiproj_theta[i]->GetNbinsX(),itgerr);
   itg = itg/TMath::Pi();
   itgerr = itgerr/TMath::Pi();
   itg = itg/norm;
   itgerr = itgerr/norm;
   printf("[%.2f:%.2f][[%.2f:%.2f]\n",h1D_THETA->GetBinContent(i+1),h1D_THETA->GetBinError(i+1),itg,itgerr);
 }
}
