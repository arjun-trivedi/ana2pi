#ifndef MOMCORR
#define MOMCORR

#include <TLorentzVector.h>

#define NSect 6

/* Theta Binning for Theta correction */
#define ThetaC_n 144
Double_t ThetaC_min = 0;
Double_t ThetaC_max = 144;
Double_t ThetaC_wid = 1.;

/* Theta Binning for Momentum correction */
#define MomC_T_n 48
Double_t MomC_T_min = 1;
Double_t MomC_T_max = 145;
Double_t MomC_T_wid = 3.;

/* Size of Vectors with parameters for phi dependence of Momentum Correction
   phi between -30 and 30 deg */
#define Npar 4

Int_t GetSector(Double_t phi) {
	//phi between -180 and 180
	//phi2 between 0 and 360
	Double_t phi2 = phi;
	if (phi<0) phi2 = phi2+360;
	Int_t sect = 1 + ((int)phi2+30)/60;
	if (sect == 7) sect = 1;
	return sect;
}


/* CLASS DEFINITION */


class MomCorr_e1f {
public:
  MomCorr_e1f(char* ddir /*= "MomCorr"*/);


private:
  char* datadir;
  /* Vectors with parameters for phi dependence of Theta Correction 
     phi between -30 and 30 deg */
  Double_t c0_theta[ThetaC_n][NSect];
  Double_t c1_theta[ThetaC_n][NSect];
  Double_t c2_theta[ThetaC_n][NSect];

  /* Vectors with parameters for phi dependence of Electron Momentum Correction
     phi between -30 and 30 deg */
  Double_t c0_mom[MomC_T_n][NSect][Npar];
  Double_t c1_mom[MomC_T_n][NSect][Npar];

  /* Vectors with parameters for phi dependence of pi+ Momentum Correction
     phi between -30 and 30 deg */
  Double_t d0_mom[MomC_T_n][NSect][Npar];
  Double_t d1_mom[MomC_T_n][NSect][Npar];


  /* REad angle correction parameters */
  void read_theta_par();
  /* Read momentum correction parameters for electrons */
  void read_mom_par();
  /* Read momentum correction parameters for pi+ */
  void read_mom_pip_par();

  /* Angle correction */
  Double_t theta_corr(Double_t , Double_t , Int_t );

  /* momentum correction for electrons */
  Double_t mom_corr(Double_t , Double_t , Double_t , Int_t );

  /* momentum correction for hadrons */
  Double_t mom_corr_pip(Double_t , Double_t , Double_t , Int_t);



public:
  TLorentzVector PcorN(TLorentzVector Pin, Int_t charge, Int_t ipart);


};



/********************************************/
/* Class Constructor */


MomCorr_e1f::MomCorr_e1f(char* ddir = "MomCorr/") {

  datadir = ddir;
  /* Reading angle correction parameters */
  fprintf(stdout, "reading mom cor parameters \n");

  read_theta_par();
  read_mom_par();
  read_mom_pip_par();

}


/**********************************************/
/* CLASS METHODES */
void MomCorr_e1f::read_theta_par() {

  /* Reading angle correction parameters */

  memset(&c0_theta[0][0], 0, ThetaC_n*NSect*sizeof(Double_t));
  memset(&c1_theta[0][0], 0, ThetaC_n*NSect*sizeof(Double_t));
  memset(&c2_theta[0][0], 0, ThetaC_n*NSect*sizeof(Double_t));

  char file[100];
  for (Int_t s=1; s<=NSect; s++) {
    sprintf(file, "%s/angles_s%d.out", datadir, s);
    FILE *fp = fopen(&file[0], "r");
    if (fp) {
      fprintf(stdout, "Theta Correction sector %d\n", s);
      for (Int_t j=0; j<ThetaC_n; j++) {
	Double_t x, y;
	Double_t a1,a2,a3;
	char str[200];

	while (fgets(str, sizeof(str), fp)) {
	  sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &x, &a1, &y, &a2, &y, &a3, &y);
	  Int_t bin = (Int_t)((x-ThetaC_min)/ThetaC_wid);
	  c0_theta[bin][s-1] = a1;
	  c1_theta[bin][s-1] = a2; 
	  c2_theta[bin][s-1] = a3;
	  fprintf(stdout, "%6.2f     %12.8f  %12.8f  %12.8f \n", x, c0_theta[bin][s-1], c1_theta[bin][s-1], c2_theta[bin][s-1]);
	}

      }

      fclose(fp);
    }
    else 
      fprintf(stdout, "===>>> WARNING: Cannot read angle correction from file: %s\n", file);

  }
}

void MomCorr_e1f::read_mom_par() {

  /* Reading momentum correction parameters for electrons */

  memset(&c0_mom[0][0][0], 0, MomC_T_n*NSect*Npar*sizeof(Double_t));
  memset(&c1_mom[0][0][0], 0, MomC_T_n*NSect*Npar*sizeof(Double_t));

  char file[100];
  for (Int_t s=1; s<=NSect; s++) {
    fprintf(stdout, "-------------------------------------- \n");
    fprintf(stdout, "electron Momentum Correction sector %d\n", s);

    for (Int_t k=0; k<2; k++) {

      sprintf(file, "%s/momentum2_s%d_c%d.out", datadir, s, k);
      FILE *fp = fopen(&file[0], "r");
      if (fp) {
	fprintf(stdout, "  Reading Parameter c%d\n", k);
	for (Int_t j=0; j<MomC_T_n; j++) {
	  Double_t x1, y;
	  Double_t a1,a2,a3,a4;
	  char str[250];
	  
	  while (fgets(str, sizeof(str), fp)) {
	    sscanf(str, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &x1, &a1, &y, &a2, &y, &a3, &y, &a4, &y);
	    
	    Int_t bin = (Int_t)((x1-MomC_T_min)/MomC_T_wid);
	    if (k==0) {
	      c0_mom[bin][s-1][0] = a1;
	      c0_mom[bin][s-1][1] = a2; 
	      c0_mom[bin][s-1][2] = a3;
	      c0_mom[bin][s-1][3] = a4;
	      fprintf(stdout, "%6.2f     %12.8f  %12.8f  %12.8f  %12.8f \n", x1, c0_mom[bin][s-1][0], c0_mom[bin][s-1][1], c0_mom[bin][s-1][2], c0_mom[bin][s-1][3]);
	    }
	    else if (k==1) {
	      c1_mom[bin][s-1][0] = a1;
	      c1_mom[bin][s-1][1] = a2; 
	      c1_mom[bin][s-1][2] = a3;
	      c1_mom[bin][s-1][3] = a4;
	      fprintf(stdout, "%6.2f     %12.8f  %12.8f  %12.8f  %12.8f \n", x1, c1_mom[bin][s-1][0], c1_mom[bin][s-1][1], c1_mom[bin][s-1][2], c1_mom[bin][s-1][3]);
	    }
	  }

	}
	fclose(fp);
      }
      else 
	fprintf(stdout, "===>>> WARNING: Cannot read electron momentum correction from file: %s\n", file);

    }


  }

}

void MomCorr_e1f::read_mom_pip_par() {

  /* Reading momentum correction parameters for pi+ */

  memset(&d0_mom[0][0][0], 0, MomC_T_n*NSect*Npar*sizeof(Double_t));
  memset(&d1_mom[0][0][0], 0, MomC_T_n*NSect*Npar*sizeof(Double_t));

  char file[100];
  for (Int_t s=1; s<=NSect; s++) {
    fprintf(stdout, "-------------------------------------- \n");
    fprintf(stdout, "pi+ Momentum Correction sector %d\n", s);

    for (Int_t k=0; k<2; k++) {

      sprintf(file, "%s/momentum3_s%d_c%d.out", datadir, s, k);
      FILE *fp = fopen(&file[0], "r");
      if (fp) {
	fprintf(stdout, "  Reading Parameter d%d\n", k);
	for (Int_t j=0; j<MomC_T_n; j++) {
	  Double_t x1, y;
	  Double_t a1,a2,a3;
	  char str[250];
	  
	  while (fgets(str, sizeof(str), fp)) {
	    sscanf(str, "%lf %lf %lf %lf %lf %lf %lf", &x1, &a1, &y, &a2, &y, &a3, &y);
	    
	    Int_t bin = (Int_t)((x1-MomC_T_min)/MomC_T_wid);
	    if (k==0) {
	      d0_mom[bin][s-1][0] = a1;
	      d0_mom[bin][s-1][1] = a2; 
	      d0_mom[bin][s-1][2] = a3;
	      fprintf(stdout, "%6.2f     %12.8f  %12.8f  %12.8f  %12.8f \n", x1, d0_mom[bin][s-1][0], d0_mom[bin][s-1][1], d0_mom[bin][s-1][2], d0_mom[bin][s-1][3]);
	    }
	    else if (k==1) {
	      d1_mom[bin][s-1][0] = a1;
	      d1_mom[bin][s-1][1] = a2; 
	      d1_mom[bin][s-1][2] = a3;
	      fprintf(stdout, "%6.2f     %12.8f  %12.8f  %12.8f  %12.8f \n", x1, d1_mom[bin][s-1][0], d1_mom[bin][s-1][1], d1_mom[bin][s-1][2], d1_mom[bin][s-1][3]);
	    }
	  }

	}
	fclose(fp);
      }
      else 
	fprintf(stdout, "===>>> WARNING: Cannot read pi+ momentum correction from file: %s\n", file);

    }


  }



}

/* ================================================================= */
Double_t MomCorr_e1f::theta_corr(Double_t ThetaM, Double_t PhiM, Int_t Sector) {

  /* input: phi between 0 and 360 deg */
  
  Double_t phis = PhiM - 60*(Sector-1);
  if (phis > 330.) phis = phis - 360;

  Int_t bin = (Int_t)((ThetaM - ThetaC_min)/ThetaC_wid);
  Double_t ThetaBin = ThetaC_min + ThetaC_wid * (bin+0.5);

  Int_t bin2 = bin + 1;
  if (ThetaM < ThetaBin) bin2 = bin - 1;
  if (bin2 < 0) bin2 = 0;
  if (bin2 >= ThetaC_n) bin2 = ThetaC_n - 1;
  Double_t ThetaBin2 = ThetaC_min + ThetaC_wid * (bin2+0.5);

  //fprintf(stdout, " ----------- \n");
  //fprintf(stdout, "T=%lf  Tb1=%lf  Tb2=%lf     b1=%d  b2=%d\n", ThetaM, ThetaBin, ThetaBin2, bin, bin2);

  //fprintf(stdout, "S=%d  th=%f  bin=%d  c0=%lf \n", Sector, ThetaM, bin, c0_theta[bin][Sector-1]);

  Double_t ThetaC = ThetaM - (c0_theta[bin][Sector-1] + c1_theta[bin][Sector-1]*phis + c2_theta[bin][Sector-1]*phis*phis);
  Double_t ThetaC2 = ThetaM - (c0_theta[bin2][Sector-1] + c1_theta[bin2][Sector-1]*phis + c2_theta[bin2][Sector-1]*phis*phis);
  Double_t fw = TMath::Abs(ThetaM - ThetaBin) / ThetaC_wid;
  Double_t fw2 = 1. - fw;
  Double_t Theta = fw2*ThetaC + fw*ThetaC2;
  //fprintf(stdout, " d1=%lf  w1=%lf      d2=%lf  w2=%lf     dm=%lf\n", ThetaC, fw, ThetaC2, fw2, Theta);

  return Theta;

}

Double_t MomCorr_e1f::mom_corr(Double_t MomM, Double_t ThetaM, Double_t PhiM, Int_t Sector) 
{
  //fprintf(stdout, " ----------- \n");


  /* input: phi between 0 and 360 deg */
  Double_t phis = PhiM - 60*(Sector-1);
  if (phis > 330.) phis = phis - 360;


  Int_t bin = (Int_t)((ThetaM - MomC_T_min)/MomC_T_wid);
  Double_t ThetaBin = MomC_T_min + MomC_T_wid*(bin+0.5);
  Double_t fw = TMath::Abs(ThetaM - ThetaBin) / MomC_T_wid;

  Double_t A0 = 0.;
  Double_t A1 = 0.;
  for (Int_t j=0; j<Npar; j++) {
    A0 = A0 + c0_mom[bin][Sector-1][j] *TMath::Power(MomM, j);
    A1 = A1 + c1_mom[bin][Sector-1][j] *TMath::Power(MomM, j);
  }
  Double_t MomC = MomM - (A0 + A1 * phis);


  Int_t bin2 = bin + 1;
  if (ThetaM < ThetaBin) bin2 = bin - 1;
  if (bin2 < 0) bin2 = 0;
  if (bin2 >= MomC_T_n) bin2 = MomC_T_n - 1;
  Double_t ThetaBin2 = MomC_T_min + MomC_T_wid*(bin2+0.5);
  Double_t fw2 = 1. - fw;


  A0 = A1 = 0.;
  for (Int_t j=0; j<Npar; j++) {
    A0 = A0 + c0_mom[bin2][Sector-1][j] *TMath::Power(MomM, j);
    A1 = A1 + c1_mom[bin2][Sector-1][j] *TMath::Power(MomM, j);
  }
  Double_t MomC2 = MomM - (A0 + A1 * phis);

  Double_t Mom = fw2*MomC + fw*MomC2;

  //fprintf(stdout, "T=%lf  Tb1=%lf  Tb2=%lf     b1=%d  b2=%d\n", ThetaM, ThetaBin, ThetaBin2, bin, bin2);
  //fprintf(stdout, "P=%lf    d1=%lf  w1=%lf      d2=%lf  w2=%lf     dm=%lf\n", MomM, MomC, fw, MomC2, fw2, Mom);

  return Mom;


}

Double_t MomCorr_e1f::mom_corr_pip(Double_t MomM, Double_t ThetaM, Double_t PhiM, Int_t Sector) 
{
  /* input: phi between 0 and 360 deg */
  Double_t phis = PhiM - 60*(Sector-1);
  if (phis > 330.) phis = phis - 360;


  Int_t bin = (Int_t)((ThetaM - MomC_T_min)/MomC_T_wid);
  Double_t ThetaBin = MomC_T_min + MomC_T_wid*(bin+0.5);
  Double_t fw = TMath::Abs(ThetaM - ThetaBin) / MomC_T_wid;

  Double_t A0 = 0.;
  Double_t A1 = 0.;
  for (Int_t j=0; j<Npar; j++) {
    A0 = A0 + d0_mom[bin][Sector-1][j] *TMath::Power(MomM, j);
    A1 = A1 + d1_mom[bin][Sector-1][j] *TMath::Power(MomM, j);
  }
  Double_t MomC = MomM - (A0 + A1 * phis);


  Int_t bin2 = bin + 1;
  if (ThetaM < ThetaBin) bin2 = bin - 1;
  if (bin2 < 0) bin2 = 0;
  if (bin2 >= MomC_T_n) bin2 = MomC_T_n - 1;
  Double_t ThetaBin2 = MomC_T_min + MomC_T_wid*(bin2+0.5);
  Double_t fw2 = 1. - fw;


  A0 = A1 = 0.;
  for (Int_t j=0; j<Npar; j++) {
    A0 = A0 + d0_mom[bin2][Sector-1][j] *TMath::Power(MomM, j);
    A1 = A1 + d1_mom[bin2][Sector-1][j] *TMath::Power(MomM, j);
  }
  Double_t MomC2 = MomM - (A0 + A1 * phis);

  Double_t Mom = fw2*MomC + fw*MomC2;

  //fprintf(stdout, " ----------- \n");
  //fprintf(stdout, "T=%lf  Tb1=%lf  Tb2=%lf     b1=%d  b2=%d\n", ThetaM, ThetaBin, ThetaBin2, bin, bin2);

  //fprintf(stdout, "P=%lf    d1=%lf  w1=%lf      d2=%lf  w2=%lf     dm=%lf\n", MomM, MomC, fw, MomC2, fw2, Mom);

  return Mom;


}

/* ================================================================= */
TLorentzVector MomCorr_e1f::PcorN(TLorentzVector Pin, Int_t charge, Int_t ipart) {
  TLorentzVector Pc;

  Double_t mass = Pin.M();

  Double_t theta = Pin.Theta() * 180./TMath::Pi();
  Double_t phi = Pin.Phi() * 180./TMath::Pi();
  if (phi < 0) phi = phi + 360.;
  Int_t sect = GetSector(phi);

  Double_t theta_c = theta_corr(theta, phi, sect);
  Double_t mom_c = Pin.P();
  if (ipart == 11) {/*   ELETTRONE */
    mom_c = mom_corr(Pin.P(), theta, phi, sect);
  }
  else if (ipart == -211) {/* PI- */
    mom_c = mom_corr(Pin.P(), theta, phi, sect);
  }
  else if (ipart == -321) {/* K- */
    mom_c = mom_corr(Pin.P(), theta, phi, sect);
  }
  else if (ipart == 211) {/* PI+ */
    mom_c = mom_corr_pip(Pin.P(), theta, phi, sect);
  }
  else if (ipart == 2212) {/* PROTONE */
    mom_c = mom_corr_pip(Pin.P(), theta, phi, sect);
  }
  else if (ipart == 321) {/* K+ */
    mom_c = mom_corr_pip(Pin.P(), theta, phi, sect);
  }
  Double_t px = mom_c * TMath::Sin(theta_c*TMath::Pi()/180.) * TMath::Cos(Pin.Phi());
  Double_t py = mom_c * TMath::Sin(theta_c*TMath::Pi()/180.) * TMath::Sin(Pin.Phi());
  Double_t pz = mom_c * TMath::Cos(theta_c*TMath::Pi()/180.);

  Pc.SetXYZM(px, py, pz, mass);

  return Pc;

}

#endif /* MOMCORR */
