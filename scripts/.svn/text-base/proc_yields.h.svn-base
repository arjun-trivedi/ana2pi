#ifndef PROCYIELDS_H
#define PROCYIELDS_H

#include <TROOT.h>
#include <THnSparse.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>
#include <THStack.h>
#include <TCanvas.h>

#include "constants.h"
#include "CanvasManager.h"

class ProcYields{
public: 
	TFile* _fin;
	TFile* _fin_th;
	TFile* _fout;
	TFile* _fout_sim;
	
	bool _issim;
	bool _isexp;
	bool _usehel;
		 
	TString _topList;

	//! Q2 and W binning: "intrinsic" and user input !//
	struct Q2Wbng{
		Int_t nQ2bins;
		Double_t Q2min;
		Double_t Q2max;
		Double_t Q2binw;
	
		Int_t nWbins;
		Double_t Wmin;
		Double_t Wmax;
		Double_t Wbinw;
	
		Int_t nQ2Wbins; 
	} _intrinsic,_user; 
	
	//! [Q2,W] bin information !//
	char _Q2Wdirname[100];
	Int_t _iQ2Wbin; 
	Double_t _Q2low;
	Double_t _Q2high;
	Double_t _Wlow;
	Double_t _Whigh;
  
	//! Objects read from input file !//
	THnSparseF* _hY8D[nVARSET][2]; //[2] = {TH, RECO}
		
	//! Objects written to output file !//
	TDirectory* _dirQ2W;
	
	//!hY5D !//
	THnSparse* _hY5D[nVARSET][nSEQ];
		
	//! Display set per [Q2,W] bin
	TH2D* _hY2D[nVARSET][nVAR][nSEQ];
	TH1D* _hY1D[nVARSET][nVAR][nSEQ];
	TH1D* _hY1D_phi_pos[nVARSET][nSEQ];
	TH1D* _hY1D_phi_neg[nVARSET][nSEQ];
	THStack*  _hs_hY2D[nSEQ];
	THStack*  _hs_hY1D[nSEQ];
	THStack*  _hs_hAsym[nSEQ];
	CanvasManager* _cm_hY2D;
	CanvasManager* _cm_hY1D;
	CanvasManager* _cm_hAsym;
	TString* _cname_hY2D[nSEQ];
	TString* _ctitle_hY2D[nSEQ];
	TString* _cname_hY1D[nSEQ];
	TString* _ctitle_hY1D[nSEQ];
	TString* _cname_hAsym[nSEQ];
	TString* _ctitle_hAsym[nSEQ];
	
	//! Canvas Manager settings !//
	static const Int_t _c_ww = 450;
	static const Int_t _c_wh = 350;

	//! misc vars !//
	char _hname[100];
	char _htitle[100];
	char _name[100];
	char _title[100];

	ProcYields(char* fyield, TString topList, 
		       Int_t nQ2bins, Double_t Q2min, Double_t Q2max, 
		       Int_t nWbins, Double_t Wmin, Double_t Wmax,
		       bool usehel=false);
	~ProcYields();

	/* *************************************************** */
	/*       ---------------------------------------       */
	/*        Methods for obtaining Observable Sets        */
	/*       ---------------------------------------       */
	/* *************************************************** */
  
  
	void Proc();
	void Proc_hYW_Dir(); 
	void Proc_q2w(); 
	void Proc_hY5D(hel_t hel=UNPOL);
  	void Proc_hY1D(hel_t hel=UNPOL);
  	void Proc_hPhi(hel_t hel=UNPOL);
  	void Proc_asym();
  	void Proc_hYW();
  	
  	Float_t calcNorm(THnSparse* h5D_exp_ACC_CORR, THnSparse* h5D_sim_ACC_CORR);
		
	void setCanvasStackNameTitle(hel_t hel=UNPOL); //for all sets of observables in a given [Q2,W] bin: I, II & III
    
	void setM1M2axisrange(TH1* h, Int_t iVarset, var_t var /* =, M1/M2 */);
	Float_t getIntegral(THnSparse* hN);
  
	void setHistNameTitleAxes_hYW(TH1D* h, Int_t iVarset, Int_t iSeq); //atrivedi: 04-28-13
	void setHistNameTitleAxes_hY5D(THnSparse* h, Int_t iVarset, seq_t seq, hel_t hel=UNPOL);
	void setHistNameTitleAxes_hY2D(TH2D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel=UNPOL);
	void setHistNameTitleAxes_hY1D(TH1D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel=UNPOL);
	void setHistNameTitleAxes_hAsym(TH1D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel=UNPOL);

	/************** for Luminosity & vgflux ********************/
	double nu(double w, double q2) {
    	return (w*w-MP*MP+q2)/(2*MP);
	}

	double epsilon(double w, double q2) {
  		double n = nu(w,q2);
  		double e0 = E1F_E0;
  		double e1 = e0-n;
  		double epsInv = 1+2*(q2+n*n)/(4*e0*e1-q2);
    	return 1.0/epsInv;
	}

	double getvgflux(double w, double q2, double e0 = E1F_E0) {
  		double eps = epsilon(w,q2);
  		return A*w*(w*w-MP*MP)/(4*PI*e0*e0*MP*MP*q2*(1-eps));
	}
	/**************************/
};

#endif //PROCYIELDS_H
