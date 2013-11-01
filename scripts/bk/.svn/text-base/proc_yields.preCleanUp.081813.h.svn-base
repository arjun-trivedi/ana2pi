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
	bool _useHel;
		 
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
	
	//! [Q2,W] bin information
	char _Q2Wdirname[100];
	Int_t _iQ2Wbin; 
	Double_t _Q2low;
	Double_t _Q2high;
	Double_t _Wlow;
	Double_t _Whigh;
  
	//! Objects read from input file !//
	THnSparseF* hY8D[nVarset][2]; //[2] = {TH, RECO}
		
	//! Objects written to output file !//
	TDirectory* _dirQ2W;
	TDirectory* _dirVarset;
	TDirectory* _dirObs;
	
	//! hYW 
	//TH1D* hYW[nVarset][nSeq];
	
	//! Following objects made in bins of [Q2,W] !//
	TH2D* hrecoQ2vW[nVarset];
	TH2D* hthQ2vW[nVarset];  
	
	//!hY5D 
	THnSparse* hY5D[nVarset][nSeq];
		
	//! Display set per [Q2,W] bin
	TH2D* hY2D[nVarset][nVar][nSeq];
	TH1D* hY1D[nVarset][nVar][nSeq];
	TH1D* hY1D_phi_pos[nVarset][nSeq];
	TH1D* hY1D_phi_neg[nVarset][nSeq];
	THStack*  hs_hY2D[nSeq];
	THStack*  hs_hY1D[nSeq];
	THStack*  hs_hAsym[nSeq];
	CanvasManager* cm_hY2D;
	CanvasManager* cm_hY1D;
	CanvasManager* cm_hAsym;
	TString* cname_hY2D[nSeq];
	TString* ctitle_hY2D[nSeq];
	TString* cname_hY1D[nSeq];
	TString* ctitle_hY1D[nSeq];
	TString* cname_hAsym[nSeq];
	TString* ctitle_hAsym[nSeq];
	
	//! Canvas Manager settings
	static const Int_t _c_ww = 450;
	static const Int_t _c_wh = 350;

	//! misc vars
	char _hname[100];
	char _htitle[100];
	char _name[100];
	char _title[100];

	ProcYields(char* fyield, TString topList, 
		       Int_t nQ2bins, Double_t Q2min, Double_t Q2max, 
		       Int_t nWbins, Double_t Wmin, Double_t Wmax,
		       bool useHel=false);
	~ProcYields();

	/* *************************************************** */
	/*       ---------------------------------------       */
	/*        Methods for obtaining Observable Sets        */
	/*       ---------------------------------------       */
	/* *************************************************** */
  
  
	void proc();
	void proc_hYW_Dir(); 
	void proc_q2w(); 
	void proc_hY5D(hel_t hel=UNPOL);
  	void proc_hY1D(hel_t hel=UNPOL);
  	void proc_asym();
  	void proc_hYW();
  	
  	Float_t calcNorm(THnSparse* h5D_exp_ACC_CORR, THnSparse* h5D_sim_ACC_CORR);
		
	void setCanvasStackNameTitle(hel_t hel=UNPOL); //for all sets of observables in a given [Q2,W] bin: I, II & III
    
	void setM1M2axisrange(TH1* h, Int_t iVarset, var_t var /* =, M1/M2 */);
	Float_t getIntegral(THnSparse* hN);
  
	void setHistNameTitleAxes_hYW(TH1D* h, Int_t iVarset, Int_t iSeq); //atrivedi: 04-28-13
	void setHistNameTitleAxes_hY5D(THnSparse* h, Int_t iVarset, seq_t seq, hel_t hel=UNPOL);
	void setHistNameTitleAxes_hY2D(TH2D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel=UNPOL);
	void setHistNameTitleAxes_hY1D(TH1D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel=UNPOL);
};

#endif //PROCYIELDS_H
