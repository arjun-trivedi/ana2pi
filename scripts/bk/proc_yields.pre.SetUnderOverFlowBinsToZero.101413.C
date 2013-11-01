#include "proc_yields.h"
#include <TLegend.h>
#include <TString.h>
#include <TKey.h>
#include <TStyle.h>

#include <iostream>
#include <vector>

using namespace std;

ProcYields::ProcYields(char* fyield, TString topList, 
					   Int_t nQ2bins, Double_t Q2min, Double_t Q2max, 
					   Int_t nWbins, Double_t Wmin, Double_t Wmax,
					   bool usehel/*=false*/){
  
  _fin      = NULL;
  _fin_th   = NULL;
  _fout     = NULL;
  _fout_sim = NULL;
  
  _issim  = false;
  _isexp  = false;
  _usehel = usehel;
    
  _dirQ2W    = NULL;
    
  for (Int_t iVarset=0;iVarset<nVARSET;iVarset++){
	_hY8D[iVarset][RECO] = NULL;
	_hY8D[iVarset][TH]   = NULL;
  }
  
  //! Set up _intrinsic Q2W binning !//
  _intrinsic.nQ2bins = 100;
  _intrinsic.Q2min   = 0;
  _intrinsic.Q2max   = 5;
  _intrinsic.Q2binw  = (_intrinsic.Q2max-_intrinsic.Q2min)/_intrinsic.nQ2bins;
  _intrinsic.nWbins  = 400;
  _intrinsic.Wmin    = 1;
  _intrinsic.Wmax    = 3;
  _intrinsic.Wbinw   = (_intrinsic.Wmax-_intrinsic.Wmin)/_intrinsic.nWbins;
  
  _intrinsic.nQ2Wbins = _intrinsic.nQ2bins*_intrinsic.nWbins;
  
  //! open fyield & determine if data _issim OR _isexp !//
  _fin = TFile::Open(fyield);
  
  if (_fin->Get("mom") == NULL) _issim=true;
  else _isexp = true;
  
  if (_issim) _fin_th = TFile::Open("mcyield.root", "READ");
  
  //! Determine Top !//
  _topList = topList;
  vector<Int_t> vTopList;
  TCollection* cTopList = _topList.Tokenize(":");
  TIter iter(cTopList);
  while(TObjString *objStr = (TObjString*)iter.Next()) {
	  TString str = objStr->GetString();
	  vTopList.push_back(atoi(str.Data()));
  }
  Info("ProcYields()", "Going to process yields from the following top selections:");
  for (vector<Int_t>::iterator it = vTopList.begin(); it != vTopList.end(); ++it)
    printf(" %d", *it);
  printf("\n \n");
    
  //! get user supplied Q2,W binning
  _user.nQ2bins = nQ2bins;
  _user.Q2min   = Q2min;
  _user.Q2max   = Q2max;
  _user.Q2binw  = (_user.Q2max-_user.Q2min)/_user.nQ2bins;
  _user.nWbins  = nWbins;
  _user.Wmin    = Wmin;
  _user.Wmax    = Wmax;
  _user.Wbinw   = (_user.Wmax-_user.Wmin)/_user.nWbins;
  
  _user.nQ2Wbins = _user.nQ2bins*_user.nWbins;
  
  Info("ProcYields()", "[nQ2bins, Q2min, Q2max][nWbins, Wmin, Wmax] = [%d, %4.3f, %4.3f][%d, %4.3f, %4.3f]\n",
  						_user.nQ2bins, _user.Q2min, _user.Q2max,_user.nWbins,  _user.Wmin,  _user.Wmax);
  Info("ProcYields()", "[Q2binw][Wbinw] = [%4.3f][%4.3f]\n", _user.Q2binw, _user.Wbinw);
  
  //! create output file !//
  char foutname[100];
  if     (_issim) {
  	if (_usehel) sprintf(foutname, "%s__%d-%4.3f-%4.3f__%d-%4.3f-%4.3f__pol__sim.root", _topList.Data(), _user.nQ2bins, _user.Q2min, _user.Q2max, _user.nWbins, _user.Wmin, _user.Wmax);
  	else         sprintf(foutname, "%s__%d-%4.3f-%4.3f__%d-%4.3f-%4.3f__sim.root", _topList.Data(), _user.nQ2bins, _user.Q2min, _user.Q2max, _user.nWbins, _user.Wmin, _user.Wmax); 
  }
  else if(_isexp){
  	if (_usehel){sprintf(foutname, "%s__%d-%4.3f-%4.3f__%d-%4.3f-%4.3f__pol__exp.root", _topList.Data(), _user.nQ2bins, _user.Q2min, _user.Q2max, _user.nWbins, _user.Wmin, _user.Wmax);}	
  	else        {sprintf(foutname, "%s__%d-%4.3f-%4.3f__%d-%4.3f-%4.3f__exp.root",      _topList.Data(), _user.nQ2bins, _user.Q2min, _user.Q2max, _user.nWbins, _user.Wmin, _user.Wmax);}
  }else{
	  Info("ProcYields()", "input file is neither mc, sim or exp. Exiting...\n");
	  return; 
  }
  _fout = new TFile(foutname,"RECREATE");
  
  //! if processing experimental data, get corresponding _sim.root file for acceptances !//
  if (_isexp){
	TString tmp(_fout->GetName());
	tmp.ReplaceAll("__pol", "");
	tmp.ReplaceAll("exp", "sim");
	TString fname = TString::Format("simdir/%s", tmp.Data());
	printf("Going to get acceptance parameters from %s\n", fname.Data());
	if(TFile::Open(fname.Data()) != NULL)	{
		_fout_sim = TFile::Open(fname.Data());
	}else{
		Info("ProcYields()", "%s DOES NOT EXIST!\n",tmp.Data());
		return;
	}
  }else{_fout_sim = NULL;}
  
  Info("ProcYields()", "Input files: %s %s\n", _fin->GetName(), _fin_th != NULL?_fin_th->GetName():"");
  Info("ProcYields()", "Output file: %s\n", _fout->GetName());
 
  //! get yields !//
  char hreco[100];
  char hth[100];
  for (vector<Int_t>::iterator it = vTopList.begin(); it != vTopList.end(); ++it){//loop over top
    Int_t top = *it;
    for (Int_t iVarset=0;iVarset<nVARSET;iVarset++){ //loop over Varsets
	  sprintf(hreco, "top/top%d/yield_varset%d", top, iVarset+1);
	  sprintf(hth,   "top/mc/yield_varset%d", iVarset+1);
	  if(_hY8D[iVarset][RECO]==NULL){
	    _hY8D[iVarset][RECO] = (THnSparseF*)_fin->Get(hreco);
	    if(_issim) _hY8D[iVarset][TH] = (THnSparseF*)_fin_th->Get(hth);
	  }else{
	    _hY8D[iVarset][RECO]->Add((THnSparseF*)_fin->Get(hreco)); 	  
	  }
	}
  }
  _fin->Close();
  
}

ProcYields::~ProcYields(){
  if (_fout_sim != NULL){_fout_sim->Close();}
  if (_fin_th != NULL){_fin_th->Close();}
}

void ProcYields::Proc(){
  Info("Proc()", "");
     
  //!Before obtaining hY5Ds in Q2W bins, for a given Q2, get differential cross section vs W 'Directly'
  if (!_usehel)Proc_hYW_Dir();
  
  //!Now loop over [Q2,W] bin and Proc_5D()
  Proc_q2w();
  
  //!hYW
  if (!_usehel)Proc_hYW();
  
  Info("Proc()", "Going to Write output file");
  _fout->Write();
  Info("Proc()", "output file written");
  Info("Proc()", "done\n");
  
  return;
}

void ProcYields::Proc_q2w(){
	Info("Proc_q2w()", "");
	_iQ2Wbin = 0;
	for(Int_t iQ2bin = 0; iQ2bin < _user.nQ2bins; iQ2bin++){
		for(Int_t iWbin = 0; iWbin < _user.nWbins; iWbin++){
			printf("--------------------\n");
			printf("[iQ2bin:iWbin] = [%d:%d]\n", iQ2bin, iWbin);
      
			_Q2low  = _user.Q2min + (iQ2bin*_user.Q2binw);
			_Q2high = _user.Q2min + (iQ2bin*_user.Q2binw) + _user.Q2binw;
			_Wlow   = _user.Wmin  + (iWbin*_user.Wbinw);
			_Whigh  = _user.Wmin  + (iWbin*_user.Wbinw)  + _user.Wbinw;
      
			sprintf(_Q2Wdirname, "[%4.3f,%4.3f]_[%4.3f,%4.3f]", _Q2low,  _Q2high, _Wlow, _Whigh);
			printf("[Q2low,Q2high]_[Wlow,Whigh] = %s\n", _Q2Wdirname);
      
			_dirQ2W = _fout->mkdir(_Q2Wdirname);
			_dirQ2W->cd();
      
			//! Loop nVARSET & set hY8D bin ranges as per [Q2,W] bin
			TH2D* hrecoQ2vW[nVARSET];
			TH2D* hthQ2vW[nVARSET]; 
			for(Int_t iVarset=0;iVarset<nVARSET;iVarset++){
				Int_t Q2binl, Q2binh, Wbinl, Wbinh =0;
				if(_issim || _isexp){
					//! set bin ranges for axes Q2 & W !//
					Q2binl =  _hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_Q2low  + _intrinsic.Q2binw/2);
					Q2binh =  _hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_Q2high - _intrinsic.Q2binw/2); 
					Wbinl  =  _hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_Wlow  + _intrinsic.Wbinw/2);
					Wbinh  =  _hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_Whigh - _intrinsic.Wbinw/2); 
					printf("[Q2binl,Q2binh]_[Wbinl,Wbinh] = [%d,%d]_[%d,%d]\n",Q2binl, Q2binh, Wbinl, Wbinh);
						
					_hY8D[iVarset][RECO]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
					_hY8D[iVarset][RECO]->GetAxis(W)->SetRange(Wbinl, Wbinh); 
						
					hrecoQ2vW[iVarset] = (TH2D*)_hY8D[iVarset][RECO]->Projection(1,2,"colz");
					sprintf(_hname, "recoQ2vW_Varset%d", iVarset+1);
					hrecoQ2vW[iVarset]->SetNameTitle(_hname, (TString(_hname) + TString(_Q2Wdirname)).Data());
				}
		
				if(_issim){
					//! set bin ranges for axes Q2 & W !//
					Q2binl =  _hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_Q2low  + _intrinsic.Q2binw/2);
					Q2binh =  _hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_Q2high - _intrinsic.Q2binw/2); 
					Wbinl  =  _hY8D[iVarset][TH]->GetAxis(W)->FindBin(_Wlow  + _intrinsic.Wbinw/2);
					Wbinh  =  _hY8D[iVarset][TH]->GetAxis(W)->FindBin(_Whigh - _intrinsic.Wbinw/2); 
					printf("[Q2binl,Q2binh]_[Wbinl,Wbinh] = [%d,%d]_[%d,%d]\n",Q2binl, Q2binh, Wbinl, Wbinh);
						
					_hY8D[iVarset][TH]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
					_hY8D[iVarset][TH]->GetAxis(W)->SetRange(Wbinl, Wbinh);
						
					hthQ2vW[iVarset] = (TH2D*)_hY8D[iVarset][TH]->Projection(1,2,"colz");
					sprintf(_hname, "thQ2vW_Varset%d", iVarset+1);
					hthQ2vW[iVarset]->SetNameTitle(_hname, (TString(_hname) + TString(_Q2Wdirname)).Data());
				}
			}
	  
			if (_usehel){
				//! 1. hY5D+
				//! 2. hY5D-
				//! 3. hY5D+/- -> (hY5D+ - hY5D-) -> hAsym
				Proc_hY5D(); //UNPOL 
				Proc_hY5D(POS); //1
				Proc_hY5D(NEG);	//2
				Proc_hPhi(); //UNPOL
				Proc_hPhi(POS);
				Proc_hPhi(NEG);
				//Proc_asym();	//3
	  		}else{
	  			setCanvasStackNameTitle(UNPOL);
	  			//! 1. hY5D
	  			//! 2. hY5D -> hY2D -> hY1D 
	  			Proc_hY5D(); //1
				Proc_hY1D(); //2
			}
	  	  	  	  
			_iQ2Wbin += 1;
		}
	}
	Info("Proc_q2w()", "Done\n");
}

void ProcYields::Proc_hY5D(hel_t hel/*=UNPOL*/){
	Info(TString::Format("Proc_hY5D(%s)",helTitle[hel].Data()), "");
	
	TDirectory* dirY5D=NULL;
	if      (hel==UNPOL)            dirY5D = _dirQ2W->mkdir("hY5D");
	else if (hel==POS||hel==NEG)    dirY5D = _dirQ2W->mkdir( TString::Format("hY5D_%s",helTitle[hel].Data()) );
	
	TDirectory* dirVarset=NULL;	
	for(Int_t iVarset=0;iVarset<nVARSET;iVarset++){
		dirVarset = dirY5D->mkdir(TString::Format("Varset%d", iVarset+1));
		dirVarset->cd();
		
		//!RECO
		if(_issim || _isexp){
			if      (hel==UNPOL)         _hY8D[iVarset][RECO]->GetAxis(H)->SetRange();
			else if (hel==POS||hel==NEG) _hY8D[iVarset][RECO]->GetAxis(H)->SetRange(hel,hel);
			_hY5D[iVarset][RECO] = (THnSparse*)_hY8D[iVarset][RECO]->Projection(5,PROJDIMS,"E");
			setHistNameTitleAxes_hY5D(_hY5D[iVarset][RECO], iVarset, RECO, hel);
			_hY5D[iVarset][RECO]->Write();
	    }
		//! TH, ACC
		if (_issim){
			if      (hel==UNPOL)         _hY8D[iVarset][TH]->GetAxis(H)->SetRange();
			else if (hel==POS||hel==NEG) _hY8D[iVarset][TH]->GetAxis(H)->SetRange(hel,hel);
			_hY5D[iVarset][TH] = (THnSparse*)_hY8D[iVarset][TH]->Projection(5,PROJDIMS,"E");
			setHistNameTitleAxes_hY5D(_hY5D[iVarset][TH], iVarset, TH, hel);
			_hY5D[iVarset][TH]->Write();
						
			//!calculate acceptance
			_hY5D[iVarset][ACC] = (THnSparse*)_hY5D[iVarset][RECO]->Clone();
			_hY5D[iVarset][ACC]->Divide(_hY5D[iVarset][TH]);
			setHistNameTitleAxes_hY5D(_hY5D[iVarset][ACC], iVarset, ACC, hel);
			_hY5D[iVarset][ACC]->Write();
		}else if(_isexp){//if _isexp, then open _fout_sim to get ACC
			sprintf(_hname, "%s/%s/Varset%d/hY5D_%s", _Q2Wdirname, "hY5D", iVarset+1, seqTitle[ACC].Data()); 
			_hY5D[iVarset][ACC]=(THnSparse*)_fout_sim->Get(_hname);
			_hY5D[iVarset][ACC]->Write();
		}
		//! ACC_CORR
		if(_issim || _isexp){
			_hY5D[iVarset][ACC_CORR] = (THnSparse*)_hY5D[iVarset][RECO]->Clone();
			_hY5D[iVarset][ACC_CORR]->Divide(_hY5D[iVarset][ACC]);
			setHistNameTitleAxes_hY5D(_hY5D[iVarset][ACC_CORR], iVarset, ACC_CORR, hel);
			_hY5D[iVarset][ACC_CORR]->Write();
		}

		/*if (hel==POS || hel==NEG){
			Info(TString::Format("Proc_hY5D(%s)",helTitle[hel].Data()), "done\n");
			return;
		}*/

		//! HOLE
		if(_issim){
			_hY5D[iVarset][HOLE] = (THnSparse*)_hY5D[iVarset][TH]->Clone();
			_hY5D[iVarset][HOLE]->Add(_hY5D[iVarset][ACC_CORR], -1);
			setHistNameTitleAxes_hY5D(_hY5D[iVarset][HOLE], iVarset, HOLE, hel);
			_hY5D[iVarset][HOLE]->Write();
		}else if (_isexp){//if _isexp, then open _fout_sim to get sim-HOLE, sim-ACC_CORR 
			//!first copy sim-HOLE into exp-HOLE
			sprintf(_hname, "%s/%s/Varset%d/hY5D_%s", _Q2Wdirname, "hY5D", iVarset+1, seqTitle[HOLE].Data());
			_hY5D[iVarset][HOLE]=(THnSparse*)_fout_sim->Get(_hname);
			
			//!now get sim-ACC_CORR for obtaining normalization factor
			sprintf(_hname, "%s/%s/Varset%d/hY5D_%s", _Q2Wdirname, "hY5D", iVarset+1, seqTitle[ACC_CORR].Data());
			THnSparse* hY5D_sim_ACC_CORR = (THnSparse*)_fout_sim->Get(_hname);
			
			//!calculate normalization factor
			//Float_t norm = _hY5D[iVarset][ACC_CORR]->GetSumw()/hY5D_sim_ACC_CORR->GetSumw();
			Float_t norm = calcNorm(_hY5D[iVarset][ACC_CORR], hY5D_sim_ACC_CORR);
			cout << "norm = " << norm << endl;
			
			//!normalize exp-HOLE
			_hY5D[iVarset][HOLE]->Scale(norm);
			
			_hY5D[iVarset][HOLE]->Write();
		}
		//! FULL
		if (_issim || _isexp){
			_hY5D[iVarset][FULL] = (THnSparse*)_hY5D[iVarset][ACC_CORR]->Clone();
			_hY5D[iVarset][FULL]->Add(_hY5D[iVarset][HOLE],1 );
			setHistNameTitleAxes_hY5D(_hY5D[iVarset][FULL], iVarset, FULL, hel);
			_hY5D[iVarset][FULL]->Write();
		}
	}
	
	Info(TString::Format("Proc_hY5D(%s)",helTitle[hel].Data()), "done\n");
}

void ProcYields::Proc_hY1D(hel_t hel/*=UNPOL*/){
	Info(TString::Format("Proc_hY1D(%s)",helTitle[hel].Data()), "");
	
	//! Create Canvas Managers !//
	_cm_hY2D = new CanvasManager(_c_ww, _c_wh);
    _cm_hY1D = new CanvasManager(_c_ww, _c_wh);
	//! Create Stacks for TH, RECO, ACC_CORR, FULL !//
	for(Int_t iSeq=0;iSeq<nSEQ;iSeq++){
		seq_t seq = static_cast<seq_t>(iSeq);
		//if (seq==ACC || seq==HOLE) continue;
		if (seq==ACC) continue;
		if (_isexp && seq==TH) continue;
		_hs_hY2D[iSeq] = new THStack(*_cname_hY2D[iSeq], *_ctitle_hY2D[iSeq]);
		_hs_hY1D[iSeq] = new THStack(*_cname_hY1D[iSeq], *_ctitle_hY1D[iSeq]);	
	}
		
	//! hY5D -> hY2D -> hY1D !//
	TDirectory* dirY2D=NULL;
	TDirectory* dirY1D=NULL;
	if (hel==UNPOL)	{
		dirY2D = _dirQ2W->mkdir("hY2D");
		dirY1D = _dirQ2W->mkdir("hY1D");
	}else if(hel==POS||hel==NEG){
		dirY2D = _dirQ2W->mkdir(TString::Format("hY2D_%s",helTitle[hel].Data()));
		dirY1D = _dirQ2W->mkdir(TString::Format("hY1D_%s",helTitle[hel].Data()));
	}

	TDirectory* dirY2D_Varset;
	TDirectory* dirY1D_Varset;
	for (Int_t iVarset=0;iVarset<nVARSET;iVarset++){
		dirY2D_Varset = dirY2D->mkdir(TString::Format("Varset%d", iVarset+1));
		dirY1D_Varset = dirY1D->mkdir(TString::Format("Varset%d", iVarset+1));
		
		for(Int_t iVar=0; iVar<nVAR; iVar++){
			var_t var = static_cast<var_t>(iVar);
			if (var == ALPHA) continue; 
			for(Int_t iSeq=0;iSeq<nSEQ;iSeq++){
				seq_t seq = static_cast<seq_t>(iSeq);
				//if (seq==ACC || seq==HOLE) continue;
				if (seq==ACC) continue;
				if (_isexp && seq==TH) continue;
								
				//! hY2D !//
				dirY2D_Varset->cd();
				_hY2D[iVarset][iVar][iSeq] = (TH2D*)_hY5D[iVarset][iSeq]->Projection(PHI, iVar, "E");//0819
				setHistNameTitleAxes_hY2D(_hY2D[iVarset][iVar][iSeq], iVarset, iVar, iSeq);
						
				//! hY1D !//
				dirY1D_Varset->cd();
				_hY1D[iVarset][iVar][iSeq] = (TH1D*)_hY2D[iVarset][iVar][iSeq]->ProjectionX("_px",0,-1,"e")->Clone();//0819
				_hY1D[iVarset][iVar][iSeq]->SetMinimum(0);
				setHistNameTitleAxes_hY1D(_hY1D[iVarset][iVar][iSeq], iVarset, iVar, iSeq);
				            			
				if (var == M1 or var == M2){
					setM1M2axisrange(_hY2D[iVarset][iVar][iSeq], iVarset, var);
					setM1M2axisrange(_hY1D[iVarset][iVar][iSeq], iVarset, var);
				}
			
				//! add to stack !//
				_hs_hY2D[iSeq]->Add(_hY2D[iVarset][iVar][iSeq], "colz");
				_hs_hY1D[iSeq]->Add(_hY1D[iVarset][iVar][iSeq], "e1");
        
			}//end nSEQ loop
		}//end nVAR loop
	}//end nVARSET loop
	
	if(_isexp || _issim){
		_cm_hY2D->add(_hs_hY2D[RECO]);
		_cm_hY1D->add(_hs_hY1D[RECO]);
		_cm_hY2D->add(_hs_hY2D[ACC_CORR]);
		_cm_hY1D->add(_hs_hY1D[ACC_CORR]);
		_cm_hY2D->add(_hs_hY2D[HOLE]);
		_cm_hY1D->add(_hs_hY1D[HOLE]);
		_cm_hY2D->add(_hs_hY2D[FULL]);
		_cm_hY1D->add(_hs_hY1D[FULL]);
	}
	if(_issim){
		_cm_hY2D->add(_hs_hY2D[TH]);
		_cm_hY1D->add(_hs_hY1D[TH]);
	}
	
	//! Save Canvases to file !//
	_dirQ2W->cd();
	_cm_hY2D->draw();
	_cm_hY2D->save();
	_cm_hY2D->close();
	
	_cm_hY1D->draw();
	_cm_hY1D->save();
	_cm_hY1D->close();
	
	Info(TString::Format("Proc_hY5D(%s)",helTitle[hel].Data()), "done\n");
}

void ProcYields::Proc_hPhi(hel_t hel/*=UNPOL*/){
	Info(TString::Format("Proc_hPhi(%s)",helTitle[hel].Data()), "");

	TDirectory* dir_phi=NULL;
	if (hel==UNPOL)	{
		dir_phi = _dirQ2W->mkdir("hPhi");
	}else if(hel==POS||hel==NEG){
		dir_phi = _dirQ2W->mkdir(TString::Format("hPhi_%s",helTitle[hel].Data()));
	}

	TDirectory* dir_varset=NULL;
	for (Int_t iVarset=0;iVarset<nVARSET;iVarset++){
		dir_varset = dir_phi->mkdir(TString::Format("Varset%d", iVarset+1));
		dir_varset->cd();
		
		//!Get relevent h5D
		THnSparse* hY5D = NULL;
		if (hel==UNPOL)	{
			hY5D = (THnSparseF*)_dirQ2W->Get(TString::Format("hY5D/Varset%d/hY5D_FULL",iVarset+1));
			//! To make hphi from TH evts//hY5D = (THnSparseF*)_dirQ2W->Get(TString::Format("hY5D/Varset%d/hY5D_TH",iVarset+1));
		}else if (hel==POS||hel==NEG){
			hY5D = (THnSparseF*)_dirQ2W->Get(TString::Format("hY5D_%s/Varset%d/hY5D_FULL",helTitle[hel].Data(),iVarset+1));
			//! To make hphi from TH evts//hY5D = (THnSparseF*)_dirQ2W->Get(TString::Format("hY5D_%s/Varset%d/hY5D_TH",helTitle[hel].Data(),iVarset+1));
		}

		int nphibins = hY5D->GetAxis(PHI)->GetNbins();
		float phimin = hY5D->GetAxis(PHI)->GetXmin();
		float phimax = hY5D->GetAxis(PHI)->GetXmax();
		TH1D* hsinphi = new TH1D("hsinphi", "hsinphi", nphibins, phimin, phimax);
		for (int ibin = 0; ibin < hsinphi->GetNbinsX(); ibin++){
			Double_t phi = hsinphi->GetBinLowEdge(ibin+1) * TMath::DegToRad();
			hsinphi->SetBinContent(ibin+1, TMath::Sin(phi));
			hsinphi->SetBinError(ibin+1, 0);
		}

		//! Make hRvVar for each Var
		TDirectory* dir_var=NULL;
		for(Int_t iVar=0; iVar<nVAR; iVar++){
			if (iVar==ALPHA || iVar==PHI) continue; 
			dir_var = dir_varset->mkdir(varName[iVar].Data());
			dir_var->cd();

			//! Create hRvVar
			int nvarbins = hY5D->GetAxis(iVar)->GetNbins();
			float varmin = hY5D->GetAxis(iVar)->GetXmin();
			float varmax = hY5D->GetAxis(iVar)->GetXmax();
			//printf("numphibins:numvarbins = %d, %d\n", nphibins,nvarbins);
			TString xtitle = TString::Format("%s", varTitle[iVarset][iVar].Data());
			TString ytitle = TString::Format("hPR");
			TString title = TString::Format("%s vs %s [h=%s] %s", ytitle.Data(), xtitle.Data(), helTitle[hel].Data(), _Q2Wdirname);
			TH1F* hRvVar = new TH1F("hRvVar",title, nvarbins, varmin, varmax);
			hRvVar->SetXTitle(xtitle);
			hRvVar->SetYTitle(ytitle);

			//! Loop over number of bins in Var and make phi projections for each
			if (iVar==M1 || iVar==M2){
				for (int ivarbin=0; ivarbin<5; ivarbin++){
					hRvVar->SetBinContent(ivarbin+1, 20);
				}
			}else{
				for (int ivarbin=0; ivarbin<nvarbins; ivarbin++){
					Float_t varbin_lowedge = hY5D->GetAxis(iVar)->GetBinLowEdge(ivarbin+1);
					Float_t varbin_highedge = varbin_lowedge + hY5D->GetAxis(iVar)->GetBinWidth(ivarbin+1);
					//! Create Phi projection histogram
					hY5D->GetAxis(iVar)->SetRange(ivarbin+1, ivarbin+1);
					TH1D* hphi = (TH1D*)hY5D->Projection(PHI,"E");
					hphi->Sumw2();
					TString name = TString::Format("hphi_proj_%02d",ivarbin+1);
					TString title = TString::Format("#phi projection for %s = [%.2f,%.2f) | h=%s | top=%s | q2w = %s", 
						                            varTitle[iVarset][iVar].Data(), varbin_lowedge, varbin_highedge,
						                            helTitle[hel].Data(), _topList.Data(), _Q2Wdirname);
					hphi->SetName(name);
					hphi->SetTitle(title);


					//! Obtain Lum and VgFlux normalized Phi distribution
					float vgflux = getvgflux(_Wlow,_Q2low);
      				float factor = 1000000000;
      				float norm = LUM*vgflux*_user.Q2binw*_user.Wbinw*factor;
      				printf("[Wlow:Q2low:Q2binw:Wbinw]= %f:%f:%f:%f\n",_Wlow,_Q2low,_user.Q2binw,_user.Wbinw);
      				printf("norm = %f\n", norm);
      				TH1D* hphinorm = (TH1D*)hphi->Clone(name+"_norm");
      				hphinorm->Scale(1/norm);

					//! Apply Method 2. and obtain R
					/*TH1D* hphi_sinphi = (TH1D*)hphi->Clone(TString::Format("hphi_sinphi%d", ivarbin+1));
					hphi_sinphi->Multiply(hsinphi);*/
					TH1D* hphinorm_sinphi = (TH1D*)hphinorm->Clone(TString::Format("hphinorm_sinphi%d", ivarbin+1));
					hphinorm_sinphi->Multiply(hsinphi);


					Double_t integerr = 0.0;
					//! Note how Underflow and Overflow bins are explicity avoided when calling TH1::IntegralAndError()
					//Double_t integ = hphi_sinphi->IntegralAndError(1, hphi_sinphi->GetNbinsX(), integerr);
					Double_t integ = hphinorm_sinphi->IntegralAndError(1, hphinorm_sinphi->GetNbinsX(), integerr);
					hRvVar->SetBinContent(ivarbin+1, integ/TMath::Pi());
					hRvVar->SetBinError(ivarbin+1, integerr);
				}
			}
			
		}//end nVAR loop
	}//end nVARSET loop
	Info(TString::Format("Proc_hPhi(%s)",helTitle[hel].Data()), "done\n");
}

void ProcYields::Proc_asym()
{
	Info("Proc_asym", "");
	
	TDirectory* dirAsym=NULL;
	if(_dirQ2W->GetDirectory("hAsym") == NULL) dirAsym = _dirQ2W->mkdir("hAsym");

	TDirectory* dirVarset=NULL;
	for (Int_t iVarset=0;iVarset<nVARSET;iVarset++){//!VARSET loop!//
		dirVarset = dirAsym->mkdir(TString::Format("Varset%d", iVarset+1));
		dirVarset->cd();

		THnSparse* hY5D_pos = (THnSparseF*)_dirQ2W->Get(TString::Format("hY5D_POS/Varset%d/hY5D_ACC_CORR",iVarset+1));
		THnSparse* hY5D_neg = (THnSparseF*)_dirQ2W->Get(TString::Format("hY5D_NEG/Varset%d/hY5D_ACC_CORR",iVarset+1));
		THnSparse* hY5D_asym = (THnSparseF*)hY5D_pos->Clone();
		setHistNameTitleAxes_hY5D(hY5D_asym, iVarset, ACC_CORR, ASYM);
		hY5D_asym->Add(hY5D_neg,-1); 
		hY5D_asym->Write();
		TH1D* hAsym = (TH1D*)hY5D_asym->Projection(PHI,"E");
		setHistNameTitleAxes_hAsym(hAsym, iVarset, PHI, ACC_CORR, ASYM);
	}
	Info("Proc_asym", "done\n");
}

void ProcYields::Proc_hYW(){
	Info("Proc_hYW()", "");
		
	TDirectory* dirhYW = _fout->mkdir("hYW");
	TH1F* hYW[nVARSET];
	
	TDirectory* dirVarset=NULL;
	for(Int_t iVarset=0;iVarset<nVARSET;iVarset++){
		Info("Proc_hYW()","Varset = Varset%d", iVarset+1);
		dirVarset=dirhYW->mkdir(TString::Format("Varset%d",iVarset+1));
		dirVarset->cd();
		hYW[iVarset] = new TH1F("hYW","hYW", _user.nWbins, _user.Wmin, _user.Wmax);
		hYW[iVarset]->SetXTitle("W[GeV]");
		
		//!Loop over Q2W dirs, get h5Ds and their yields
		TIter nextkey(_fout->GetListOfKeys());
		TKey *key;
		while (key = (TKey*)nextkey()) {
			TString Q2Wdirname = key->GetName();
			if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
			Info("Proc_hYW()","Q2Wdir = %s", Q2Wdirname.Data());
			TString wrange = Q2Wdirname.Tokenize("_")->At(1)->GetName();
			TString wlow = wrange.Tokenize(",")->At(0)->GetName();
			wlow.Remove(0,1); //remove "["
			//Float_t w = wlow.Atof();
			Double_t w = wlow.Atof();
									
			sprintf(_hname, "%s/hY5D/Varset%d/hY5D_FULL", Q2Wdirname.Data(),iVarset+1);
			THnSparse* hY5D_FULL = (THnSparse*)_fout->Get(_hname);
			if (hY5D_FULL == NULL) cout <<"could not get h5D" << endl;
			//Float_t yield = getIntegral(hY5D_FULL);
			Double_t yield = getIntegral(hY5D_FULL);
			//hYW[iVarset]->Fill(w, yield);
			hYW[iVarset]->SetBinContent(hYW[iVarset]->FindBin(w+_intrinsic.Wbinw), yield);
			Info("Proc_hYW()","W = %f, bin# = %d, yield = %f\n", w, hYW[iVarset]->FindBin(w+_intrinsic.Wbinw), yield);
		}
	}
	Info("Proc_hYW()", "done\n");
}

void ProcYields::Proc_hYW_Dir(){
	Info("Proc_hYW_Dir()", "");
	
	TDirectory* dirhYW_Dir = _fout->mkdir("hYW_Dir");
	TH1D* hYW[nVARSET][nSEQ];
	
	TDirectory* dirVarset;
	for(Int_t iVarset=0;iVarset<nVARSET;iVarset++){
		dirVarset = dirhYW_Dir->mkdir(TString::Format("Varset%d", iVarset+1));
		dirVarset->cd();
		
		//! RECO		
		if(_issim || _isexp){
			//! set bin ranges for axes Q2 & W !//
			Int_t Q2binl =  _hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_user.Q2min  + _intrinsic.Q2binw/2);
			Int_t Q2binh =  _hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_user.Q2max - _intrinsic.Q2binw/2); 
			Int_t Wbinl  =  _hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_user.Wmin  + _intrinsic.Wbinw/2);
			Int_t Wbinh  =  _hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_user.Wmax - _intrinsic.Wbinw/2);
			
			_hY8D[iVarset][RECO]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
			_hY8D[iVarset][RECO]->GetAxis(W)->SetRange(Wbinl, Wbinh); 
			hYW[iVarset][RECO] = (TH1D*)_hY8D[iVarset][RECO]->Projection(2,"E");//0819
			hYW[iVarset][RECO]->Rebin(5);//060313: Since Evegeny Wbinw = 25 Mev and _instrinsic.Wbinw (used in DataAna) is 5MeV
			setHistNameTitleAxes_hYW(hYW[iVarset][RECO], iVarset, RECO);
		}
		//! TH, ACC
		if(_issim){
			//! set bin ranges for axes Q2 & W !//
			Int_t Q2binl =  _hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_user.Q2min  + _intrinsic.Q2binw/2);
			Int_t Q2binh =  _hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_user.Q2max - _intrinsic.Q2binw/2); 
			Int_t Wbinl  =  _hY8D[iVarset][TH]->GetAxis(W)->FindBin(_user.Wmin  + _intrinsic.Wbinw/2);
			Int_t Wbinh  =  _hY8D[iVarset][TH]->GetAxis(W)->FindBin(_user.Wmax - _intrinsic.Wbinw/2);
			
			_hY8D[iVarset][TH]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
			_hY8D[iVarset][TH]->GetAxis(W)->SetRange(Wbinl, Wbinh); 
			hYW[iVarset][TH] = (TH1D*)_hY8D[iVarset][TH]->Projection(2,"E"); //0819
			hYW[iVarset][TH]->Rebin(5);//060313: Since Evegeny Wbinw = 25 Mev and _instrinsic.Wbinw (used in DataAna) is 5MeV
			setHistNameTitleAxes_hYW(hYW[iVarset][TH], iVarset, TH);
								
			//! calculate acceptance !//
			hYW[iVarset][ACC] = (TH1D*)hYW[iVarset][RECO]->Clone();
			hYW[iVarset][ACC]->Divide(hYW[iVarset][TH]);
			setHistNameTitleAxes_hYW(hYW[iVarset][ACC], iVarset, ACC);
		}else if(_isexp){ //if _isexp, then open _fout_sim to get acceptances
			sprintf(_hname, "hYW_Dir/Varset%d/hYW_%s", iVarset+1, seqTitle[ACC].Data()); 
			hYW[iVarset][ACC]=(TH1D*)_fout_sim->Get(_hname);
		}		
		//! ACC_CORR
		if(_issim || _isexp){
			hYW[iVarset][ACC_CORR] = (TH1D*)hYW[iVarset][RECO]->Clone();
			hYW[iVarset][ACC_CORR]->Divide(hYW[iVarset][ACC]);
			/*hYW[iVarset][ACC_CORR] = (TH1D*)hYW[iVarset][RECO]->Clone();
			hYW[iVarset][ACC_CORR]->Divide(hYW[iVarset][ACC]);
			hYW[iVarset][ACC_CORR]->Rebin(5); //060313: Since Evegeny Wbinw = 25 Mev and _instrinsic.Wbinw (used in DataAna) is 5MeV*/
			setHistNameTitleAxes_hYW(hYW[iVarset][ACC_CORR], iVarset, ACC_CORR);
		}
	}
	
	//Reset Q2, W ranges
	for(Int_t iVarset=0;iVarset<nVARSET;iVarset++){
		if(_issim || _isexp){
			_hY8D[iVarset][RECO]->GetAxis(Q2)->SetRange();
			_hY8D[iVarset][RECO]->GetAxis(W)->SetRange(); 
		}
		if(_issim){
			_hY8D[iVarset][TH]->GetAxis(Q2)->SetRange();
			_hY8D[iVarset][TH]->GetAxis(W)->SetRange();
		}
	}
	
	Info("Proc_hYW_Dir()", "done\n");
	
}

Float_t ProcYields::calcNorm(THnSparse* h5D_expACCCORR, THnSparse* h5D_simACCCORR){
	Float_t norm = 0;
	Float_t nExpEvts = 0;
	Float_t nSimEvts = 0;
	
	Int_t expBinCoord[5];
	Int_t nExpBins = h5D_expACCCORR->GetNbins();
	
	for(Int_t iExpBin=0 ; iExpBin<nExpBins ; iExpBin++){
		nExpEvts += h5D_expACCCORR->GetBinContent(iExpBin,expBinCoord);
		
		Int_t iSimBin = h5D_simACCCORR->GetBin(expBinCoord);
		nSimEvts += h5D_simACCCORR->GetBinContent(iSimBin);
	}
	norm = nExpEvts/nSimEvts;
	return norm;
}

void ProcYields::setCanvasStackNameTitle(hel_t hel/*=UNPOL*/){
	char tmp_name[200];
	char tmp_title[200];
	
	for (Int_t iSeq=0;iSeq<nSEQ;iSeq++){
		/* ***hY2D*** */
		sprintf(tmp_name,  "hY2D_%s_%s",  seqTitle[iSeq].Data(), helTitle[hel].Data());
		sprintf(tmp_title, "hY2D_top%s_%s_%s_%s", _Q2Wdirname, _topList.Data(), seqTitle[iSeq].Data(), helTitle[hel].Data());
		_cname_hY2D[iSeq]  = new TString(tmp_name);
		_ctitle_hY2D[iSeq] = new TString(tmp_title);
				
		/* ***hY1D*** */
		sprintf(tmp_name,  "hY1D_%s_%s", seqTitle[iSeq].Data(), helTitle[hel].Data());
		sprintf(tmp_title, "hY1D_top%s_%s_%s_%s", _Q2Wdirname, _topList.Data(), seqTitle[iSeq].Data(), helTitle[hel].Data());
		_cname_hY1D[iSeq]  = new TString(tmp_name);
		_ctitle_hY1D[iSeq] = new TString(tmp_title);

		/* ***hAsym*** */
		sprintf(tmp_name,  "hAsym_%s_%s", seqTitle[iSeq].Data(), helTitle[hel].Data());
		sprintf(tmp_title, "hAsym_top%s_%s_%s_%s", _Q2Wdirname, _topList.Data(), seqTitle[iSeq].Data(), helTitle[hel].Data());
		_cname_hAsym[iSeq]  = new TString(tmp_name);
		_ctitle_hAsym[iSeq] = new TString(tmp_title);
	}
}

void ProcYields::setM1M2axisrange(TH1* h, Int_t iVarset, var_t var /* =, M1/M2 */){
	if (varTitle[iVarset][var] == "M_{p#pi^{+}}" || varTitle[iVarset][var] == "M_{p#pi^{-}}") {
		h->GetXaxis()->SetRangeUser(1.1000, _Whigh-0.14);
	}
	if (varTitle[iVarset][var] == "M_{#pi^{+}#pi^{-}}") {
		h->GetXaxis()->SetRangeUser(0.2780, _Whigh-0.938);
	}
}

Float_t ProcYields::getIntegral(THnSparse* hN){
  Float_t integral=0;

  Int_t nbins = hN->GetNbins();
  //cout << "number of filled bins = " << nbins << endl;
  for (Int_t ibin=0;ibin<nbins;ibin++){
    integral += hN->GetBinContent(ibin);
  }
  //cout << "with underflow and overlow" << integral + hN->GetBinContent(0) + hN->GetBinContent(nbins) << endl;
  return integral;
}

void ProcYields::setHistNameTitleAxes_hYW(TH1D* h, Int_t iVarset, Int_t iSeq){
	sprintf(_hname, "hYW_%s", seqTitle[iSeq].Data());
	sprintf(_htitle, "hYW_%s", seqTitle[iSeq].Data());
		
    h->SetNameTitle(_hname, _htitle);
	h->SetXTitle("W (GeV)");
}

void ProcYields::setHistNameTitleAxes_hY5D(THnSparse* h, Int_t iVarset, seq_t seq, hel_t hel/*=UNPOL*/){
	sprintf(_hname, "hY5D_%s", seqTitle[seq].Data());
	//! Make _title more descriptive.
	if (hel==UNPOL) sprintf(_htitle, "hY5D(M1, M2, %s, %s, %s)",                          varTitle[iVarset][THETA].Data(), varTitle[iVarset][PHI].Data(),varTitle[iVarset][ALPHA].Data());
	else            sprintf(_htitle, "hY5D_%s(M1, M2, %s, %s, %s)", helTitle[hel].Data(), varTitle[iVarset][THETA].Data(), varTitle[iVarset][PHI].Data(),varTitle[iVarset][ALPHA].Data());
	                    
	h->SetNameTitle(_hname, _htitle);
}

void ProcYields::setHistNameTitleAxes_hY2D(TH2D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel/*=UNPOL*/){
	sprintf(_hname, "hY2D_%s_phiV%s", seqTitle[iSeq].Data(), varName[iVar].Data());
	sprintf(_htitle, "%s %s : %s, Varset = %s", seqTitle[iSeq].Data(), (varTitle[iVarset][PHI] + TString("V") + varTitle[iVarset][iVar]).Data(), _Q2Wdirname, asgnmtTitle[iVarset].Data());
	
    h->SetNameTitle(_hname, _htitle);
	h->SetXTitle((varTitle[iVarset][iVar]+varUnitName[iVar]).Data());
	h->SetYTitle((varTitle[iVarset][PHI]+varUnitName[PHI]).Data());
}

void ProcYields::setHistNameTitleAxes_hY1D(TH1D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel/*=UNPOL*/){
	sprintf(_hname, "hY1D_%s_%s",  seqTitle[iSeq].Data(), varName[iVar].Data());
	//make _title more descriptive
	if (hel==UNPOL) sprintf(_htitle, "%s %s : %s",                          seqTitle[iSeq].Data(), varTitle[iVarset][iVar].Data(), _Q2Wdirname);
	else            sprintf(_htitle, "%s %s %s : %s", helTitle[hel].Data(), seqTitle[iSeq].Data(), varTitle[iVarset][iVar].Data(), _Q2Wdirname);	
    h->SetNameTitle(_hname, _htitle);
	h->SetXTitle((varTitle[iVarset][iVar]+varUnitName[iVar]).Data());
}

void ProcYields::setHistNameTitleAxes_hAsym(TH1D* h, Int_t iVarset, Int_t iVar, Int_t iSeq, hel_t hel/*=UNPOL*/){
	sprintf(_hname, "hAsym_%s_%s",  seqTitle[iSeq].Data(), varName[iVar].Data());
	//make _title more descriptive
	if (hel==UNPOL) sprintf(_htitle, "%s %s : %s",                          seqTitle[iSeq].Data(), varTitle[iVarset][iVar].Data(), _Q2Wdirname);
	else            sprintf(_htitle, "%s %s %s : %s", helTitle[hel].Data(), seqTitle[iSeq].Data(), varTitle[iVarset][iVar].Data(), _Q2Wdirname);	
    h->SetNameTitle(_hname, _htitle);
	h->SetXTitle((varTitle[iVarset][iVar]+varUnitName[iVar]).Data());
}
