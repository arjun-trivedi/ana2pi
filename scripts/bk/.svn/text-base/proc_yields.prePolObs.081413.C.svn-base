#include "proc_yields.h"
#include <TLegend.h>
#include <TString.h>
#include <TKey.h>

#include <iostream>
#include <vector>

using namespace std;

ProcYields::ProcYields(char* fyield, TString topList, 
					   Int_t nQ2bins, Double_t Q2min, Double_t Q2max, 
					   Int_t nWbins, Double_t Wmin, Double_t Wmax){
  
  _fin = NULL;
  _fin_th = NULL;
  _fout = NULL;
  _fout_sim = NULL;
  
  _issim = false;
  _isexp = false;
    
  _dirQ2W = NULL;
  _dirVarset = NULL;
  _dirObs = NULL;
  
  for (Int_t iVarset=0;iVarset<nVarset;iVarset++){
	hY8D[iVarset][RECO] = NULL;
	hY8D[iVarset][TH]   = NULL;
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
  
  Info("ProcYields()", "[nQ2bins, Q2min, Q2max][nWbins, Wmin, Wmax] = [%d, %4.3f, %4.3f][%d, %4.3f, %4.3f]\n",_user.nQ2bins, _user.Q2min, _user.Q2max,
                                                                                                _user.nWbins,  _user.Wmin,  _user.Wmax);
  Info("ProcYields()", "[Q2binw][Wbinw] = [%4.3f][%4.3f]\n", _user.Q2binw, _user.Wbinw);
  
  //! create output file !//
  char foutname[100];
  if(_issim)        {sprintf(foutname, "%s__%d-%4.3f-%4.3f__%d-%4.3f-%4.3f__sim.root",       _topList.Data(), _user.nQ2bins, _user.Q2min, _user.Q2max, _user.nWbins, _user.Wmin, _user.Wmax);}
  else if(_isexp)   {sprintf(foutname, "%s__%d-%4.3f-%4.3f__%d-%4.3f-%4.3f__exp.root",       _topList.Data(), _user.nQ2bins, _user.Q2min, _user.Q2max, _user.nWbins, _user.Wmin, _user.Wmax);}
  else{
	  Info("ProcYields()", "input file is neither mc, sim or exp. Exiting...\n");
	  return; 
   }
  _fout = new TFile(foutname,"RECREATE");
  
  //! if processing experimental data, get corresponding _sim.root file for acceptances !//
  if (_isexp){
	TString tmp(_fout->GetName());
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
    for (Int_t iVarset=0;iVarset<nVarset;iVarset++){ //loop over Varsets
	  sprintf(hreco, "top/top%d/yield_varset%d", top, iVarset+1);
	  sprintf(hth,   "top/mc/yield_varset%d", iVarset+1);
	  if(hY8D[iVarset][RECO]==NULL){
	    hY8D[iVarset][RECO] = (THnSparseF*)_fin->Get(hreco);
	    if(_issim) hY8D[iVarset][TH] = (THnSparseF*)_fin_th->Get(hth);
	  }else{
	    hY8D[iVarset][RECO]->Add((THnSparseF*)_fin->Get(hreco)); 	  
	  }
	}
  }
  _fin->Close();
  
}

ProcYields::~ProcYields(){
  if (_fout_sim != NULL){_fout_sim->Close();}
  if (_fin_th != NULL){_fin_th->Close();}
}

void ProcYields::proc(){
  Info("proc()", "");
     
  //!Before obtaining hY5Ds in Q2W bins, for a given Q2, get differential cross section vs W 'Directly'
  proc_hYW_Dir();
  
  //!Now loop over [Q2,W] bin and proc_5D()
  proc_q2w();
  
  //!hYW
  proc_hYW();
  
  Info("proc()", "Going to Write output file");
  _fout->Write();
  Info("proc()", "output file written");
  Info("proc()", "done\n");
  
  return;
}

void ProcYields::proc_q2w(){
	Info("proc_q2w()", "");
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
      
			//! Loop nVarset & set hy8D bin ranges as per [Q2,W] bin
			for(Int_t iVarset=0;iVarset<nVarset;iVarset++){
				Int_t Q2binl, Q2binh, Wbinl, Wbinh =0;
				if(_issim || _isexp){
					//! set bin ranges for axes Q2 & W !//
					Q2binl =  hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_Q2low  + _intrinsic.Q2binw/2);
					Q2binh =  hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_Q2high - _intrinsic.Q2binw/2); 
					Wbinl  =  hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_Wlow  + _intrinsic.Wbinw/2);
					Wbinh  =  hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_Whigh - _intrinsic.Wbinw/2); 
					printf("[Q2binl,Q2binh]_[Wbinl,Wbinh] = [%d,%d]_[%d,%d]\n",Q2binl, Q2binh, Wbinl, Wbinh);
						
					hY8D[iVarset][RECO]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
					hY8D[iVarset][RECO]->GetAxis(W)->SetRange(Wbinl, Wbinh); 
						
					hrecoQ2vW[iVarset] = (TH2D*)hY8D[iVarset][RECO]->Projection(1,2,"colz");
					sprintf(_hname, "recoQ2vW_Varset%d", iVarset+1);
					hrecoQ2vW[iVarset]->SetNameTitle(_hname, (TString(_hname) + TString(_Q2Wdirname)).Data());
				}
		
				if(_issim){
					//! set bin ranges for axes Q2 & W !//
					Q2binl =  hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_Q2low  + _intrinsic.Q2binw/2);
					Q2binh =  hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_Q2high - _intrinsic.Q2binw/2); 
					Wbinl  =  hY8D[iVarset][TH]->GetAxis(W)->FindBin(_Wlow  + _intrinsic.Wbinw/2);
					Wbinh  =  hY8D[iVarset][TH]->GetAxis(W)->FindBin(_Whigh - _intrinsic.Wbinw/2); 
					printf("[Q2binl,Q2binh]_[Wbinl,Wbinh] = [%d,%d]_[%d,%d]\n",Q2binl, Q2binh, Wbinl, Wbinh);
						
					hY8D[iVarset][TH]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
					hY8D[iVarset][TH]->GetAxis(W)->SetRange(Wbinl, Wbinh);
						
					hthQ2vW[iVarset] = (TH2D*)hY8D[iVarset][TH]->Projection(1,2,"colz");
					sprintf(_hname, "thQ2vW_Varset%d", iVarset+1);
					hthQ2vW[iVarset]->SetNameTitle(_hname, (TString(_hname) + TString(_Q2Wdirname)).Data());
				}
			}
	  
			setCanvasStackNameTitle();
	  
			//! 1. hY5D
			proc_hY5D();
			//! 2. hY5D -> hY2D -> hY1D 
			proc_hY1D();
	  	  	  	  
			_iQ2Wbin += 1;
		}
	}
	Info("proc_q2w()", "Done\n");
}

void ProcYields::proc_hY5D(){
	Info("proc_hY5D()", "");
			
	_dirObs = _dirQ2W->mkdir("hY5D");
	
	for(Int_t iVarset=0;iVarset<nVarset;iVarset++){
		_dirVarset = _dirObs->mkdir(TString::Format("Varset%d", iVarset+1));
		_dirVarset->cd();
		
		//!RECO
		if(_issim || _isexp){
			hY5D[iVarset][RECO] = (THnSparse*)hY8D[iVarset][RECO]->Projection(5,dim,"E");
			setHistNameTitleAxes_hY5D(hY5D[iVarset][RECO], iVarset, RECO);
			hY5D[iVarset][RECO]->Write();
	    }
		//! TH, ACC
		if (_issim){
			hY5D[iVarset][TH] = (THnSparse*)hY8D[iVarset][TH]->Projection(5,dim,"E");
			setHistNameTitleAxes_hY5D(hY5D[iVarset][TH], iVarset, TH);
			hY5D[iVarset][TH]->Write();
						
			//!calculate acceptance
			hY5D[iVarset][ACC] = (THnSparse*)hY5D[iVarset][RECO]->Clone();
			hY5D[iVarset][ACC]->Divide(hY5D[iVarset][TH]);
			setHistNameTitleAxes_hY5D(hY5D[iVarset][ACC], iVarset, ACC);
			hY5D[iVarset][ACC]->Write();
		}else if(_isexp){//if _isexp, then open _fout_sim to get ACC
			sprintf(_hname, "%s/%s/Varset%d/hY5D_%s", _Q2Wdirname, "hY5D", iVarset+1, seqTitle[ACC].Data()); 
			hY5D[iVarset][ACC]=(THnSparse*)_fout_sim->Get(_hname);
			hY5D[iVarset][ACC]->Write();
		}
		//! ACC_CORR
		if(_issim || _isexp){
			hY5D[iVarset][ACC_CORR] = (THnSparse*)hY5D[iVarset][RECO]->Clone();
			hY5D[iVarset][ACC_CORR]->Divide(hY5D[iVarset][ACC]);
			setHistNameTitleAxes_hY5D(hY5D[iVarset][ACC_CORR], iVarset, ACC_CORR);
			hY5D[iVarset][ACC_CORR]->Write();
		}
		//! HOLE
		if(_issim){
			hY5D[iVarset][HOLE] = (THnSparse*)hY5D[iVarset][TH]->Clone();
			hY5D[iVarset][HOLE]->Add(hY5D[iVarset][ACC_CORR], -1);
			setHistNameTitleAxes_hY5D(hY5D[iVarset][HOLE], iVarset, HOLE);
			hY5D[iVarset][HOLE]->Write();
		}else if (_isexp){//if _isexp, then open _fout_sim to get sim-HOLE, sim-ACC_CORR 
			//!first copy sim-HOLE into exp-HOLE
			sprintf(_hname, "%s/%s/Varset%d/hY5D_%s", _Q2Wdirname, "hY5D", iVarset+1, seqTitle[HOLE].Data());
			hY5D[iVarset][HOLE]=(THnSparse*)_fout_sim->Get(_hname);
			
			//!now get sim-ACC_CORR for obtaining 
			sprintf(_hname, "%s/%s/Varset%d/hY5D_%s", _Q2Wdirname, "hY5D", iVarset+1, seqTitle[ACC_CORR].Data());
			THnSparse* hY5D_sim_ACC_CORR = (THnSparse*)_fout_sim->Get(_hname);
			
			//!calculate normalization factor
			//Float_t norm = hY5D[iVarset][ACC_CORR]->GetSumw()/hY5D_sim_ACC_CORR->GetSumw();
			Float_t norm = calcNorm(hY5D[iVarset][ACC_CORR], hY5D_sim_ACC_CORR);
			cout << "norm = " << norm << endl;
			
			//!normalize exp-HOLE
			hY5D[iVarset][HOLE]->Scale(norm);
			
			hY5D[iVarset][HOLE]->Write();
		}
		//! FULL
		if (_issim || _isexp){
			hY5D[iVarset][FULL] = (THnSparse*)hY5D[iVarset][ACC_CORR]->Clone();
			hY5D[iVarset][FULL]->Add(hY5D[iVarset][HOLE],1 );
			setHistNameTitleAxes_hY5D(hY5D[iVarset][FULL], iVarset, FULL);
			hY5D[iVarset][FULL]->Write();
		}
	}
	
	Info("proc_hY5D()", "done\n");
}

void ProcYields::proc_hY1D(){
	Info("proc_hY1D()", "");
	
	//! Create Canvas Managers !//
	cm_hY2D = new CanvasManager(_c_ww, _c_wh);
    cm_hY1D = new CanvasManager(_c_ww, _c_wh);
	//! Create Stacks for TH, RECO, ACC_CORR, FULL !//
	for(Int_t iSeq=0;iSeq<nSeq;iSeq++){
		seq_t seq = static_cast<seq_t>(iSeq);
		if (seq==ACC || seq==HOLE) continue;
		if (_isexp && seq==TH) continue;
		hs_hY2D[iSeq] = new THStack(*cname_hY2D[iSeq], *ctitle_hY2D[iSeq]);
		hs_hY1D[iSeq] = new THStack(*cname_hY1D[iSeq], *ctitle_hY1D[iSeq]);	
	}
		
	//! hY5D -> hY2D -> hY1D !//
	if(_dirQ2W->GetDirectory("hY2D") == NULL) _dirQ2W->mkdir("hY2D");
	if(_dirQ2W->GetDirectory("hY1D") == NULL) _dirQ2W->mkdir("hY1D");
	for (Int_t iVarset=0;iVarset<nVarset;iVarset++){
		_dirQ2W->cd("hY2D");
		if(gDirectory->GetDirectory(TString::Format("Varset%d", iVarset+1)) == NULL) gDirectory->mkdir(TString::Format("Varset%d", iVarset+1));
		_dirQ2W->cd("hY1D");
		if(gDirectory->GetDirectory(TString::Format("Varset%d", iVarset+1)) == NULL) gDirectory->mkdir(TString::Format("Varset%d", iVarset+1));
		for(Int_t iVar=0; iVar<nVar; iVar++){
			var_t var = static_cast<var_t>(iVar);
			if (var == ALPHA) continue; 
			for(Int_t iSeq=0;iSeq<nSeq;iSeq++){
				seq_t seq = static_cast<seq_t>(iSeq);
				if (seq==ACC || seq==HOLE) continue;
				if (_isexp && seq==TH) continue;
				
				TDirectory* dir;
				//! hY2D !//
				_dirQ2W->cd(TString::Format("hY2D/Varset%d", iVarset+1));
				hY2D[iVarset][iVar][iSeq] = (TH2D*)hY5D[iVarset][iSeq]->Projection(PHI, iVar, "");
				setHistNameTitleAxes_hY2D(hY2D[iVarset][iVar][iSeq], iVarset, iVar, iSeq);
						
				//! hY1D !//
				_dirQ2W->cd(TString::Format("hY1D/Varset%d", iVarset+1));
				hY1D[iVarset][iVar][iSeq] = (TH1D*)hY2D[iVarset][iVar][iSeq]->ProjectionX()->Clone();
				hY1D[iVarset][iVar][iSeq]->SetMinimum(0);
				setHistNameTitleAxes_hY1D(hY1D[iVarset][iVar][iSeq], iVarset, iVar, iSeq);
            			
				if (var == M1 or var == M2){
					setM1M2axisrange(hY2D[iVarset][iVar][iSeq], iVarset, var);
					setM1M2axisrange(hY1D[iVarset][iVar][iSeq], iVarset, var);
				}
			
				//! add to stack !//
				hs_hY2D[iSeq]->Add(hY2D[iVarset][iVar][iSeq], "colz");
				hs_hY1D[iSeq]->Add(hY1D[iVarset][iVar][iSeq], "e1");
        
			}//end nSeq loop
		}//end nVar loop
	}//end nVarset loop
	
	if(_issim || _isexp){
		cm_hY2D->add(hs_hY2D[RECO]);
		cm_hY1D->add(hs_hY1D[RECO]);
	}
	if (_issim){
		cm_hY2D->add(hs_hY2D[TH]);
		cm_hY1D->add(hs_hY1D[TH]);
	}
	if(_issim || _isexp){
		cm_hY2D->add(hs_hY2D[ACC_CORR]);
		cm_hY1D->add(hs_hY1D[ACC_CORR]);
	}
	if(_issim || _isexp){
		cm_hY2D->add(hs_hY2D[FULL]);
		cm_hY1D->add(hs_hY1D[FULL]);
	}
	
	//! Save Canvases to file !//
	_dirQ2W->cd();
	cm_hY2D->draw();
	cm_hY2D->save();
	cm_hY2D->close();
	
	cm_hY1D->draw();
	cm_hY1D->save();
	cm_hY1D->close();
	
	Info("proc_hY1D()", "done\n");
}

void ProcYields::proc_hYW(){
	Info("proc_hYW()", "");
		
	_dirObs = _fout->mkdir("hYW");
	TH1F* hYW[nVarset];
	
	for(Int_t iVarset=0;iVarset<nVarset;iVarset++){
		Info("proc_hYW()","Varset = Varset%d", iVarset+1);
		_dirVarset=_dirObs->mkdir(TString::Format("Varset%d",iVarset+1));
		_dirVarset->cd();
		hYW[iVarset] = new TH1F("hYW","hYW", _user.nWbins, _user.Wmin, _user.Wmax);
		hYW[iVarset]->SetXTitle("W[GeV]");
		
		//!Loop over Q2W dirs, get h5Ds and their yields
		TIter nextkey(_fout->GetListOfKeys());
		TKey *key;
		while (key = (TKey*)nextkey()) {
			TString Q2Wdirname = key->GetName();
			if(Q2Wdirname.EqualTo("hYW_Dir") || Q2Wdirname.EqualTo("hYW"))continue;
			Info("proc_hYW()","Q2Wdir = %s", Q2Wdirname.Data());
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
			Info("proc_hYW()","W = %f, bin# = %d, yield = %f\n", w, hYW[iVarset]->FindBin(w+_intrinsic.Wbinw), yield);
		}
	}
	Info("proc_hYW()", "done\n");
}

void ProcYields::proc_hYW_Dir(){
	Info("proc_hYW_Dir()", "");
	
	_dirObs = _fout->mkdir("hYW_Dir");
	TH1D* hYW[nVarset][nSeq];
	
	for(Int_t iVarset=0;iVarset<nVarset;iVarset++){
		_dirVarset = _dirObs->mkdir(TString::Format("Varset%d", iVarset+1));
		_dirVarset->cd();
		
		//! RECO		
		if(_issim || _isexp){
			//! set bin ranges for axes Q2 & W !//
			Int_t Q2binl =  hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_user.Q2min  + _intrinsic.Q2binw/2);
			Int_t Q2binh =  hY8D[iVarset][RECO]->GetAxis(Q2)->FindBin(_user.Q2max - _intrinsic.Q2binw/2); 
			Int_t Wbinl  =  hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_user.Wmin  + _intrinsic.Wbinw/2);
			Int_t Wbinh  =  hY8D[iVarset][RECO]->GetAxis(W)->FindBin(_user.Wmax - _intrinsic.Wbinw/2);
			
			hY8D[iVarset][RECO]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
			hY8D[iVarset][RECO]->GetAxis(W)->SetRange(Wbinl, Wbinh); 
			hYW[iVarset][RECO] = (TH1D*)hY8D[iVarset][RECO]->Projection(2);
			hYW[iVarset][RECO]->Rebin(5);//060313: Since Evegeny Wbinw = 25 Mev and _instrinsic.Wbinw (used in DataAna) is 5MeV
			setHistNameTitleAxes_hYW(hYW[iVarset][RECO], iVarset, RECO);
		}
		//! TH, ACC
		if(_issim){
			//! set bin ranges for axes Q2 & W !//
			Int_t Q2binl =  hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_user.Q2min  + _intrinsic.Q2binw/2);
			Int_t Q2binh =  hY8D[iVarset][TH]->GetAxis(Q2)->FindBin(_user.Q2max - _intrinsic.Q2binw/2); 
			Int_t Wbinl  =  hY8D[iVarset][TH]->GetAxis(W)->FindBin(_user.Wmin  + _intrinsic.Wbinw/2);
			Int_t Wbinh  =  hY8D[iVarset][TH]->GetAxis(W)->FindBin(_user.Wmax - _intrinsic.Wbinw/2);
			
			hY8D[iVarset][TH]->GetAxis(Q2)->SetRange(Q2binl, Q2binh);
			hY8D[iVarset][TH]->GetAxis(W)->SetRange(Wbinl, Wbinh); 
			hYW[iVarset][TH] = (TH1D*)hY8D[iVarset][TH]->Projection(2);
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
	for(Int_t iVarset=0;iVarset<nVarset;iVarset++){
		if(_issim || _isexp){
			hY8D[iVarset][RECO]->GetAxis(Q2)->SetRange();
			hY8D[iVarset][RECO]->GetAxis(W)->SetRange(); 
		}
		if(_issim){
			hY8D[iVarset][TH]->GetAxis(Q2)->SetRange();
			hY8D[iVarset][TH]->GetAxis(W)->SetRange();
		}
	}
	
	Info("proc_hYW_Dir()", "done\n");
	
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

void ProcYields::setCanvasStackNameTitle(){
	char tmp_name[200];
	char tmp_title[200];
	
	for (Int_t iSeq=0;iSeq<nSeq;iSeq++){
		/* ***hY2D*** */
		sprintf(tmp_name,  "hY2Dab_%s",  seqTitle[iSeq].Data());
		sprintf(tmp_title, "hY2Dab_top%s_%s_%s",   _topList.Data(), seqTitle[iSeq].Data(), _Q2Wdirname);
		cname_hY2D[iSeq]  = new TString(tmp_name);
		ctitle_hY2D[iSeq] = new TString(tmp_title);
				
		/* ***hY1D*** */
		sprintf(tmp_name,  "hY1Da_%s", seqTitle[iSeq].Data());
		sprintf(tmp_title, "hY1Da_top%s_%s_%s", _topList.Data(), seqTitle[iSeq].Data(), _Q2Wdirname);
		cname_hY1D[iSeq]  = new TString(tmp_name);
		ctitle_hY1D[iSeq] = new TString(tmp_title);
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

void ProcYields::setHistNameTitleAxes_hY5D(THnSparse* h, Int_t iVarset, seq_t seq){
	sprintf(_hname, "hY5D_%s", seqTitle[seq].Data());
	sprintf(_htitle, "M1, M2, %s, %s, %s", varTitle[iVarset][THETA].Data(), varTitle[iVarset][PHI].Data(),varTitle[iVarset][ALPHA].Data());
	                    
	h->SetNameTitle(_hname, _htitle);
}

void ProcYields::setHistNameTitleAxes_hY2D(TH2D* h, Int_t iVarset, Int_t iVar, Int_t iSeq){
	sprintf(_hname, "hY2D_%s_phiV%s", seqTitle[iSeq].Data(), varName[iVar].Data());
	sprintf(_htitle, "%s %s : %s, Varset = %s", seqTitle[iSeq].Data(), (varTitle[iVarset][PHI] + TString("V") + varTitle[iVarset][iVar]).Data(), _Q2Wdirname, asgnmtTitle[iVarset].Data());
	
    h->SetNameTitle(_hname, _htitle);
	h->SetXTitle((varTitle[iVarset][iVar]+varUnitName[iVar]).Data());
	h->SetYTitle((varTitle[iVarset][PHI]+varUnitName[PHI]).Data());
}

void ProcYields::setHistNameTitleAxes_hY1D(TH1D* h, Int_t iVarset, Int_t iVar, Int_t iSeq){
	sprintf(_hname, "hY1D_%s_%s", seqTitle[iSeq].Data(), varName[iVar].Data());
	sprintf(_htitle, "%s %s : %s", seqTitle[iSeq].Data(), varTitle[iVarset][iVar].Data(), _Q2Wdirname);
		
    h->SetNameTitle(_hname, _htitle);
	h->SetXTitle((varTitle[iVarset][iVar]+varUnitName[iVar]).Data());
}
