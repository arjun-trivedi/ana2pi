#ifndef PROCD2PI_H
#define PROCD2PI_H

#include "ep_processor.h" // Base class: EpProcessor
#include "cuts.h"
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TRandom.h> //atrivedi[06-13-14]: for _rand: tmp to generate random values for alpha
#include <math.h>
#include "proc_eid.h"
#include "proc_pid.h"

using namespace TMath;
using namespace ParticleConstants;

/******************************************************
[10-28-15]
+ Note
	+ For Reconstructed data, ProcD2pi must follow a determined "chain of processors" (set up in scripts)
		+ While there are several reasons that are encapsulated into "chain of processors", one critical reason
		  is that since [10-18-15], ProcD2pi does not do its own EID and PID ; it relies on ProcEid to pass 
		  only events with identified e- and uses h10's gpart indices of identified hadrons by ProcPidNew 
		  (dAna->pidnew.h10IdxP/Pip/Pim) to set up 2pi kinematics. Note that if ProcPid is used, then the code has to be, for now, manually changed
		  to use h10's gpart indices of identified hadrons setup up ProcPid (dAna->pid.h10IdxP/Pip/Pim)
	+ For Thrown data, it is directly called (setup in scripts)
+ For Reconstruction,Pass Events only if:
	+ t1: p,pip,pim are found AND MMppippim satisfied
	OR	
	+ t2: p,pip are found AND pim is not found AND MMppip satisfied 
	OR
	+ t3: p,pim are found AND pip is not found  AND MMppim satisfied
	OR
	+ t4: pip,pim are found AND p is not found AND MMpippim satisfied

[11-16-15]
+ ProcD2pi() now takes in 'study_BG' as a parameter which is used for 'study_BG'. When study_BG=kTRUE,
  then missing mass cuts are not applied and therefore the TTree made contains the entire MM spectrum,
  which can be studied.
*******************************************************/

class ProcD2pi : public EpProcessor {

public:
	ProcD2pi(TDirectory *td,DataH10* dataH10,DataAna* dataAna,
			 bool procT, bool procR, bool make_tree=kFALSE, bool study_BG=kFALSE);
	~ProcD2pi();
	void handle();
	//at-h8 void write();
protected:
	void McKin();
	
	void UpdateEkin(Bool_t useMc = kFALSE);
	
	void ResetLvs();
	void SetLabFrmHadronLvs();
	void UpdateD2pi_Q2_W_MM(Bool_t ismc  = kFALSE);
	void UpdateD2pi(Bool_t ismc = kFALSE);

	//void AddBranches(TTree* t, Bool_t ismc=kFALSE);
	
	Float_t getTheta(TLorentzVector lv); //angle in degrees between lv and _lvQCMS
	Float_t getPhi(TLorentzVector lv);   //spherical phi angle in degrees for lv 
	Float_t getAlpha(TVector3 uv_Gf,TVector3 uv_Gp,TVector3 uv_Bf,TVector3 uv_Bp);
	Float_t invTan(Float_t y, Float_t x); //returns angle in radians [0, 2pi]; uses ATan which returns angle in radians [-pi/2, pi/2]
	
	TRandom* _rand; //atrivedi[06-13-14]: for _rand
	ProcEid* _proc_eid;
	ProcPid* _proc_pid;
	bool _procT,_procR;
	bool _make_tree;
	bool _study_BG;
	TObjArray* _hists_ana_MM;    
	TObjArray** _h8R[NTOPS];
	TObjArray* _hists_MM_R[NTOPS];
	TObjArray* _hists_ekin_R[NTOPS];
	TObjArray** _h8T;
	TObjArray* _hists_MM_T;
	TObjArray* _hists_ekin_T;
	//TLorentzVector _lvE0_ST;
	TLorentzVector _lvE0;
	TLorentzVector _lvP0;
	TLorentzVector _lvQ;
	TLorentzVector _lvW;
	TLorentzVector _lvE;
	TLorentzVector _lvP;
	TLorentzVector _lvPip;
	TLorentzVector _lvPim;
	TLorentzVector _lvMM[NTOPS];
	TLorentzVector _lvQCMS;
	TLorentzVector _lvP0CMS;
	TLorentzVector _lvPCMS;
	TLorentzVector _lvPipCMS;
	TLorentzVector _lvPimCMS;
	TTree* _tT;
	TTree* _tR;
protected:
	static const Int_t NUM_EVTCUTS=14;//6
	enum { EVT_NULL, EVT, 
		   EVT_E,
		   EVT_P_FOUND, EVT_PIP_FOUND, EVT_PIM_FOUND,
		   EVT_PPIPPIM_EX, EVT_PPIP_EX, EVT_PPIM_EX, EVT_PIPPIM_EX,
		   EVT_T1, EVT_T2, EVT_T3, EVT_T4, EVT_OTHER };
};

ProcD2pi::ProcD2pi(TDirectory *td,DataH10* dataH10,DataAna* dataAna,
				   bool procT, bool procR, bool make_tree/*=kFALSE*/,bool study_BG/*=kFALSE*/)
				 :EpProcessor(td, dataH10, dataAna) {
	_rand=new TRandom(); //atrivedi[06-13-14]: for _rand
	_proc_eid=new ProcEid(dataH10,dataAna);
	_proc_pid=new ProcPid(dataH10,dataAna);
	_procT=procT;
	_procR=procR;
	_make_tree=make_tree;
	_study_BG=study_BG;
	_lvE0 = dH10->lvE0;
	_lvP0 = dH10->lvP0;
		
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_E,"e^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_FOUND,"p found");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIP_FOUND,"#pi^{+} found");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIM_FOUND,"#pi^{-} found");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM_EX,"p + #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP_EX,   "p + #pi^{+} + (pi^{-}_{m})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM_EX,   "p + (#pi^{+}_{m}) + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM_EX, "(p_{m})+ #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T1,"T1:p + #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T2,"T2:p + #pi^{+} + (#pi^{-}_{m})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T3,"T3:p + #pi^{+}_{m} + (#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T4,"T4:(p_{m}) + #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");
	hevtsum->SetMinimum(0);

	if (_procT) {
		TDirectory* subdir=NULL;
		if (_procR) subdir=dirout->mkdir("T");
		else subdir=dirout;

		subdir->cd();
		_h8T    =dAna->makeYields();
		_hists_MM_T  =dAna->makeHistsMM();
		_hists_ekin_T=dAna->makeHistsEkin();
				
		if (_make_tree){
			_tT = new TTree("tT","Tree containing Thrown data for 2pi events");
			//AddBranches(_tT,kTRUE);
			dAna->addBranches_Data2pi(_tT,kTRUE);
		}
	}
	if (_procR) {
		TDirectory* subdir=NULL;
		if (_procT) subdir=dirout->mkdir("R");
		else subdir=dirout;

		subdir->cd();
		_hists_ana_MM = dAna->makeHistsMM();
		for (int iTop = 0; iTop < NTOPS; ++iTop)
		{
			subdir->mkdir(TString::Format("top%d",iTop+1))->cd();
			_h8R[iTop]    =dAna->makeYields();
			_hists_MM_R[iTop]  =dAna->makeHistsMM();
			_hists_ekin_R[iTop]=dAna->makeHistsEkin();
		}
		
		subdir->cd();
		if (_make_tree){
			_tR = new TTree("tR","TTree containing Reconstructed data for 2pi events");
			//AddBranches(_tR);
			dAna->addBranches_DataEid(_tR);
			dAna->addBranches_DataEkin(_tR);
			dAna->addBranches_DataPidNew(_tR);
			dAna->addBranches_Data2pi(_tR);
		}
	}
}

ProcD2pi::~ProcD2pi() {
	delete hevtsum;
	delete _hists_ana_MM;
	for (int iTop = 0; iTop < NTOPS; ++iTop)
	{
		delete _h8R[iTop];
		delete _hists_MM_R[iTop];
		delete _hists_ekin_R[iTop];
	}
	delete _h8T;
	delete _hists_MM_T;
	delete _hists_ekin_T;
	delete _rand; //atrivedi[06-13-14]: for _rand
	delete _proc_eid;
	delete _proc_pid;
}

void ProcD2pi::handle() {
	//Info("In ProcD2pi::handle()");
	pass = kFALSE;
	
	/*const TLorentzVector lvE0 = dH10->lvE0;
	const TLorentzVector lvP0 = dH10->lvP0;*/
	TLorentzVector lvE0=_lvE0;
	TLorentzVector lvP0=_lvP0;

	ResetLvs();
	
	hevtsum->Fill(EVT);
	
	if (_procT && !_procR){
		McKin();
		dAna->fillYields(_h8T, dAna->d2pi_mc.W, kTRUE);
		dAna->fillHistsMM(_hists_MM_T, kTRUE);
		dAna->fillHistsEkin(_hists_ekin_T, kTRUE);
		
		if (_make_tree){
			_tT->Fill();
		}
		EpProcessor::handle(); 
		return;
	}
	
	//! ResetLvs() before _procR
	//! (In case of _procT, Lvs will have been set to Thrown values)
	ResetLvs(); 

	//! + Set particles' h10 indices as per ProcEid,ProcPidNew (or ProcPid, in which case manual changes required!)
	//! + hevtsum filling for EID and PID is redundant, but helps in debugging in ProcD2pi correctly followed
	//!	  EID and PID processors

	//! e: Since events have passed ProcEid, h10 index of e=0
	dAna->h10idxE=0;
	//! p,pip,pim
	//! pidnew
	dAna->h10idxP=dAna->pidnew.h10IdxP;
	dAna->h10idxPip=dAna->pidnew.h10IdxPip;
	dAna->h10idxPim=dAna->pidnew.h10IdxPim;
	if (dAna->pidnew.h10IdxP>0)   hevtsum->Fill(EVT_P_FOUND);
	if (dAna->pidnew.h10IdxPip>0) hevtsum->Fill(EVT_PIP_FOUND);
	if (dAna->pidnew.h10IdxPim>0) hevtsum->Fill(EVT_PIM_FOUND);
	//! pidold
	/*dAna->h10idxP=dAna->pid.h10IdxP;
	dAna->h10idxPip=dAna->pid.h10IdxPip;
	dAna->h10idxPim=dAna->pid.h10IdxPim;
	if (dAna->pid.h10IdxP>0) hevtsum->Fill(EVT_P_FOUND);
	if (dAna->pid.h10IdxPip>0) hevtsum->Fill(EVT_PIP_FOUND);
	if (dAna->pid.h10IdxPim>0) hevtsum->Fill(EVT_PIM_FOUND);*/

	//! Determine if identified particles match requirement of tops
	Bool_t t1a=dAna->h10idxP>0 && dAna->h10idxPip>0 && dAna->h10idxPim>0;
	Bool_t t2a=dAna->h10idxP>0 && dAna->h10idxPip>0 && dAna->h10idxPim==0;
	Bool_t t3a=dAna->h10idxP>0 && dAna->h10idxPip==0 && dAna->h10idxPim>0;
	Bool_t t4a=dAna->h10idxP==0 && dAna->h10idxPip>0 && dAna->h10idxPim>0;
	
	hevtsum->Fill(EVT_E);
	if (t1a) hevtsum->Fill(EVT_PPIPPIM_EX);
	if (t2a) hevtsum->Fill(EVT_PPIP_EX);
	if (t3a) hevtsum->Fill(EVT_PPIM_EX);
	if (t4a) hevtsum->Fill(EVT_PIPPIM_EX);
	
	//! Now move on to event selection
	if ( t1a || t2a || t3a || t4a ) { //top:particle selection
		//! Q2,W
		Double_t mom = dH10->p[dAna->h10idxE];
		Double_t px = mom*dH10->cx[dAna->h10idxE];
		Double_t py = mom*dH10->cy[dAna->h10idxE];
		Double_t pz = mom*dH10->cz[dAna->h10idxE];
		Double_t energy = Sqrt(mom*mom+MASS_E*MASS_E);
		_lvE.SetPxPyPzE(px,py,pz,energy);
		_lvQ = lvE0-_lvE;
		_lvW = _lvQ+lvP0;
						
		//! _lvP, _lvPip, _lvPim, _lvMM[TOP1/2/3/4] 
		//! lvs for only particular top is set; rest are 0
		SetLabFrmHadronLvs();
				  
		//! 
		//! MM for only particular top is set; rest are 0  
		Float_t mm2ppippim = _lvMM[TOP1].Mag2();
		Float_t mm2ppip    = _lvMM[TOP2].Mag2();
		Float_t mm2ppim    = _lvMM[TOP3].Mag2();
		Float_t mm2pippim  = _lvMM[TOP4].Mag2();

		//used to study and determine top MM cut
		UpdateD2pi_Q2_W_MM();
		dAna->fillHistsMM(_hists_ana_MM); //method "knows" the particular top's MM hist to fill
		//! [07-17-15] The following should have already been called 
		//! by the previous processors: ProcEid,ProcPidElast 
		//_proc_eid->updateEkin();
		//_proc_eid->updateEid();
		//_proc_pid->updatePid();
						
		//! Make final top selection cut
		Bool_t t1b=TMath::Abs(mm2ppippim) < 0.0005;
		Bool_t t2b=mm2ppip>-0.04 && mm2ppip<0.06;//AT:(0.00,0.04),EI:(-0.04,0.06)
		Bool_t t3b=mm2ppim>0 && mm2ppim<0.04;//0.16;
		Bool_t t4b=mm2pippim>0.8 && mm2pippim<1.00;//mm2pippim>0.0 && mm2pippim<1.25;

		Bool_t t1=t1a && t1b;
		Bool_t t2=t2a && t2b;
		Bool_t t3=t3a && t3b;
		Bool_t t4=t4a && t4b;

		if (_study_BG){//! then over-ride above selection with that based purely on 'a' selection
			t1=t1a;
			t2=t2a;
			t3=t3a;
			t4=t4a;
		}
			
		if ( t1 || t2 || t3 || t4) { //top:particle selection + MM cut
			pass = kTRUE;
				
			if (t1) {
				hevtsum->Fill(EVT_T1);
				dAna->d2pi.top = 1;
			}else if (t2) {
				hevtsum->Fill(EVT_T2);
				dAna->d2pi.top = 2;
				_lvPim = _lvMM[TOP2];
			}else if (t3) {
				hevtsum->Fill(EVT_T3);
				dAna->d2pi.top = 3;
				_lvPip = _lvMM[TOP3];
			}else if (t4) {
				hevtsum->Fill(EVT_T4);
				dAna->d2pi.top = 4;
				_lvP = _lvMM[TOP4];
			}

			UpdateD2pi(); //MM part of d2pi already updated
			dAna->fillYields(_h8R[dAna->d2pi.top-1],dAna->d2pi.W);
			dAna->fillHistsMM(_hists_MM_R[dAna->d2pi.top-1]);
			dAna->fillHistsEkin(_hists_ekin_R[dAna->d2pi.top-1]);
			if (_make_tree){
				_tR->Fill();
			}
			if (_procT){
				ResetLvs();
				McKin();
				dAna->fillYields(_h8T, dAna->d2pi_mc.W, kTRUE);
				dAna->fillHistsMM(_hists_MM_T, kTRUE);
				dAna->fillHistsEkin(_hists_ekin_T, kTRUE);
				if(_make_tree){
					_tT->Fill();
				}
			}
		} else (hevtsum->Fill(EVT_OTHER));
	}

	if (pass) {
		EpProcessor::handle();
	}
}

void ProcD2pi::McKin() {
	const TLorentzVector lvE0 = dH10->lvE0;
	const TLorentzVector lvP0 = dH10->lvP0;
	/*TLorentzVector lvE0 = _lvE0_ST;
	TLorentzVector lvP0 = _lvP0;*/
	
	//printf("num mc = %i\n",dAna->h10.mcnentr);
	for (Int_t idx = 0; idx < dH10->mcnentr; idx++) {
		Int_t _id = dH10->mcid[idx];
		Float_t _p = dH10->mcp[idx];
		Float_t _theta = dH10->mctheta[idx]*DegToRad();
		Float_t _phi = dH10->mcphi[idx]*DegToRad();
		Float_t _mass = dH10->mcm[idx];
		Float_t _energy = Sqrt(_p*_p+_mass*_mass);
		Float_t _pz = _p*Cos(_theta);
		Float_t _py = _p*Sin(_phi)*Sin(_theta);
		Float_t _px = _p*Cos(_phi)*Sin(_theta);
		Int_t _sector = (dH10->mcphi[idx]+30)/60.0+1;
		while (_sector<0) _sector+=6;
		while (_sector>6) _sector-=6;

		Float_t phiD = dH10->mcphi[idx];
		Float_t thetaD = dH10->mctheta[idx];
		TLorentzVector lvtmp;
		TLorentzVector lvtmp2;
		//printf("id = %i\n");
		switch(_id) {
		case ELECTRON:
			_lvE.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PROTON:
			_lvP.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIP:
			_lvPip.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIM:
			_lvPim.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		default:
			break;
		}
	}

	_lvQ = lvE0-_lvE;
	_lvW = _lvQ+lvP0;
	_lvMM[TOP1] = (_lvW-(_lvP+_lvPip+_lvPim));
	_lvMM[TOP2] = (_lvW-(_lvP+_lvPip));
	_lvMM[TOP3] = (_lvW-(_lvP+_lvPim));  
	_lvMM[TOP4] = (_lvW-(_lvPip+_lvPim)); 
	//at-h8
	/*dAna->d2pi_mc.Q2 = -1*(_lvQ.Mag2());
	dAna->d2pi_mc.W = _lvW.Mag();*/
	
	//UpdateEkin(kTRUE);
	_proc_eid->updateEkin(kTRUE);
	UpdateD2pi_Q2_W_MM(kTRUE);
	UpdateD2pi(kTRUE);
}

void ProcD2pi::UpdateEkin(Bool_t useMc /*= kFALSE*/) {
	/*const TLorentzVector lvE0 = dH10->lvE0;
	const TLorentzVector lvP0 = dH10->lvP0;*/
	TLorentzVector lvE0;
	lvE0=_lvE0;
	//if (useMc) lvE0 =_lvE0_ST;
	TLorentzVector lvP0 = _lvP0;
	TLorentzVector lvE1;
		
	DataEkin *ekin = &dAna->eKin;
	if (useMc) ekin = &dAna->eKin_mc;
	
	Float_t mom;
	Float_t px;
	Float_t py;
	Float_t pz;
	if (!useMc){
		mom = dH10->p[0];
		px = mom*dH10->cx[0];
		py = mom*dH10->cy[0];
		pz = mom*dH10->cz[0];
		lvE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
	} else{
		for (Int_t idx = 0; idx < dH10->mcnentr; idx++) {
			Int_t _id = dH10->mcid[idx];
			if (_id != ELECTRON) continue;
			mom = dH10->mcp[idx];
			Float_t _theta = dH10->mctheta[idx]*DegToRad();
			Float_t _phi = dH10->mcphi[idx]*DegToRad();
			px = mom*Cos(_phi)*Sin(_theta);
			py = mom*Sin(_phi)*Sin(_theta);
			pz = mom*Cos(_theta);
			lvE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
		}
	}
	
	TLorentzVector lvQ = lvE0-lvE1;
	ekin->W = (lvQ+lvP0).Mag();
	ekin->Q2 = -1*lvQ.Mag2();
	ekin->nu = lvQ.E();
	ekin->xb = ekin->Q2/(2*MASS_P*ekin->nu);
	ekin->E1 = lvE1.E();
	ekin->theta1 = lvE1.Theta()*RadToDeg();
	Double_t phitmp = lvE1.Phi()*RadToDeg();
	ekin->phi1 = phitmp<-30 ? phitmp+360 : phitmp;
	ekin->theta = lvQ.Theta()*RadToDeg();
	phitmp = lvQ.Phi()*RadToDeg();
	ekin->phi = phitmp<-30 ? phitmp+360 : phitmp;
	//if (!useMc) {ekin->eSector = dH10->.sc_sect[dH10->.sc[0]-1];}
}

void ProcD2pi::ResetLvs(){
	_lvQ.SetXYZT(0,0,0,0);
	_lvW.SetXYZT(0,0,0,0);
	_lvE.SetXYZT(0,0,0,0);
	_lvP.SetXYZT(0,0,0,0);
	_lvPip.SetXYZT(0,0,0,0);
	_lvPim.SetXYZT(0,0,0,0);
	_lvMM[TOP1].SetXYZT(0,0,0,0);
	_lvMM[TOP2].SetXYZT(0,0,0,0);
	_lvMM[TOP3].SetXYZT(0,0,0,0);
	_lvMM[TOP4].SetXYZT(0,0,0,0);
	_lvQCMS.SetXYZT(0,0,0,0);
	_lvP0CMS.SetXYZT(0,0,0,0);
	_lvPCMS.SetXYZT(0,0,0,0);
	_lvPipCMS.SetXYZT(0,0,0,0);
	_lvPimCMS.SetXYZT(0,0,0,0);
}

void ProcD2pi::SetLabFrmHadronLvs(){
	//Track Identified as hadrons
	if (dAna->h10idxP>0) {
		Double_t mom = dH10->p[dAna->h10idxP];
		Double_t px = mom*dH10->cx[dAna->h10idxP];
		Double_t py = mom*dH10->cy[dAna->h10idxP];
		Double_t pz = mom*dH10->cz[dAna->h10idxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		_lvP.SetPxPyPzE(px,py,pz,energy);
	}
	if (dAna->h10idxPip>0) {
		Double_t mom = dH10->p[dAna->h10idxPip];
		Double_t px = mom*dH10->cx[dAna->h10idxPip];
		Double_t py = mom*dH10->cy[dAna->h10idxPip];
		Double_t pz = mom*dH10->cz[dAna->h10idxPip];
		Double_t energy = Sqrt(mom*mom+MASS_PIP*MASS_PIP);
		_lvPip.SetPxPyPzE(px,py,pz,energy);
	}
	if (dAna->h10idxPim>0) {
		Double_t mom = dH10->p[dAna->h10idxPim];
		Double_t px = mom*dH10->cx[dAna->h10idxPim];
		Double_t py = mom*dH10->cy[dAna->h10idxPim];
		Double_t pz = mom*dH10->cz[dAna->h10idxPim];
		Double_t energy = Sqrt(mom*mom+MASS_PIM*MASS_PIM);
		_lvPim.SetPxPyPzE(px,py,pz,energy);
	}
	
	
	if (dAna->h10idxP>0 && dAna->h10idxPip>0 && dAna->h10idxPim>0){
		_lvMM[TOP1] = (_lvW-(_lvP+_lvPip+_lvPim));
	}
	if (dAna->h10idxP>0 && dAna->h10idxPip>0 && dAna->h10idxPim==0){//dAna->h10idxPim==-1
		_lvMM[TOP2] = (_lvW-(_lvP+_lvPip));
	}
	if (dAna->h10idxP>0 && dAna->h10idxPim>0 && dAna->h10idxPip==0){//dAna->h10idxPip==-1
		_lvMM[TOP3] = (_lvW-(_lvP+_lvPim));  
	}
	if (dAna->h10idxPip>0 && dAna->h10idxPim>0 && dAna->h10idxP==0){//dAna->h10idxP==-1
		_lvMM[TOP4] = (_lvW-(_lvPip+_lvPim)); 
	}
}

void ProcD2pi::UpdateD2pi_Q2_W_MM(Bool_t ismc /* = kFALSE */) {
	Data2pi *tp = &(dAna->d2pi);
	if (ismc) tp = &(dAna->d2pi_mc);

	//!Q2,W
	tp->Q2 = -1*(_lvQ.Mag2());
	tp->W = _lvW.Mag();

	//! {MMs}
	tp->mm2ppippim = _lvMM[TOP1].Mag2();
	tp->mmppippim  = _lvMM[TOP1].Mag();
	tp->mm2ppip    = _lvMM[TOP2].Mag2();
	tp->mmppip     = _lvMM[TOP2].Mag();
	tp->mm2ppim    = _lvMM[TOP3].Mag2();
	tp->mmppim     = _lvMM[TOP3].Mag();
	tp->mm2pippim  = _lvMM[TOP4].Mag2();
	tp->mmpippim   = _lvMM[TOP4].Mag();
}

void ProcD2pi::UpdateD2pi(Bool_t ismc /* = kFALSE */){
	// const TLorentzVector lvE0 = dH10->lvE0;
	// const TLorentzVector lvP0 = dH10->lvP0;
	TLorentzVector lvE0;
	lvE0=_lvE0;
	//if (ismc) lvE0 = _lvE0_ST;
	TLorentzVector lvP0 = _lvP0;
	
	Data2pi *tp = &(dAna->d2pi);
	if (ismc) tp = &(dAna->d2pi_mc);

	//! Initial Beam Energy
	tp->p_e0=lvE0.E();
	//! From DC: sector information for detected particles
	if (!ismc){
		tp->sector_e=dH10->sc_sect[dH10->sc[dAna->h10idxE]-1];
		tp->sector_p=dH10->sc_sect[dH10->sc[dAna->h10idxP]-1];
		tp->sector_pip=dH10->sc_sect[dH10->sc[dAna->h10idxPip]-1];
		tp->sector_pim=dH10->sc_sect[dH10->sc[dAna->h10idxPim]-1];
	}
	//! Reconstructed Kinematics 
	//! for e',p',p,pip,pim at e' vertex
	tp->p_e=_lvE.P();
	tp->p_p=_lvP.P();
	tp->p_pip=_lvPip.P();
	tp->p_pim=_lvPim.P();
	tp->theta_e=_lvE.Theta()*RadToDeg();
	tp->theta_p=_lvP.Theta()*RadToDeg();
	tp->theta_pip=_lvPip.Theta()*RadToDeg();
	tp->theta_pim=_lvPim.Theta()*RadToDeg();
	//! [05-06-16] phi angles in lab frame between [-30,330]
	//! + Changed when doing substudies/study_efid/e16 (to correspond with CLAS sectors)
	/*tp->phi_e=getPhi(_lvE);//_lvE.Phi()*RadToDeg();
	tp->phi_p=getPhi(_lvP);//_lvP.Phi()*RadToDeg();
	tp->phi_pip=getPhi(_lvPip);//_lvPip.Phi()*RadToDeg();
	tp->phi_pim=getPhi(_lvPim);//_lvPim.Phi()*RadToDeg();*/
	//! First get values between [-180,180]
	float phi_e=_lvE.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_p=_lvP.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_pip=_lvPip.Phi()*TMath::RadToDeg();// [-180,180]
	float phi_pim=_lvPim.Phi()*TMath::RadToDeg();// [-180,180]
	//! phi: now convert to [-30,330]
	tp->phi_e=phi_e<-30?phi_e+360:phi_e; // [-30,330]
	tp->phi_p=phi_p<-30?phi_p+360:phi_p; // [-30,330]
	tp->phi_pip=phi_pip<-30?phi_pip+360:phi_pip; // [-30,330]
	tp->phi_pim=phi_pim<-30?phi_pim+360:phi_pim; // [-30,330]

	//! Reconstructed e' Vertex
	if (!ismc){
		tp->vx_e=dH10->vx[dAna->h10idxE];
		tp->vx_p=dH10->vx[dAna->h10idxP];
		tp->vx_pip=dH10->vx[dAna->h10idxPip];
		tp->vx_pim=dH10->vx[dAna->h10idxPim];
		tp->vy_e=dH10->vy[dAna->h10idxE];
		tp->vy_p=dH10->vy[dAna->h10idxP];
		tp->vy_pip=dH10->vy[dAna->h10idxPip];
		tp->vy_pim=dH10->vy[dAna->h10idxPim];
		tp->vz_e=dH10->vz[dAna->h10idxE];
		tp->vz_p=dH10->vz[dAna->h10idxP];
		tp->vz_pip=dH10->vz[dAna->h10idxPip];
		tp->vz_pim=dH10->vz[dAna->h10idxPim];
	}else{
		Float_t vx_x_el=dH10->mcvx_x_el;
		Float_t vx_y_el=dH10->mcvx_y_el;
		Float_t vx_z_el=dH10->mcvx_z_el;
		//! In MCTK banks, there is only an electron
		//! vertex defined, unlike in SEB banks, where
		//! each particle has a vertex position.
		//! For all other particles, in MCTK, I assign them to the
		//! electron vertex
		tp->vx_e=vx_x_el;
		tp->vy_e=vx_y_el;
		tp->vz_e=vx_z_el;
		tp->vx_p=vx_x_el;
		tp->vy_p=vx_y_el;
		tp->vz_p=vx_z_el;
		tp->vx_pip=vx_x_el;
		tp->vy_pip=vx_y_el;
		tp->vz_pip=vx_z_el;
		tp->vx_pim=vx_x_el;
		tp->vy_pim=vx_y_el;
		tp->vz_pim=vx_z_el;
		/*for (Int_t idx = 0; idx < dH10->mcnentr; idx++) {
		Int_t _id = dH10->mcid[idx];
		Float_t tmp_vx=dH10->mcvx[idx];
		Float_t tmp_vy=dH10->mcvy[idx];
		Float_t tmp_vz=dH10->mcvz[idx];
		
		switch(_id) {
			case ELECTRON:
				tp->vx_e=tmp_vx;
				tp->vy_e=tmp_vy;
				tp->vz_e=tmp_vz;
				break;
			case PROTON:
				tp->vx_p=tmp_vx;
				tp->vy_p=tmp_vy;
				tp->vz_p=tmp_vz;
				break;
			case PIP:
				tp->vx_pip=tmp_vx;
				tp->vy_pip=tmp_vy;
				tp->vz_pip=tmp_vz;
				break;
			case PIM:
				tp->vx_pim=tmp_vx;
				tp->vy_pim=tmp_vy;
				tp->vz_pim=tmp_vz;
				break;
			default:
				break;
			}
		}*/
	}

	//!Q2,W should be updated already by UpdateD2pi_Q2_W_MM()
	/*tp->Q2 = -1*(_lvQ.Mag2());
	tp->W = _lvW.Mag();*/

	//Helicity
	if(!ismc) {tp->h = dH10->evthel;} //e1f

	//!{MMs} should be updated already by UpdateD2pi_Q2_W_MM()

	//! Varsets
	//! Calculate rotation: taken from Evan's phys-ana-omega on 08-05-13
	TVector3 uz = _lvQ.Vect().Unit();
	TVector3 ux = (lvE0.Vect().Cross(_lvE.Vect())).Unit();
	ux.Rotate(-TMath::Pi()/2,uz);
	TRotation r3;// = new TRotation();
	r3.SetZAxis(uz,ux).Invert();
	//! _w and _q are in z-direction
	TVector3 boost(-1*_lvW.BoostVector());
	TLorentzRotation r4(r3); //*_boost);
	r4 *= boost; //*_3rot;
	
	_lvQCMS   = _lvQ;
	_lvP0CMS  = lvP0;
	_lvPCMS   = _lvP;
	_lvPipCMS = _lvPip;
	_lvPimCMS = _lvPim;
	/*_lvQCMS.Boost(-1*_lvW.BoostVector());
	_lvP0CMS.Boost(-1*_lvW.BoostVector());
	_lvPCMS.Boost(-1*_lvW.BoostVector());
	_lvPipCMS.Boost(-1*_lvW.BoostVector());
	_lvPimCMS.Boost(-1*_lvW.BoostVector());*/
	_lvQCMS.Transform(r4);
	_lvP0CMS.Transform(r4);
	_lvPCMS.Transform(r4);
	_lvPipCMS.Transform(r4);
	_lvPimCMS.Transform(r4);
	
	Float_t Mppip = (_lvPCMS + _lvPipCMS).Mag();
	Float_t Mppim = (_lvPCMS + _lvPimCMS).Mag();
	Float_t Mpippim = (_lvPipCMS + _lvPimCMS).Mag();
	
	tp->M_ppip=Mppip;
	tp->M_ppim=Mppim;
	tp->M_pippim=Mpippim;

	//! Should be directly be able to get the Theta() angle since
	//! virtual photon defines the coordinate system
	tp->theta_cms_p=_lvPCMS.Theta()*RadToDeg();
	tp->theta_cms_pip=_lvPipCMS.Theta()*RadToDeg();
	tp->theta_cms_pim=_lvPimCMS.Theta()*RadToDeg();
	/*tp->theta_cms_p=getTheta(_lvPCMS);
	tp->theta_cms_pip=getTheta(_lvPipCMS);
	tp->theta_cms_pim=getTheta(_lvPimCMS);*/

	tp->phi_cms_p=getPhi(_lvPCMS);
	tp->phi_cms_pip=getPhi(_lvPipCMS);
	tp->phi_cms_pim=getPhi(_lvPimCMS);

	//! alpha angle
	//! Following vectors are used in the calculation of alpha angle; see doc. for getAlpha()
	TVector3 G_f(0,0,0);
	TVector3 G_p(0,0,0);
	TVector3 B_f(0,0,0);
	TVector3 B_p(0,0,0);


	G_f=uz*-1; //! Gamma vector used in calculation of alpha always "follows" -uz

	//! alpha[p',pip][p,pim]
	G_p=_lvPimCMS.Vect().Unit();
	B_f=_lvPipCMS.Vect().Unit();
	B_p=G_p;
	tp->alpha_1=getAlpha(G_f,G_p,B_f,B_p);
	//! alpha[pip,pim][p,p']
	G_p=_lvPCMS.Vect().Unit();
	B_f=_lvPipCMS.Vect().Unit();
	B_p=G_p;
	tp->alpha_2=getAlpha(G_f,G_p,B_f,B_p);
	//! alpha[p',pim][p,pip]
	G_p=_lvPipCMS.Vect().Unit();
	B_f=_lvPimCMS.Vect().Unit();//_lvPCMS.Vect().Unit();
	B_p=G_p;
	tp->alpha_3=getAlpha(G_f,G_p,B_f,B_p);


	/*tp->alpha_1=_rand->Uniform(0,1)*360;//180;
	tp->alpha_2=_rand->Uniform(0,1)*360;//180;
	tp->alpha_3=_rand->Uniform(0,1)*360;//180;*/
	//cout<<"alphas="<<tp->alpha_1<<":"<<tp->alpha_2<<":"<<tp->alpha_3<<endl;

	/*tp->varset1.M1 = Mppip;
	tp->varset1.M2 = Mpippim;
	tp->varset1.theta = getTheta(_lvPimCMS);
	tp->varset1.phi = getPhi(_lvPimCMS);
	//tp->varset1.alpha = getAlpha(_lvPimCMS);
	tp->varset1.alpha = 180;
	
	
	tp->varset2.M1 = Mppip;
	tp->varset2.M2 = Mpippim;
	tp->varset2.theta = getTheta(_lvPCMS);
	tp->varset2.phi = getPhi(_lvPCMS);
	//tp->varset2.alpha = getAlpha(_lvPCMS);
	tp->varset2.alpha = 180;
	
	tp->varset3.M1 = Mppip;
	tp->varset3.M2 = Mppim;
	tp->varset3.theta = getTheta(_lvPipCMS);
	tp->varset3.phi = getPhi(_lvPipCMS);
	//tp->varset1.alpha = getAlpha(_lvPipCMS);
	tp->varset3.alpha = 180;*/

	//! [01-27-16] Added to study scpd
	if(!ismc){
		tp->sc_pd_e=  100*dH10->dc_sect[dH10->dc[dAna->h10idxE]-1]  +dH10->sc_pd[dH10->sc[dAna->h10idxE]-1];
		tp->sc_pd_p=  100*dH10->dc_sect[dH10->dc[dAna->h10idxP]-1]  +dH10->sc_pd[dH10->sc[dAna->h10idxP]-1];
		tp->sc_pd_pip=100*dH10->dc_sect[dH10->dc[dAna->h10idxPip]-1]+dH10->sc_pd[dH10->sc[dAna->h10idxPip]-1];
		tp->sc_pd_pim=100*dH10->dc_sect[dH10->dc[dAna->h10idxPim]-1]+dH10->sc_pd[dH10->sc[dAna->h10idxPim]-1];
	}
}

/*void ProcD2pi::AddBranches(TTree* t, Bool_t ismc/*=kFALSE*///){
	/*Data2pi *tp = &(dAna->d2pi);
	if (ismc) tp = &(dAna->d2pi_mc);

	//! Initial Beam Energy
	t->Branch("p_e0",&tp->p_e0);
	//! Reconstructed Kinematics 
	//! for e',p',p,pip,pim at e' vertex
	t->Branch("p_e",&tp->p_e);
	t->Branch("p_p",&tp->p_p);
	t->Branch("p_pip",&tp->p_pip);
	t->Branch("p_pim",&tp->p_pim);
	t->Branch("theta_e",&tp->theta_e);
	t->Branch("theta_p",&tp->theta_p);
	t->Branch("theta_pip",&tp->theta_pip);
	t->Branch("theta_pim",&tp->theta_pim);
	t->Branch("phi_e",&tp->phi_e);
	t->Branch("phi_p",&tp->phi_p);
	t->Branch("phi_pip",&tp->phi_pip);
	t->Branch("phi_pim",&tp->phi_pim);
	//! Reconstructed e' Vertex
	t->Branch("vx_e",&tp->vx_e);
	t->Branch("vx_p",&tp->vx_p);
	t->Branch("vx_pip",&tp->vx_pip);
	t->Branch("vx_pim",&tp->vx_pim);
	t->Branch("vy_e",&tp->vy_e);
	t->Branch("vy_p",&tp->vy_p);
	t->Branch("vy_pip",&tp->vy_pip);
	t->Branch("vy_pim",&tp->vy_pim);
	t->Branch("vz_e",&tp->vz_e);
	t->Branch("vz_p",&tp->vz_p);
	t->Branch("vz_pip",&tp->vz_pip);
	t->Branch("vz_pim",&tp->vz_pim);
	//! Q2, W
	t->Branch("Q2",&tp->Q2);
	t->Branch("W",&tp->W);
	//! Helicity
	t->Branch("h",&tp->h);
	//! {MMs}
	t->Branch("mm2ppippim",&tp->mm2ppippim);
	t->Branch("mmppippim",&tp->mmppippim);
	t->Branch("mm2ppip",&tp->mm2ppip);
	t->Branch("mmppip",&tp->mmppip);
	t->Branch("mm2ppim",&tp->mm2ppim);
	t->Branch("mmppim",&tp->mmppim);
	t->Branch("mm2pippim",&tp->mm2pippim);
	t->Branch("mmpippim",&tp->mmpippim);
	//! Topology
	t->Branch("top",&tp->top);
	//! Varsets
	t->Branch("M_ppip",&tp->M_ppip);
	t->Branch("M_ppim",&tp->M_ppim);
	t->Branch("M_pippim",&tp->M_pippim);
	t->Branch("theta_cms_p",&tp->theta_cms_p);
	t->Branch("theta_cms_pip",&tp->theta_cms_pip);
	t->Branch("theta_cms_pim",&tp->theta_cms_pim);
	t->Branch("phi_cms_p",&tp->phi_cms_p);
	t->Branch("phi_cms_pip",&tp->phi_cms_pip);
	t->Branch("phi_cms_pim",&tp->phi_cms_pim);
	t->Branch("alpha_1",&tp->alpha_1);
	t->Branch("alpha_2",&tp->alpha_2);
	t->Branch("alpha_3",&tp->alpha_3);
}*/

Float_t ProcD2pi::getTheta(TLorentzVector lv){
	Float_t retVal = 0;
	/* ***atrivedi 09-13-12*** */
	TVector3 lv_3vec = lv.Vect();
	TVector3 mLvQCMS_3vec = _lvQCMS.Vect();
	//retVal = RadToDeg()*ACos( (lv.Dot(_lvQCMS))/(lv.P()*_lvQCMS.P()) );
	retVal = RadToDeg()*ACos( (lv_3vec.Dot(mLvQCMS_3vec))/(lv_3vec.Mag()*mLvQCMS_3vec.Mag()) );
	/***************************/
	return retVal;
}

Float_t ProcD2pi::getPhi(TLorentzVector lv) {
	Float_t retVal = 0;
	retVal = RadToDeg()*invTan(lv.Py(), lv.Px());
	return retVal;
}

//	+Computes the alpha angle between the two hadron planes as described by Iulia in :
//  https://clasweb.jlab.org/wiki/index.php/2%CF%80_kinematics
//
//	+ The input arguments are:
//		+ uv_Gf: unit vector that Gamma "follows"
//		+ uv_Gp: unit vector that Gamma is perpendicular to
//		+ uv_Bf: unit vector that Beta "follows"
//		+ uv_Bp: unit vector that Beta is perpendicular to
// 
//	+ The computations proceeds as:
//		1. v_G=aG*uv_Gf + bG*uv_Gp (look up code/website/hand-written notes to see how aG,bG are computed)
//		2. v_B=aB*uv_Bf + bB*uv_Bp (look up code/website/hand-written notes to see how aB,bB are computed)
//		3. v_D=v_G Cross v_B
//		4. Calculate v_D_angle_uv_Gp. 
//			(Since the angle should either 0 or 180 degrees, matter of precision related to computing is tackled
//			by using 'floor()' when v_D_angle_uv_Gp<90 and 'ceil()' when v_D_angle_uv_Gp>90)
//			+ If v_D is colinear with uv_Gp (i.e. angle=0):
//				alpha=acos(v_G Dot v_B)
//			+ Else If v_D is anti-colinear with uv_Gp (i.e. angle=180):
//				alpha=2*pi - acos(v_G Dot v_B)
//			+ Else:
//				Error! v_D should either be colinear or anticolinear only! alpha is set to=-9999.w .
//		5. return alpha
Float_t ProcD2pi::getAlpha(TVector3 uv_Gf,TVector3 uv_Gp,TVector3 uv_Bf,TVector3 uv_Bp){
	//! step 1.
	Float_t aG=TMath::Sqrt( 1/( 1-TMath::Power(uv_Gf.Dot(uv_Gp),2) ) );
	Float_t bG=-uv_Gf.Dot(uv_Gp)*aG;
	TVector3 v_G=aG*uv_Gf + bG*uv_Gp;
	/*Info("ProcD2pi::getAlpha()","Magnitude of G=%f",TMath::Sqrt(v_G.Dot(v_G)));
	Info("ProcD2pi::getAlpha()","Angle of G with uv_Gp=%f",ACos(v_G.Dot(uv_Gp))*RadToDeg());*/

	//! step 2.
	Float_t aB=TMath::Sqrt( 1/( 1-TMath::Power(uv_Bf.Dot(uv_Bp),2) ) );
	Float_t bB=-uv_Bf.Dot(uv_Bp)*aB;
	TVector3 v_B=aB*uv_Bf + bB*uv_Bp;
	/*Info("ProcD2pi::getAlpha()","Magnitude of B=%f",TMath::Sqrt(v_B.Dot(v_B)));
	Info("ProcD2pi::getAlpha()","Angle of B with uv_Bp=%f",ACos(v_B.Dot(uv_Bp))*RadToDeg());*/

	//! step 3.
	TVector3 v_D=v_G.Cross(v_B);

	//! step 4.
	Float_t alpha=-9999;
	Float_t v_D_angle_uv_Gp=ACos(v_D.Unit().Dot(uv_Gp))*RadToDeg();
	//Info("ProcD2pi::getAlpha()","angle=%f",v_D_angle_uv_Gp);
	(v_D_angle_uv_Gp<90?v_D_angle_uv_Gp=floor(v_D_angle_uv_Gp):v_D_angle_uv_Gp=ceil(v_D_angle_uv_Gp));
	//Info("ProcD2pi::getAlpha()","angle after floor/ceil=%d",int(v_D_angle_uv_Gp));
	if (int(v_D_angle_uv_Gp)==0){
		alpha=ACos(v_G.Dot(v_B))*RadToDeg();
	}else if (int(v_D_angle_uv_Gp)==180){
		alpha=( 2*Pi()-ACos(v_G.Dot(v_B)) )*RadToDeg() ;
	}else{
		Info("ProcD2pi::getAlpha()","Error: v_D is niether colinear nor anticolinear(angle=%f) with uv_Gp! Setting alpha=-9999",v_D_angle_uv_Gp);
	}

	//! step 5.
	return alpha;
}

Float_t ProcD2pi::invTan(Float_t y, Float_t x){
	Float_t retVal = 0; 
	if (x > 0 && y > 0)  retVal = ATan(y/x);          //1st Quad.
	if (x < 0 && y > 0)  retVal = ATan(y/x) + Pi();   //2nd Quad
	if (x < 0 && y < 0)  retVal = ATan(y/x) + Pi();   //3rd Quad
	if (x > 0 && y < 0)  retVal = ATan(y/x) + 2*Pi(); //4th Quad 
	if (x == 0 && y > 0) retVal = Pi()/2;
	if (x == 0 && y < 0) retVal = 3*Pi()/2; 
	return retVal;  
}


//at-h8
/*void ProcD2pi::write(){
	Info("ProcD2pi::write()", "");
	
	if (_hists_ana_MM != NULL) {
		dirout->cd();
		hevtsum->Write();
		_hists_ana_MM->Write();
		for (int iTop = 0; iTop < NTOPS; ++iTop)
		{
			dirout->cd(TString::Format("top%d",iTop+1));
			_h8R[iTop]->Write();
			_hists_MM_R[iTop]->Write();
			_hists_ekin_R[iTop]->Write();
		}
	}
	
	if(_procT){
		dirout->cd("mc");
		_h8T->Write();
		_hists_MM_T->Write();
		_hists_ekin_T->Write();
	}
}*/

#endif // PROCD2PI_H
