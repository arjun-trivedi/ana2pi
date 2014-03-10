#ifndef PROCD2PI_H
#define PROCD2PI_H

#include "ep_processor.h" // Base class: EpProcessor
#include "cuts.h"
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcD2pi : public EpProcessor {

public:
	ProcD2pi(TDirectory *td,DataH10* dataH10,DataAna* dataAna,
			 bool procT, bool procR);
	~ProcD2pi();
	void handle();
	void write();
protected:
	void McKin();
	
	void UpdateEkin(Bool_t useMc = kFALSE);
	
	void ResetLvs();
	void SetLabFrmHadronLvs();
	void UpdateD2pi_MM(Bool_t ismc  = kFALSE);
	void UpdateD2pi(Bool_t ismc = kFALSE);

	void AddBranches(TTree* t, Bool_t ismc=kFALSE);
	
	Float_t getTheta(TLorentzVector lv); //angle in degrees between lv and _lvQCMS
	Float_t getPhi(TLorentzVector lv);   //spherical phi angle in degrees for lv 
    Float_t invTan(Float_t y, Float_t x); //returns angle in radians [0, 2pi]; uses ATan which returns angle in radians [-pi/2, pi/2]
    
    bool _procT; 
    bool _procR;
    TObjArray* _hists_ana_MM;    
    TObjArray* _yields_R[NTOPS];
	TObjArray* _hists_MM_R[NTOPS];
	TObjArray* _hists_ekin_R[NTOPS];
	TObjArray* _yields_T;
	TObjArray* _hists_MM_T;
	TObjArray* _hists_ekin_T;
	TLorentzVector _lvE0_ST;
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
	static const Int_t NUM_EVTCUTS = 16;
	enum { EVT_NULL,  EVT,       EVT_GPART0,  EVT_GPARTEQ1, EVT_GPARTEQ2, EVT_GPARTEQ3, EVT_GPARTEQ4, EVT_GPART4,
	       EVT_GOODE, EVT_GOODP, EVT_GOODPIP, EVT_GOODPIM,
	       EVT_T1,    EVT_T2,    EVT_T3,      EVT_T4,       EVT_OTHER
	     };
};

ProcD2pi::ProcD2pi(TDirectory *td,DataH10* dataH10,DataAna* dataAna,
				   bool procT, bool procR)
				 :EpProcessor(td, dataH10, dataAna) {
	_procT=procT;
	_procR=procR;
	_lvE0_ST.SetPxPyPzE(0,0,5.497,TMath::Sqrt(5.497*5.497+MASS_E*MASS_E));
	_lvE0 = dH10->lvE0;
	_lvP0 = dH10->lvP0;
		
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPART0,"gpart > 0");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ1,"gpart = 1");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ2,"gpart = 2");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ3,"gpart = 3");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPARTEQ4,"gpart = 4");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GPART4,"gpart > 4");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODE,"good e");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODP,"good p");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODPIP,"good pi^{+}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_GOODPIM,"good pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T1,"Type1(p#pi^{+}#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T2,"Type2(p#pi^{+})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T3,"Type3(p#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T4,"Type4(#pi^{+}#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");

	if (_procT) {
		TDirectory* subdir=NULL;
		if (_procR) subdir=dirout->mkdir("T");
		else subdir=dirout;

		subdir->cd();
		_yields_T    =dAna->makeYields();
		_hists_MM_T  =dAna->makeHistsMM();
		_hists_ekin_T=dAna->makeHistsEkin();
		
		_tT = new TTree("tT","Tree containing Thrown data for 2pi events");
		AddBranches(_tT,kTRUE);
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
			_yields_R[iTop]    =dAna->makeYields();
			_hists_MM_R[iTop]  =dAna->makeHistsMM();
			_hists_ekin_R[iTop]=dAna->makeHistsEkin();
		}
		subdir->cd();
		_tR = new TTree("tR","TTree containing Reconstructed data for 2pi events");
		AddBranches(_tR);
	}
}

ProcD2pi::~ProcD2pi() {
	delete hevtsum;
	delete _hists_ana_MM;
	for (int iTop = 0; iTop < NTOPS; ++iTop)
    {
    	delete _yields_R[iTop];
    	delete _hists_MM_R[iTop];
    	delete _hists_ekin_R[iTop];
    }
    delete _yields_T;
	delete _hists_MM_T;
	delete _hists_ekin_T;
}

void ProcD2pi::handle() {
	//Info("In ProcD2pi::handle()");
	pass = kFALSE;
	
	/*const TLorentzVector lvE0 = dH10->lvE0;
	const TLorentzVector lvP0 = dH10->lvP0;*/
	TLorentzVector lvE0=_lvE0;
	TLorentzVector lvP0 = _lvP0;

	ResetLvs();
	
	hevtsum->Fill(EVT);
	
	if (_procT && !_procR){
		McKin();
		dAna->fillYields(_yields_T, kTRUE);
		dAna->fillHistsMM(_hists_MM_T, kTRUE);
		dAna->fillHistsEkin(_hists_ekin_T, kTRUE);
		_tT->Fill();
		EpProcessor::handle(); 
		return;
	}
	/*if(!_procR){
		EpProcessor::handle(); 
		return;
	}*/

	//! ResetLvs() before _procR
	//! (In case of _procT, Lvs will have been set to Thrown values)
	ResetLvs(); 

	if (dH10->gpart>0) {hevtsum->Fill(EVT_GPART0);}
	if (dH10->gpart==1) {hevtsum->Fill(EVT_GPARTEQ1);}
	if (dH10->gpart==2) {hevtsum->Fill(EVT_GPARTEQ2);}
	if (dH10->gpart==3) {hevtsum->Fill(EVT_GPARTEQ3);}
	if (dH10->gpart==4) {hevtsum->Fill(EVT_GPARTEQ4);}
	if (dH10->gpart>4) {hevtsum->Fill(EVT_GPART4);}
	
	Bool_t gE, gP, gPip, gPim;
	gE = gP = gPip = gPim = kFALSE;
	for (Int_t i = 0; i < dH10->gpart; i++) {
		switch (dH10->id[i]) {
		case ELECTRON:
			//if (i == 0){
				dAna->h10idxE = i;
				gE = kTRUE;
				hevtsum->Fill(EVT_GOODE); 
			//}
			break;
		case PROTON:
			//if (i > 0 && dH10->q[i] == 1 && dH10->dc[i] > 0 && dH10->sc[i] > 0) {
				dAna->h10idxP = i;
				gP = kTRUE;
				hevtsum->Fill(EVT_GOODP); 
			//}
			break;
		case PIP:
			//if (i > 0 && dH10->q[i] == 1) {
				dAna->h10idxPip = i;
				gPip = kTRUE;
				hevtsum->Fill(EVT_GOODPIP);
			//}
			break;
		case PIM:
		    //if (i > 0 && dH10->q[i] == -1) {
				dAna->h10idxPim = i;
				gPim = kTRUE;
				hevtsum->Fill(EVT_GOODPIM);
			//}
			break;
		default:
			break;
		}
	}
	if (gE) { 
		if ( (gP && gPip && gPim) || (gP && gPip) || (gP && gPim) || (gPip && gPim) ) {
			/* *** Q2, W *** */
			Double_t mom = dH10->p[dAna->h10idxE];
	        Double_t px = mom*dH10->cx[dAna->h10idxE];
	        Double_t py = mom*dH10->cy[dAna->h10idxE];
	        Double_t pz = mom*dH10->cz[dAna->h10idxE];
	        Double_t energy = Sqrt(mom*mom+MASS_E*MASS_E);
	        _lvE.SetPxPyPzE(px,py,pz,energy);
	        _lvQ = lvE0-_lvE;
	        _lvW = _lvQ+lvP0;
	        /*dAna->d2pi.Q2 = -1*(_lvQ.Mag2());
	        dAna->d2pi.W = _lvW.Mag();*/
	        
	        /* *** _lvP, _lvPip, _lvPim, _lvMM[TOP1/2/3/4] *** */
	        SetLabFrmHadronLvs();
	              
	        /* *** MM *** */   
	        Float_t mm2ppippim = _lvMM[TOP1].Mag2();
			Float_t mm2ppip    = _lvMM[TOP2].Mag2();
			Float_t mm2ppim    = _lvMM[TOP3].Mag2();
			Float_t mm2pippim  = _lvMM[TOP4].Mag2();
						
			Bool_t t1a = gP && gPip && gPim;// && dH10->gpart == 4; //atrivedi: 071613
			Bool_t t1b = TMath::Abs(mm2ppippim) < 0.0005;
			Bool_t t2a = gP && gPip && !gPim;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t2b = mm2ppip>0 && mm2ppip<0.04;// && dAna->d2pi.W < 1.1; //atrivedi tmp
			Bool_t t3a = gP && gPim && !gPip;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t3b = mm2ppim>0 && mm2ppim<0.04;// && dAna->d2pi.top.W < 1.1; //atrivedi tmp
			Bool_t t4a = gPip && gPim && !gP;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t4b = mm2pippim>0.8 && mm2pippim<1;
			
			Bool_t t1 = t1a && t1b;
			Bool_t t2 = t2a && t2b;
			Bool_t t3 = t3a && t3b;
			Bool_t t4 = t4a && t4b;
			
			if (t1a || t2a || t3a || t4a) { //used to determine top MM cut
				UpdateD2pi_MM();
				dAna->fillHistsMM(_hists_ana_MM);
				UpdateEkin();
			}
			if ( t1 || t2 || t3 || t4) { //final top selection post MM cut
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
				dAna->fillYields(_yields_R[dAna->d2pi.top-1]);
				dAna->fillHistsMM(_hists_MM_R[dAna->d2pi.top-1]);
				dAna->fillHistsEkin(_hists_ekin_R[dAna->d2pi.top-1]);
				_tR->Fill();
				if (_procT){
					ResetLvs();
					McKin();
					dAna->fillYields(_yields_T, kTRUE);
					dAna->fillHistsMM(_hists_MM_T, kTRUE);
					dAna->fillHistsEkin(_hists_ekin_T, kTRUE);
					_tT->Fill();
				}
			} else (hevtsum->Fill(EVT_OTHER));
		}
	}
	if (pass) {
		EpProcessor::handle();
	}
}

void ProcD2pi::McKin() {
	/*const TLorentzVector lvE0 = dH10->lvE0;
	const TLorentzVector lvP0 = dH10->lvP0;*/
	TLorentzVector lvE0 = _lvE0_ST;
	TLorentzVector lvP0 = _lvP0;
	
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
	dAna->d2pi_mc.Q2 = -1*(_lvQ.Mag2());
	dAna->d2pi_mc.W = _lvW.Mag();
	
	UpdateEkin(kTRUE);
	UpdateD2pi_MM(kTRUE);
	UpdateD2pi(kTRUE);
}

void ProcD2pi::UpdateEkin(Bool_t useMc /*= kFALSE*/) {
	/*const TLorentzVector lvE0 = dH10->lvE0;
	const TLorentzVector lvP0 = dH10->lvP0;*/
	TLorentzVector lvE0;
	lvE0=_lvE0;
	if (useMc) lvE0 =_lvE0_ST;
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
	
	_lvMM[TOP1] = (_lvW-(_lvP+_lvPip+_lvPim));
	_lvMM[TOP2] = (_lvW-(_lvP+_lvPip));
    _lvMM[TOP3] = (_lvW-(_lvP+_lvPim));  
	_lvMM[TOP4] = (_lvW-(_lvPip+_lvPim)); 
}

void ProcD2pi::UpdateD2pi_MM(Bool_t ismc /* = kFALSE */) {
	Data2pi *tp = &(dAna->d2pi);
	if (ismc) tp = &(dAna->d2pi_mc);
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
	if (ismc) lvE0 = _lvE0_ST;
	TLorentzVector lvP0 = _lvP0;
	
	Data2pi *tp = &(dAna->d2pi);
	if (ismc) tp = &(dAna->d2pi_mc);

	//! Initial Beam Energy
	tp->p_e0=lvE0.P();

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
	tp->phi_e=_lvE.Phi()*RadToDeg();
	tp->phi_p=_lvP.Phi()*RadToDeg();
	tp->phi_pip=_lvPip.Phi()*RadToDeg();
	tp->phi_pim=_lvPim.Phi()*RadToDeg();
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

	//!Q2,W
	tp->Q2 = -1*(_lvQ.Mag2());
	tp->W = _lvW.Mag();

	//Helicity
	if(!ismc) {tp->h = dH10->evthel;} //e1f

	//!{MMs} should be updated already by UpdateD2pi_MM()

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

	tp->theta_cms_p=getTheta(_lvPCMS);
	tp->theta_cms_pip=getTheta(_lvPipCMS);
	tp->theta_cms_pim=getTheta(_lvPimCMS);

	tp->phi_cms_p=getPhi(_lvPCMS);
	tp->phi_cms_pip=getPhi(_lvPipCMS);
	tp->phi_cms_pim=getPhi(_lvPimCMS);

	tp->alpha_1=180;
	tp->alpha_2=180;
	tp->alpha_3=180;

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
}

void ProcD2pi::AddBranches(TTree* t, Bool_t ismc/*=kFALSE*/){
	Data2pi *tp = &(dAna->d2pi);
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

	//! to demo proposal::eid 03-09-13
	t->Branch("p",&dH10->p[0]); //momentum is redundant, but kept for integrity check
	t->Branch("etot",&dH10->etot[0]); 
	t->Branch("ec_eo",&dH10->ec_eo[0]);
	t->Branch("ec_ei",&dH10->ec_ei[0]);
	t->Branch("nphe",&dH10->nphe[0]);
}

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


void ProcD2pi::write(){
	Info("ProcD2pi::write()", "");
	
	if (_hists_ana_MM != NULL) {
		dirout->cd();
		hevtsum->Write();
		_hists_ana_MM->Write();
		for (int iTop = 0; iTop < NTOPS; ++iTop)
		{
			dirout->cd(TString::Format("top%d",iTop+1));
			_yields_R[iTop]->Write();
			_hists_MM_R[iTop]->Write();
			_hists_ekin_R[iTop]->Write();
		}
	}
	
	if(_procT){
		dirout->cd("mc");
		_yields_T->Write();
		_hists_MM_T->Write();
		_hists_ekin_T->Write();
	}
}

#endif // PROCD2PI_H
