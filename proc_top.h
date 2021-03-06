#ifndef PROCTOP_H
#define PROCTOP_H

#include "ep_processor.h" // Base class: EpProcessor
#include "cuts.h"
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcTop : public EpProcessor {

public:
	ProcTop(TDirectory *td, DataAna* dataAna, TString h10type);
	~ProcTop();
	void handle(DataH10* dH10);
	void write();
protected:
	void McKin(DataH10 *dH10);
	
	void UpdateEkin(DataH10* dH10, Bool_t useMc = kFALSE);
	
	void setLabFrame4vecHadrons(DataH10* dH10);
	void UpdateMM(Bool_t ismc  = kFALSE);
	void UpdateVarsets(DataH10* dH10, Bool_t ismc = kFALSE);
	
	Float_t getTheta(TLorentzVector lv); //angle in degrees between lv and mLvQCMS
	Float_t getPhi(TLorentzVector lv);   //spherical phi angle in degrees for lv 
    Float_t invTan(Float_t y, Float_t x); //returns angle in radians [0, 2pi]; uses ATan which returns angle in radians [-pi/2, pi/2]
    
    Bool_t mMakeOnlyMCYields; //atrivedi 051013
    TObjArray* mHistsAnaMM;    // for 't1a || t2a || t3a || t4a'. For technical reasons, "hists" was used in place of histsMM
    TObjArray* mYields[nTOP];
	TObjArray* mHistsMM[nTOP];
	TObjArray* mHistseKin[nTOP];
	TObjArray* mYields_Th;
	TObjArray* mHistsMM_Th;
	TObjArray* mHistseKin_Th;
	TLorentzVector mLvQ;
	TLorentzVector mLvW;
	TLorentzVector mLvE;
	TLorentzVector mLvP;
	TLorentzVector mLvPip;
	TLorentzVector mLvPim;
	TLorentzVector mLvMM[nTOP];
	TLorentzVector mLvQCMS;
	TLorentzVector mLvP0CMS;
	TLorentzVector mLvPCMS;
	TLorentzVector mLvPipCMS;
	TLorentzVector mLvPimCMS;
protected:
	static const Int_t NUM_EVTCUTS = 16;
	enum { EVT_NULL,  EVT,       EVT_GPART0,  EVT_GPARTEQ1, EVT_GPARTEQ2, EVT_GPARTEQ3, EVT_GPARTEQ4, EVT_GPART4,
	       EVT_GOODE, EVT_GOODP, EVT_GOODPIP, EVT_GOODPIM,
	       EVT_T1,    EVT_T2,    EVT_T3,      EVT_T4,       EVT_OTHER
	     };
};

ProcTop::ProcTop(TDirectory *td, DataAna* dAna, TString h10type) : EpProcessor(td, dAna, h10type) {
	//atrivedi 051013 mMakeOnlyMCYields if: 
	//if 'is_h10sim' & 'procorder = top'
	mMakeOnlyMCYields = kFALSE;
	if (is_h10sim && EpProcessor::isFirstProc()) mMakeOnlyMCYields = kTRUE;
	
	mHistsAnaMM = NULL; //for techincal reasons, hists was used
	for (int iTop = 0; iTop < nTOP; ++iTop)
	{
		mYields[iTop]   =NULL;
		mHistsMM[iTop]  =NULL;
		mHistseKin[iTop]=NULL;
	}
	mYields_Th   =NULL;
	mHistsMM_Th  =NULL;
	mHistseKin_Th=NULL; 
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
}

ProcTop::~ProcTop() {
	delete hevtsum;
	delete mHistsAnaMM;
	for (int iTop = 0; iTop < nTOP; ++iTop)
    {
    	delete mYields[iTop];
    	delete mHistsMM[iTop];
    	delete mHistseKin[iTop];
    }
    delete mYields_Th;
	delete mHistsMM_Th;
	delete mHistseKin_Th;
}

void ProcTop::handle(DataH10* dH10) {
	//Info("In ProcTop::handle()");
	pass = kFALSE;
	
	const TLorentzVector _4vE0 = dH10->_4vE0;
	const TLorentzVector _4vP0 = dH10->_4vP0;
	
	hevtsum->Fill(EVT);
	
	mLvQ.SetXYZT(0,0,0,0);
	mLvW.SetXYZT(0,0,0,0);
	mLvE.SetXYZT(0,0,0,0);
	mLvP.SetXYZT(0,0,0,0);
	mLvPip.SetXYZT(0,0,0,0);
	mLvPim.SetXYZT(0,0,0,0);
	mLvMM[iTOP1].SetXYZT(0,0,0,0);
	mLvMM[iTOP2].SetXYZT(0,0,0,0);
	mLvMM[iTOP3].SetXYZT(0,0,0,0);
	mLvMM[iTOP4].SetXYZT(0,0,0,0);
	mLvQCMS.SetXYZT(0,0,0,0);
	mLvP0CMS.SetXYZT(0,0,0,0);
	mLvPCMS.SetXYZT(0,0,0,0);
	mLvPipCMS.SetXYZT(0,0,0,0);
	mLvPimCMS.SetXYZT(0,0,0,0);
	
	if (mHistsAnaMM==NULL && !mMakeOnlyMCYields) {
		dirout->cd();
		mHistsAnaMM = dAna->makeHistsMM();
		for (int iTop = 0; iTop < nTOP; ++iTop)
		{
			dirout->mkdir(TString::Format("top%d",iTop+1))->cd();
			mYields[iTop]   =dAna->makeYields();
			mHistsMM[iTop]  =dAna->makeHistsMM();
			mHistseKin[iTop]=dAna->makeHistsEkin();
		}
	}else if (mHistsMM_Th==NULL && mMakeOnlyMCYields){
		dirout->cd();
		dirout->mkdir("mc")->cd();
		mYields_Th   =dAna->makeYields();
		mHistsMM_Th  =dAna->makeHistsMM();
		mHistseKin_Th=dAna->makeHistsEkin();
	}
	if (mMakeOnlyMCYields) {
		McKin(dH10);
		dAna->fillYields(mYields_Th, kTRUE);
		dAna->fillHistsMM(mHistsMM_Th, kTRUE);
		dAna->fillHistsEkin(mHistseKin_Th, kTRUE);
		EpProcessor::handle(dH10); //atrivedi 120113
		return; //atrivedi 051013
	}
	
	
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
	        mLvE.SetPxPyPzE(px,py,pz,energy);
	        mLvQ = _4vE0-mLvE;
	        mLvW = mLvQ+_4vP0;
	        dAna->dTop.Q2 = -1*(mLvQ.Mag2());
	        dAna->dTop.W = mLvW.Mag();
	        
	        /* *** mLvP, mLvPip, mLvPim, mLvMM[iTOP1/2/3/4] *** */
	        setLabFrame4vecHadrons(dH10);
	              
	        /* *** MM *** */   
	        Float_t mm2ppippim = mLvMM[iTOP1].Mag2();
			Float_t mm2ppip    = mLvMM[iTOP2].Mag2();
			Float_t mm2ppim    = mLvMM[iTOP3].Mag2();
			Float_t mm2pippim  = mLvMM[iTOP4].Mag2();
						
			Bool_t t1a = gP && gPip && gPim;// && dH10->gpart == 4; //atrivedi: 071613
			Bool_t t1b = TMath::Abs(mm2ppippim) < 0.0005;
			Bool_t t2a = gP && gPip && !gPim;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t2b = mm2ppip>0 && mm2ppip<0.04;// && dAna->dTop.W < 1.1; //atrivedi tmp
			Bool_t t3a = gP && gPim && !gPip;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t3b = mm2ppim>0 && mm2ppim<0.04;// && dAna->top.W < 1.1; //atrivedi tmp
			Bool_t t4a = gPip && gPim && !gP;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t4b = mm2pippim>0.8 && mm2pippim<1;
			
			Bool_t t1 = t1a && t1b;
			Bool_t t2 = t2a && t2b;
			Bool_t t3 = t3a && t3b;
			Bool_t t4 = t4a && t4b;
			
			if (t1a || t2a || t3a || t4a) { //used to determine top MM cut
				UpdateMM();
				dAna->fillHistsMM(mHistsAnaMM);
								
				UpdateEkin(dH10);
				
			}
			if ( t1 || t2 || t3 || t4) { //final top selection post MM cut
				pass = kTRUE;
				
				if (t1) {
					hevtsum->Fill(EVT_T1);
					dAna->top = 1;
				}else if (t2) {
					hevtsum->Fill(EVT_T2);
					dAna->top = 2;
					mLvPim = mLvMM[iTOP2];
				}else if (t3) {
					hevtsum->Fill(EVT_T3);
					dAna->top = 3;
					mLvPip = mLvMM[iTOP3];
				}else if (t4) {
					hevtsum->Fill(EVT_T4);
					dAna->top = 4;
					mLvP = mLvMM[iTOP4];
				}
				UpdateVarsets(dH10);
				dAna->fillYields(mYields[dAna->top-1]);
				dAna->fillHistsMM(mHistsMM[dAna->top-1]);
				dAna->fillHistsEkin(mHistseKin[dAna->top-1]);
			} else (hevtsum->Fill(EVT_OTHER));
		}
	}
	if (pass) {
		EpProcessor::handle(dH10);
	}
}

void ProcTop::McKin(DataH10 *dH10) {
	const TLorentzVector _4vE0 = dH10->_4vE0;
	const TLorentzVector _4vP0 = dH10->_4vP0;
	
	mLvQ.SetXYZT(0,0,0,0);
	mLvW.SetXYZT(0,0,0,0);
	mLvE.SetXYZT(0,0,0,0);
	mLvP.SetXYZT(0,0,0,0);
	mLvPip.SetXYZT(0,0,0,0);
	mLvPim.SetXYZT(0,0,0,0);
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
			mLvE.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PROTON:
			mLvP.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIP:
			mLvPip.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIM:
			mLvPim.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		default:
			break;
		}
	}

	mLvQ = _4vE0-mLvE;
	mLvW = mLvQ+_4vP0;
	dAna->dTop_mc.Q2 = -1*(mLvQ.Mag2());
	dAna->dTop_mc.W = mLvW.Mag();
	
	UpdateEkin(dH10, kTRUE);
	UpdateMM(kTRUE);
	UpdateVarsets(dH10, kTRUE);
}

void ProcTop::UpdateEkin(DataH10* dH10, Bool_t useMc /*= kFALSE*/) {
	const TLorentzVector _4vE0 = dH10->_4vE0;
	const TLorentzVector _4vP0 = dH10->_4vP0;
	TLorentzVector _4vE1;
		
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
		_4vE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
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
			_4vE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
		}
	}
	
	TLorentzVector _4vQ = _4vE0-_4vE1;
	ekin->W = (_4vQ+_4vP0).Mag();
	ekin->Q2 = -1*_4vQ.Mag2();
	ekin->nu = _4vQ.E();
	ekin->xb = ekin->Q2/(2*MASS_P*ekin->nu);
	ekin->E1 = _4vE1.E();
	ekin->theta1 = _4vE1.Theta()*RadToDeg();
	Double_t phitmp = _4vE1.Phi()*RadToDeg();
	ekin->phi1 = phitmp<-30 ? phitmp+360 : phitmp;
	ekin->theta = _4vQ.Theta()*RadToDeg();
	phitmp = _4vQ.Phi()*RadToDeg();
	ekin->phi = phitmp<-30 ? phitmp+360 : phitmp;
	//if (!useMc) {ekin->eSector = dH10->.sc_sect[dH10->.sc[0]-1];}
}

void ProcTop::setLabFrame4vecHadrons(DataH10* dH10){
	//Track Identified as hadrons
	if (dAna->h10idxP>0) {
		Double_t mom = dH10->p[dAna->h10idxP];
		Double_t px = mom*dH10->cx[dAna->h10idxP];
		Double_t py = mom*dH10->cy[dAna->h10idxP];
		Double_t pz = mom*dH10->cz[dAna->h10idxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		mLvP.SetPxPyPzE(px,py,pz,energy);
	}
	if (dAna->h10idxPip>0) {
		Double_t mom = dH10->p[dAna->h10idxPip];
		Double_t px = mom*dH10->cx[dAna->h10idxPip];
		Double_t py = mom*dH10->cy[dAna->h10idxPip];
		Double_t pz = mom*dH10->cz[dAna->h10idxPip];
		Double_t energy = Sqrt(mom*mom+MASS_PIP*MASS_PIP);
		mLvPip.SetPxPyPzE(px,py,pz,energy);
	}
	if (dAna->h10idxPim>0) {
		Double_t mom = dH10->p[dAna->h10idxPim];
		Double_t px = mom*dH10->cx[dAna->h10idxPim];
		Double_t py = mom*dH10->cy[dAna->h10idxPim];
		Double_t pz = mom*dH10->cz[dAna->h10idxPim];
		Double_t energy = Sqrt(mom*mom+MASS_PIM*MASS_PIM);
		mLvPim.SetPxPyPzE(px,py,pz,energy);
	}
	
	mLvMM[iTOP1] = (mLvW-(mLvP+mLvPip+mLvPim));
	mLvMM[iTOP2] = (mLvW-(mLvP+mLvPip));
    mLvMM[iTOP3] = (mLvW-(mLvP+mLvPim));  
	mLvMM[iTOP4] = (mLvW-(mLvPip+mLvPim)); 
}

void ProcTop::UpdateMM(Bool_t ismc /* = kFALSE */) {
	DataTop *tp = &(dAna->dTop);
	if (ismc) tp = &(dAna->dTop_mc);
	tp->mm2ppippim = mLvMM[iTOP1].Mag2();
	tp->mmppippim  = mLvMM[iTOP1].Mag();
	tp->mm2ppip    = mLvMM[iTOP2].Mag2();
	tp->mmppip     = mLvMM[iTOP2].Mag();
	tp->mm2ppim    = mLvMM[iTOP3].Mag2();
	tp->mmppim     = mLvMM[iTOP3].Mag();
	tp->mm2pippim  = mLvMM[iTOP4].Mag2();
	tp->mmpippim   = mLvMM[iTOP4].Mag();
}

void ProcTop::UpdateVarsets(DataH10* dH10, Bool_t ismc /* = kFALSE */){
	const TLorentzVector _4vE0 = dH10->_4vE0;
	const TLorentzVector _4vP0 = dH10->_4vP0;
	
	DataTop *tp = &(dAna->dTop);
	if (ismc) tp = &(dAna->dTop_mc);

	//Calculate rotation: taken from Evan's phys-ana-omega on 08-05-13
	TVector3 uz = mLvQ.Vect().Unit();
    TVector3 ux = (_4vE0.Vect().Cross(mLvE.Vect())).Unit();
    ux.Rotate(-TMath::Pi()/2,uz);
    TRotation r3;// = new TRotation();
    r3.SetZAxis(uz,ux).Invert();
    //_w and _q are in z-direction
    TVector3 boost(-1*mLvW.BoostVector());
    TLorentzRotation r4(r3); //*_boost);
    r4 *= boost; //*_3rot;
	
	mLvQCMS   = mLvQ;
	mLvP0CMS  = _4vP0;
	mLvPCMS   = mLvP;
	mLvPipCMS = mLvPip;
	mLvPimCMS = mLvPim;
	/*mLvQCMS.Boost(-1*mLvW.BoostVector());
	mLvP0CMS.Boost(-1*mLvW.BoostVector());
	mLvPCMS.Boost(-1*mLvW.BoostVector());
	mLvPipCMS.Boost(-1*mLvW.BoostVector());
	mLvPimCMS.Boost(-1*mLvW.BoostVector());*/
	mLvQCMS.Transform(r4);
	mLvP0CMS.Transform(r4);
	mLvPCMS.Transform(r4);
	mLvPipCMS.Transform(r4);
	mLvPimCMS.Transform(r4);
	
	Float_t Mppip = (mLvPCMS + mLvPipCMS).Mag();
	Float_t Mppim = (mLvPCMS + mLvPimCMS).Mag();
	Float_t Mpippim = (mLvPipCMS + mLvPimCMS).Mag();
	
	tp->varset1.M1 = Mppip;
	tp->varset1.M2 = Mpippim;
	tp->varset1.theta = getTheta(mLvPimCMS);
	tp->varset1.phi = getPhi(mLvPimCMS);
	//tp->varset1.alpha = getAlpha(mLvPimCMS);
	tp->varset1.alpha = 180;
	
	
	tp->varset2.M1 = Mppip;
	tp->varset2.M2 = Mpippim;
	tp->varset2.theta = getTheta(mLvPCMS);
	tp->varset2.phi = getPhi(mLvPCMS);
	//tp->varset2.alpha = getAlpha(mLvPCMS);
	tp->varset2.alpha = 180;
	
	tp->varset3.M1 = Mppip;
	tp->varset3.M2 = Mppim;
	tp->varset3.theta = getTheta(mLvPipCMS);
	tp->varset3.phi = getPhi(mLvPipCMS);
	//tp->varset1.alpha = getAlpha(mLvPipCMS);
	tp->varset3.alpha = 180;
	
	//helicity
	if(!ismc) {tp->h = dH10->evthel;} //e1f
	
}

Float_t ProcTop::getTheta(TLorentzVector lv){
	Float_t retVal = 0;
	/* ***atrivedi 09-13-12*** */
	TVector3 lv_3vec = lv.Vect();
	TVector3 mLvQCMS_3vec = mLvQCMS.Vect();
	//retVal = RadToDeg()*ACos( (lv.Dot(mLvQCMS))/(lv.P()*mLvQCMS.P()) );
	retVal = RadToDeg()*ACos( (lv_3vec.Dot(mLvQCMS_3vec))/(lv_3vec.Mag()*mLvQCMS_3vec.Mag()) );
	/***************************/
	return retVal;
}

Float_t ProcTop::getPhi(TLorentzVector lv) {
	Float_t retVal = 0;
	retVal = RadToDeg()*invTan(lv.Py(), lv.Px());
	return retVal;
}

Float_t ProcTop::invTan(Float_t y, Float_t x){
	Float_t retVal = 0; 
	if (x > 0 && y > 0)  retVal = ATan(y/x);          //1st Quad.
	if (x < 0 && y > 0)  retVal = ATan(y/x) + Pi();   //2nd Quad
	if (x < 0 && y < 0)  retVal = ATan(y/x) + Pi();   //3rd Quad
	if (x > 0 && y < 0)  retVal = ATan(y/x) + 2*Pi(); //4th Quad 
	if (x == 0 && y > 0) retVal = Pi()/2;
	if (x == 0 && y < 0) retVal = 3*Pi()/2; 
	return retVal;  
}


void ProcTop::write(){
	Info("ProcTop::write()", "");
	
	if (mHistsAnaMM != NULL) {
		dirout->cd();
		hevtsum->Write();
		mHistsAnaMM->Write();
		for (int iTop = 0; iTop < nTOP; ++iTop)
		{
			dirout->cd(TString::Format("top%d",iTop+1));
			mYields[iTop]->Write();
			mHistsMM[iTop]->Write();
			mHistseKin[iTop]->Write();
		}
	}
	
	if(mMakeOnlyMCYields){
		dirout->cd("mc");
		mYields_Th->Write();
		mHistsMM_Th->Write();
		mHistseKin_Th->Write();
	}
}

#endif // PROCTOP_H
