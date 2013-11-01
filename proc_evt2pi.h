#ifndef PROCTOP_H
#define PROCTOP_H

#include "ep_processor.h" // Base class: EpProcessor
#include "cuts.h"
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcTop : public EpProcessor {

public:
	ProcTop(TDirectory *td, DataAna* dataAna, TString h10type);
	~ProcTop();
	void handle(DataH10* dH10);
	void write();
protected:
	void mcKin(DataH10 *dH10);
	
	void updateEkin(DataH10* dH10, Bool_t useMc = kFALSE);
	
	void set4vecHadrons(DataH10* dH10);
	void updateTop_MM(Bool_t ismc  = kFALSE);
	void updateTop_VarSets(DataH10* dH10, Bool_t ismc = kFALSE);
	
	Float_t getTheta(TLorentzVector lv); //angle in degrees between lv and lvQCMS
	Float_t getPhi(TLorentzVector lv);   //spherical phi angle in degrees for lv 
    Float_t invTan(Float_t y, Float_t x); //returns angle in radians [0, 2pi]; uses ATan which returns angle in radians [-pi/2, pi/2]
    
    Bool_t makeOnlyMCYields; //atrivedi 051013
    TObjArray *histsMM;    // for 't1a || t2a || t3a || t4a'. For technical reasons, "hists" was used in place of histsMM
    TObjArray *histseKin;  // for 't1a || t2a || t3a || t4a'
    TObjArray *yields_top0;// for 't1 || t2 || t3 || t4'
    TObjArray *yields_top1;
	TObjArray *yields_top2;
	TObjArray *yields_top3;
	TObjArray *yields_top4;
	TObjArray *histsMM_top0;// for 't1 || t2 || t3 || t4'
	TObjArray *histsMM_top1;
	TObjArray *histsMM_top2;
	TObjArray *histsMM_top3;
	TObjArray *histsMM_top4;
	TObjArray *histseKin_top0;// for 't1 || t2 || t3 || t4'
	TObjArray *histseKin_top1;
	TObjArray *histseKin_top2;
	TObjArray *histseKin_top3;
	TObjArray *histseKin_top4; 
	TObjArray *yields_mc;
	TObjArray *histsMM_mc;
	TObjArray *histseKin_mc;
	TLorentzVector lvQ;
	TLorentzVector lvW;
	TLorentzVector lvE;
	TLorentzVector lvP;
	TLorentzVector lvPip;
	TLorentzVector lvPim;
	TLorentzVector lvMMtop1;
	TLorentzVector lvMMtop2;
	TLorentzVector lvMMtop3;
	TLorentzVector lvMMtop4;
	TLorentzVector lvQCMS;
	TLorentzVector lvP0CMS;
	TLorentzVector lvPCMS;
	TLorentzVector lvPipCMS;
	TLorentzVector lvPimCMS;
protected:
	//TH1F *hevtsum;
	static const Int_t NUM_EVTCUTS = 17;
	enum { EVT_NULL, EVT, EVT_GPART0, EVT_GPARTEQ1, EVT_GPARTEQ2, EVT_GPARTEQ3, EVT_GPARTEQ4, EVT_GPART4,
	       EVT_GOODE, EVT_GOODP, EVT_GOODPIP, EVT_GOODPIM,
	       EVT_T0, EVT_T1, EVT_T2, EVT_T3, EVT_T4, EVT_OTHER
	     };
};

ProcTop::ProcTop(TDirectory *td, DataAna* dAna, TString h10type) : EpProcessor(td, dAna, h10type) {
	//atrivedi 051013 makeOnlyMCYields if: 
	//if 'is_h10sim' & 'procorder = top'
	makeOnlyMCYields = kFALSE;
	if (is_h10sim && EpProcessor::isFirstProc()) makeOnlyMCYields = kTRUE;
	
	histsMM = NULL; //for techincal reasons, hists was used
	histseKin = NULL;
	yields_top0 = NULL;
	yields_top1 = NULL;
	yields_top2 = NULL;
	yields_top3 = NULL;
	yields_top4 = NULL;
	histsMM_top0 = NULL;
	histsMM_top1 = NULL;
	histsMM_top2 = NULL;
	histsMM_top3 = NULL;
	histsMM_top4 = NULL;
	histseKin_top0 = NULL;
	histseKin_top1 = NULL;
	histseKin_top2 = NULL;
	histseKin_top3 = NULL;
	histseKin_top4 = NULL;  
	yields_mc = NULL;
	histsMM_mc = NULL;
	histseKin_mc = NULL; 
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
	hevtsum->GetXaxis()->SetBinLabel(EVT_T0,"Type0"); //TYPEO = PPIPPIM||PPIP||PPIM||PIPPIM
	hevtsum->GetXaxis()->SetBinLabel(EVT_T1,"Type1(p#pi^{+}#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T2,"Type2(p#pi^{+})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T3,"Type3(p#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_T4,"Type4(#pi^{+}#pi^{-})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");
	//hevtsum->SetDirectory(dirout);
}

ProcTop::~ProcTop() {
	delete hevtsum;
	delete histsMM;
	delete histseKin;
    delete yields_top0;
    delete yields_top1;
	delete yields_top2;
	delete yields_top3;
	delete yields_top4;
	delete histsMM_top0;
	delete histsMM_top1;
	delete histsMM_top2;
	delete histsMM_top3;
	delete histsMM_top4;
	delete histseKin_top0;
	delete histseKin_top1;
	delete histseKin_top2;
	delete histseKin_top3;
	delete histseKin_top4; 
	delete yields_mc;
	delete histsMM_mc;
	delete histseKin_mc;
}

void ProcTop::handle(DataH10* dH10) {
	//Info("In ProcTop::handle()");
	pass = kFALSE;
	
	const TLorentzVector _4vE0 = dH10->_4vE0;
	const TLorentzVector _4vP0 = dH10->_4vP0;
	
	hevtsum->Fill(EVT);
	
	lvQ.SetXYZT(0,0,0,0);
	lvW.SetXYZT(0,0,0,0);
	lvE.SetXYZT(0,0,0,0);
	lvP.SetXYZT(0,0,0,0);
	lvPip.SetXYZT(0,0,0,0);
	lvPim.SetXYZT(0,0,0,0);
	lvMMtop1.SetXYZT(0,0,0,0);
	lvMMtop2.SetXYZT(0,0,0,0);
	lvMMtop3.SetXYZT(0,0,0,0);
	lvMMtop4.SetXYZT(0,0,0,0);
	if (histsMM==NULL && !makeOnlyMCYields) {
		dirout->cd();
		histsMM = dAna->makeHistsMM();
		histseKin = dAna->makeHistsEkin();
		dirout->mkdir("top0")->cd();
		yields_top0 = dAna->makeYields();
		histsMM_top0  = dAna->makeHistsMM();
		histseKin_top0 = dAna->makeHistsEkin();
		dirout->mkdir("top1")->cd();
		yields_top1 = dAna->makeYields();
		histsMM_top1  = dAna->makeHistsMM();
		histseKin_top1 = dAna->makeHistsEkin();
		dirout->mkdir("top2")->cd();
		yields_top2 = dAna->makeYields();
		histsMM_top2  = dAna->makeHistsMM();
		histseKin_top2 = dAna->makeHistsEkin();
		dirout->mkdir("top3")->cd();
		yields_top3 = dAna->makeYields();
		histsMM_top3  = dAna->makeHistsMM();
		histseKin_top3 = dAna->makeHistsEkin();
		dirout->mkdir("top4")->cd();
		yields_top4 = dAna->makeYields();
		histsMM_top4  = dAna->makeHistsMM();
		histseKin_top4 = dAna->makeHistsEkin();
	}else if (histsMM_mc==NULL && makeOnlyMCYields){
		dirout->cd();
		dirout->mkdir("mc")->cd();
		yields_mc = dAna->makeYields();
		histsMM_mc  = dAna->makeHistsMM();
		histseKin_mc = dAna->makeHistsEkin();
	}
	
	//atrivedi 051013 
	if (makeOnlyMCYields) {
		mcKin(dH10);
		dAna->fillYields(yields_mc, kTRUE);
		dAna->fillHistsMM(histsMM_mc, kTRUE);
		dAna->fillHistsEkin(histseKin_mc, kTRUE);
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
	        lvE.SetPxPyPzE(px,py,pz,energy);
	        lvQ = _4vE0-lvE;
	        lvW = lvQ+_4vP0;
	        dAna->dTop.Q2 = -1*(lvQ.Mag2());
	        dAna->dTop.W = lvW.Mag();
	        
	        /* *** lvP, lvPip, lvPim, lvMMtop1(2,3,4) *** */
	        set4vecHadrons(dH10);
	              
	        /* *** MM *** */   
	        Float_t mm2ppippim = lvMMtop1.Mag2();
			Float_t mm2ppip    = lvMMtop2.Mag2();
			Float_t mm2ppim    = lvMMtop3.Mag2();
			Float_t mm2pippim  = lvMMtop4.Mag2();
						
			Bool_t t1a = gP && gPip && gPim;// && dH10->gpart == 4; //atrivedi: 071613
			Bool_t t1b = TMath::Abs(mm2ppippim) < 0.0005;
			Bool_t t2a = gP && gPip && !gPim;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t2b = mm2ppip>0 && mm2ppip<0.04;// && dAna->dTop.W < 1.1; //atrivedi tmp
			Bool_t t3a = gP && gPim && !gPip;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t3b = mm2ppim>0 && mm2ppim<0.04;// && dAna->dTop.W < 1.1; //atrivedi tmp
			Bool_t t4a = gPip && gPim && !gP;// && dH10->gpart == 3; //atrivedi: 071613
			Bool_t t4b = mm2pippim>0.8 && mm2pippim<1;
			
			Bool_t t1 = t1a && t1b;
			Bool_t t2 = t2a && t2b;
			Bool_t t3 = t3a && t3b;
			Bool_t t4 = t4a && t4b;
			
			if (t1a || t2a || t3a || t4a) { //used to determine top MM cut
				updateTop_MM();
				dAna->fillHistsMM(histsMM);
								
				updateEkin(dH10);
				dAna->fillHistsEkin(histseKin);
			}
			if ( t1 || t2 || t3 || t4) { //final top selection post MM cut
				pass = kTRUE;
				
				if (t1) {
					hevtsum->Fill(EVT_T1);
					dAna->top = 1;
					
					updateTop_VarSets(dH10);
					dAna->fillYields(yields_top1);
					
					dAna->fillHistsMM(histsMM_top1);
					dAna->fillHistsEkin(histseKin_top1);
				}
				if (t2) {
					hevtsum->Fill(EVT_T2);
					dAna->top = 2;
					
					lvPim = lvMMtop2;
					updateTop_VarSets(dH10);
					dAna->fillYields(yields_top2);
					
					dAna->fillHistsMM(histsMM_top2);
					dAna->fillHistsEkin(histseKin_top2);
				}
				if (t3) {
					hevtsum->Fill(EVT_T3);
					dAna->top = 3;
					
					lvPip = lvMMtop3;
					updateTop_VarSets(dH10);
					dAna->fillHistsMM(histsMM_top3);
					
					dAna->fillHistsEkin(histseKin_top3);
					dAna->fillYields(yields_top3);
				}
				if (t4) {
					hevtsum->Fill(EVT_T4);
					dAna->top = 4;
					
					lvP = lvMMtop4;
					updateTop_VarSets(dH10);
					dAna->fillHistsMM(histsMM_top4);
					
					dAna->fillHistsEkin(histseKin_top4);
					dAna->fillYields(yields_top4);
				}
				
				hevtsum->Fill(EVT_T0);
				dAna->fillHistsMM(histsMM_top0);
				dAna->fillHistsEkin(histseKin_top0);
				//atrivedi:080213
				//dAna->fillYields(yields_top0);
				
			} else (hevtsum->Fill(EVT_OTHER));
		}
	}
	if (pass) {
		EpProcessor::handle(dH10);
	}
}

void ProcTop::mcKin(DataH10 *dH10) {
	const TLorentzVector _4vE0 = dH10->_4vE0;
	const TLorentzVector _4vP0 = dH10->_4vP0;
	
	lvQ.SetXYZT(0,0,0,0);
	lvW.SetXYZT(0,0,0,0);
	lvE.SetXYZT(0,0,0,0);
	lvP.SetXYZT(0,0,0,0);
	lvPip.SetXYZT(0,0,0,0);
	lvPim.SetXYZT(0,0,0,0);
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
			lvE.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PROTON:
			lvP.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIP:
			lvPip.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		case PIM:
			lvPim.SetPxPyPzE(_px,_py,_pz,_energy);
			break;
		default:
			break;
		}
	}

	lvQ = _4vE0-lvE;
	lvW = lvQ+_4vP0;
	dAna->dTop_mc.Q2 = -1*(lvQ.Mag2());
	dAna->dTop_mc.W = lvW.Mag();
	
	updateEkin(dH10, kTRUE);
	updateTop_MM(kTRUE);
	updateTop_VarSets(dH10, kTRUE);
}

void ProcTop::updateEkin(DataH10* dH10, Bool_t useMc /*= kFALSE*/) {
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

void ProcTop::set4vecHadrons(DataH10* dH10){
	//Track Identified as hadrons
	if (dAna->h10idxP>0) {
		Double_t mom = dH10->p[dAna->h10idxP];
		Double_t px = mom*dH10->cx[dAna->h10idxP];
		Double_t py = mom*dH10->cy[dAna->h10idxP];
		Double_t pz = mom*dH10->cz[dAna->h10idxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		lvP.SetPxPyPzE(px,py,pz,energy);
	}
	if (dAna->h10idxPip>0) {
		Double_t mom = dH10->p[dAna->h10idxPip];
		Double_t px = mom*dH10->cx[dAna->h10idxPip];
		Double_t py = mom*dH10->cy[dAna->h10idxPip];
		Double_t pz = mom*dH10->cz[dAna->h10idxPip];
		Double_t energy = Sqrt(mom*mom+MASS_PIP*MASS_PIP);
		lvPip.SetPxPyPzE(px,py,pz,energy);
	}
	if (dAna->h10idxPim>0) {
		Double_t mom = dH10->p[dAna->h10idxPim];
		Double_t px = mom*dH10->cx[dAna->h10idxPim];
		Double_t py = mom*dH10->cy[dAna->h10idxPim];
		Double_t pz = mom*dH10->cz[dAna->h10idxPim];
		Double_t energy = Sqrt(mom*mom+MASS_PIM*MASS_PIM);
		lvPim.SetPxPyPzE(px,py,pz,energy);
	}
	
	lvMMtop1 = (lvW-(lvP+lvPip+lvPim));
	lvMMtop2 = (lvW-(lvP+lvPip));
    lvMMtop3 = (lvW-(lvP+lvPim));  
	lvMMtop4 = (lvW-(lvPip+lvPim)); 
}

void ProcTop::updateTop_MM(Bool_t ismc /* = kFALSE */) {
	DataTop *tp = &(dAna->dTop);
	if (ismc) tp = &(dAna->dTop_mc);
	tp->mm2ppippim = lvMMtop1.Mag2();
	tp->mmppippim  = lvMMtop1.Mag();
	tp->mm2ppip    = lvMMtop2.Mag2();
	tp->mmppip     = lvMMtop2.Mag();
	tp->mm2ppim    = lvMMtop3.Mag2();
	tp->mmppim     = lvMMtop3.Mag();
	tp->mm2pippim  = lvMMtop4.Mag2();
	tp->mmpippim   = lvMMtop4.Mag();
}

void ProcTop::updateTop_VarSets(DataH10* dH10, Bool_t ismc /* = kFALSE */){
	const TLorentzVector _4vE0 = dH10->_4vE0;
	const TLorentzVector _4vP0 = dH10->_4vP0;
	
	DataTop *tp = &(dAna->dTop);
	if (ismc) tp = &(dAna->dTop_mc);
	
	lvQCMS   = lvQ;
	lvP0CMS  = _4vP0;
	lvPCMS   = lvP;
	lvPipCMS = lvPip;
	lvPimCMS = lvPim;
	lvQCMS.Boost(-1*lvW.BoostVector());
	lvP0CMS.Boost(-1*lvW.BoostVector());
	lvPCMS.Boost(-1*lvW.BoostVector());
	lvPipCMS.Boost(-1*lvW.BoostVector());
	lvPimCMS.Boost(-1*lvW.BoostVector());
	
	Float_t Mppip = (lvPCMS + lvPipCMS).Mag();
	Float_t Mppim = (lvPCMS + lvPimCMS).Mag();
	Float_t Mpippim = (lvPipCMS + lvPimCMS).Mag();
	
	tp->varset1.M1 = Mppip;
	tp->varset1.M2 = Mpippim;
	tp->varset1.theta = getTheta(lvPimCMS);
	tp->varset1.phi = getPhi(lvPimCMS);
	//tp->varset1.alpha = getAlpha(lvPimCMS);
	tp->varset1.alpha = 180;
	
	
	tp->varset2.M1 = Mppip;
	tp->varset2.M2 = Mpippim;
	tp->varset2.theta = getTheta(lvPCMS);
	tp->varset2.phi = getPhi(lvPCMS);
	//tp->varset2.alpha = getAlpha(lvPCMS);
	tp->varset2.alpha = 180;
	
	tp->varset3.M1 = Mppip;
	tp->varset3.M2 = Mppim;
	tp->varset3.theta = getTheta(lvPipCMS);
	tp->varset3.phi = getPhi(lvPipCMS);
	//tp->varset1.alpha = getAlpha(lvPipCMS);
	tp->varset3.alpha = 180;
	
	//helicity
	if(!ismc) {tp->h = dH10->evthel;} //e1f
	
}

Float_t ProcTop::getTheta(TLorentzVector lv){
	Float_t retVal = 0;
	/* ***atrivedi 09-13-12*** */
	TVector3 lv_3vec = lv.Vect();
	TVector3 lvQCMS_3vec = lvQCMS.Vect();
	//retVal = RadToDeg()*ACos( (lv.Dot(lvQCMS))/(lv.P()*lvQCMS.P()) );
	retVal = RadToDeg()*ACos( (lv_3vec.Dot(lvQCMS_3vec))/(lv_3vec.Mag()*lvQCMS_3vec.Mag()) );
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
	
	if (histsMM != NULL) {
		dirout->cd();
		hevtsum->Write();
		histsMM->Write();
		histseKin->Write();
		dirout->cd("top0");
		yields_top0->Write();
		histsMM_top0->Write();
		histseKin_top0->Write();
		dirout->cd("top1");
		yields_top1->Write();
		histsMM_top1->Write();
		histseKin_top1->Write();
		dirout->cd("top2");
		yields_top2->Write();
		histsMM_top2->Write();
		histseKin_top2->Write();
		dirout->cd("top3");
		yields_top3->Write();
		histsMM_top3->Write();
		histseKin_top3->Write();
		dirout->cd("top4");
		yields_top4->Write();
		histsMM_top4->Write();
		histseKin_top4->Write();
	}
	
	if(dirout->cd("mc")){
		yields_mc->Write();
		histsMM_mc->Write();
		histseKin_mc->Write();
	}
}

#endif // PROCTOP_H
