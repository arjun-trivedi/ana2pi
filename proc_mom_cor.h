#ifndef PROCMOMCOR_H
#define PROCMOMCOR_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "mom_corr.cpp"
#include "wrpr_e_corr_sub.h" //! For e16:pcorr:e
#include "wrpr_pi_corr_sub.h" //! For e16:pcorr:pip
#include "wrpr_pr.h" //! For e16:pcorr:p (actually eloss, not pcorr)
#include "particle_constants.h"
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <stdlib.h> 

/*
[01-08-16]
+ Added capability to apply pcorr:{e,pip,p}. Note that for p, it is actually 'eloss' 
  rather than 'pcorr':
  	+ pcorr=Correction of reconstructed momenta due to "B-field" irregularities
  	+ eloss=Correction for classical energy lost by particle traversing matter
+ Note that as a result, for e16, and for the sake of integriy for e1f too,
  this processor must be called after 'pidnew' (if 'pid' is used, "manual intervention" will be required;
  for more on the detail of this "manual intervention" see notes for proc_d2pi.h written on [10-28-15]) 
  since h10-indices set up by pidnew are used to obtain the kinematic information for pip and p.
*/

using namespace ParticleConstants;

class ProcMomCor : public EpProcessor
{
public:
	ProcMomCor(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcMomCor();
	
	void handle();
	//virtual void write(); 
	
protected:
	void updateMomCor();
	void updateMomCor_hadron(TString hdrn_name);
	void updateEkin(Bool_t useMc = kFALSE);
	
	class MomCorr_e1f *_pcorr;
	//TH1F *hdcx, *hdcy, *hdcz, *hdp;
	//TH2F *hdpVp;
};

ProcMomCor::ProcMomCor(TDirectory *td,DataH10* dataH10,DataAna* dataAna)
						 :EpProcessor(td, dataH10, dataAna)
{
	TString path;
  	path=getenv("WORKSPACE");
  	TString fpcorr=TString::Format("%s/ana2pi/MomCorr",path.Data());
  	_pcorr = new MomCorr_e1f((char *)fpcorr.Data());
	//_pcorr = new MomCorr_e1f("/home/trivedia/CLAS/workspace/ana2pi/MomCorr");
}

ProcMomCor::~ProcMomCor()
{
	
}

void ProcMomCor::handle()
{
	//Info("ProcMomCor::handle()", "");
	pass = kFALSE;
	
	if (dAna->d2pi.top==0 && hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
		TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
		dAna->makeHistsMomCor(hists[MONMODE][EVTINC], dirmon);
		dAna->makeHistsEkin(histsEkin[MONMODE][EVTINC], dirmon);
	}else if(dAna->d2pi.top!=0 && hists[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
		for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor%d",iTop));
			dAna->makeHistsMomCor(hists[MONMODE][iTop], dirmon);
			dAna->makeHistsEkin(histsEkin[MONMODE][iTop], dirout);
		}
	}
	
	if ( dH10->id[0] == ELECTRON && dH10->p[0] > 0 ) {
		updateMomCor();		
		updateEkin();
		if (dAna->d2pi.top == 0) { //i.e inclusive event
			dAna->fillHistsMomCor(hists[MONMODE][EVTINC]);
			dAna->fillHistsEkin(histsEkin[MONMODE][EVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsMomCor(hists[MONMODE][dAna->d2pi.top-1]);
			dAna->fillHistsEkin(histsEkin[MONMODE][dAna->d2pi.top-1]);
		}
	}
	pass = kTRUE; //atrivedi
	EpProcessor::handle();
}

void ProcMomCor::updateMomCor() {
		
	if (dH10->expt=="e1f"){
		TLorentzVector *pin = new TLorentzVector();
		Float_t cx = dH10->cx[0];
		Float_t cy = dH10->cy[0];
		Float_t cz = dH10->cz[0];
		Float_t p = dH10->p[0];
		Int_t q = dH10->q[0];
		Int_t id = dH10->id[0];
		Float_t px = p*cx;
		Float_t py = p*cy;
		Float_t pz = p*cz;
		Float_t mass = GetPdgMass(id);
		pin->SetXYZM(px,py,pz,mass);
		TLorentzVector cp = _pcorr->PcorN(*pin,q,id);
		dH10->p[0] = cp.P();
		dH10->cx[0] = cp.Px()/cp.P();
		dH10->cy[0] = cp.Py()/cp.P();
		dH10->cz[0] = cp.Pz()/cp.P();
		dAna->mom.sector = dH10->sc_sect[dH10->sc[0]-1];
		dAna->mom.p = p;
		dAna->mom.dcx = dH10->cx[0]-cx;
		dAna->mom.dcy = dH10->cy[0]-cy;
		dAna->mom.dcz = dH10->cz[0]-cz;
		dAna->mom.dp = dH10->p[0]-p;
		delete pin;
	}else if (dH10->expt=="e16"){
		//! Get all required quanties for e_corr_sub
		float p = dH10->p[0];
		float cx = dH10->cx[0];
		float cy = dH10->cy[0];
		float cz = dH10->cz[0];
		float px = p*cx;
		float py = p*cy;
		float pz = p*cz;
		TLorentzVector lvE1(px,py,pz,Sqrt(p*p+MASS_E*MASS_E));
		float theta=lvE1.Theta()*RadToDeg();
		float phi=lvE1.Phi()*RadToDeg(); //! [-180,180] = valid phi range for e16 .f subroutines =[-180,180]
		float torcur=3375;
		int sector=dH10->sc_sect[dH10->sc[0]-1];
		//! Now call e_corr_sub and pass by reference theta_corr & mom_corr
		float theta_corr=0;
		float p_corr=0;
		e_corr_sub(theta,phi,p,torcur,sector,theta_corr,p_corr);
		//! Now form updated TLorentzVector
		float px_corr=p_corr*Sin(theta_corr*DegToRad())*Cos(phi*DegToRad());
		float py_corr=p_corr*Sin(theta_corr*DegToRad())*Sin(phi*DegToRad());
		float pz_corr=p_corr*Cos(theta_corr*DegToRad());
		TLorentzVector lvE1_corr(px_corr,py_corr,pz_corr,Sqrt(p_corr*p_corr+MASS_E*MASS_E));		
		//! Now updated dH10
		dH10->p[0] = lvE1_corr.P();
		dH10->cx[0] = lvE1_corr.Px()/lvE1_corr.P();
		dH10->cy[0] = lvE1_corr.Py()/lvE1_corr.P();
		dH10->cz[0] = lvE1_corr.Pz()/lvE1_corr.P();
		//! Update DataMom
		dAna->mom.sector = sector;
		dAna->mom.p = p;
		dAna->mom.dcx = dH10->cx[0]-cx;
		dAna->mom.dcy = dH10->cy[0]-cy;
		dAna->mom.dcz = dH10->cz[0]-cz;
		dAna->mom.dp = dH10->p[0]-p;

		if (dH10->rctn=="2pi"){
			//! pcorr of pip and p
			updateMomCor_hadron("pip");
			updateMomCor_hadron("p");
		}
	}
}

void ProcMomCor::updateMomCor_hadron(TString hdrn_name){
	int h10idx;
	float mass_hdrn;
	if (hdrn_name=="pip"){
		h10idx=dAna->pidnew.h10IdxPip;
		mass_hdrn=MASS_PIP;
	}else if (hdrn_name=="p"){
		h10idx=dAna->pidnew.h10IdxP;
		mass_hdrn=MASS_P;
	}else{
		Info("ProcMomCor::updateMomCor_hadron()", "Not implemented for hdrn_name=%s",hdrn_name.Data());
		return;
	}

	//! Get all required quanties 
	float p = dH10->p[h10idx];
	float cx = dH10->cx[h10idx];
	float cy = dH10->cy[h10idx];
	float cz = dH10->cz[h10idx];
	float px = p*cx;
	float py = p*cy;
	float pz = p*cz;
	TLorentzVector lv(px,py,pz,Sqrt(p*p+mass_hdrn*mass_hdrn));
	float theta=lv.Theta()*RadToDeg();
	float phi=lv.Phi()*RadToDeg(); //! [-180,180] = valid phi range for e16 .f subroutines =[-180,180]
	float torcur=3375;
	int sector=dH10->sc_sect[dH10->sc[h10idx]-1];

	//! + Prepare to call .f routines
	//! + Note for p 'eloss' is done and therefore theta_corr=theta
	//!   i.e. theta is not corrected
	float theta_corr=0;
	float p_corr=0;
	if (hdrn_name=="pip"){
		pi_corr_sub(theta,phi,p,torcur,sector,theta_corr,p_corr);
	}else if (hdrn_name=="p"){
		preloss(p,theta,p_corr);
		theta_corr=theta;
	}
	//! Now form updated TLorentzVector
	float px_corr=p_corr*Sin(theta_corr*DegToRad())*Cos(phi*DegToRad());
	float py_corr=p_corr*Sin(theta_corr*DegToRad())*Sin(phi*DegToRad());
	float pz_corr=p_corr*Cos(theta_corr*DegToRad());
	TLorentzVector lv_corr(px_corr,py_corr,pz_corr,Sqrt(p_corr*p_corr+mass_hdrn*mass_hdrn));		
	//! Now updated dH10
	dH10->p[h10idx] = lv_corr.P();
	dH10->cx[h10idx] = lv_corr.Px()/lv_corr.P();
	dH10->cy[h10idx] = lv_corr.Py()/lv_corr.P();
	dH10->cz[h10idx] = lv_corr.Pz()/lv_corr.P();
	//! Update DataMom
	if (hdrn_name=="pip"){
		dAna->mom.sector_pip = sector;
		dAna->mom.p_pip = p;
		dAna->mom.dcx_pip = dH10->cx[h10idx]-cx;
		dAna->mom.dcy_pip = dH10->cy[h10idx]-cy;
		dAna->mom.dcz_pip = dH10->cz[h10idx]-cz;
		dAna->mom.dp_pip = dH10->p[h10idx]-p;
	}else if (hdrn_name=="p"){
		dAna->mom.sector_p = sector;
		dAna->mom.p_p = p;
		dAna->mom.dcx_p = dH10->cx[h10idx]-cx;
		dAna->mom.dcy_p = dH10->cy[h10idx]-cy;
		dAna->mom.dcz_p = dH10->cz[h10idx]-cz;
		dAna->mom.dp_p = dH10->p[h10idx]-p;
	}
	return;
}

void ProcMomCor::updateEkin(Bool_t useMc /*= kFALSE*/) {
	const TLorentzVector lvE0 = dH10->lvE0;
	const TLorentzVector lvP0 = dH10->lvP0;
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
	if (!useMc) {ekin->sector = dH10->sc_sect[dH10->sc[0]-1];}
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
}

#endif // PROCMOMCOR_H
