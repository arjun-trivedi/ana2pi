#ifndef PROC_DNSTR_DIS_H
#define PROC_DNSTR_DIS_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "cuts.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcDnstrDIS : public EpProcessor {

public:
	ProcDnstrDIS(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcDnstrDIS();
	void handle();
	void write();
protected:
	Bool_t _handle_ST;
	Int_t _ne;
	Int_t _np;
	Int_t _h10idx_p;
	TLorentzVector _lvE0;
	TLorentzVector _lvP0;
	DataNstrDIS* _dNstrDIS;
	TH2F* _hmmepvW;
	TTree* _t;
	TTree* _t_ST;
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT,EVT_PASS};//, EVT_GPART, EVT_NE, EVT_NP};

	void UpdateDataNstrDIS();
	//** Following 2 functions used to apply efid to ST-electrons
	bool inFid(); 
	float GetPhi(TLorentzVector lv);
	int GetSector(TLorentzVector lv); 
};

ProcDnstrDIS::ProcDnstrDIS(TDirectory *td,DataH10* dataH10,DataAna* dataAna) : EpProcessor(td, dataH10, dataAna) {
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"nevts");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"passed");

	_hmmepvW = new TH2F("hmmepvW","hmmepvW",100,0,4,100,-1,1);

		
	_lvE0 = dH10->lvE0;
	_lvP0 = dH10->lvP0;
	_dNstrDIS = &dAna->dNstrDIS;
	
	 /*_t = new TTree("t","TTree containing data for Elastic events");
	 _t->Branch("d","DataNstrDIS",_dNstrDIS,32000,0);
	 if(_handle_ST){
	 	_t_ST = new TTree("t_ST","TTree containing ST data for Elastic events");
	 	_t_ST->Branch("d","DataNstrDIS",_dNstrDIS_ST,32000,0);
	 }*/

}

ProcDnstrDIS::~ProcDnstrDIS() {
	//delete hevtsum;
	
}

void ProcDnstrDIS::handle() {
	//Info("In ProcDnstrDIS::handle()");
	
	pass = kFALSE;
	hevtsum->Fill(EVT);

	//! Nstrs, DIS events selection
	//! Basically need to remove DVCS events: ep->e'p'gamma
	_ne=0;
	_np=0;
	_h10idx_p=0;
	bool det_e=kFALSE;
	bool det_p=kFALSE;
		
	//if (dH10->gpart>4) return;
	for (int igpart = 0; igpart < dH10->gpart; igpart++)	{
		if (igpart==0 && dH10->id[igpart]==ELECTRON){
			_ne+=1;	
		} 
		if (dH10->id[igpart]==PROTON){
			_np+=1; 
			_h10idx_p=igpart;
		}   
	}
	if(_ne==1)det_e=kTRUE;
	if(_np>=1)det_p=kTRUE;
	if (det_e && !det_p){ //directly pass, since DVCS requires p'
		pass=kTRUE;
	}else if (det_e && det_p){//apply MM cut before passing
		UpdateDataNstrDIS();
		_hmmepvW->Fill(_dNstrDIS->W,_dNstrDIS->MMep);
	}

	if (pass) {
		//Info("In ProcDnstrDIS::handle()","pass\n");
		hevtsum->Fill(EVT_PASS);
		UpdateDataNstrDIS();
		_t->Fill();
		EpProcessor::handle();
	}
}

void ProcDnstrDIS::UpdateDataNstrDIS(){
	/*_dNstrDIS->gpart=dH10->gpart;
	_dNstrDIS->ne=_ne;
	_dNstrDIS->np=_np;*/

	//_dNstrDIS->sector=dH10->sc_sect[dH10->sc[0]-1];
	//_dNstrDIS->sector=dH10->sc_part;

	//create vxE,vxP,lvE,lvP
	//TVector3 vxE(dH10->vx[0],dH10->vy[0],dH10->vz[0]);
	Float_t p_e = dH10->p[0];
	Float_t px_e = p_e*dH10->cx[0];
	Float_t py_e = p_e*dH10->cy[0];
	Float_t pz_e = p_e*dH10->cz[0];
	TLorentzVector lvE(px_e,py_e,pz_e,TMath::Sqrt(p_e*p_e+MASS_E*MASS_E));

	//TVector3 vxP(dH10->vx[_h10idx_p],dH10->vy[_h10idx_p],dH10->vz[_h10idx_p]);
	Float_t p_p = dH10->p[_h10idx_p];
	Float_t px_p = p_p*dH10->cx[_h10idx_p];
	Float_t py_p = p_p*dH10->cy[_h10idx_p];
	Float_t pz_p = p_p*dH10->cz[_h10idx_p];
	TLorentzVector lvP(px_p,py_p,pz_p,TMath::Sqrt(p_p*p_p+MASS_P*MASS_P));

	//_dNstrDIS->vxE = vxE;
	//_dNstrDIS->vxP = vxP;
	//_dNstrDIS->lvE = lvE;
	//_dNstrDIS->lvP = lvP;

	//Q2,W
	TLorentzVector lvQ = _lvE0-lvE;
	TLorentzVector lvW = lvQ+_lvP0;
	_dNstrDIS->Q2 = -1*(lvQ.Mag2());
	_dNstrDIS->W = lvW.Mag();

	//MMp
	TLorentzVector lvMMep = lvW-(lvE+lvP);
	_dNstrDIS->MMep = lvMMep.Mag();
}

bool ProcDnstrDIS::inFid() {
	TLorentzVector lvE(0,0,0,0);
	for (Int_t inprt=0; inprt<dH10->nprt;inprt++)	{
		if (dH10->pidpart[inprt]==3){
			Float_t px=dH10->pxpart[inprt];
			Float_t py=dH10->pypart[inprt];
			Float_t pz=dH10->pzpart[inprt];
			Float_t e=dH10->epart[inprt];
			lvE.SetPxPyPzE(px,py,pz,e);
		}
	}
	Int_t id = 11;
	Double_t p = lvE.P();
	Float_t theta = lvE.Theta()*RadToDeg();
	Float_t phi = lvE.Phi()*RadToDeg(); //phitmp = [-180,180]
	Int_t sector=GetSector(lvE);
	Int_t paddle=0;
	Bool_t inFid = Cuts::Fiducial(id,p,theta,phi,sector,paddle);
	return inFid;
}

int ProcDnstrDIS::GetSector(TLorentzVector lv){
	float phi = GetPhi(lv);

	int sector = -9999;
	if( (330 < phi && phi <= 360) || (0 <= phi && phi <= 30) ){
		sector = 1;
	}
	if(30 < phi && phi <= 90){
		sector = 2;
	}
	if(90 < phi && phi <= 150){
		sector = 3;
	}
	if(150 < phi && phi <= 210){
		sector = 4;
	}
	if(210 < phi && phi <= 270){
		sector = 5;
	}
	if(270 < phi && phi <= 330){
		sector = 6;
	}
	return sector;
}

float ProcDnstrDIS::GetPhi(TLorentzVector lv){
	double phi = -9999;
	/* TLorentzVector Phi() returns angle between -pi and pi using ATan2(y,x)
	   The following transforms [-pi,pi]-->[0,2pi] */
	if(lv.Phi()*(180/TMath::Pi()) < 0){
		phi = lv.Phi()*(180/TMath::Pi()) + 360;
	} else {
		phi = lv.Phi()*(180/TMath::Pi());
	}
	return phi;
}


void ProcDnstrDIS::write(){

}

#endif // PROC_DNSTR_DIS_H
