#ifndef PROCDELAST_H
#define PROCDELAST_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcDelast : public EpProcessor {

public:
	ProcDelast(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcDelast();
	void handle();
	void write();
protected:
	Int_t _ne;
	Int_t _np;
	Int_t _h10idx_p;
	TLorentzVector _lvE0;
	TLorentzVector _lvP0;
	DataElastic* _dElast;
	TTree* _t;
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT,EVT_PASS};//, EVT_GPART, EVT_NE, EVT_NP};

	void UpdateDataElastic();
};

ProcDelast::ProcDelast(TDirectory *td,DataH10* dataH10,DataAna* dataAna) : EpProcessor(td, dataH10, dataAna) {
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"nevts");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"passed");
	
	_lvE0 = dH10->lvE0;
	_lvP0 = dH10->lvP0;
	_dElast = &dAna->dElast;

	 _t = new TTree("t","TTree containing data for Elastic events");
	 _t->Branch("d","DataElastic",_dElast,32000,0);
}

ProcDelast::~ProcDelast() {
	//delete hevtsum;
	
}

void ProcDelast::handle() {
	//Info("In ProcDelast::handle()");
	pass = kFALSE;

	hevtsum->Fill(EVT);

	//elastic event selection
	_ne=0;
	_np=0;
	_h10idx_p=0;
	bool det_e=kFALSE;
	bool det_p=kFALSE;
		
	if (dH10->gpart>4) return;
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
	if (det_e && det_p)pass=kTRUE;

	if (pass) {
		//Info("In ProcDelast::handle()","pass\n");
		hevtsum->Fill(EVT_PASS);
		UpdateDataElastic();
		_t->Fill();
		EpProcessor::handle();
	}/*else{
		Info("In ProcDelast::handle()","not pass\n");
		UpdateDataElastic();
		_t->Fill();
	}*/
}

void ProcDelast::UpdateDataElastic(){
	_dElast->gpart=dH10->gpart;
	_dElast->ne=_ne;
	_dElast->np=_np;

	//create lvE,lvP
	Float_t p_e = dH10->p[0];
	Float_t px_e = p_e*dH10->cx[0];
	Float_t py_e = p_e*dH10->cy[0];
	Float_t pz_e = p_e*dH10->cz[0];
	TLorentzVector lvE(px_e,py_e,pz_e,TMath::Sqrt(p_e*p_e+MASS_E*MASS_E));

	Float_t p_p = dH10->p[_h10idx_p];
	Float_t px_p = p_p*dH10->cx[_h10idx_p];
	Float_t py_p = p_p*dH10->cy[_h10idx_p];
	Float_t pz_p = p_p*dH10->cz[_h10idx_p];
	TLorentzVector lvP(px_p,py_p,pz_p,TMath::Sqrt(p_p*p_p+MASS_P*MASS_P));

	_dElast->lvE = lvE;
	_dElast->lvP = lvP;

	//Q2,W
	TLorentzVector lvQ = _lvE0-lvE;
	TLorentzVector lvW = lvQ+_lvP0;
	_dElast->Q2 = -1*(lvQ.Mag2());
	_dElast->W = lvW.Mag();

	//MMp
	TLorentzVector lvMMp = lvQ+_lvP0-lvP;
	_dElast->MMp = lvMMp.Mag();
}


void ProcDelast::write(){

}

#endif // PROCDELAST_H
