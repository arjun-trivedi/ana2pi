#ifndef PROC_DINCLSV_H
#define PROC_DINCLSV_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "cuts.h"
#include "particle_constants.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

using namespace TMath;
using namespace ParticleConstants;

class ProcDinclsv : public EpProcessor {

public:
	ProcDinclsv(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcDinclsv();
	void handle();
	void write();
protected:
	TLorentzVector _lvE0;
	TLorentzVector _lvP0;
	DataInclsv* _dInclsv;
	TTree* _t;
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT,EVT_PASS};//, EVT_GPART, EVT_NE, EVT_NP};

	void UpdateDataInclsv();
	void AddBranches(TTree* t);
};

ProcDinclsv::ProcDinclsv(TDirectory *td,DataH10* dataH10,DataAna* dataAna) : EpProcessor(td, dataH10, dataAna) {
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"nevts");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"passed");

	_lvE0 = dH10->lvE0;
	_lvP0 = dH10->lvP0;
		
	_t = new TTree("t","TTree containing data for inclusive events");
	AddBranches(_t);

}

ProcDinclsv::~ProcDinclsv() {
	//delete hevtsum;
	
}

void ProcDinclsv::handle() {
	//Info("In ProcDinclsv::handle()");
	
	pass = kFALSE;
	hevtsum->Fill(EVT);

	UpdateDataInclsv();
	_t->Fill();
	pass = kTRUE;
	
	if (pass) {
		hevtsum->Fill(EVT_PASS);
		EpProcessor::handle();
	}
}

void ProcDinclsv::UpdateDataInclsv(){
	DataInclsv *d = &(dAna->dInclsv);

	d->gpart = dH10->gpart;
	//!for eid
	d->nphe = dH10->nphe[dH10->cc[0]-1];
	d->etot = dH10->etot[dH10->ec[0]-1];
	d->ec_ei = dH10->ec_ei[dH10->ec[0]-1];
	d->ec_eo = dH10->ec_eo[dH10->ec[0]-1];
	//! e- kinmatics
	Float_t p = dH10->p[0];
	Float_t px = p*dH10->cx[0];
	Float_t py = p*dH10->cy[0];
	Float_t pz = p*dH10->cz[0];
	TLorentzVector lvE(px,py,pz,TMath::Sqrt(p*p+MASS_E*MASS_E));
	TLorentzVector lvQ = _lvE0-lvE;
	TLorentzVector lvW = lvQ+_lvP0;
	d->p = lvE.P();
	d->theta = lvE.Theta()*RadToDeg();
	d->phi = lvE.Phi()*RadToDeg();
	d->Q2 = -1*(lvQ.Mag2());
	d->nu = lvQ.E();
	d->q = lvQ.P();
	d->x = (d->Q2)/(2*MASS_P*d->nu);
	d->W = lvW.Mag();	
}

void ProcDinclsv::AddBranches(TTree* t){
	DataInclsv *d = &(dAna->dInclsv);
	
	t->Branch("gpart",&d->gpart);
	//! for eid
	t->Branch("nphe",&d->nphe);
	t->Branch("etot",&d->etot);
	t->Branch("ec_ei",&d->ec_ei);
	t->Branch("ec_eo",&d->ec_eo);
	//! e- kinematics
	t->Branch("p",&d->p);
	t->Branch("theta",&d->theta);
	t->Branch("phi",&d->phi);
	t->Branch("Q2",&d->Q2);
	t->Branch("nu",&d->nu);
	t->Branch("q",&d->q);
	t->Branch("x",&d->x);
	t->Branch("W",&d->W);
}

void ProcDinclsv::write(){

}

#endif // PROC_DINCLSV_H
