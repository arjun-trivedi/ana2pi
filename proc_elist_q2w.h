#ifndef PROCELISTQ2W_H
#define PROCELISTQ2W_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TH1.h>
#include <TEntryList.h>

using namespace ParticleConstants;

class ProcEListQ2W : public EpProcessor
{

public:
	ProcEListQ2W(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcEListQ2W();
	
	void handle();

	TEntryList* _el[kQ2W_NCrsBins];
	
protected:
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT, EVT_PASS };
	
	void updateEkin(Bool_t useMc = kFALSE);
	bool _useMc;
};

ProcEListQ2W::ProcEListQ2W(TDirectory *td,DataH10* dataH10,DataAna* dataAna)
						 : EpProcessor(td, dataH10, dataAna)
{
	
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total Events");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"Events in ana-region");

	//! Create Q2W TEntryLists
	for(int i=0;i<kQ2W_NCrsBins;i++){
		_el[i]=new TEntryList(TString::Format("elist_q2w_%d",i+1),TString::Format("elist_q2w_%d",i+1));
	}

	//! atrivedi [06-13-14]: The following logic is based on
	//! the assumption that when making subsets of 'q2wFull',
	//!   - for simulation, Q2W filtering will be done at ST level
	//!   - for experiment, Q2W filtering will be done at ER level
	_useMc=kFALSE;
	if (dH10->dtyp=="sim") {
		_useMc=kTRUE;
	}
}

ProcEListQ2W::~ProcEListQ2W()
{
	delete _el;
	delete hevtsum;
	delete hists;
	//delete _hists_ekin;
}

void ProcEListQ2W::handle()
{
	// Info("ProcEListQ2W::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);
	
	updateEkin(_useMc);
	// Info("ProcEListQ2W::handle()", "updated Ekin");	
	DataEkin *ekin = &dAna->eKin;
	if (_useMc) ekin = &dAna->eKin_mc;
	for(int i=0;i<kQ2_NCrsBins;i++){
		for(int j=0;j<kW_NCrsBins;j++){
			if (ekin->Q2>=kQ2_CrsBin[i].xmin && ekin->Q2<kQ2_CrsBin[i].xmax && 
				ekin->W>=kW_CrsBin[j].xmin && ekin->W<kW_CrsBin[j].xmax){
				Int_t bin=(j+1)+(i*7);
				Int_t bin_idx=bin-1;
				// Info("ProcEListQ2W::handle()", "i,j,bin,bin_idx=%d,%d,%d,%d",i,j,bin,bin_idx);
				_el[bin_idx]->Enter(dH10->get_ientry_h10chain(),dH10->h10chain);
				// Info("ProcEListQ2W::handle()", "update _el");

				hevtsum->Fill(EVT_PASS);
				pass = kTRUE;
				EpProcessor::handle();
			}
		}
	}
}

void ProcEListQ2W::updateEkin(Bool_t useMc /*= kFALSE*/) {
	const TLorentzVector _4vE0 = dH10->lvE0;
	const TLorentzVector _4vP0 = dH10->lvP0;
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
	if (!useMc) {ekin->sector = dH10->sc_sect[dH10->sc[0]-1];}
	ekin->W = (_4vQ+_4vP0).Mag();
	ekin->Q2 = -1*_4vQ.Mag2();
	// ekin->nu = _4vQ.E();
	// ekin->xb = ekin->Q2/(2*MASS_P*ekin->nu);
	// ekin->E1 = _4vE1.E();
	// ekin->theta1 = _4vE1.Theta()*RadToDeg();
	// Double_t phitmp = _4vE1.Phi()*RadToDeg(); //phitmp = [-180,180]
	// ekin->phi1 = phitmp<-30 ? phitmp+360 : phitmp;
	// ekin->theta = _4vQ.Theta()*RadToDeg();
	// phitmp = _4vQ.Phi()*RadToDeg(); //phitmp = [-180,180]
	// ekin->phi = phitmp<-30 ? phitmp+360 : phitmp;
}

#endif // PROCELISTQ2W_H
