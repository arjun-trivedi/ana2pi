#ifndef PROCSKIMQ2W_H
#define PROCSKIMQ2W_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TH1.h>

using namespace ParticleConstants;

class ProcSkimQ2W : public EpProcessor
{

public:
	ProcSkimQ2W(TDirectory *td, DataAna* dataAna, TString h10type);
	~ProcSkimQ2W();
	
	void handle(DataH10* dH10);
	
protected:
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT, EVT_PASS };
	
	void updateEkin(DataH10* dH10, Bool_t useMc = kFALSE);
};

ProcSkimQ2W::ProcSkimQ2W(TDirectory *td, DataAna* dataAna, TString h10type) : EpProcessor(td, dataAna, h10type)
{
	dirout->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->SetMinimum(0);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"pass Q2W skim");
}

ProcSkimQ2W::~ProcSkimQ2W()
{
	delete hevtsum;
	delete hists;
}

void ProcSkimQ2W::handle(DataH10* dH10)
{
	//Info("ProcSkimQ2W::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);
	
	if (dAna->top==0 && histsEkin[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
		TDirectory* dirmon = dirout->mkdir(TString::Format("mon"));
		dAna->makeHistsEkin(histsEkin[MONMODE][EVTINC], dirmon);
	}else if(dAna->top!=0 && histsEkin[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
		for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
			TDirectory* dirmon = dirout->mkdir(TString::Format("mon%d",iTop));
			dAna->makeHistsEkin(histsEkin[MONMODE][iTop], dirmon);
		}
	}
	
	updateEkin(dH10);
	
	if (dAna->top == 0) { //i.e inclusive event
		dAna->fillHistsEkin(histsEkin[MONMODE][EVTINC]);
	}else{ //i.e 2pi event
		dAna->fillHistsEkin(histsEkin[MONMODE][dAna->top]);
	}
	
	
	//if ( (dAna->eKin.Q2 >= 1.3 && dAna->eKin.Q2 <= 1.6 && dAna->eKin.W >=1.5 && dAna->eKin.W <= 1.9) ) { //!q2w range 1
	if ( (dAna->eKin.Q2 >= 1.9 && dAna->eKin.Q2 <= 2.5 && dAna->eKin.W >=1.3 && dAna->eKin.W <= 1.9) ) { //!q2w range 2; Evgeny Isupov	 
		if (histsEkin[CUTMODE][EVTINC][SECTOR0]==NULL) {
			TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
			dAna->makeHistsEkin(histsEkin[CUTMODE][EVTINC], dircut);
		}
		dAna->fillHistsEkin(histsEkin[CUTMODE][EVTINC]);
		
		hevtsum->Fill(EVT_PASS);
		pass = kTRUE;
		EpProcessor::handle(dH10);
	}
	
}

void ProcSkimQ2W::updateEkin(DataH10* dH10, Bool_t useMc /*= kFALSE*/) {
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
	if (!useMc) {ekin->sector = dH10->sc_sect[dH10->sc[0]-1];}
	ekin->W = (_4vQ+_4vP0).Mag();
	ekin->Q2 = -1*_4vQ.Mag2();
	ekin->nu = _4vQ.E();
	ekin->xb = ekin->Q2/(2*MASS_P*ekin->nu);
	ekin->E1 = _4vE1.E();
	ekin->theta1 = _4vE1.Theta()*RadToDeg();
	Double_t phitmp = _4vE1.Phi()*RadToDeg(); //phitmp = [-180,180]
	ekin->phi1 = phitmp<-30 ? phitmp+360 : phitmp;
	ekin->theta = _4vQ.Theta()*RadToDeg();
	phitmp = _4vQ.Phi()*RadToDeg(); //phitmp = [-180,180]
	ekin->phi = phitmp<-30 ? phitmp+360 : phitmp;
}

#endif // PROCSKIMQ2W_H
