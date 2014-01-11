#ifndef PROCMOMCOR_H
#define PROCMOMCOR_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "mom_corr.cpp"
#include "particle_constants.h"
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>

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
	void updateEkin(Bool_t useMc = kFALSE);
	
	class MomCorr_e1f *_pcorr;
	//TH1F *hdcx, *hdcy, *hdcz, *hdp;
	//TH2F *hdpVp;
};

ProcMomCor::ProcMomCor(TDirectory *td,DataH10* dataH10,DataAna* dataAna)
						 :EpProcessor(td, dataH10, dataAna)
{
	_pcorr = new MomCorr_e1f("/home/trivedia/CLAS/workspace/ana2pi/MomCorr");
}

ProcMomCor::~ProcMomCor()
{
	
}

void ProcMomCor::handle()
{
	//Info("ProcMomCor::handle()", "");
	pass = kFALSE;
	
	if (dAna->top==0 && hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
		TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
		dAna->makeHistsMomCor(hists[MONMODE][EVTINC], dirmon);
		dAna->makeHistsEkin(histsEkin[MONMODE][EVTINC], dirmon);
	}else if(dAna->top!=0 && hists[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
		for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor%d",iTop));
			dAna->makeHistsMomCor(hists[MONMODE][iTop], dirmon);
			dAna->makeHistsEkin(histsEkin[MONMODE][iTop], dirout);
		}
	}
	
	if ( dH10->id[0] == ELECTRON && dH10->p[0] > 0 ) {
		updateMomCor();		
		updateEkin();
		if (dAna->top == 0) { //i.e inclusive event
			dAna->fillHistsMomCor(hists[MONMODE][EVTINC]);
			dAna->fillHistsEkin(histsEkin[MONMODE][EVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsMomCor(hists[MONMODE][dAna->top-1]);
			dAna->fillHistsEkin(histsEkin[MONMODE][dAna->top-1]);
		}
	}
	pass = kTRUE; //atrivedi
	EpProcessor::handle();
}

void ProcMomCor::updateMomCor() {
		
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
