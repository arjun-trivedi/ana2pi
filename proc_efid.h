#ifndef PROCEFID_H
#define PROCEFID_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "cuts.h"
#include <TMath.h>

using namespace TMath;

class ProcEFid : public EpProcessor {

public:
	ProcEFid(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE);
	~ProcEFid();
	
	void handle();
	//void write();
	
protected:	
	Bool_t inFid();
	void updateEFid();
};

ProcEFid::ProcEFid(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                   Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/)
                   :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	
}

ProcEFid::~ProcEFid()
{
	
}

void ProcEFid::handle()
{
	//Info("ProcEFid::handle()", "");
	pass = kFALSE;
	
	if (mon||mononly)
	{
		if (dAna->top==0 && hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
			dAna->makeHistsEFid(hists[MONMODE][EVTINC], dirmon);
		}else if(dAna->top!=0 && hists[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
			for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
				TDirectory* dirmon = dirout->mkdir(TString::Format("monitor%d",iTop));
				dAna->makeHistsEFid(hists[MONMODE][iTop], dirmon);
			}
		}

		updateEFid();
		if (dAna->top == 0) { //i.e inclusive event
			dAna->fillHistsEFid(hists[MONMODE][EVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsEFid(hists[MONMODE][dAna->top]);
		}
	}
	
		
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}
	
	dAna->efid.fidE = inFid();
	if (dAna->efid.fidE)
	{
		if (mon)
		{
			if (hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsEFid(hists[CUTMODE][EVTINC], dircut);
			}
			dAna->fillHistsEFid(hists[CUTMODE][EVTINC]);/* code */
		}

		pass = kTRUE;
		EpProcessor::handle();
	}
}

Bool_t ProcEFid::inFid() {
		
	Int_t id = dH10->id[0];
	Double_t p = dH10->p[0];
	Float_t theta = RadToDeg()*ACos(dH10->cz[0]);
	Float_t phi = RadToDeg()*ATan2(dH10->cy[0],dH10->cx[0]); //phitmp = [-180,180]
	Int_t scidx = dH10->sc[0]-1;
	Int_t sector = 1;
	Int_t paddle = 0;
	if (scidx >= 0) {
		sector = dH10->sc_sect[scidx];
		paddle = dH10->sc_pd[scidx];
	}
	Bool_t inFid = Cuts::Fiducial(id, p, theta, phi, sector, paddle);
	return inFid;
}

void ProcEFid::updateEFid() {
	TLorentzVector _4vE1;
		
	Float_t mom;
	Float_t px;
	Float_t py;
	Float_t pz;
	mom = dH10->p[0];
	px = mom*dH10->cx[0];
	py = mom*dH10->cy[0];
	pz = mom*dH10->cz[0];
	_4vE1.SetPxPyPzE(px,py,pz,Sqrt(mom*mom+MASS_E*MASS_E));
	
	Int_t sector = dH10->sc_sect[dH10->sc[0]-1];
	Double_t phi = _4vE1.Phi()*RadToDeg(); //phi = [-180,180]
	//ensure that phi is between -30 and 330
    phi = phi < -30 ? phi += 360 : phi;
    		
	dAna->efid.sector = sector;	
	dAna->efid.phi = phi;
	dAna->efid.theta = _4vE1.Theta()*RadToDeg();
	dAna->efid.dc_xsc = dH10->dc_xsc[dH10->dc[0]-1];
	dAna->efid.dc_ysc = dH10->dc_ysc[dH10->dc[0]-1];
	dAna->efid.ech_x = dH10->ech_x[dH10->ec[0]-1];
	dAna->efid.ech_y = dH10->ech_y[dH10->ec[0]-1];
}

#endif // PROCEFID_H
