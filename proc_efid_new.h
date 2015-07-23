#ifndef PROCEFIDNEW_H
#define PROCEFIDNEW_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "fidfuncs.C"
#include "wrpr_cut_fid_e16.h"
#include <TMath.h>

using namespace TMath;

class ProcEFidNew : public EpProcessor {

public:
	ProcEFidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcEFidNew();
	
	void handle();
	//void write();
	
protected:	
	//Bool_t inFid();
	void updateEFid();
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT, EVT_PASS};

	TDirectory* _dirmon;
	TDirectory* _dircut;
};

ProcEFidNew::ProcEFidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna)
                   :EpProcessor(td, dataH10, dataAna)
{	
	if (dH10->expt=="e1f") {
    	Info("ProcEFidNew::ProcEFidNew()", "dH10.expt==E1F. Will use E1F Fiducial cuts"); 
    }else if (dH10->expt=="e16") {
    	Info("ProcEFidNew::ProcEFidNew()", "dH10.expt==E16. Will use E16 Fiducial cuts"); 
    }else{
    	Info("ProcEFidNew::ProcEFidNew()", "Could not determine dH10.expt! Will use E1F Fiducial cuts");
    }

	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"pass");

	//! Monitor output objects
	_dirmon = dirout->mkdir(TString::Format("monitor"));
	dAna->makeHistsEFid(hists[MONMODE][EVTINC], _dirmon);

	//! Cut output objects
	_dircut = dirout->mkdir(TString::Format("cut"));
	dAna->makeHistsEFid(hists[CUTMODE][EVTINC], _dircut);
}

ProcEFidNew::~ProcEFidNew()
{
	delete _dirmon;
	delete _dircut;
}

void ProcEFidNew::handle()
{
	//Info("ProcEFidNew::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);
	
	
	//! Monitoring mode
	updateEFid();
	dAna->fillHistsEFid(hists[MONMODE][EVTINC]);
		
	//! Cut mode
	//dAna->efid.fidE = inFid();
	dAna->efid.fidE = inFid(dAna->efid.p,dAna->efid.theta,dAna->efid.phi,dAna->efid.sector);
	if (dAna->efid.fidE)
	{
		hevtsum->Fill(EVT_PASS);
		dAna->fillHistsEFid(hists[CUTMODE][EVTINC]);/* code */
		
		pass = kTRUE;
		EpProcessor::handle();
	}
}

/*Bool_t ProcEFidNew::inFid() {
		
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

	Bool_t inFid=false;
	if (dH10->expt=="e1f"){
		inFid = Cuts::Fiducial(id, p, theta, phi, sector, paddle);
	}else if (dH10->expt=="e16"){
        inFid = Fiducial_e16_elctrn(id,p, theta, phi);
    }else{//! [05-17-15] For now, if expt is not determined, use E1F Fiducial cuts
    	inFid = Cuts::Fiducial(id, p, theta, phi, sector, paddle);
    }
	return inFid;
}*/

void ProcEFidNew::updateEFid() {
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
	dAna->efid.p = mom;
	dAna->efid.phi = phi;
	dAna->efid.theta = _4vE1.Theta()*RadToDeg();
	dAna->efid.dc_xsc = dH10->dc_xsc[dH10->dc[0]-1];
	dAna->efid.dc_ysc = dH10->dc_ysc[dH10->dc[0]-1];
	dAna->efid.ech_x = dH10->ech_x[dH10->ec[0]-1];
	dAna->efid.ech_y = dH10->ech_y[dH10->ec[0]-1];
}

#endif // PROCEFIDNEW_H
