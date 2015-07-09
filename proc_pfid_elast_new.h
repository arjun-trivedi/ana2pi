#ifndef PROCPFIDELASTNEW_H
#define PROCPFIDELASTNEW_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "wrpr_cut_fid_e16.h"
#include <TMath.h>

using namespace TMath;

/*
 [06-28-15]
 + Processor must follow ProcSkimQElast:ProcPidElastNew

 + Pass Events only if:
	+ dAna->pid_elastnew.h10IdxP>0 && each of p is in fiducial region of detector
*/
class ProcPFidElastNew : public EpProcessor {

public:
	ProcPFidElastNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE);
	~ProcPFidElastNew();
	
	void handle();
	//void write();
	
protected:	
	static const Int_t NUM_EVTCUTS=3;
	enum { EVT_NULL, EVT, 
		   EVT_P_EX,
		   EVT_P_EX_PASS};

	void updatePFid();
};

ProcPFidElastNew::ProcPFidElastNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                   Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/)
                   :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_EX,"p");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_EX_PASS,"p pass");
	hevtsum->SetMinimum(0.);
}

ProcPFidElastNew::~ProcPFidElastNew()
{
	
}

void ProcPFidElastNew::handle()
{
	//Info("ProcPFidElastNew::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);
	
	if (hists[MONMODE][EVTINC][SECTOR0]==NULL) { 
		TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
		dAna->makeHistsPFidElast(hists[MONMODE][EVTINC], dirmon);
	}
		
	updatePFid();
	dAna->fillHistsPFidElastNew(hists[MONMODE][EVTINC]);
	
	
		
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}
	

	bool in_fid_rgn=kTRUE;
	if(dAna->skimq_elast.isEVT_EQGT_1POS){
		if (dAna->pid_elastnew.h10IdxP>0) {
			hevtsum->Fill(EVT_P_EX);
			bool in_fid_rgn_p=  Fiducial_e16_hdrn(dAna->pfid_elast.theta_p,  dAna->pfid_elast.phi_p,  dAna->pfid_elast.sector_p);
			in_fid_rgn=in_fid_rgn_p;
			if (in_fid_rgn){
				hevtsum->Fill(EVT_P_EX_PASS);
			}
		}
	}
	
	if (in_fid_rgn)
	{
		if (hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
			TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
			dAna->makeHistsPFidElast(hists[CUTMODE][EVTINC], dircut);
		}
		dAna->fillHistsPFidElastNew(hists[CUTMODE][EVTINC]);
		
		pass = kTRUE;
		EpProcessor::handle();
	}
}

void ProcPFidElastNew::updatePFid() {
	TLorentzVector lvP,lvPip,lvPim;
	if (dAna->pid_elastnew.h10IdxP>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxP];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxP];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxP];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		lvP.SetPxPyPzE(px,py,pz,energy);

		dAna->pfid_elast.sector_p=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxP]-1];
		dAna->pfid_elast.theta_p=lvP.Theta()*RadToDeg();
		dAna->pfid_elast.phi_p=lvP.Phi()*RadToDeg(); //phi = [-180,180]
		dAna->pfid_elast.p_p=lvP.P();
	}
}

#endif // PROCPFIDELASTNEW_H
