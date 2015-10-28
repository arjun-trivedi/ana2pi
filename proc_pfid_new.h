#ifndef PROCPFIDNEW_H
#define PROCPFIDNEW_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "wrpr_cut_fid_e16.h"
#include <TMath.h>

using namespace TMath;

/*
 [11-28-15]
 + Pass Events only if:
	+ t1: dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0
	      && each of p,pip,pim are in fiducial region of detector
	OR	
	+ t2: dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim==0
	      && each of p,pip are in fiducial region of detector
	OR
	+ t3: dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip==0 && dAna->pidnew.h10IdxPim>0
	  && each of p,pim are in fiducial region of detector
	OR
	+ t4: dAna->pidnew.h10IdxP==0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0
	  && each of pip,pim are in fiducial region of detector

*/

//! CLAS sector phi angle coverage
const Int_t PHI_FULL_SECTOR[6][2]={
	{-30,30},
	{30,90},
	{90,150},
	{150,210},
	{210,270},
	{270,330}
};

class ProcPFidNew : public EpProcessor {

public:
	ProcPFidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE);
	~ProcPFidNew();
	
	void handle();
	//void write();
	
protected:	
	static const Int_t NUM_EVTCUTS=9;
	enum { EVT_NULL, EVT, 
		   EVT_PPIPPIM_EX,      EVT_PPIP_EX,      EVT_PPIM_EX,      EVT_PIPPIM_EX,
		   EVT_PPIPPIM_EX_PASS, EVT_PPIP_EX_PASS, EVT_PPIM_EX_PASS, EVT_PIPPIM_EX_PASS};

	void updatePFid();
	int get_sector(float phi);
};

ProcPFidNew::ProcPFidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                   Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/)
                   :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM_EX,"p + #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP_EX,   "p + #pi^{+} + (pi^{-}_{m})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM_EX,   "p + (pi^{+}_{m}) + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM_EX, "(p_{m})+ #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM_EX_PASS,"p + #pi^{+} + #pi^{-} infid");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP_EX_PASS,   "p + #pi^{+} + (pi^{-}_{m}) infid");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM_EX_PASS,   "p + (pi^{+}_{m}) + #pi^{-} infid");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM_EX_PASS, "(p_{m})+ #pi^{+} + #pi^{-} infid");
	hevtsum->SetMinimum(0.);
}

ProcPFidNew::~ProcPFidNew()
{
	
}

void ProcPFidNew::handle()
{
	//Info("ProcPFidNew::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);
	
	//if (mon||mononly)
	//{
		if (dAna->d2pi.top==0 && hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
			dAna->makeHistsPFid(hists[MONMODE][EVTINC], dirmon);
		}else if(dAna->d2pi.top!=0 && hists[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
			for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
				TDirectory* dirmon = dirout->mkdir(TString::Format("monitor_top%d",iTop+1));
				dAna->makeHistsPFid(hists[MONMODE][iTop], dirmon);
			}
		}

		updatePFid();
		if (dAna->d2pi.top == 0) { //i.e inclusive event
			dAna->fillHistsPFidNew(hists[MONMODE][EVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsPFidNew(hists[MONMODE][dAna->d2pi.top-1]);
		}
	//}
	
		
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}
	

	bool in_fid_rgn=kFALSE;
	if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0) {
		hevtsum->Fill(EVT_PPIPPIM_EX);
		bool in_fid_rgn_p=  Fiducial_e16_hdrn(dAna->pfid.theta_p,  dAna->pfid.phi_p,  dAna->pfid.sector_p);
		bool in_fid_rgn_pip=Fiducial_e16_hdrn(dAna->pfid.theta_pip,dAna->pfid.phi_pip,dAna->pfid.sector_pip);
		bool in_fid_rgn_pim=Fiducial_e16_hdrn(dAna->pfid.theta_pim,dAna->pfid.phi_pim,dAna->pfid.sector_pim);
		in_fid_rgn=in_fid_rgn_p && in_fid_rgn_pip && in_fid_rgn_pim;
		if (in_fid_rgn){
			hevtsum->Fill(EVT_PPIPPIM_EX_PASS);
		}
	}else if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim==0) {
		hevtsum->Fill(EVT_PPIP_EX);
		bool in_fid_rgn_p=  Fiducial_e16_hdrn(dAna->pfid.theta_p,  dAna->pfid.phi_p,  dAna->pfid.sector_p);
		bool in_fid_rgn_pip=Fiducial_e16_hdrn(dAna->pfid.theta_pip,dAna->pfid.phi_pip,dAna->pfid.sector_pip);
		in_fid_rgn=in_fid_rgn_p && in_fid_rgn_pip;
		if (in_fid_rgn){
			hevtsum->Fill(EVT_PPIP_EX_PASS);
		}
	}else if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip==0 && dAna->pidnew.h10IdxPim>0) {
		hevtsum->Fill(EVT_PPIM_EX);
		bool in_fid_rgn_p=  Fiducial_e16_hdrn(dAna->pfid.theta_p,  dAna->pfid.phi_p,  dAna->pfid.sector_p);
		bool in_fid_rgn_pim=Fiducial_e16_hdrn(dAna->pfid.theta_pim,dAna->pfid.phi_pim,dAna->pfid.sector_pim);
		in_fid_rgn=in_fid_rgn_p && in_fid_rgn_pim;
		if (in_fid_rgn){
			hevtsum->Fill(EVT_PPIM_EX_PASS);
		}
	}else if (dAna->pidnew.h10IdxP==0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0) {
		hevtsum->Fill(EVT_PIPPIM_EX);
		bool in_fid_rgn_pip=Fiducial_e16_hdrn(dAna->pfid.theta_pip,dAna->pfid.phi_pip,dAna->pfid.sector_pip);
		bool in_fid_rgn_pim=Fiducial_e16_hdrn(dAna->pfid.theta_pim,dAna->pfid.phi_pim,dAna->pfid.sector_pim);
		in_fid_rgn=in_fid_rgn_pip && in_fid_rgn_pim;
		if (in_fid_rgn){
			hevtsum->Fill(EVT_PIPPIM_EX_PASS);
		}
	}
	
	if (in_fid_rgn)
	{
		//if (mon)
		//{
			if (dAna->d2pi.top==0 && hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsPFid(hists[CUTMODE][EVTINC], dircut);
			}else if(dAna->d2pi.top!=0 && hists[CUTMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
				for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
					TDirectory* dircut = dirout->mkdir(TString::Format("cut_top%d",iTop+1));
					dAna->makeHistsPFid(hists[CUTMODE][iTop], dircut);
				}
			}

			if (dAna->d2pi.top == 0) { //i.e inclusive event
				dAna->fillHistsPFidNew(hists[CUTMODE][EVTINC]);
			}else { //i.e 2pi event
				dAna->fillHistsPFidNew(hists[CUTMODE][dAna->d2pi.top-1]);
			}
		//}

		pass = kTRUE;
		EpProcessor::handle();
	}
}

void ProcPFidNew::updatePFid() {
	TLorentzVector lvP,lvPip,lvPim;
	if (dAna->pidnew.h10IdxP>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxP];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxP];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxP];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		lvP.SetPxPyPzE(px,py,pz,energy);

		//dAna->pfid.sector_p=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxP]-1];
		dAna->pfid.theta_p=lvP.Theta()*RadToDeg();
		dAna->pfid.phi_p=lvP.Phi()*RadToDeg(); //phi = [-180,180]
		dAna->pfid.phi_p=dAna->pfid.phi_p<-30?dAna->pfid.phi_p+360:dAna->pfid.phi_p; // [-30,330]
		dAna->pfid.p_p=lvP.P();
		dAna->pfid.sector_p=get_sector(dAna->pfid.phi_p);
	}
	if (dAna->pidnew.h10IdxPip>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxPip];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxPip];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxPip];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxPip];
		Double_t energy = Sqrt(mom*mom+MASS_PIP*MASS_PIP);
		lvPip.SetPxPyPzE(px,py,pz,energy);

		//dAna->pfid.sector_pip=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxPip]-1];
		dAna->pfid.theta_pip=lvPip.Theta()*RadToDeg();
		dAna->pfid.phi_pip=lvPip.Phi()*RadToDeg(); //phi = [-180,180]
		dAna->pfid.phi_pip=dAna->pfid.phi_pip<-30?dAna->pfid.phi_pip+360:dAna->pfid.phi_pip; // [-30,330]
		dAna->pfid.p_pip=lvPip.P();
		dAna->pfid.sector_pip=get_sector(dAna->pfid.phi_pip);
	}
	if (dAna->pidnew.h10IdxPim>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxPim];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxPim];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxPim];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxPim];
		Double_t energy = Sqrt(mom*mom+MASS_PIM*MASS_PIM);
		lvPim.SetPxPyPzE(px,py,pz,energy);

		//dAna->pfid.sector_pim=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxPim]-1];
		dAna->pfid.theta_pim=lvPim.Theta()*RadToDeg();
		dAna->pfid.phi_pim=lvPim.Phi()*RadToDeg(); //phi = [-180,180]
		dAna->pfid.phi_pim=dAna->pfid.phi_pim<-30?dAna->pfid.phi_pim+360:dAna->pfid.phi_pim; // [-30,330]
		dAna->pfid.p_pim=lvPim.P();
		dAna->pfid.sector_pim=get_sector(dAna->pfid.phi_pim);
	}
}

/*
phi has to be [-30,330] for the following function to work.
+ [10-25-15] Copied from elast_lite/h10looper_e1f.cpp
*/
int ProcPFidNew::get_sector(float phi){
	int sector=-1;
	if ( (phi>=PHI_FULL_SECTOR[0][0]  && phi<PHI_FULL_SECTOR[0][1]) || phi==PHI_FULL_SECTOR[5][1]){
		sector=1;
	}else if (phi>=PHI_FULL_SECTOR[1][0] && phi<PHI_FULL_SECTOR[1][1]){
		sector=2;
	}else if (phi>=PHI_FULL_SECTOR[2][0] && phi<PHI_FULL_SECTOR[2][1]){
		sector=3;
	}else if (phi>=PHI_FULL_SECTOR[3][0] && phi<PHI_FULL_SECTOR[3][1]){
		sector=4;
	}else if (phi>=PHI_FULL_SECTOR[4][0] && phi<PHI_FULL_SECTOR[4][1]){
		sector=5;
	}else if (phi>=PHI_FULL_SECTOR[5][0] && phi<PHI_FULL_SECTOR[5][1]){
		sector=6;
	}else{
		Info("ProcPFidNew::get_sector","sector not found for phi=%.2f. Returning sector=-1",phi);
	}
	return sector;
}

#endif // PROCPFIDNEW_H
