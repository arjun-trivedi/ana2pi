#ifndef PROCPEFFNEW_H
#define PROCPEFFNEW_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "cut_eff.h"
#include <TMath.h>

using namespace TMath;

/*
 [03-20-15]
 + Processor must follow ProcSkimQ:ProcPidNew

 + Pass Events only if:
	+ dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0 && each of p,pip,pim 
	  are in efficient region of detector
	OR	
	+ dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && each of p,pip
	  are in efficient region of detector
	OR
	+ dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPim>0 && each of p,pim
      are in efficient region of detector
	OR
	+ dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0 && each of pip,pim
      are in efficient region of detector

 + Since cuts were developed for 2pi and that too mostly relying on top1,
   the monitor histograms will not look as expected since usually this processor
   will be called before d2piR and therefore the sample of electrons will be different.
     + However, to test functionality of this processor, call it after d2piR (with all
       other event selections prior to d2piR applied)
    
*/
class ProcPEffNew : public EpProcessor {

public:
	ProcPEffNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE);
	~ProcPEffNew();
	
	void handle();
	//void write();
	
protected:	
	CutEff* _cut_eff_tool;
		
	static const Int_t NUM_EVTCUTS=9;
	enum { EVT_NULL, EVT, 
		   EVT_PPIPPIM_EX,      EVT_PPIP_EX,      EVT_PPIM_EX,      EVT_PIPPIM_EX,
		   EVT_PPIPPIM_EX_PASS, EVT_PPIP_EX_PASS, EVT_PPIM_EX_PASS, EVT_PIPPIM_EX_PASS};

	void updatePEff();
};

ProcPEffNew::ProcPEffNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                   Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/)
                   :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	_cut_eff_tool=new CutEff();
	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM_EX,"p#pi^{+}#pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP_EX,"p#pi^{+}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM_EX,"p#pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM_EX,"#pi^{+}#pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM_EX_PASS,"p#pi^{+}#pi^{-} pass");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP_EX_PASS,"p#pi^{+} pass");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM_EX_PASS,"p#pi^{-} pass");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM_EX_PASS,"#pi^{+}#pi^{-} pass");
	hevtsum->SetMinimum(0.);
}

ProcPEffNew::~ProcPEffNew()
{
	delete _cut_eff_tool;
}

void ProcPEffNew::handle()
{
	//Info("ProcPEffNew::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);
	
	//if (mon||mononly)
	//{
		if (dAna->d2pi.top==0 && hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
			dAna->makeHistsPEff(hists[MONMODE][EVTINC], dirmon);
		}else if(dAna->d2pi.top!=0 && hists[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
			for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
				TDirectory* dirmon = dirout->mkdir(TString::Format("monitor_top%d",iTop+1));
				dAna->makeHistsPEff(hists[MONMODE][iTop], dirmon);
			}
		}

		updatePEff();
		if (dAna->d2pi.top == 0) { //i.e inclusive event
			dAna->fillHistsPEffNew(hists[MONMODE][EVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsPEffNew(hists[MONMODE][dAna->d2pi.top-1]);
		}
	//}
	
		
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}
	

	bool in_eff_rgn=kTRUE;
	if(dAna->skimq.isEVT_ETGT_2POS_ETGT_1NEG){
		if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0) {
			hevtsum->Fill(EVT_PPIPPIM_EX);
			bool in_eff_rgn_p=  _cut_eff_tool->InEfficientRegion(PROTON,dAna->peff.sector_p,  dAna->peff.theta_p,  dAna->peff.p_p);
			bool in_eff_rgn_pip=_cut_eff_tool->InEfficientRegion(PIP,   dAna->peff.sector_pip,dAna->peff.theta_pip,dAna->peff.p_pip);
			bool in_eff_rgn_pim=_cut_eff_tool->InEfficientRegion(PIM,   dAna->peff.sector_pim,dAna->peff.theta_pim,dAna->peff.p_pim);
			in_eff_rgn=in_eff_rgn_p && in_eff_rgn_pip && in_eff_rgn_pim;
			if (in_eff_rgn){
				hevtsum->Fill(EVT_PPIPPIM_EX_PASS);
			}
		}
	}
	if(dAna->skimq.isEVT_ETGT_2POS){
		if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPip>0){
			hevtsum->Fill(EVT_PPIP_EX);
			bool in_eff_rgn_p=  _cut_eff_tool->InEfficientRegion(PROTON,dAna->peff.sector_p,  dAna->peff.theta_p,  dAna->peff.p_p);
			bool in_eff_rgn_pip=_cut_eff_tool->InEfficientRegion(PIP,   dAna->peff.sector_pip,dAna->peff.theta_pip,dAna->peff.p_pip);
			in_eff_rgn=in_eff_rgn_p && in_eff_rgn_pip;
			if (in_eff_rgn){
				hevtsum->Fill(EVT_PPIP_EX_PASS);
			}
		}
	}
	if(dAna->skimq.isEVT_ETGT_1POS_ETGT_1NEG){
		if (dAna->pidnew.h10IdxP>0 && dAna->pidnew.h10IdxPim>0){
			hevtsum->Fill(EVT_PPIM_EX);
			bool in_eff_rgn_p=  _cut_eff_tool->InEfficientRegion(PROTON,dAna->peff.sector_p,  dAna->peff.theta_p,  dAna->peff.p_p);
			bool in_eff_rgn_pim=_cut_eff_tool->InEfficientRegion(PIM,   dAna->peff.sector_pim,dAna->peff.theta_pim,dAna->peff.p_pim);
			in_eff_rgn=in_eff_rgn_p && in_eff_rgn_pim;
			if (in_eff_rgn){
				hevtsum->Fill(EVT_PPIM_EX_PASS);
			}
		}
		if(dAna->pidnew.h10IdxPip>0 && dAna->pidnew.h10IdxPim>0){
			hevtsum->Fill(EVT_PIPPIM_EX);
			bool in_eff_rgn_pip=_cut_eff_tool->InEfficientRegion(PIP,   dAna->peff.sector_pip,dAna->peff.theta_pip,dAna->peff.p_pip);
			bool in_eff_rgn_pim=_cut_eff_tool->InEfficientRegion(PIM,   dAna->peff.sector_pim,dAna->peff.theta_pim,dAna->peff.p_pim);
			in_eff_rgn=in_eff_rgn_pip && in_eff_rgn_pim;
			if (in_eff_rgn){
				hevtsum->Fill(EVT_PIPPIM_EX_PASS);
			}
		}
	}
	//bool in_efficient_region=_cut_eff_tool->InEfficientRegion(ELECTRON,dAna->eeff.sector,dAna->eeff.theta,dAna->eeff.p);
	if (in_eff_rgn)
	{
		//if (mon)
		//{
			if (dAna->d2pi.top==0 && hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsPEff(hists[CUTMODE][EVTINC], dircut);
			}else if(dAna->d2pi.top!=0 && hists[CUTMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
				for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
					TDirectory* dircut = dirout->mkdir(TString::Format("cut_top%d",iTop+1));
					dAna->makeHistsPEff(hists[CUTMODE][iTop], dircut);
				}
			}

			if (dAna->d2pi.top == 0) { //i.e inclusive event
				dAna->fillHistsPEffNew(hists[CUTMODE][EVTINC]);
			}else { //i.e 2pi event
				dAna->fillHistsPEffNew(hists[CUTMODE][dAna->d2pi.top-1]);
			}
		//}

		pass = kTRUE;
		EpProcessor::handle();
	}
}

void ProcPEffNew::updatePEff() {
	TLorentzVector lvP,lvPip,lvPim;
	if (dAna->pidnew.h10IdxP>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxP];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxP];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxP];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxP];
		Double_t energy = Sqrt(mom*mom+MASS_P*MASS_P);
		lvP.SetPxPyPzE(px,py,pz,energy);

		dAna->peff.sector_p=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxP]-1];
		dAna->peff.theta_p=lvP.Theta()*RadToDeg();
		dAna->peff.p_p=lvP.P();
	}
	if (dAna->pidnew.h10IdxPip>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxPip];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxPip];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxPip];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxPip];
		Double_t energy = Sqrt(mom*mom+MASS_PIP*MASS_PIP);
		lvPip.SetPxPyPzE(px,py,pz,energy);

		dAna->peff.sector_pip=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxPip]-1];
		dAna->peff.theta_pip=lvPip.Theta()*RadToDeg();
		dAna->peff.p_pip=lvPip.P();
	}
	if (dAna->pidnew.h10IdxPim>0) {
		Double_t mom = dH10->p[dAna->pidnew.h10IdxPim];
		Double_t px = mom*dH10->cx[dAna->pidnew.h10IdxPim];
		Double_t py = mom*dH10->cy[dAna->pidnew.h10IdxPim];
		Double_t pz = mom*dH10->cz[dAna->pidnew.h10IdxPim];
		Double_t energy = Sqrt(mom*mom+MASS_PIM*MASS_PIM);
		lvPim.SetPxPyPzE(px,py,pz,energy);

		dAna->peff.sector_pim=dH10->sc_sect[dH10->sc[dAna->pidnew.h10IdxPim]-1];
		dAna->peff.theta_pim=lvPim.Theta()*RadToDeg();
		dAna->peff.p_pim=lvPim.P();
	}
}

#endif // PROCPEFFNEW_H
