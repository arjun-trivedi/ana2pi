#ifndef PROCEEFF_H
#define PROCEEFF_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "cut_eff.h"
#include <TMath.h>

using namespace TMath;

/*
 [03-20-15]
 + Event is passed if electron is within efficient region of the detector

 + Since cuts were developed for 2pi and that too mostly relying on top1,
   the monitor histograms will not look as expected since usually this processor
   will be called before d2piR and therefore the sample of electrons will be different.
     + However, to test functionality of this processor, call it after d2piR (with all
       other event selections prior to d2piR applied)
    
*/
class ProcEEff : public EpProcessor {

public:
	ProcEEff(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE);
	~ProcEEff();
	
	void handle();
	//void write();
	
protected:	
	CutEff* _cut_eff_tool;
		
	static const Int_t NUM_EVTCUTS = 2;
	enum { EVT_NULL, EVT, EVT_PASS};

	void updateEEff();
};

ProcEEff::ProcEEff(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                   Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/)
                   :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	_cut_eff_tool=new CutEff();
	td->cd();
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PASS,"pass");
	hevtsum->SetMinimum(0.);
}

ProcEEff::~ProcEEff()
{
	delete _cut_eff_tool;
}

void ProcEEff::handle()
{
	//Info("ProcEEff::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);
	
	//if (mon||mononly)
	//{
		if (dAna->d2pi.top==0 && hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
			dAna->makeHistsEEff(hists[MONMODE][EVTINC], dirmon);
		}else if(dAna->d2pi.top!=0 && hists[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
			for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
				TDirectory* dirmon = dirout->mkdir(TString::Format("monitor_top%d",iTop+1));
				dAna->makeHistsEEff(hists[MONMODE][iTop], dirmon);
			}
		}

		updateEEff();
		if (dAna->d2pi.top == 0) { //i.e inclusive event
			dAna->fillHistsEEff(hists[MONMODE][EVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsEEff(hists[MONMODE][dAna->d2pi.top-1]);
		}
	//}
	
		
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}
	
	bool in_efficient_region=_cut_eff_tool->InEfficientRegion(ELECTRON,dAna->eeff.sector,dAna->eeff.theta,dAna->eeff.p);
	if (in_efficient_region)
	{
		hevtsum->Fill(EVT_PASS);
		//if (mon)
		//{
			if (dAna->d2pi.top==0 && hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsEEff(hists[CUTMODE][EVTINC], dircut);
			}else if(dAna->d2pi.top!=0 && hists[CUTMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
				for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
					TDirectory* dircut = dirout->mkdir(TString::Format("cut_top%d",iTop+1));
					dAna->makeHistsEEff(hists[CUTMODE][iTop], dircut);
				}
			}

			if (dAna->d2pi.top == 0) { //i.e inclusive event
				dAna->fillHistsEEff(hists[CUTMODE][EVTINC]);
			}else { //i.e 2pi event
				dAna->fillHistsEEff(hists[CUTMODE][dAna->d2pi.top-1]);
			}
		//}

		pass = kTRUE;
		EpProcessor::handle();
	}
}

void ProcEEff::updateEEff() {
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
	    		
	dAna->eeff.sector=sector;	
	dAna->eeff.theta=_4vE1.Theta()*RadToDeg();
	dAna->eeff.p=_4vE1.P();
}

#endif // PROCEEFF_H
