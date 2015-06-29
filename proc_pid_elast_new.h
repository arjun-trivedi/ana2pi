#ifndef PROCPIDELASTNEW_H
#define PROCPIDELASTNEW_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include "pid.h"

using namespace ParticleConstants;

/******************************************************
[06-28-15]
+ Note, ProcPidElastNew must follow ProcSkimQElast
+ Pass Events only if:
	+ dAna->skimq_elast.isEVT_EQGT_1POS AND proton is identified
*******************************************************/

class ProcPidElastNew : public EpProcessor
{

public:
	ProcPidElastNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		    Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE);
	ProcPidElastNew(DataH10* dataH10,DataAna* dataAna);
	~ProcPidElastNew();

	void handle();
	//void write();
	void updatePidNew();
				
protected:
	Pid* _pid_tool;

	static const Int_t NUM_EVTCUTS=3;
	enum { EVT_NULL, EVT, EVT_P_EX, EVT_OTHER};   
	     
	//void updatePid();
	Float_t getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc);
};

ProcPidElastNew::ProcPidElastNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                 Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/)
                 :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	td->cd();	
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_EX,"p");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");

	_pid_tool=new Pid(dH10->dtyp);
}

ProcPidElastNew::ProcPidElastNew(DataH10* dataH10,DataAna* dataAna)
                 :EpProcessor(dataH10, dataAna)
{
	delete _pid_tool;
}

ProcPidElastNew::~ProcPidElastNew()
{
	
}

void ProcPidElastNew::handle()
{
	//Info("ProcPidElastNew::handle()", "");
	pass = kFALSE;
	hevtsum->Fill(EVT);

	//! 1. Obtain DataPidNew. NOTE, id[ntrk] is untouched!
	//! + At this point all the necessary data to perform PID is obtained and id[ntrk] is untouched i.e.
	//!   it continues to be equal to 0
	//! + PID is done later and there id[ntrk] & h10IdxP(Pip,Pim) are updated
	updatePidNew();
	
	if (mon||mononly)
	{
		if (hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
			dAna->makeHistsPidElastMon(hists[MONMODE][EVTINC], dirmon);
		}	
	
		if (dAna->d2pi.top == 0) { //i.e inclusive event
			dAna->fillHistsPidElastMon(hists[MONMODE][EVTINC]);
		}
	}
	
		
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}	

	//! 2. Apply PID cut using DataPidElastNew
	//! Here id[ntrk] is set along with h10IdxP
	//! Additionally, dH10->id[] is also set
	DataPidElastNew* dpid = &dAna->pid_elastnew;
	for (int itrk=0;itrk<dpid->ntrk;itrk++){//! loop over ntrk in the event
		if (dpid->q[itrk]==1){//! Select +ve trks
			if (_pid_tool->is_proton(dpid->dt_p[itrk],dpid->p[itrk])){
				dpid->id[itrk]=PROTON;
				dpid->h10IdxP=dpid->h10_idx[itrk];
				dH10->id[dpid->h10_idx[itrk]]==PROTON;
			}
		}
	}

	//! 3. Finally decide if event passes selection criterion based on qskim and PID 
	if(dAna->skimq_elast.isEVT_EQGT_1POS){
	 	if (dpid->h10IdxP>0) {
			hevtsum->Fill(EVT_P_EX);
			pass = kTRUE;
		}else{
			hevtsum->Fill(EVT_OTHER);
		}
	}
	
	if (pass) {
		if (mon)
		{
			if (hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsPidElastCut(hists[CUTMODE][EVTINC], dircut);
			}
			dAna->fillHistsPidElastCut(hists[CUTMODE][EVTINC]);
		}
		EpProcessor::handle();
	}
}

void ProcPidElastNew::updatePidNew()
{
	DataPidElastNew* dpid = &dAna->pid_elastnew;

	//! Directly measured electron quantities
	Float_t l_e=dH10->sc_r[dH10->sc[0]-1];
	Float_t t_e=dH10->sc_t[dH10->sc[0]-1];
	Float_t t_off=t_e-(l_e/SOL);

	//! Copy data into Branch variables
	dpid->l_e=l_e;
	dpid->t_e=t_e;
	dpid->t_off=t_off;
	
	for (Int_t i=1;i<dH10->gpart;i++) {
		if (dH10->q[i]==1){
			//! Directly measured quantities
			Int_t q=dH10->q[i];
			Int_t dc=dH10->dc[i];
			Int_t sc=dH10->sc[i];
			Float_t p=dH10->p[i];
			Float_t l=dH10->sc_r[dH10->sc[i]-1];
			Float_t t=dH10->sc_t[dH10->sc[i]-1];
			Float_t b=( l/(t-t_off) )/SOL; 
			Int_t sector=dH10->sc_sect[dH10->sc[i]-1];
			Int_t id=dH10->id[i];

			//! Quantities under particle assumption
			Float_t b_p  =TMath::Sqrt((p*p)/(MASS_P*MASS_P+p*p));
			Float_t b_pip=TMath::Sqrt((p*p)/(MASS_PIP*MASS_PIP+p*p));
			Float_t b_pim=TMath::Sqrt((p*p)/(MASS_PIM*MASS_PIM+p*p));

			Float_t dt_p  =l/(b_p*SOL)+t_off- t;
			Float_t dt_pip=l/(b_pip*SOL)+t_off-t;
			Float_t dt_pim=l/(b_pim*SOL)+t_off-t;

			//! Copy data into Branch variables
			dpid->ntrk+=1;
			dpid->q[dpid->ntrk-1]=q;
			dpid->dc[dpid->ntrk-1]=dc;
			dpid->sc[dpid->ntrk-1]=sc;
			dpid->q[dpid->ntrk-1]=q;
			dpid->p[dpid->ntrk-1]=p;
			dpid->l[dpid->ntrk-1]=l;
			dpid->t[dpid->ntrk-1]=t;
			dpid->b[dpid->ntrk-1]=b;
			dpid->sector[dpid->ntrk-1]=sector;
			dpid->h10_idx[dpid->ntrk-1]=i;
			//dpid->id[dpid->ntrk-1]=0;//! Deliberately not set here! Set while performing PID
			dpid->b_p[dpid->ntrk-1]=b_p;
			dpid->b_pip[dpid->ntrk-1]=b_pip;
			dpid->b_pim[dpid->ntrk-1]=b_pim;
			dpid->dt_p[dpid->ntrk-1]=dt_p;
			dpid->dt_pip[dpid->ntrk-1]=dt_pip;
			dpid->dt_pim[dpid->ntrk-1]=dt_pim;
		}
	}
}

Float_t ProcPidElastNew::getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc){
	//Define CC plane equation: Ax + By + Cz + D = 0
	Float_t A = -0.000785;
	Float_t B = 0;
	Float_t C = -0.00168;
	Float_t D = 1;
	Float_t magS = TMath::Sqrt(A*A + B*B + C*C);

	//length of line perpendicular to CC plane from hit in SC plane
	Float_t h = 0;

	//cosine of angle between h & t vector
	Float_t cAlpha = 0;

	//magnitude of vector t
	Float_t t = 0;

	h = ((A*x_sc + B*y_sc + C*z_sc) + D)/magS;

	cAlpha = (cx_sc*x_sc + cy_sc*y_sc + cz_sc*z_sc)/magS;

	t = h/cAlpha;

	Float_t x_cc = x_sc - t*cx_sc;
	Float_t y_cc = y_sc - t*cy_sc;
	Float_t z_cc = z_sc - t*cz_sc;

	return TMath::ACos(z_cc/TMath::Sqrt(x_cc*x_cc + y_cc*y_cc + z_cc*z_cc)) * (180/TMath::Pi());
}

#endif // PROCPIDELASTNEW_H
