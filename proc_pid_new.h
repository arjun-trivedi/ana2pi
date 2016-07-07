#ifndef PROCPIDNEW_H
#define PROCPIDNEW_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include "pid.h"

using namespace ParticleConstants;

/******************************************************
[06-09-15] New effort to develop PID for ana2pi

[11-28-15] Documentation
+ Pass Events only if:
	+ t1: p,pip,pim are found
	OR	
	+ t2: p,pip are found AND pim is not found
	OR
	+ t3: p,pim are found AND pip is not found
	OR
	+ t4: pip,pim are found AND p is not found
+ Import notes:
	1, The 1st identified hadron in gpart is assumed to be the true hadron and
  once it is found, the remaining tracks in gpart are not considered.
	2. A track can exclusively be identified as either p or pip or pim
	  i.e. it cannot be identified as both, at it can happen, for example,
	  in the high momentum region.
+ Output
	+ hevtsum: showing stats for number of events processed and number passing selection
	+ In 'mon' and 'mononly' mode: PID monitoring hists
	+ In 'make_tree' mode: Tree containing DataPidNew  (see DataPidNew for Tree's structure)

[06-27-16] Added 'gpart' and 'stat[ntrk]' to output Tree in order to debug SS-bands due to
           gpart-pid*stat-pid
*******************************************************/

class ProcPidNew : public EpProcessor
{

public:
	ProcPidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		    Bool_t monitor=kFALSE,Bool_t monitorOnly=kFALSE, Bool_t make_tree=kFALSE);
	ProcPidNew(DataH10* dataH10,DataAna* dataAna);
	~ProcPidNew();

	void handle();
	void updatePidNew();

protected:
	Bool_t _make_tree;
	Pid* _pid_tool;
	ProcEid* _proc_eid;
	
	TTree* _tree;

	static const Int_t NUM_EVTCUTS=9;
	enum { EVT_NULL, EVT, EVT_P_FOUND, EVT_PIP_FOUND, EVT_PIM_FOUND,
		   EVT_PPIPPIM_EX, EVT_PPIP_EX, EVT_PPIM_EX, EVT_PIPPIM_EX, EVT_OTHER};   
	     
	Float_t getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc);
};

ProcPidNew::ProcPidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
                 Bool_t monitor/* = kFALSE*/,Bool_t monitorOnly /*= kFALSE*/,Bool_t make_tree/*=kFALSE*/)
                 :EpProcessor(td, dataH10, dataAna, monitor, monitorOnly)
{
	td->cd();	
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_P_FOUND,"p found");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIP_FOUND,"#pi^{+} found");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIM_FOUND,"#pi^{-} found");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM_EX,"p + #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP_EX,   "p + #pi^{+} + (pi^{-}_{m})");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM_EX,   "p + (pi^{+}_{m}) + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM_EX, "(p_{m})+ #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");
	hevtsum->SetMinimum(0);

	_make_tree=make_tree;
	_pid_tool=new Pid(dH10->dtyp);
	_proc_eid=new ProcEid(dataH10,dataAna);
	
	if(_make_tree){
		_tree=new TTree("tree","TTree containing data for PID");
		//! Add PID data as branches to _t
		dAna->addBranches_DataPidNew(_tree);
		dAna->addBranches_DataEid(_tree);
	}
}

ProcPidNew::ProcPidNew(DataH10* dataH10,DataAna* dataAna)
                 :EpProcessor(dataH10, dataAna)
{
	
}

ProcPidNew::~ProcPidNew()
{
	delete _tree;
	delete _pid_tool;
	delete _proc_eid;
}

void ProcPidNew::handle()
{
	//Info("ProcPidNew::handle()", "");
	pass=kFALSE;
	hevtsum->Fill(EVT);

	//! 1. Obtain DataPidNew. NOTE, id[ntrk] and h10IdxP(Pip,Pim) are untouched!
	//! + At this point all the necessary data to perform PID is obtained and id[ntrk] is untouched i.e.
	//!   it continues to be equal to 0
	//! + PID is done later and there id[ntrk] & h10IdxP(Pip,Pim) are updated
	updatePidNew();

	//! If user wants to DataPidNew to study PID params, DataPidNew along with Ekin data is put in a Tree
	if (_make_tree){
		_proc_eid->updateEkin();
		_tree->Fill();

		pass=kTRUE;
		EpProcessor::handle();
		return;
	}

	//! If user want to monitor what PID histograms look like before applying PID cut
	if (mon||mononly)
	{
		if (dAna->d2pi.top==0 && hists[MONMODE][EVTINC][SECTOR0]==NULL) { //i.e. inclusive event
			TDirectory* dirmon = dirout->mkdir(TString::Format("monitor"));
			dAna->makeHistsPidMon(hists[MONMODE][EVTINC], dirmon);
		}else if(dAna->d2pi.top!=0 && hists[MONMODE][TOP1][SECTOR0]==NULL){ //i.e. 2pi event
			for(Int_t iTop=TOP1;iTop<NTOPS;iTop++){
				TDirectory* dirmon = dirout->mkdir(TString::Format("monitor%d",iTop));
				dAna->makeHistsPidMon(hists[MONMODE][iTop], dirmon);
			}
		}
	
		if (dAna->d2pi.top == 0) { //i.e inclusive event
			dAna->fillHistsPidMon(hists[MONMODE][EVTINC]);
		}else { //i.e 2pi event
			dAna->fillHistsPidMon(hists[MONMODE][dAna->d2pi.top-1]);
		}
	}
	
	if (mononly){
		pass = kTRUE;
		EpProcessor::handle();
		return;
	}

	
	//! 2. Apply PID cut using DataPidNew
	//! + Here id[ntrk] is set along with h10IdxP(Pip,Pim)
	//! + Additionally, dH10->id[] is also set
	//!
	//! + After finding first hadron, exit gpart loop.
	bool found_p,found_pip,found_pim;
	found_p=found_pip=found_pim=kFALSE;
	DataPidNew* dpid = &dAna->pidnew;
	for (int itrk=0;itrk<dpid->ntrk;itrk++){//! loop over ntrk in the event
		if (dpid->dc[itrk]>0 && dpid->sc[itrk]>0){
			if (dpid->q[itrk]==1){//! Select +ve trks
				if (!found_p && _pid_tool->is_proton(dpid->dt_p[itrk],dpid->p[itrk])){
					found_p=kTRUE;
					hevtsum->Fill(EVT_P_FOUND);
					dpid->id[itrk]=PROTON;
					dpid->h10IdxP=dpid->h10_idx[itrk];
					dH10->id[dpid->h10_idx[itrk]]=PROTON;
				}else if (!found_pip && _pid_tool->is_pip(dpid->dt_pip[itrk],dpid->p[itrk])){
					found_pip=kTRUE;
					hevtsum->Fill(EVT_PIP_FOUND);
					dpid->id[itrk]=PIP;
					dpid->h10IdxPip=dpid->h10_idx[itrk];
					dH10->id[dpid->h10_idx[itrk]]=PIP;
				}
			}
			if (dpid->q[itrk]==-1){//! Select -ve trks
				if (!found_pim && _pid_tool->is_pim(dpid->dt_pim[itrk],dpid->p[itrk])){
					found_pim=kTRUE;
					hevtsum->Fill(EVT_PIM_FOUND);
					dpid->id[itrk]=PIM;
					dpid->h10IdxPim=dpid->h10_idx[itrk];
					dH10->id[dpid->h10_idx[itrk]]=PIM;
				}
			}
		}
	}

	if (dpid->h10IdxP>0 && dpid->h10IdxPip>0 && dpid->h10IdxPim>0) {
		hevtsum->Fill(EVT_PPIPPIM_EX);
		pass = kTRUE;
	}else if (dpid->h10IdxP>0 && dpid->h10IdxPip>0 && dpid->h10IdxPim==0){
		hevtsum->Fill(EVT_PPIP_EX);
		pass = kTRUE;
	}else if (dpid->h10IdxP>0 && dpid->h10IdxPip==0 && dpid->h10IdxPim>0){
		hevtsum->Fill(EVT_PPIM_EX);
		pass = kTRUE;
	}else if (dpid->h10IdxP==0 && dpid->h10IdxPip>0 && dpid->h10IdxPim>0){
		hevtsum->Fill(EVT_PIPPIM_EX);
		pass = kTRUE;
	}else{
		hevtsum->Fill(EVT_OTHER);
	}

	
	//! If event passes, the call next processor
	if (pass) {
		//! If in mon mode fill PID cut monitoring histograms
		if (mon)
		{
			if (hists[CUTMODE][EVTINC][SECTOR0]==NULL) {
				TDirectory* dircut = dirout->mkdir(TString::Format("cut"));
				dAna->makeHistsPidCut(hists[CUTMODE][EVTINC], dircut);
			}
			dAna->fillHistsPidCut(hists[CUTMODE][EVTINC]);
		}
		
		EpProcessor::handle();
	}
}

void ProcPidNew::updatePidNew()
{
	DataPidNew* dpid = &dAna->pidnew;

	//! Directly measured electron quantities
	Float_t l_e=dH10->sc_r[dH10->sc[0]-1];
	Float_t t_e=dH10->sc_t[dH10->sc[0]-1];
	Float_t t_off=t_e-(l_e/SOL);

	//! Copy data into Branch variables
	//! directly measured quantities for e-
	dpid->l_e=l_e;
	dpid->t_e=t_e;
	dpid->t_off=t_off;
	//!gpart (should be equal to ntrk)
	dpid->gpart=dH10->gpart;
	//! update ntrk and measured quantities for ntrk
	for (Int_t i=1;i<dH10->gpart;i++) {
		if (dH10->q[i]==1 || dH10->q[i]==-1){
			//! Directly measured quantities
			Int_t q=dH10->q[i];
			Int_t dc=dH10->dc[i];
			Int_t sc=dH10->sc[i];
			Int_t stat=dH10->stat[i];
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
			dpid->stat[dpid->ntrk-1]=stat;
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

Float_t ProcPidNew::getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc){
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

#endif // PROCPIDNEW_H
