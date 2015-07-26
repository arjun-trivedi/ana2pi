#ifndef PROCPIDNEW_H
#define PROCPIDNEW_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include "pid.h"

using namespace ParticleConstants;

/******************************************************
[06-09-15] New effort to develop PID for ana2pi

[06-24-15] Documentation
+ Note, ProcPid must follow ProcSkimQ
+ Pass Events only if:
	+ dAna->skimq.isEVT_ETGT_2POS_ETGT_1NEG AND proton,pip,pim are identified
	OR	+ 
	+ dAna->skimq.isEVT_ETGT_2POS     AND proton,pip are identified
	OR
	+ dAna->skimq.isEVT_ETGT_1POS_ETGT_1NEG AND (p,pim identified OR pip,pim identfied)
+ Output
	+ hevtsum: showing stats for number of events processed and number passing selection
	+ In 'mon' and 'mononly' mode: PID monitoring hists
	+ In 'make_tree' mode: Tree containing DataPidNew  (see DataPidNew for Tree's structure)
*******************************************************/

class ProcPidNew : public EpProcessor
{

public:
	ProcPidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna, 
		     Bool_t make_tree=kFALSE);
	ProcPidNew(DataH10* dataH10,DataAna* dataAna);
	~ProcPidNew();

	void handle();
	void updatePidNew();

protected:
	Bool_t _make_tree;
	Pid* _pid_tool;
	ProcEid* _proc_eid;
	
	//TTree* _tree;
	TDirectory* _dirmon;
	TDirectory* _dircut;
	TTree* _t[NPROCMODES];

	static const Int_t NUM_EVTCUTS=6;
	enum { EVT_NULL, EVT, EVT_PPIPPIM_EX, EVT_PPIP_EX, EVT_PPIM_EX, EVT_PIPPIM_EX, EVT_OTHER};   
	     
	Float_t getCCtheta(Float_t x_sc, Float_t y_sc, Float_t z_sc, Float_t cx_sc, Float_t cy_sc, Float_t cz_sc);
};

ProcPidNew::ProcPidNew(TDirectory *td,DataH10* dataH10,DataAna* dataAna,Bool_t make_tree/*=kFALSE*/)
                 :EpProcessor(td, dataH10, dataAna)
{
	td->cd();	
	hevtsum = new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
	hevtsum->GetXaxis()->SetBinLabel(EVT,"Total");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIPPIM_EX,"p + #pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIP_EX,"p + #pi^{+}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PPIM_EX,"p + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_PIPPIM_EX,"#pi^{+} + #pi^{-}");
	hevtsum->GetXaxis()->SetBinLabel(EVT_OTHER,"other");

	_make_tree=make_tree;
	_dirmon=NULL;
	_dircut=NULL;
	_pid_tool=new Pid(dH10->dtyp);
	_proc_eid=new ProcEid(dataH10,dataAna);

	//! Create monitor objects
	_dirmon = dirout->mkdir(TString::Format("monitor"));
	dAna->makeHistsPidMon(hists[MONMODE][EVTINC], _dirmon);
			
	if(_make_tree){
		_dirmon->cd();
		_t[MONMODE]=new TTree("t","TTree containing data for PID");
		//! Add PID data as branches to _t
		dAna->addBranches_DataPidNew(_t[MONMODE]);
		//dAna->addBranches_DataEid(_t[MONMODE]);
	}

	//! Create cut objects
	_dircut = dirout->mkdir(TString::Format("cut"));
	dAna->makeHistsPidMon(hists[CUTMODE][EVTINC], _dircut);
			
	if(_make_tree){
		_dircut->cd();
		_t[CUTMODE]=new TTree("t","TTree containing data for PID");
		//! Add PID data as branches to _t
		dAna->addBranches_DataPidNew(_t[CUTMODE]);
		//dAna->addBranches_DataEid(_t[CUTMODE]);
	}
}

ProcPidNew::ProcPidNew(DataH10* dataH10,DataAna* dataAna)
                 :EpProcessor(dataH10, dataAna)
{
	
}

ProcPidNew::~ProcPidNew()
{
	//delete _tree;
	delete _t;
	delete _dirmon;
	delete _dircut;
	delete _pid_tool;
	delete _proc_eid;
}

void ProcPidNew::handle()
{
	//Info("ProcPidNew::handle()", "");
	pass=kFALSE;
	hevtsum->Fill(EVT);

	//! Monitor
	//! 1. Obtain DataPidNew. NOTE, id[ntrk],np,npip,npim are untouched!
	//! + At this point all the necessary data to perform PID is obtained and id[ntrk] is untouched i.e.
	//!   it continues to be equal to 0
	//! + PID is done later and there id[ntrk] & h10IdxP(Pip,Pim) are updated
	updatePidNew();

	dAna->fillHistsPidMon(hists[MONMODE][EVTINC]);
	if (_make_tree){
		_proc_eid->updateEkin();
		_t[MONMODE]->Fill();
	}

	
	//! Cut mode	
	//! 2. Apply PID cut using DataPidNew
	//! Here id[ntrk] is set along with h10IdxP(Pip,Pim)
	//! Additionally, dH10->id[] is also set
	DataPidNew* dpid = &dAna->pidnew;
	for (int itrk=0;itrk<dpid->ntrk;itrk++){//! loop over ntrk in the event
		if (dpid->q[itrk]==1){//! Select +ve trks
			if (_pid_tool->is_proton(dpid->dt_p[itrk],dpid->p[itrk])){
				dpid->np+=1;
				dpid->id[itrk]=PROTON;
				dpid->h10IdxP=dpid->h10_idx[itrk];
				dH10->id[dpid->h10_idx[itrk]]==PROTON;
			}else if (_pid_tool->is_pip(dpid->dt_pip[itrk],dpid->p[itrk])){
				dpid->npip+=1;
				dpid->id[itrk]=PIP;
				dpid->h10IdxPip=dpid->h10_idx[itrk];
				dH10->id[dpid->h10_idx[itrk]]==PIP;
			}
		}
		if (dpid->q[itrk]==-1){//! Select -ve trks
			if (_pid_tool->is_pim(dpid->dt_pim[itrk],dpid->p[itrk])){
				dpid->npim+=1;
				dpid->id[itrk]=PIM;
				dpid->h10IdxPim=dpid->h10_idx[itrk];
				dH10->id[dpid->h10_idx[itrk]]==PIM;
			}
		}
	}

	//! 3. Finally decide if event passes selection criterion based on qskim and PID for ana2pi
	if(dAna->skimq.isEVT_ETGT_2POS_ETGT_1NEG){
	 	if (dpid->h10IdxP>0 && dpid->h10IdxPip>0 && dpid->h10IdxPim>0) {
			hevtsum->Fill(EVT_PPIPPIM_EX);
			pass = kTRUE;
		}else{
			hevtsum->Fill(EVT_OTHER);
		}
	}
	if(dAna->skimq.isEVT_ETGT_2POS){
		if (dpid->h10IdxP>0 && dpid->h10IdxPip>0) {
			hevtsum->Fill(EVT_PPIP_EX);
			pass = kTRUE;
		}else{
			hevtsum->Fill(EVT_OTHER);
		}
	}
	if(dAna->skimq.isEVT_ETGT_1POS_ETGT_1NEG){
		if(dpid->h10IdxP>0 && dpid->h10IdxPim>0) {
			hevtsum->Fill(EVT_PPIM_EX);
			pass = kTRUE;
		}else if(dpid->h10IdxPip>0 && dpid->h10IdxPim>0) {
			hevtsum->Fill(EVT_PIPPIM_EX);
			pass = kTRUE;
		}else{
			hevtsum->Fill(EVT_OTHER);	
		}
	}
	
	//! If event passes, the call next processor
	pass=kTRUE;
	if (pass) {
		dAna->fillHistsPidMon(hists[CUTMODE][EVTINC]);
		if (_make_tree){
			_proc_eid->updateEkin();
			_t[CUTMODE]->Fill();
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
	dpid->l_e=l_e;
	dpid->t_e=t_e;
	dpid->t_off=t_off;
	dpid->gpart=dH10->gpart;
	for (Int_t i=1;i<dH10->gpart;i++) {
		if (dH10->q[i]==1 || dH10->q[i]==-1){
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
			if      (dH10->q[i]==1)  dpid->npos+=1;
			else if (dH10->q[i]==-1) dpid->nneg+=1;
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
		}else{
			dpid->nzro+=1;
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
