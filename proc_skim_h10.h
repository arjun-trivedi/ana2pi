#ifndef PROCFILLSKIM_H
#define PROCFILLSKIM_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TTree.h>

using namespace ParticleConstants;

/**********************************
/* [02-15-15]
/* This processer is used to keep on selected Branches of a h10
/* Tree depending on the need of the analysis. 
***********************************/

class ProcSkimH10 : public EpProcessor
{

public:
	ProcSkimH10(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcSkimH10();
	
	void handle();
	void write();
	
protected:
	TTree* _tH10skim;
	static const Int_t NUM_EVTCUTS=2;
	enum {EVT_NULL, EVT,EVT_PASS};
};

ProcSkimH10::ProcSkimH10(TDirectory *td,DataH10* dataH10,DataAna* dataAna)
						 : EpProcessor(td, dataH10, dataAna)
{
	_tH10skim = NULL;
	dirout->cd();
	hevtsum=new TH1D("hevtsum","Event Statistics",NUM_EVTCUTS,0.5,NUM_EVTCUTS+0.5);
}

ProcSkimH10::~ProcSkimH10()
{
	delete _tH10skim;
}

void ProcSkimH10::handle()
{
	//Info("ProcSkimH10::handle()", "");
	hevtsum->Fill(EVT);
	pass = kFALSE;
	if (_tH10skim==NULL){
		//dirout->cd();
		TDirectory* dir=(TDirectory*)dirout->GetMother();
		dir->cd();
		_tH10skim = new TTree("h10","skim");
		//! Branches needed to make d2pi
		if(dH10->expt=="e1f" && dH10->dtyp=="exp"){
			_tH10skim->Branch("evthel",&dH10->evthel,"evthel/F");
		}
		_tH10skim->Branch("gpart",  &dH10->gpart,"gpart/I");
		_tH10skim->Branch("p",      dH10->p,"p[gpart]/F");
		_tH10skim->Branch("cx",     dH10->cx,"cx[gpart]/F");
		_tH10skim->Branch("cy",     dH10->cy,"cy[gpart]/F");
		_tH10skim->Branch("cz",     dH10->cz,"cz[gpart]/F");
		_tH10skim->Branch("id",     dH10->id,"id[gpart]/S");
		_tH10skim->Branch("stat",   dH10->stat,"stat[gpart]/B"); //atrivedi: 042413
		_tH10skim->Branch("q",      dH10->q,"q[gpart]/B");
		_tH10skim->Branch("dc",     dH10->dc,"dc[gpart]/b"); //atrivedi: 042413
		_tH10skim->Branch("cc",     dH10->cc,"cc[gpart]/b"); //atrivedi: 042413
		_tH10skim->Branch("sc",     dH10->sc,"sc[gpart]/b"); //atrivedi: 042413
		_tH10skim->Branch("ec",     dH10->ec,"ec[gpart]/b"); //atrivedi: 042413
		_tH10skim->Branch("dc_part",&dH10->dc_part,"dc_part/I"); //atrivedi: 042413
		_tH10skim->Branch("dc_stat",dH10->dc_stat,"dc_stat[dc_part]/B"); //atrivedi: 042413
		_tH10skim->Branch("sc_part",&dH10->sc_part,"sc_part/I"); //atrivedi: 042413
		_tH10skim->Branch("sc_sect",dH10->sc_sect,"sc_sect[sc_part]/b"); //atrivedi: 042413
		_tH10skim->Branch("sc_pd",  dH10->sc_pd,"sc_pd[sc_part]/b"); //atrivedi: 042413
		_tH10skim->Branch("ec_part",&dH10->ec_part,"ec_part/I"); //atrivedi: 042413
		//! Branches needed to calculate Luminosity (currently made)
		_tH10skim->Branch("q_l",    &dH10->q_l,"q_l/F");
		_tH10skim->Branch("t_l",    &dH10->t_l,"t_l/F");
		//! 
		/*_tH10skim->Branch("tr_time",&dH10->vz,"tr_time/F");
		_tH10skim->Branch("etot",   dH10->etot,"etot[ec_part]/F"); //atrivedi: 042413
		_tH10skim->Branch("ec_eo",  dH10->ec_eo,"ec_eo[ec_part]/F"); 
		_tH10skim->Branch("ec_ei",  dH10->ec_ei,"ec_ei[ec_part]/F");
		_tH10skim->Branch("vx",     dH10->vx,"vx[gpart]/F");
		_tH10skim->Branch("vy",     dH10->vy,"vy[gpart]/F");
		_tH10skim->Branch("vz",     dH10->vz,"vz[gpart]/F");*/
		
		if(dH10->expt=="e1f" && dH10->dtyp=="sim"){
			_tH10skim->Branch("mcnentr",&dH10->mcnentr,"mcnentr/I");
			_tH10skim->Branch("mcid",   dH10->mcid,"mcid[mcnentr]/I");
			_tH10skim->Branch("mcp",    dH10->mcp,"mcp[mcnentr]/F");
			_tH10skim->Branch("mctheta",dH10->mctheta,"mctheta[mcnentr]/F");
			_tH10skim->Branch("mcphi",  dH10->mcphi,"mcphi[mcnentr]/F");
			_tH10skim->Branch("mcm",    dH10->mcm,"mcm[mcnentr]/F");
		}
		
	}
		
	_tH10skim->Fill();
	
	pass = kTRUE;
	hevtsum->Fill(EVT);
	EpProcessor::handle();
	return;
}

void ProcSkimH10::write()
{
	Info("ProcSkimH10::write()", "");
	//dirout->cd();
	TDirectory* dir=(TDirectory*)dirout->GetMother();
	dir->cd();
	_tH10skim->Write();
}

#endif // PROCFILLSKIM_H
