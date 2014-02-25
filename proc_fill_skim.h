#ifndef PROCFILLSKIM_H
#define PROCFILLSKIM_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TTree.h>

using namespace ParticleConstants;

class ProcFillSkim : public EpProcessor
{

public:
	ProcFillSkim(TDirectory *td, DataH10* dataH10,DataAna* dataAna);
	~ProcFillSkim();
	
	void handle();
	void write();
	
protected:
	TTree* _tH10skim;
};

ProcFillSkim::ProcFillSkim(TDirectory *td, DataH10* dataH10,DataAna* dataAna) : EpProcessor(td, dataH10, dataAna)
{
	_tH10skim = NULL;
}

ProcFillSkim::~ProcFillSkim()
{
	delete _tH10skim;
}

void ProcFillSkim::handle()
{
	//Info("ProcFillSkim::handle()", "");
	pass = kFALSE;
	if (_tH10skim==NULL){
		dirout->cd();
		_tH10skim = new TTree("h10","skim");
		if(dH10->expt=="e1f" && dH10->dtyp=="exp")_tH10skim->Branch("evthel",&dH10->evthel,"evthel/B");
		_tH10skim->Branch("gpart",&dH10->gpart,"gpart/I");
		_tH10skim->Branch("p",dH10->p,"p[gpart]/F");
		_tH10skim->Branch("cx",dH10->cx,"cx[gpart]/F");
		_tH10skim->Branch("cy",dH10->cy,"cy[gpart]/F");
		_tH10skim->Branch("cz",dH10->cz,"cz[gpart]/F");
		_tH10skim->Branch("id",dH10->id,"id[gpart]/S");
		_tH10skim->Branch("stat",dH10->stat,"stat[gpart]/B"); //atrivedi: 042413
		_tH10skim->Branch("q",dH10->q,"q[gpart]/B");
		_tH10skim->Branch("dc",dH10->dc,"dc[gpart]/b"); //atrivedi: 042413
		_tH10skim->Branch("cc",dH10->cc,"cc[gpart]/b"); //atrivedi: 042413
		_tH10skim->Branch("sc",dH10->sc,"sc[gpart]/b"); //atrivedi: 042413
		_tH10skim->Branch("ec",dH10->ec,"ec[gpart]/b"); //atrivedi: 042413
		
		
		_tH10skim->Branch("dc_part",&dH10->dc_part,"dc_part/I"); //atrivedi: 042413
		_tH10skim->Branch("dc_stat",dH10->dc_stat,"dc_stat[dc_part]/B"); //atrivedi: 042413
		
		
		_tH10skim->Branch("sc_part",&dH10->sc_part,"sc_part/I"); //atrivedi: 042413
		_tH10skim->Branch("sc_sect",dH10->sc_sect,"sc_sect[sc_part]/b"); //atrivedi: 042413
		_tH10skim->Branch("sc_pd",dH10->sc_pd,"sc_pd[sc_part]/b"); //atrivedi: 042413
		
		
		_tH10skim->Branch("ec_part",&dH10->ec_part,"ec_part/I"); //atrivedi: 042413
		_tH10skim->Branch("etot",dH10->etot,"etot[ec_part]/F"); //atrivedi: 042413
		
		if(dH10->dtyp=="sim"){
			_tH10skim->Branch("mcnentr",&dH10->mcnentr,"mcnentr/I");
			_tH10skim->Branch("mcid",dH10->mcid,"mcid[mcnentr]/I");
			_tH10skim->Branch("mcp",dH10->mcp,"mcp[mcnentr]/F");
			_tH10skim->Branch("mctheta",dH10->mctheta,"mctheta[mcnentr]/F");
			_tH10skim->Branch("mcphi",dH10->mcphi,"mcphi[mcnentr]/F");
			_tH10skim->Branch("mcm",dH10->mcm,"mcm[mcnentr]/F");
		}
		
	}
		
	_tH10skim->Fill();
	
	pass = kTRUE;
	EpProcessor::handle();
	return;
}

void ProcFillSkim::write()
{
	dirout->cd();
	_tH10skim->Write();
}

#endif // PROCFILLSKIM_H
