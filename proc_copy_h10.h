#ifndef PROCCOPYH10_H
#define PROCCOPYH10_H

#include "ep_processor.h" // Base class: EpProcessor
#include "data_h10.h"
#include "particle_constants.h"
#include <TTree.h>

using namespace ParticleConstants;

class ProcCopyH10 : public EpProcessor
{

public:
	ProcCopyH10(TDirectory *td, DataAna* dataAna, TString h10type);
	~ProcCopyH10();
	
	void handle(DataH10* dH10);
	void write();
	
protected:
	TTree* _tH10copy;
};

ProcCopyH10::ProcCopyH10(TDirectory *td, DataAna* dataAna, TString h10type) : EpProcessor(td, dataAna, h10type)
{
	_tH10copy = NULL;
}

ProcCopyH10::~ProcCopyH10()
{
	delete _tH10copy;
}

void ProcCopyH10::handle(DataH10* dH10)
{
	//Info("ProcCopyH10::handle()", "");
	pass = kFALSE;
	if (_tH10copy==NULL){
		dirout->cd();
		_tH10copy = (TTree*) dH10->fChain->GetTree()->CloneTree(0);
	}
		
	_tH10copy->Fill();
	
	pass = kTRUE;
	EpProcessor::handle(dH10);
	return;
}

void ProcCopyH10::write()
{
	dirout->cd();
	_tH10copy->Write();
}

#endif // PROCCOPYH10_H
