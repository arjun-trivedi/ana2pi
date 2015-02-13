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
	ProcCopyH10(TDirectory *td,DataH10* dataH10,DataAna* dataAna);
	~ProcCopyH10();
	
	void handle();
	void write();
	
protected:
	TTree* _tH10copy;
};

ProcCopyH10::ProcCopyH10(TDirectory *td,DataH10* dataH10,DataAna* dataAna)
						 :EpProcessor(td, dataH10, dataAna)
{
	_tH10copy = NULL;
}

ProcCopyH10::~ProcCopyH10()
{
	delete _tH10copy;
}

void ProcCopyH10::handle()
{
	//Info("ProcCopyH10::handle()", "");
	pass = kFALSE;
	if (_tH10copy==NULL){
		TDirectory* dir=(TDirectory*)dirout->GetMother();
		dir->cd();
		//dirout->cd();
		_tH10copy = (TTree*) dH10->h10chain->GetTree()->CloneTree(0);
	}
		
	_tH10copy->Fill();
	
	pass = kTRUE;
	EpProcessor::handle();
	return;
}

void ProcCopyH10::write()
{
	//dirout->cd();
	TDirectory* dir=(TDirectory*)dirout->GetMother();
	dir->cd();
	_tH10copy->Write();
}

#endif // PROCCOPYH10_H
