#ifndef EPPROCESSOR_CXX
#define EPPROCESSOR_CXX

#include "ep_processor.h"

EpProcessor::EpProcessor(TDirectory *td, DataH10* dataH10, DataAna* dataAna, Bool_t monitor, Bool_t monitorOnly /*= kFALSE*/) {
	_is_first_proc = kFALSE;
	pass = kFALSE;
	mon = monitor;
	mononly = monitorOnly;

	for(Int_t iProcMode=0;iProcMode<NPROCMODES;iProcMode++){
		for(Int_t iEvtSel=0;iEvtSel<NEVTSELS;iEvtSel++){
			//hists[iProcMode][iEvtSel] = NULL;
			//histsEkin[iProcMode][iEvtSel] = NULL;
			hists[iProcMode][iEvtSel] = new TObjArray*[NSECTORS];
			histsEkin[iProcMode][iEvtSel] = new TObjArray*[NSECTORS];
			for(Int_t iSector=0;iSector<NSECTORS;iSector++){
				hists[iProcMode][iEvtSel][iSector] = NULL;
				histsEkin[iProcMode][iEvtSel][iSector] = NULL;
			}
		}
	}
		
	if (td != NULL) dirout = td;
	else dirout = gDirectory;
	TDirectory* dirMother = (TDirectory*)dirout->GetMother();
	if (dirMother!=NULL && dirMother->GetListOfKeys()->GetEntries() == 1) _is_first_proc=kTRUE;
	if (_is_first_proc) Info("EpProcessor::EpProcessor()", "%s isFirstProc!\n", td->GetName());
				
	dH10 = dataH10;
	if (dataAna != NULL) dAna = dataAna;
	
	/* *** */
		
	_next_proc = 0;
}

EpProcessor::EpProcessor(DataH10* dataH10, DataAna* dataAna) {
	dH10 = dataH10;
	if (dataAna != NULL) dAna = dataAna;
}

EpProcessor::~EpProcessor(){
	Info("~EpProcessor()", "");
	for(Int_t iProcMode=0;iProcMode<NPROCMODES;iProcMode++){
		for(Int_t iEvtSel=0;iEvtSel<NEVTSELS;iEvtSel++){
			//dAna->deleteHists(hists[iProcMode][iEvtSel]);
			//dAna->deleteHists(histsEkin[iProcMode][iEvtSel]);
			for(Int_t iSector=0;iSector<NSECTORS;iSector++){
				delete hists[iProcMode][iEvtSel][iSector];
				delete histsEkin[iProcMode][iEvtSel][iSector];
			}
			delete hists[iProcMode][iEvtSel];
			delete histsEkin[iProcMode][iEvtSel];
		}
	}
}

void EpProcessor::add(EpProcessor *n) {
	if (_next_proc) _next_proc->add(n);
	else _next_proc = n;
}

void EpProcessor::setNext(EpProcessor *n) {
	_next_proc = n;
}

Bool_t EpProcessor::isFirstProc(){
	return _is_first_proc;
}

void EpProcessor::handle() {
	//printf("In EpProcessor::handle\n");
	if (_next_proc) _next_proc->handle();
}

void EpProcessor::write(){
	Info("EpProcessor::write()", "");
	if(hevtsum!= NULL){
		dirout->cd();
		hevtsum->Write();
	}
	TDirectory* _dirout;
	for(Int_t iProcMode=0;iProcMode<NPROCMODES;iProcMode++){
		for(Int_t iEvtSel=0;iEvtSel<NEVTSELS;iEvtSel++){
			for(Int_t iSector=0;iSector<NSECTORS;iSector++){
				if (iProcMode == MONMODE && iEvtSel == EVTINC){
					_dirout = dirout->GetDirectory(TString::Format("monitor/sector%d", iSector));
					if (_dirout != NULL) _dirout->cd();
				}else if (iProcMode == MONMODE && iEvtSel != EVTINC){
					_dirout = dirout->GetDirectory(TString::Format("monitor%d/sector%d", iEvtSel, iSector));
					if (_dirout != NULL) _dirout->cd();
				}else if (iProcMode == CUTMODE){
					_dirout = dirout->GetDirectory(TString::Format("cut/sector%d", iSector));
					if (_dirout != NULL) _dirout->cd();
				}
				if (hists[iProcMode][iEvtSel][iSector] != NULL) {
					//Info("EpProcessor::write()", "hists[%d][%d][%d]", iProcMode, iEvtSel, iSector);
					//dAna->writeHists(hists[iProcMode][iEvtSel], _dirout);
					hists[iProcMode][iEvtSel][iSector]->Write();
					//Info("EpProcessor::write()", "done hists[%d][%d][%d]", iProcMode, iEvtSel, iSector);
				}
				if (histsEkin[iProcMode][iEvtSel][iSector] != NULL)	{
					//Info("EpProcessor::write()", "histsEkin[%d][%d]", iProcMode, iEvtSel);
					//dAna->writeHists(histsEkin[iProcMode][iEvtSel], _dirout);
					histsEkin[iProcMode][iEvtSel][iSector]->Write();
					//Info("EpProcessor::write()", "done histsEkin[%d][%d]", iProcMode, iEvtSel);
				}
			}
		}
	}
}

#endif // EPPROCESSOR_H
