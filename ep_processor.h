#ifndef EPPROCESSOR_H
#define EPPROCESSOR_H

#include <TROOT.h>
#include <TObjArray.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TLorentzVector.h> //atrivedi 020113: I am not sure why 
							//after including TLorentzVector, TObject::Info()
							//becomes accessible to EpProcessor

//using namespace TMath;
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include "data_h10.h"
#include "data_ana.h"
#include "particle_constants.h"

using namespace std;
using namespace AnalysisConstants;

class EpProcessor
{
	EpProcessor *next;
	Bool_t _isFirstProc;
public:
	static const Int_t nPROCMODE = 2;
	enum{iMODE_MON, iMODE_CUT};
	Bool_t mMon;
	Bool_t mMonOnly;
	
	Bool_t pass;
	
    TDirectory *dirout;
    TH1D* hevtsum;
	TObjArray** hists[nPROCMODE][nEVTSEL];
	TObjArray** histsEkin[nPROCMODE][nEVTSEL];
	
	DataAna* dAna;
	
	Bool_t is_h10e1f,is_h10e16,is_h10exp,is_h10sim;
	TString h10exp,h10dtype,h10skim;
	
	EpProcessor(){
		pass = kFALSE;
		next = 0;
	}
	EpProcessor(TDirectory *td, DataAna* dataAna, TString h10type, Bool_t mon = kFALSE, Bool_t monOnly = kFALSE);
	~EpProcessor();
		
	void add(EpProcessor *n);
	void setNext(EpProcessor *n);
	Bool_t isFirstProc();
	virtual void handle(DataH10* dH10);
	virtual void write();
};

//EpProcessor::EpProcessor(TDirectory *td /*= NULL*/, DataAna* dataAna /*= NULL*/, Bool_t monOnly /*= kFALSE*/) {
/*	pass = kFALSE;
	mMonOnly = monOnly;
	
	for(Int_t iProcMode=0;iProcMode<nPROCMODE;iProcMode++){
		for(Int_t iEvtSel=0;iEvtSel<nEVTSEL;iEvtSel++){
			hists[iProcMode][iEvtSel] = NULL;
			histsEkin[iProcMode][iEvtSel] = NULL;
		}
	}
		
	if (td != NULL) dirout = td;
	else dirout = gDirectory;
		
	if (dataAna != NULL) dAna = dataAna;
		
	next = 0;
}

EpProcessor::~EpProcessor(){
	for(Int_t iProcMode=0;iProcMode<nPROCMODE;iProcMode++){
		for(Int_t iEvtSel=0;iEvtSel<nEVTSEL;iEvtSel++){
			dAna->deleteHists(hists[iProcMode][iEvtSel]);
			dAna->deleteHists(histsEkin[iProcMode][iEvtSel]);
		}
	}
}

void EpProcessor::setNext(EpProcessor *n) {
	next = n;
}

void EpProcessor::add(EpProcessor *n) {
	if (next) next->add(n);
	else next = n;
}

void EpProcessor::handle(DataH10* dH10) {
	//printf("In EpProcessor::handle\n");
	if (next) next->handle(dH10);
}

void EpProcessor::write(){
	Info("EpProcessor::write()", "");
	if(hevtsum!= NULL){
		dirout->cd();
		hevtsum->Write();
	}
	TDirectory* _dirout;
	for(Int_t iProcMode=0;iProcMode<nPROCMODE;iProcMode++){
		for(Int_t iEvtSel=0;iEvtSel<nEVTSEL;iEvtSel++){
			//for(Int_t iSector=0;iSector<nSECTOR;iSector++){
				if (iProcMode == iMODE_MON) {
					if(iEvtSel == iEVTINC) _dirout = dirout->GetDirectory(TString::Format("mon"));
					else _dirout = dirout->GetDirectory(TString::Format("mon%d", iEvtSel));
				}else if (iProcMode == iMODE_CUT){
					_dirout = dirout->GetDirectory(TString::Format("cut"));
				}
			
				if (hists[iProcMode][iEvtSel] != NULL) {
					dAna->writeHists(hists[iProcMode][iEvtSel], _dirout);
					dAna->writeHists(histsEkin[iProcMode][iEvtSel], _dirout);
				}
			//}
		}
	}
}*/

#endif // EPPROCESSOR_H
