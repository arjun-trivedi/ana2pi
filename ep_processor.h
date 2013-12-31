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
	
	DataH10* dH10;
	DataAna* dAna;
	
	Bool_t is_h10e1f,is_h10e16,is_h10exp,is_h10sim;
	TString h10exp,h10dtype,h10skim;
	
	EpProcessor(){
		pass = kFALSE;
		next = 0;
	}
	EpProcessor(TDirectory *td, DataH10* dataH10, DataAna* dataAna, Bool_t mon = kFALSE, Bool_t monOnly = kFALSE);
	~EpProcessor();
		
	void add(EpProcessor *n);
	void setNext(EpProcessor *n);
	Bool_t isFirstProc();
	virtual void handle();
	virtual void write();
};

#endif // EPPROCESSOR_H
