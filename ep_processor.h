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
	EpProcessor *_next_proc;
	Bool_t _is_first_proc;
public:
	static const Int_t NPROCMODES = 2;
	enum{MONMODE, CUTMODE};
	Bool_t mon;
	Bool_t mononly;
	
	Bool_t pass;
	
    TDirectory *dirout;
    TH1D* hevtsum;
	TObjArray** hists[NPROCMODES][NEVTSELS];
	TObjArray** histsEkin[NPROCMODES][NEVTSELS];
	
	DataH10* dH10;
	DataAna* dAna;
	
	/*Bool_t is_h10e1f,is_h10e16,is_h10exp,is_h10sim;
	TString h10exp,h10dtype,h10skim;*/
	
	EpProcessor(){
		pass = kFALSE;
		_next_proc = 0;
	}
	EpProcessor(TDirectory *td, DataH10* dataH10, DataAna* dataAna, Bool_t monitor = kFALSE, Bool_t monitorOnly = kFALSE);
	~EpProcessor();
		
	void add(EpProcessor *n);
	void setNext(EpProcessor *n);
	Bool_t isFirstProc();
	virtual void handle();
	virtual void write();
};

#endif // EPPROCESSOR_H
