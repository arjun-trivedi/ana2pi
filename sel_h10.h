//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 21 05:07:54 2012 by ROOT version 5.30/04
// from TChain h10/
//////////////////////////////////////////////////////////

#ifndef SelH10_h
#define SelH10_h

#include <TROOT.h>
#include <TH1.h>
#include <TChain.h>
#include <TFile.h>
#include <TProofOutputFile.h>
#include <TDirectory.h>
#include <TSelector.h>
#include <TStopwatch.h>
#include <vector>
#include "data_h10.h"
#include "data_ana.h"
#include "ep_processor.h"

class SelH10 : public TSelector {
public :
    SelH10(TTree * /*tree*/ =0);
	virtual ~SelH10() { ; };
	virtual Int_t   Version() const { return 2;	}
	virtual void    Begin(TTree *tree);
	virtual void    SlaveBegin(TTree *tree);
	virtual void    Init(TTree *tree);
	virtual Bool_t  Notify();
	virtual Bool_t  Process(Long64_t entry);
	virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0);
	virtual void    SlaveTerminate();
	virtual void    Terminate();
	
	virtual void    SetInputList(TList *input);
	
	ClassDef(SelH10,0);

private:
	Int_t treenumber;
	Bool_t firstFileOfRun;
	Long64_t numInChain;
	Long64_t eventnum, numfilled;
	Int_t blocksize;
	TStopwatch swmain;
	TStopwatch swgroup;
	TTree *fChain; //!pointer to the analyzed TTree or TChain
	
	//DataH10
	TNamed *h10type;
	DataH10* dH10;
	Bool_t is_h10e1f,is_h10e16,is_h10exp,is_h10sim;
	TString h10exp,h10dtype,h10skim;
	
	//DataAna
	DataAna* dAna;
	
	//Procs
	vector<EpProcessor*> procs;
	EpProcessor* topProc;
	TObjString* lastProcName;
	
	TFile *fFileOut;
	TProofOutputFile *fProofFile;
	TString *fOutFileName;
	TDirectory* mkdir(const char* dirname);
	
	void SetupProcs();
};

#endif
