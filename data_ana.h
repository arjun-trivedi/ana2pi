#ifndef DATAANA_H
#define DATAANA_H


#include "data_eid.h"
#include "data_efid.h"
#include "data_skim_q.h"
#include "data_mom.h"
#include "data_pid.h"
#include "data_ekin.h"
#include "data_top.h"
#include "data_h10.h"
#include <TObjArray.h>

class DataAna
{
public:
	DataAna();
	virtual ~DataAna();
	void Clear();
	Int_t opart;
	Int_t h10idxE;
	Int_t h10idxP;
	Int_t h10idxPip;
	Int_t h10idxPim;
	Int_t top;
	
	DataEid eid;
	DataEFid efid;
	DataSkimQ skimq;
	DataPid pid;
	DataMom mom;
	DataEkin eKin;
	DataEkin eKin_mc;
	DataTop dTop;
	DataTop dTop_mc;
	
	struct h8Dbng{
		Int_t bins;
		Double_t xmin;
		Double_t xmax;
	} bngQ2, bngW, bngMppip, bngMppim, bngMpippim, bngTheta, bngPhi, bngAlpha;
	
	void makeHistsEid(TObjArray** hists, TDirectory* dirout);
	void makeHistsEFid(TObjArray** hists, TDirectory* dirout);
	void makeHistsMomCor(TObjArray** hists, TDirectory* dirout);
	void makeHistsPid(TObjArray** hists, TDirectory* dirout);
	void makeHistsEkin(TObjArray** hists, TDirectory* dirout);
	
	TObjArray* makeHistsEkin();
	TObjArray* makeHistsMM();
	TObjArray* makeYields();
	
	void writeHists(TObjArray** hists, TDirectory *dirout);
	
	void deleteHists(TObjArray** hists);
	
	void fillHistsEid(TObjArray** hists, Bool_t useMc = kFALSE);
    void fillHistsEFid(TObjArray** hists, Bool_t useMc = kFALSE);
    void fillHistsMomCor(TObjArray** hists, Bool_t useMc = kFALSE);
	void fillHistsPid(TObjArray** hists, Bool_t useMc = kFALSE);
	void fillHistsEkin(TObjArray** hists, Bool_t useMc = kFALSE);
	
	void fillHistsEkin(TObjArray* hists, Bool_t useMc = kFALSE);
	void fillHistsMM(TObjArray *hists, Bool_t useMc = kFALSE);
	void fillYields(TObjArray *hists, Bool_t useMc = kFALSE);
};

#endif // DATAANA_H
