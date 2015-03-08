#ifndef DATAANA_H
#define DATAANA_H


#include "data_eid.h"
#include "data_efid.h"
#include "data_skim_q.h"
#include "data_skim_q_elast.h"
#include "data_mom.h"
#include "data_pid.h"
#include "data_pid_elast.h"
#include "data_ekin.h"
#include "data_2pi.h"
#include "data_h10.h"
#include "data_elastic.h"
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
	DataSkimQElast skimq_elast;
	DataPid pid;
	DataPidElast pid_elast;
	DataMom mom;
	DataEkin eKin;
	DataEkin eKin_mc;
	Data2pi d2pi;
	Data2pi d2pi_mc;
	DataElastic dElast;
	DataElastic dElast_ST;
	
	struct h8Dbng{
		Int_t bins;
		Double_t xmin;
		Double_t xmax;
	} bngQ2, bngW, bngMppip, bngMppim, bngMpippim, bngTheta, bngPhi, bngAlpha;
	
	void makeHistsEid(TObjArray** hists, TDirectory* dirout);
	void makeHistsEFid(TObjArray** hists, TDirectory* dirout);
	void makeHistsMomCor(TObjArray** hists, TDirectory* dirout);
	void makeHistsPid(TObjArray** hists, TDirectory* dirout);
	void makeHistsPidElast(TObjArray** hists, TDirectory* dirout);
	void makeHistsEkin(TObjArray** hists, TDirectory* dirout);
	
	TObjArray* makeHistsEkin();
	TObjArray* makeHistsMM();
	TObjArray** makeYields();

	TObjArray* makeHistsMMElastic();
	TObjArray* makeYieldsElastic();
	
	void writeHists(TObjArray** hists, TDirectory *dirout);
	
	void deleteHists(TObjArray** hists);
	
	void fillHistsEid(TObjArray** hists, Bool_t useMc = kFALSE);
    void fillHistsEFid(TObjArray** hists, Bool_t useMc = kFALSE);
    void fillHistsMomCor(TObjArray** hists, Bool_t useMc = kFALSE);
	void fillHistsPid(TObjArray** hists, Bool_t useMc = kFALSE);
	void fillHistsPidElast(TObjArray** hists, Bool_t useMc = kFALSE);
	void fillHistsEkin(TObjArray** hists, Bool_t useMc = kFALSE);
	
	void fillHistsEkin(TObjArray* hists, Bool_t useMc = kFALSE);
	void fillHistsMM(TObjArray *hists, Bool_t useMc = kFALSE);
	void fillYields(TObjArray** hists, Float_t w, Bool_t useMc = kFALSE);

	void fillHistsMMElastic(TObjArray *hists, Bool_t useMc = kFALSE);
	void fillYieldsElastic(TObjArray* hists, Bool_t useMc = kFALSE);

	void addBranches_DataElastic(TTree* t, Bool_t useMc=kFALSE);
	void addBranches_DataEid(TTree* t);
	void addBranches_DataPid(TTree* t);
	void addBranches_DataPidElast(TTree* t);
};

#endif // DATAANA_H
