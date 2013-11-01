#ifndef PROCTOPS_H
#define PROCTOPS_H

#include "ep_processor.h" // Base class: EpProcessor
#include <TDirectory.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TH1.h>

class ProcTops : public EpProcessor {

public:
	ProcTops(TDirectory *td);
	virtual ~ProcTops();
	virtual void handle(DataOmega* d);
	virtual void write();
	void McKin(DataOmega *d);
	void set4vecHadrons(DataOmega* d);
	void calcMM(DataOmega* d, Bool_t ismc  = kFALSE);
	void calcVars(DataOmega* d, Bool_t ismc = kFALSE);
	Float_t getTheta(TLorentzVector lv); //angle in degrees between lv and lvQCMS
	Float_t getPhi(TLorentzVector lv);   //spherical phi angle in degrees for lv 
    Float_t invTan(Float_t y, Float_t x); //returns angle in radians [0, 2pi]
    //TObjArray *histsMM;  // for 't1a || t2a || t3a || t4a'. For technical reasons, "hists" was used in place of histsMM
    TObjArray *histseKin;  // for 't1a || t2a || t3a || t4a'
    TObjArray *yields_top0;// for 't1 || t2 || t3 || t4'
    TObjArray *yields_top1;
	TObjArray *yields_top2;
	TObjArray *yields_top3;
	TObjArray *yields_top4;
	TObjArray *histsMM_top0;// for 't1 || t2 || t3 || t4'
	TObjArray *histsMM_top1;
	TObjArray *histsMM_top2;
	TObjArray *histsMM_top3;
	TObjArray *histsMM_top4;
	TObjArray *histseKin_top0;// for 't1 || t2 || t3 || t4'
	TObjArray *histseKin_top1;
	TObjArray *histseKin_top2;
	TObjArray *histseKin_top3;
	TObjArray *histseKin_top4; 
	TObjArray *yields_mc;
	TObjArray *histsMM_mc;
	TObjArray *histseKin_mc;
	TLorentzVector lvQ;
	TLorentzVector lvW;
	TLorentzVector lvE;
	TLorentzVector lvP;
	TLorentzVector lvPip;
	TLorentzVector lvPim;
	TLorentzVector lvMMtop1;
	TLorentzVector lvMMtop2;
	TLorentzVector lvMMtop3;
	TLorentzVector lvMMtop4;
	TLorentzVector lvQCMS;
	TLorentzVector lvP0CMS;
	TLorentzVector lvPCMS;
	TLorentzVector lvPipCMS;
	TLorentzVector lvPimCMS;
protected:
	//TH1F *hevtsum;
	static const Int_t NUM_EVTCUTS = 17;
	enum { EVT_NULL, EVT, EVT_GPART0, EVT_GPARTEQ1, EVT_GPARTEQ2, EVT_GPARTEQ3, EVT_GPARTEQ4, EVT_GPART4,
	       EVT_GOODE, EVT_GOODP, EVT_GOODPIP, EVT_GOODPIM,
	       EVT_TOP0, EVT_PPIPPIM, EVT_PPIP, EVT_PPIM, EVT_PIPPIM, EVT_OTHER
	     };
};

#endif // PROCTOPS_H
