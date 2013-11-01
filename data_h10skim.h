#ifndef DATAH10SKIM
#define DATAH10SKIM

#include <TROOT.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include "data_h10.h"

class DataH10skim : public DataH10
{

public:
	DataH10skim();
	~DataH10skim();
	void Bind(TTree* fChain);
	void Reconcile();
	void Clear();
	
	bool is_e1f;
	bool is_e16;
	const TLorentzVector _4vE0;
	const TLorentzVector _4vP0;

	//static const int _MAX_PARTS = 24; //atrivedi 031313: works only for exp. data
	static const int _MAX_PARTS = 40;   //atrivedi 031313: Changed since: gpart(sc_part etc) distribution = [0,40] for e1f sim. data 
	
	Char_t evthel;
	Int_t gpart;
	Short_t id[_MAX_PARTS]; //[gpart]
	Float_t p[_MAX_PARTS]; //[gpart]
	Float_t cx[_MAX_PARTS]; //[gpart]
	Float_t cy[_MAX_PARTS]; //[gpart]
	Float_t cz[_MAX_PARTS]; //[gpart]
		
	
	/* ***** EXP SPECIFIC BRANCHES ***** */
	TBranch *b_evthel;

	/* ***** COMMON BRANCHES ****** */
	TBranch *b_gpart; //!
	TBranch *b_id; //!
	TBranch *b_p; //!
	TBranch *b_cx; //!
	TBranch *b_cy; //!
	TBranch *b_cz; //!
};

#endif // DATAH10SKIM
