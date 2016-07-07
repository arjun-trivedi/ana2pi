#ifndef ELECTRONID_H
#define ELECTRONID_H

#include <fstream>
#include <string>
#include <TF1.h>
#include "epconfig.h"

using namespace std;

class Eid {

	static Eid* ms_instance;

public:
	Eid(char* eidParFileName);
	virtual ~Eid();
	//static Eid* Instance();
	static void Release();

	void Print();
	TObjArray* getFuncsCuts();
	Bool_t PassThreshold(Float_t p);
	Bool_t PassECFid(float uvw[3]);
	Bool_t PassSF(Int_t sector, Float_t p, Float_t sf);
	Bool_t Pass(Int_t sector, Float_t p, Float_t sf);
	TF1* Mean(Int_t sector) { return (TF1*)_fpol3Mean[sector-1]->Clone(); }
	TF1* Low(Int_t sector) { return (TF1*)_fpol3Low[sector-1]->Clone(); }
	TF1* High(Int_t sector) { return (TF1*)_fpol3High[sector-1]->Clone(); }
	//! *** Following added on [07-06-16] ***
	//! + For E16 only!
	//! + Taken from elast_lite/constants.h::namespace E16
	//! + Note that the values for these cuts are the same for exp and sim
	//! + Taken directly from h10looper_e1f.h/.cpp
	bool pass_ECin_min(int sctr, float ec_ei);
	bool pass_zvtx(int sctr, float zvtx);
	//! ******

	void DrawSFcuts(int sector);
   
    Bool_t eidParFileFound;

private:
	string _fname;
	fstream _f;
	epconfig::data _cutParMap;
	Double_t _ecThreshold;
	Double_t _Umin,_Umax,_Vmin,_Vmax,_Wmin,_Wmax;
	Double_t _mean[6][4], _sigma[6][4];
	TF1 *_fpol3Mean[6], *_fpol3Low[6], *_fpol3High[6];
	//! *** Following added on [07-06-16] ***
	//! + For E16 only!
	//! + Taken from elast_lite/constants.h::namespace E16
	//! + Note that the values for these cuts are the same for exp and sim
	Double_t _ECin_min[6];
	Double_t _zvtx_min[6];
	Double_t _zvtx_max[6];
	//! ******
	
};

#endif // ELECTRONID_H
