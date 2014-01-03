#ifndef DATAELASTIC_H
#define DATAELASTIC_H

#include <TROOT.h>
#include <TObject.h>
#include <TLorentzVector.h>
class DataElastic:public TObject
{

public:
	DataElastic();
	virtual ~DataElastic();
	void Clear();
	
	Int_t gpart;
	Int_t ne;
	Int_t np;
	Float_t Q2;
	Float_t W;
	Float_t MMp;
	TLorentzVector lvE;
	TLorentzVector lvP;
	/*Float_t p[2];
	Float_t px[2];
	Float_t py[2];
	Float_t pz[2];*/
};

#endif // DATAELASTIC_H
