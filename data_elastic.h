#ifndef DATAELASTIC_H
#define DATAELASTIC_H

#include <TROOT.h>
#include <TObject.h>
#include <TLorentzVector.h>
//#include <TVector3.h>
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
	TVector3 vxE;
	TVector3 vxP;
	TLorentzVector lvE;
	TLorentzVector lvP;
};

#endif // DATAELASTIC_H
