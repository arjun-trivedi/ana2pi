#ifndef DATANSTRDIS_H
#define DATANSTRDIS_H

#include <TROOT.h>
#include <TObject.h>
#include <TLorentzVector.h>
//#include <TVector3.h>
class DataNstrDIS:public TObject
{

public:
	DataNstrDIS();
	virtual ~DataNstrDIS();
	void Clear();
	
	/*Int_t gpart;
	Int_t ne;
	Int_t np;*/
	//Int_t sector;
	Float_t Q2;
	Float_t W;
	Float_t MMep;
	/*TVector3 vxE;
	TVector3 vxP;
	TLorentzVector lvE;
	TLorentzVector lvP;*/
};

#endif // DATANSTRDIS_H
