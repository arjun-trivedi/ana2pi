#ifndef DATAINCLSV_H
#define DATAINCLSV_H

#include <TROOT.h>
#include <TObject.h>
#include <TLorentzVector.h>
//#include <TVector3.h>
class DataInclsv:public TObject
{

public:
	DataInclsv();
	virtual ~DataInclsv();
	void Clear();
	
	Int_t gpart;
	//!for eid
	Int_t nphe;
	Float_t etot;
	Float_t ec_ei;
	Float_t ec_eo;
	//! e- kinematics
	Float_t p;
	Float_t theta;
	Float_t phi;
	Float_t Q2;
	Float_t nu;
	Float_t q;
	Float_t x; //Bjorken-x
	Float_t W;
};

#endif // DATAINCLSV_H
