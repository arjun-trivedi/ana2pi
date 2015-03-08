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
	
	//! gpart
	Int_t gpart;
	//!Q2,W
	Float_t Q2;
	Float_t W;
	//! MMs for Event Selection
	Float_t MMep;
	Float_t MM2ep;
	//! Reconstructed Kinematics
	//! for e',p' at e' vertex
	Float_t p_e;
	Float_t p_p;
	Float_t theta_e;
	Float_t theta_p;
	Float_t phi_e;
	Float_t phi_p;
	//! Reconstructed e' vertex
	Float_t vx_e;
	Float_t vx_p;
	Float_t vy_e;
	Float_t vy_p;
	Float_t vz_e;
	Float_t vz_p;
};

#endif // DATAELASTIC_H
