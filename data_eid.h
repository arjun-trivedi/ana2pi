#ifndef DATAEID_H
#define DATAEID_H

#include <TROOT.h>

class DataEid
{
public:
	DataEid();
	virtual ~DataEid();
	void Clear();
	//! From cooking
	Int_t id;
	//! From DC
	Float_t p;
	Float_t dc_xsc;
    Float_t dc_ysc;
    Float_t dc_zsc;
    Float_t vz;
	//! From SC
	Int_t sector;
    Float_t b;
    Float_t b_e;
    Float_t dt_e;
    //from EC
	Float_t ec_ei;
	Float_t ec_eo;
	Float_t etot;
	Float_t ech_x;
	Float_t ech_y;
	Float_t ech_z;
	Float_t ecU;
	Float_t ecV;
	Float_t ecW;
	//from CC
	Int_t cc;
	Int_t nphe;
	Int_t cc_segm;
	Int_t pmt;
	Float_t cc_theta;
};

#endif // DATAEID_H
