#ifndef DATA2PI_H
#define DATA2PI_H

#include <TROOT.h>
#include <TObject.h>

class Data2pi:public TObject
{

public:
	Data2pi();
	ClassDef(Data2pi, 1)
	virtual ~Data2pi();
	void Clear();
	
	//! Initial Beam Energy
	Float_t p_e0;
	//! Reconstructed Kinematics 
	//! for e',p',p,pip,pim at e' vertex
	Float_t p_e;
	Float_t p_p;
	Float_t p_pip;
	Float_t p_pim;
	Float_t theta_e;
	Float_t theta_p;
	Float_t theta_pip;
	Float_t theta_pim;
	Float_t phi_e;
	Float_t phi_p;
	Float_t phi_pip;
	Float_t phi_pim;
	//! Reconstructed e' Vertex
	Float_t vx_e;
	Float_t vx_p;
	Float_t vx_pip;
	Float_t vx_pim;
	Float_t vy_e;
	Float_t vy_p;
	Float_t vy_pip;
	Float_t vy_pim;
	Float_t vz_e;
	Float_t vz_p;
	Float_t vz_pip;
	Float_t vz_pim;	
	//! Q2, W
	Float_t Q2; 
	Float_t W;
	//! Helicity
	Char_t  h;
	//! {MMs}
	Float_t mm2ppippim;
    Float_t mmppippim;
	Float_t mm2ppip;
    Float_t mmppip;
	Float_t mm2ppim;
    Float_t mmppim;
	Float_t mm2pippim;
    Float_t mmpippim; 
    //! Topology
	Int_t top;
    //! Varsets
    Float_t M_ppip;
    Float_t M_ppim;
    Float_t M_pippim;
    Float_t theta_cms_p;
	Float_t theta_cms_pip;
	Float_t theta_cms_pim;
	Float_t phi_cms_p;
	Float_t phi_cms_pip;
	Float_t phi_cms_pim;
	Float_t alpha_1;
	Float_t alpha_2;
	Float_t alpha_3;

    /*struct vars{
		Float_t M1;
	    Float_t M2;
	    Float_t theta;
	    Float_t phi;
	    Float_t alpha;
    } varset1, varset2, varset3;*/
};

#endif // DATA2PI_H
