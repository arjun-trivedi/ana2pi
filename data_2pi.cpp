#include "data_2pi.h"

ClassImp(Data2pi)

Data2pi::Data2pi()
{
}

Data2pi::~Data2pi()
{
}

void Data2pi::Clear()
{
	//! Initial Beam Energy
	p_e0 = 0;
	//! From DC: sector information for detected particles
	sector_e = sector_p = sector_pip = sector_pim=0;
	//! Reconstructed Kinematics 
	//! for e',p',p,pip,pim at e' vertex
	p_e     = p_p     = p_pip     = p_pim     = 0;
	theta_e = theta_p = theta_pip = theta_pim = 0;
	phi_e =   phi_p   = phi_pip   = phi_pim   = 0;
	//! Reconstructed e' Vertex
	vx_e = vx_p = vx_pip = vx_pim = 0;
	vy_e = vy_p = vy_pip = vy_pim = 0;
	vz_e = vz_p = vz_pip = vz_pim = 0;
	//! Q2, W
	Q2 = W = 0;
	//! Helicity
	h = 0;
	//! {MMs}
	mm2ppippim = mmppippim  = mm2ppip = mmppip= mm2ppim = mmppim = mm2pippim =  mmpippim = 0;
	//! Topology
	top=0;
	//! Varsets
	M_ppip      = M_ppim        = M_pippim      = 0;
	theta_cms_p = theta_cms_pip = theta_cms_pim = 0;
	phi_cms_p   = phi_cms_pip   = phi_cms_pim   = 0;
	alpha_1     = alpha_2       = alpha_3       = 0;

	/*varset1.M1 = varset1.M2 = varset1.theta_cms = varset1.phi_ = varset1.alpha = 0;
	varset2.M1 = varset2.M2 = varset2.theta_cms = varset2.phi_ = varset2.alpha = 0;
	varset3.M1 = varset3.M2 = varset3.theta_cms = varset3.phi_ = varset3.alpha = 0;*/
}
