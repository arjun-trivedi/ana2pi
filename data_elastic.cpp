#include "data_elastic.h"

DataElastic::DataElastic()
{
}

DataElastic::~DataElastic()
{
}

void DataElastic::Clear()
{
	//! gpart
	gpart=-1;
	//!Q2,W
	Q2=0;
	W=0;
	//! MMs for Event Selection
	MMep=0;
	MM2ep=0;
	//! Reconstructed Kinematics
	//! for e',p' at e' vertex
	p_e=0;
	p_p=0;
	theta_e=0;
	theta_p=0;
	phi_e=0;
	phi_p=0;
	//! Reconstructed e' vertex
	vx_e=0;
	vx_p=0;
	vy_e=0;
	vy_p=0;
	vz_e=0;
	vz_p=0;
	//! EID
	nphe=0;
}
