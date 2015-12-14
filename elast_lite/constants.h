#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const Float_t SOL = 29.9792458;

//particle codes, usually PDG codes, but always those used in BOS
const Int_t PROTON = 2212;
const Int_t NEUTRON = 2112;
const Int_t PIP = 211;
const Int_t PIM = -211;
const Int_t PI0 = 111;
const Int_t KP = 321;
const Int_t KM = -321;
const Int_t PHOTON = 22;
const Int_t ELECTRON = 11;

//PDG particle masses in GeV/c2
const Float_t MASS_P = 0.93827203;
const Float_t MASS_N = 0.93956556;
const Float_t MASS_E = 0.000511;
const Float_t MASS_PIP = 0.13957018;
const Float_t MASS_PIM = 0.13957018;
const Float_t MASS_PI0 = 0.1349766;
const Float_t MASS_KP = 0.493677;
const Float_t MASS_KM = 0.493677;
const Float_t MASS_G = 0.0;
const Float_t MASS_OMEGA = 0.78265; //

//! CLAS sector phi angle coverage
const Int_t PHI_FULL_SECTOR[6][2]={
	{-30,30},
	{30,90},
	{90,150},
	{150,210},
	{210,270},
	{270,330}
};

const Int_t PHI_CENTRAL_SECTOR[6][2]={
	{0,2},
	{60,62},
	{120,122},
	{182,184},
	{242,244},
	{302,304}
};


namespace E1F{
	const Float_t E0_P = 5.499;

	//! eid-cut pars
	//! minimum momentum that satisfies EC threshold
	const Float_t P_MIN_ECTH = 0.64;

	//! ECin cut exp (taken from $EPCODE)
	const Float_t ECIN_MIN_EXP[6]={0.058,0.064,0.060,0.056,0.058,0.056};  
	//! ECin cut sim (taken from $EPCODE)
	const Float_t ECIN_MIN_SIM[6]={0.063,0.063,0.063,0.063,0.063,0.063}; 

	//! EC fid cut pars
	//!(as per MG and EP($EPCODE), and therefore as per "E1F run group"?)
	const Float_t UMIN=20; 
	const Float_t UMAX=400;
	const Float_t VMIN=0;
	const Float_t VMAX=375;
	const Float_t WMIN=0;
	const Float_t WMAX=410;

	//! z-vertex cut
	//! + (-27.5,-22.5) cut from Kijun's ana note (sec. 4.2.1)
	//! + The value above is adapted for exp sectors where there is a mismatch
	//! (Kijun and Isupov correct this mismatch, but I for now have decided to
	//! adapt the cut instead. See sub_studies/study_vrtx/plot_e_z_vtx.C)
	//!exp
	const Float_t ZVTX_MIN_EXP[6]={-28.25,-27.50,-27.50,-27.50,-28.25,-28.75};
	const Float_t ZVTX_MAX_EXP[6]={-23.00,-22.50,-22.25,-22.50,-23.00,-23.50};
	//! sim
	const Float_t ZVTX_MIN_SIM[6]={-27.50,-27.50,-27.50,-27.50,-27.50,-27.50};
  	const Float_t ZVTX_MAX_SIM[6]={-22.50,-22.50,-22.50,-22.50,-22.50,-22.50};

	//! SF-exp
	const Float_t SF_MEAN_EXP[6][4]={
		{0.259351, 0.0290101,-0.002351,  -1.04206e-05},
		{0.281727, 0.0579125,-0.0110491,  0.000745872},
		{0.258238, 0.0469234,-0.00793056, 0.000672446},
		{0.274411, 0.0182674, 0.00222566,-0.000565052},
		{0.239176, 0.0444721,-0.00890621, 0.000707731},
		{0.248752, 0.0356904,-0.00505897, 0.000103396}
	};
	const Float_t SF_SIGMA_EXP[6][4]={
		{0.0470147,-0.0191508, 0.00515491,-0.000496982},
		{0.058316, -0.0118746,-0.000340192,0.000365358},
		{0.0477462,-0.0172543, 0.00409656,-0.000289126},
		{0.0443633,-0.0101284, 0.00287309,-0.000345008},
		{0.049301, -0.0173756, 0.00400414,-0.000309591},
		{0.04315,  -0.0104406, 0.00127932,-3.75922e-05}
	};
	//! SF-sim
	const Float_t SF_MEAN_SIM[6][4]={
		{0.220766, 0.0147447, -0.00489355, 0.000519253},
		{0.219945, 0.0180767, -0.00639825, 0.000742092},
		{0.220313, 0.0169074, -0.00548383, 0.000568149},
		{0.221288, 0.0137292, -0.00437765, 0.000424327},
		{0.219234, 0.0177388, -0.00557197, 0.000555405},
		{0.217248, 0.019857,  -0.0075526,  0.000927366}
	};
	const Float_t SF_SIGMA_SIM[6][4]={
		{0.0354434, -0.0143688, 0.00352196, -0.000335203},
		{0.0349228, -0.0145381, 0.00368543, -0.000359392},
		{0.0347707, -0.0135492, 0.00307874, -0.000260751},
		{0.036096,  -0.0153546, 0.00394917, -0.000390688},
		{0.0353282, -0.0142703, 0.00347418, -0.000323637},
		{0.0355262, -0.0138957, 0.00315626, -0.000268765}
	};
}

namespace E16{
	const Float_t E0_P = 5.754;

	//! eid-cut pars
	//! minimum momentum that satisfies EC threshold
	//! (as per EI analysis note)
	const Float_t P_MIN_ECTH = 0.70;

	//! ECin cut exp : taken from EI's analysis note and email of [12-11-15] 
	const Float_t ECIN_MIN_EXP[6]={0.06,0.06,0.06,0.06,0.06,0.06};  
	//! ECin cut sim : taken from EI's analysis note and email of [12-11-15] 
	const Float_t ECIN_MIN_SIM[6]={0.06,0.06,0.06,0.06,0.06,0.06}; 

	//! EC fid cut pars
	//!(as per EI analysis note)
	const Float_t UMIN=40; 
	const Float_t UMAX=9999; //! no cut on UMAX 
	const Float_t VMIN=0;
	const Float_t VMAX=360;
	const Float_t WMIN=0;
	const Float_t WMAX=390;

	//! z-vertex cut
	//! [12-13-15]
	//! + sim: values as per EI; his note and email of [12-11-15]
	//! + exp: currently not applied and set to same as sim, but will have to verified
	//!        after apply z-vertex corrections
	//!exp
	const Float_t ZVTX_MIN_EXP[6]={-8.0,-8.0,-8.0,-8.0,-8.0,-8.0};
	const Float_t ZVTX_MAX_EXP[6]={-0.8,-0.8,-0.8,-0.8,-0.8,-0.8};
	//! sim
	const Float_t ZVTX_MIN_SIM[6]={-8.0,-8.0,-8.0,-8.0,-8.0,-8.0};
  	const Float_t ZVTX_MAX_SIM[6]={-0.8,-0.8,-0.8,-0.8,-0.8,-0.8};
}

#endif /* CONSTANTS_H_ */
