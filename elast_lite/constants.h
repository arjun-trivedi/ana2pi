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

//! constants for the CC
const int CC_NSEGMENTS=18;
//! The structure should match in study_nphe_e16.py
const int CC_NPMTS=3;
enum {IPMT_L,IPMT_C,IPMT_R};

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
	//!(as per Dan Carman, MG and EP($EPCODE), and therefore as per "E1F run group"?)
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
  	//! [06-12-17]
  	//! + zvtx cut determined after Empty Target BG subtraction (etgt-bg-sub)
  	//! + This cut is only supposed to cut away the foil
  	//! + Value should be as coded in in $SUBSTUDIES/study_tgt_BG/e16/obtain_and_validate_R.py
  	//!   + Latest plots in $STUDY_TGT_BG_E16_DATADIR/results_R_<latest-date>
  	//!exp
	const Float_t ZVTX_MIN_EXP_ETGT_BG_SUB[6]={-8.00,-8.00,-8.00,-8.00,-8.00,-8.00};
	const Float_t ZVTX_MAX_EXP_ETGT_BG_SUB[6]={-0.75,-0.75,-0.75,-0.75,-0.75,-0.75};
	//! sim: different because zvtx for ER neq SR: bascically from ER by +0.5 offset by 0.5
	const Float_t ZVTX_MIN_SIM_ETGT_BG_SUB[6]={-7.50,-7.50,-7.50,-7.50,-7.50,-7.50};
  	const Float_t ZVTX_MAX_SIM_ETGT_BG_SUB[6]={-0.25,-0.25,-0.25,-0.25,-0.25,-0.25};

  	//! SF-exp
  	//! [01-08-16]
  	//! + Note that in comparison to E1F, SF_LOW/HIGH are directly entered
  	//! + This is because that in the new procedure to obtain SF cuts (sub_studies/study_eid/study_SF),
  	//!   fh(fl)=fmu+3*fsg(-3*fsg) are calculated and fitted. 
  	//! + In comparison, E1F pars were obtained very early on in my PhD. experience (~2011?) 
  	//!   where I was fitting fmu,fsg and using them to obtain fh(fl)=fmu+3*fsg(-3*fsg).
  	//! + While fitting fmu,fsg and using them to obtain fh,fl is ideal (it allows obtaining fh,fl with variable x*sg),
  	//!   given that 3-sigma cut is standard, using the linearity of the relationship between mu,sg and fh,fl, one can
  	//!   be obtained from the other
  	//! [02-17-16]
  	//! + pars taken from 'sub_studies/study_eid/study_SF/results_SFvp_e16/cutpars/exp_fullSF.txt'
  	//! [01-30-17] 
  	//! + ER and SR pars re-obtained for ananote (see Tomboy nb eid_SF:tag log=013017)
  	//! + (note SR obtained for 1st time; hitherto E1F:SR pars were beig used)
  	//! + pars taken from 'sub_studies/study_eid/study_SF/results_SFvp_e16/cutpars/exp(sim)_fullSF.txt'
	const Float_t SF_HIGH_EXP[6][4]={
		{0.399590, -0.040143, 0.011540, -0.001094},
		{0.390926, 0.001999, -0.004236, 0.000599},
		{0.409147, -0.048734, 0.020250, -0.002784},
		{0.379817, -0.008548, 0.000851, -0.000317},
		{0.404339, -0.031027, 0.008351, -0.000585},
		{0.394472, -0.020601, 0.009961, -0.001609}
	};
	const Float_t SF_LOW_EXP[6][4]={
		{0.087871, 0.127279, -0.036583, 0.003719},
		{0.105699, 0.098628, -0.022632, 0.001846},
		{0.091488, 0.134494, -0.036890, 0.003485},
		{0.097464, 0.124746, -0.036889, 0.003873},
		{0.064301, 0.145412, -0.043355, 0.004373},
		{0.077697, 0.144575, -0.041525, 0.004332}
	};
	const Float_t SF_HIGH_SIM[6][4]={
		{0.327104, -0.033060, 0.008391, -0.000890},
		{0.320327, -0.024806, 0.005451, -0.000528},
		{0.318247, -0.020546, 0.003594, -0.000292},
		{0.321015, -0.026799, 0.006448, -0.000732},
		{0.320849, -0.023275, 0.004487, -0.000347},
		{0.336825, -0.041208, 0.009857, -0.000854}
	};
	const Float_t SF_LOW_SIM[6][4]={
		{0.124180, 0.051826, -0.014285, 0.001427},
		{0.137053, 0.038839, -0.009485, 0.000890},
		{0.132904, 0.042775, -0.010605, 0.000984},
		{0.127734, 0.046799, -0.012346, 0.001185},
		{0.124138, 0.051577, -0.013662, 0.001339},
		{0.129977, 0.043062, -0.011006, 0.001069}
	};
	//! SF-sim (same as E1F)
}

#endif /* CONSTANTS_H_ */
