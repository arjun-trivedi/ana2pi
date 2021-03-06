#ifndef CUTS_H_
#define CUTS_H_

#include <TMath.h>
#include <TF1.h>
#include <fstream>
#include <vector>
#include "particle_constants.h"
#include "epconfig.h"

using namespace ParticleConstants;

/** Contains static functions and constants for all cuts, currently for E1F data only. All parameters are predetermined by substudies.
 * 	
 *	Cuts include the following:
 * 	<ul>
 * 		<li>momentum thresholds</li>
 * 		<li>target position</li>
 * 		<li>vertex restrictions</li>
 * 		<li>missing mass skims</li>
 * 		<li>fiducial</li>
 * 	</ul>
 * 
 * Substudy modules usually consist of a batch-mode tree-maker, which makes TTrees for on-the-fly analysis within ROOT, and an analyzer, executable from within ROOT, usually within one class to make it easier for others to copy and customize.
 * <ul>
 * 	<li>@see h10fid to aid in extraction of fiducial cut parameters.</li>
 * </ul>
 * 
 * 	Notes:
 * 	<ul>
 * 		<li>Consider removing ROOT dependencies.</li>
 * 		<li>Make parameters configurable either by config file or setter methods. This would require refactoring of ParticleConstants as well, since it contains experiment-specific values, e.g., beam energy.</li>
 * 	</ul>
 */
namespace Cuts {

/* MOMENTUM THRESHOLDS */
  static const Float_t P_MIN_E = 0.683; //!< Electron low momentum threshold. For E1F 683 MeV, corresponding to 172 mV in EC.

/* VERTEX CUT VALUES */
  static const Float_t VX_ZMIN = -27.75; //!< Empirically-determined target downstream edge z-position in cm, for E1F.
  static const Float_t VX_ZMAX = -22.25; //!< Empirically-determined target upstream edge z-position in cm, for E1F.
  static const Float_t VX_DZ = 1.6; //!< Vertex resolution in cm, for E1F.

/* MISSING MASS SKIM REGIONS */
  static const Float_t MMP_MIN = 0.652; //!< Low boundary of omega mass skim region. gamma+p0-p1.
  static const Float_t MMP_MAX = 0.912; //!< High boundary of omega mass skim region. gamma+p0-p1.
  static const Float_t MMPPI_MIN = 0.240; //!< Low boundary of two-pion mass skim region. gamma+p0-(p1+piX).
  static const Float_t MMP2PI_MIN = 0.050; //!< Low boundary of one-pion mass skim region. gamma+p0-(p1+pip+pim).
  static const Float_t MMP2PI_MAX = 0.270; //!< High boundary of one-pion mass skim region. gamma+p0-(p1+pip+pim).

/* FIDUCIAL PARAMETER FUNCTIONS */
  static TF1 *fDphi0 = NULL;	//!< First parameter, with added momentum-dependence, of <i>standard</i> fiducial parametrization.
  static TF1 *fDphi1 = NULL;	//!< Second parameter, with added momentum-dependence, of <i>standard</i> fiducial parametrization.
  static TF1 *fDphi2 = NULL;	//!< Third parameter, with added momentum-dependence, of <i>standard</i> fiducial parametrization.
  static TF1 *fTheta0 = NULL;	//!< Momentum-dependent minimum theta function.
	static TF1 *fTheta0_e = NULL;
	static TF1 *fTheta0_p = NULL;
	
	/** Electron momentum threshold cut.
	 * @param pE Uncorrected electron momentum.
	 * @return True if above threshold; false otherwise.
	 * @see P_MIN_E
	 */
  static bool Threshold_pE(Float_t pE) {
    return pE>P_MIN_E;
  }

	/** Requires vertex position of particle to be within target boundaries.
	 * @param vxZ Corrected vertex z-position.
	 * @return True if within target boundaries; false otherwise.
	 * @see VX_ZMIN
	 * @see VX_ZMAX
	 */
  static bool Vertex(Float_t vxZ) {
    return ( vxZ > VX_ZMIN && vxZ < VX_ZMAX );
  }

	/** Requires vertex position events particles to be within vertex region of electron.
	 * @param vxZ_e Corrected vertex z-position of electron.
	 * @param vxZ_p Corrected vertex z-position of other particle.
	 * @return True if within vertex resolution; false otherwise.
	 * @see VX_DZ
	 */
  static bool Vertex(Float_t vxZ_e, Float_t vxZ_p) {
    return ( TMath::Abs(vxZ_e-vxZ_p) < VX_DZ );
  }

	/** Returns a value that encodes the skims satisfied by event.
	 * @param mmp gamma+p0-p1
	 * @param mmppip gamma+p0-(p1+piX)
	 * @param mmppim gamma+p0-(p1+piX)
	 * @param mmppippim gamma+p0-(p1+pip+pim)
	 * @return Integer between 0 and 7, an encoding of passed skims:
	 * <ul>
	 * 	<li>Lowest order bit -- mmppip &gt; MMPPI_MIN</li>
	 * 	<li>Middle order bit -- mmppim &gt; MMPPI_MIN</li>
	 * 	<li>Highest order bit -- mmppippim between MMP2PI_MIN and MMP2PI_MAX</li>
	 * </ul>
	 * @see MMP_MIN
	 * @see MMP_MAX
	 * @see MMPPI_MIN
	 * @see MMP2PI_MIN
	 * @see MMP2PI_MAX
	 */
  static Int_t MM(Float_t mmp, Float_t mmppip, Float_t mmppim, Float_t mmppippim) {
    Int_t retval = 0;
    if (mmppip>MMPPI_MIN && mmp>MMP_MIN && mmp<MMP_MAX) retval+=1;
    if (mmppim>MMPPI_MIN && mmp>MMP_MIN && mmp<MMP_MAX) retval+=2;
    if (mmppippim>MMP2PI_MIN && mmppippim<MMP2PI_MAX && mmp>MMP_MIN && mmp<MMP_MAX) retval+=4;
    return retval;
  }

	/** Returns TF1 function for minimum theta as a function of momentum.
	 * 	@param id Particle ID.
	 * 
	 * to-do:
	 * <ul>
	 * 	<li>Add a tightness parameter?</li>
	 * </ul>
	 */
  static TF1* GetTheta0F(Int_t id) {
		Double_t par0a, par0b, par0c;
		//if (fTheta0_e == NULL) {
		//	fTheta0_e = new TF1("fTheta0","[0]+[1]/((x+[2])*3375.0/2250.0)",0,6);
		//}
		if (id == ELECTRON) {
			if (fTheta0_e == NULL) {
				fTheta0_e = new TF1("fTheta0_e","[0]+[1]/((x+[2])*3375.0/2250.0)",0,6);
			}
			fTheta0 = fTheta0_e;
			par0a = 11.0; //10.7
			par0b = 17.6;
			par0c = 0.045;
		} else if (id == PROTON || id == PIP) {
			if (fTheta0_p == NULL) {
				// thmin = 8.  + 20.*(1.-p_h*rat/8.)**15.
				fTheta0_p = new TF1("fTheta0_p","[0]+[1]*(1-x*3375.0/2250.0/8)**[2]",0,6);
			}
			fTheta0 = fTheta0_p;
			//if (id == PROTON)	par0a = 6;
			//else par0a = 8;
			par0a = 8;
			par0b = 20;
			par0c = 15;
			//par0a = 7;
			//par0b = 8;
			//par0c = 0.1;
		} else if (id == PIM) {
			if (fTheta0_e == NULL) {
				fTheta0_e = new TF1("fTheta0_e","[0]+[1]/((x+[2])*3375.0/2250.0)",0,6);
			}
			fTheta0 = fTheta0_e;
			par0a = 10.6067;
			par0b = 17.7586;
			par0c = 0.0527344;
		}
		fTheta0->SetParameters(par0a,par0b,par0c);
    return fTheta0;
  }

	static TF1* GetFidParm0(Int_t id) {
		Double_t p0, p1, p2, p3;
    if (fDphi0 == NULL) {
      fDphi0 = new TF1("fDphi0","[0]+[1]*tanh([2]*(x-[3]))",0,6);
    }
		if (id == ELECTRON) {
			p0 = 34.0; //28.0;
			p1 = 3.0; //5.0
			p2 = 0.85; //0.75;
			p3 = 1.2;
		} else if (id == PROTON) {
			p0 = 31;
			p1 = p2 = p3 = 0;
		} else if (id == PIP) {
			p0 = 31; //28;
			p1 = p2 = p3 = 0;
		} else if (id == PIM) {
			p0 = 28;
			p1 = p2 = p3 = 0;
		}
		fDphi0->SetParameters(p0,p1,p2,p3);
		return fDphi0;
	}

	static TF1* GetFidParm1(Int_t id) {
		Double_t p0, p1, p2;
    if (fDphi1 == NULL) {
			fDphi1 = new TF1("fDphi1","[0]+[1]*TMath::Exp(x/[2])",0,6);
		}
		if (id == ELECTRON) {
			p0 = 0.009684;
			p1 = 0.015;
			p2 = -0.823181;
		} else if (id == PROTON || id == PIP || id == PIM) {
			p0 = 0.22;
			p1 = p2 = 0;
		}
		fDphi1->SetParameters(p0,p1,p2);
		return fDphi1;
	}

	static TF1* GetFidParm2(Int_t id) {
		Double_t p0, p1, p2, p3;
		if (fDphi2 == NULL) {
			fDphi2 = new TF1("fDphi2","[0]+[1]/TMath::SinH([2]*(x-[3]))",0,6);
		}
		if (id == ELECTRON) {
			p0 = 1.946;
			p1 = 1.711;
			p2 = 0.8655;
			p3 = 0.665;
		} else if (id == PROTON || id == PIP || id == PIM) {
			p0 = 0.15;
			p1 = p2 = p3 = 0;
		}
		fDphi2->SetParameters(p0,p1,p2,p3);
		return fDphi2;
	}
	
	/** Fiducial cut function for single momentum slice, usable in TF1 object.
	 * @param x x[0] contains the theta value, the ordinate of the function.
	 * @param parms Paramaters that can be varied to minimize a fit (parms[0] and parms[3] are expected to be fixed).
	 * Parameter meanings:
	 * <ul>
	 * 	<li>0 -- first parameter of <i>standard</i> fiducial parametrization.</li>
	 * 	<li>1 -- second parameter of <i>standard</i> fiducial parametrization.</li>
	 * 	<li>2 -- second parameter of <i>standard</i> fiducial parametrization.</li>
	 * 	<li>3 -- momentum, fixed to corrected particle momentum.</li>
	 * 	<li>4 -- theta minimum, fixed according to ftheta0.</li>
	 * </ul>
	 */
  static Double_t FidFun(Double_t *x, Double_t *parms) {
	  Double_t d2r = TMath::DegToRad();
	  Double_t theta0 = parms[4];
	  Double_t par1 = parms[0];
	  Double_t par2 = parms[1];
	  Double_t par3 = parms[2];
	  Double_t p = parms[3];
	  Double_t theta = d2r*x[0];
	  theta0 = d2r*theta0;
	  //dphi describes distance from zero in degrees
	  //for one momentum and theta value
	  Double_t dphi = par1*TMath::Power(TMath::Sin(theta-theta0),par2*TMath::Power(p*3375.0/2250.0,par3));
	  return dphi;
  }
  
	/** Returns FidFun wrapped in TF1.
	 * 	@param id Particle ID.
	 * 	@param p	Particle momentum.
	 * 	@param norm Normalization, primarily used to flip between positive-phi sector edge to negative-phi sector edge.
	 * 	@param tparm Tightness parameter.
	 * 
	 * to-do:
	 * <ul>
	 * 	<li>Implement tightness parameter.</li>
	 * 	<li>Use TF1 functions for parameter values.</li>
	 * 	<li>Add functionality for other particles.</li>
	 * </ul>
	 */
  static TF1* GetFidF(Int_t id, Double_t p, Double_t norm = 1, Double_t tparm = 0) {
	  Double_t rad = 3375.0/2250.0;
	  TF1 *ft0 = GetTheta0F(id);
	  TF1 *fp0 = GetFidParm0(id);
	  TF1 *fp1 = GetFidParm1(id);
	  TF1 *fp2 = GetFidParm2(id);
	  Double_t theta0 = ft0->Eval(p);
	  TF1 *ret = new TF1("fidf",FidFun,theta0,50,5);
	  ret->SetParameters(norm*fp0->Eval(p),fp1->Eval(p),fp2->Eval(p),p,theta0);
	  return ret;
  }
  
  /** Requires particle to be within fiducial volume (theta, phi, p) and to be detected on <i>good</i> SC/TOF paddle.
   *  @param id Particle id.
   *  @param p Particle momentum.
   *  @param phi Azimuthal angle in detector frame, in degrees, -30 to 330.
   *  @param theta Polar angle in detector frame, in degrees, 0 to 180.
   *  @param sector Sector number, from 1 to 6
   *  @param paddle SC/TOF paddle number.
   *
   *  Note: Current function assumes phi symmetry of each sector.
   */
  static bool Fiducial(Int_t id, Double_t p, Double_t theta, Double_t phi,
		       Int_t sector, Int_t paddle = 0) {
		/* ****** return true if unimplemented ****** */
		if ( !(id == ELECTRON || id == PROTON || id == PIP || id == PIM) ) return true;
		
    //ensure that phi is between -30 and 330
    phi = phi > 330 ? phi -= 360 : phi;
    phi = phi < -30 ? phi += 360 : phi;
    //make phi symmetric about zero
    phi -= (sector-1)*60;
	  
	 TF1 *dphiFunc = GetFidF(id,p);
    //dphi describes distance from zero in degrees
    //for one momentum and theta value
    Double_t dphi = dphiFunc->Eval(theta);
	 delete dphiFunc;
		
    return TMath::Abs(phi)<dphi;
  }

/* *************************************************************
 * These _Get_*_CutPar functions are from Arjun's HMeid class. *
 * All require eid.out to be present.                          *
 * ************************************************************ */
	static epconfig::data _cutParMap;
	static Int_t _nPbins;
	static Int_t _cut_nPhe;
	//electron CC cut parms
	static double _a0_low[6], _a1_low[6], _a2_low[6];
	static double _a0_high[6], _a1_high[6], _a2_high[6];
	//electron SF cut function parms
	static double _mean_a0[6];
	static double _mean_a1[6];
	static double _mean_a2[6];
	static double _mean_a3[6];
	static double _sigma_a0[6];
	static double _sigma_a1[6];
	static double _sigma_a2[6];
	static double _sigma_a3[6];
	//electron EC cut parms
	static Float_t _EC_min_p;
	
	static bool _Get_CCtheta_CutPars(Bool_t sim = kFALSE) {
		bool returnVal = false;
		ifstream f;
		if (sim) f.open("eid.mc.out");
		else f.open("eid.out");
		if (f.is_open()) {
			fprintf(stdout, "Found configuration file! \n");
			f >> _cutParMap;
			//get SFcut parameters
			for (int iSector = 0; iSector < 6; iSector++) {
				char a0_low[100], a1_low[100], a2_low[100];
				char a0_high[100], a1_high[100], a2_high[100];
				sprintf(a0_low, "%d_a0_CCtheta_low", iSector+1);
				sprintf(a1_low, "%d_a1_CCtheta_low", iSector+1);
				sprintf(a2_low, "%d_a2_CCtheta_low", iSector+1);
				sprintf(a0_high, "%d_a0_CCtheta_high", iSector+1);
				sprintf(a1_high, "%d_a1_CCtheta_high", iSector+1);
				sprintf(a2_high, "%d_a2_CCtheta_high", iSector+1);
				sscanf(_cutParMap[a0_low].c_str(),"%lf",&_a0_low[iSector]);
				sscanf(_cutParMap[a1_low].c_str(),"%lf",&_a1_low[iSector]);
				sscanf(_cutParMap[a2_low].c_str(),"%lf",&_a2_low[iSector]);
				sscanf(_cutParMap[a0_high].c_str(),"%lf",&_a0_high[iSector]);
				sscanf(_cutParMap[a1_high].c_str(),"%lf",&_a1_high[iSector]);
				sscanf(_cutParMap[a2_high].c_str(),"%lf",&_a2_high[iSector]);
			}
			returnVal =  true;
		} else {
			fprintf(stdout, "Not Found file with cut pars! \n");
			returnVal = false;
		}
		return returnVal;
	}
	
	static bool _Get_CCtiming_CutPars() {
		printf("Not implemented!\n");
		return true;
	}
	static bool _Get_nPhe_CutPars() {
		_cut_nPhe = 30;
		return true;
	}
	static bool _Get_ECThreshold_CutPars() {
		_EC_min_p = 0.64;
		return true;
	}
	static bool _Get_SF_CutPars(Bool_t sim = kFALSE) {
		bool returnVal = false;
		ifstream f;
		if (sim) {
			f.open("eid.mc.out");
		} else {
			_mean_a0[0] = 0.260961;
			_mean_a1[0] = 0.033737;
			_mean_a2[0] = -0.006303;
			_mean_a3[0] = 0.000646;
			_mean_a0[1] = 0.274190;
			_mean_a1[1] = 0.066307;
			_mean_a2[1] = -0.015193;
			_mean_a3[1] = 0.001226;
			_mean_a0[2] = 0.266103;
			_mean_a1[2] = 0.048361;
			_mean_a2[2] = -0.009831;
			_mean_a3[2] = 0.000948;
			_mean_a0[3] = 0.265564;
			_mean_a1[3] = 0.038295;
			_mean_a2[3] = -0.005993;
			_mean_a3[3] = 0.000282;
			_mean_a0[4] = 0.245044;
			_mean_a1[4] = 0.050019;
			_mean_a2[4] = -0.011594;
			_mean_a3[4] = 0.000975;
			_mean_a0[5] = 0.247762;
			_mean_a1[5] = 0.042443;
			_mean_a2[5] = -0.008261;
			_mean_a3[5] = 0.000541;
			_sigma_a0[0] = 0.043446;
			_sigma_a1[0] = -0.019802;
			_sigma_a2[0] = 0.006339;
			_sigma_a3[0] = -0.000687;
			_sigma_a0[1] = 0.049097;
			_sigma_a1[1] = -0.009840;
			_sigma_a2[1] = 0.000051;
			_sigma_a3[1] = 0.000257;
			_sigma_a0[2] = 0.043670;
			_sigma_a1[2] = -0.015904;
			_sigma_a2[2] = 0.004181;
			_sigma_a3[2] = -0.000368;
			_sigma_a0[3] = 0.036371;
			_sigma_a1[3] = -0.002412;
			_sigma_a2[3] = 0.000049;
			_sigma_a3[3] = -0.000000;
			_sigma_a0[4] = 0.041725;
			_sigma_a1[4] = -0.008642;
			_sigma_a2[4] = 0.000555;
			_sigma_a3[4] = 0.000147;
			_sigma_a0[5] = 0.038737;
			_sigma_a1[5] = -0.011680;
			_sigma_a2[5] = 0.003047;
			_sigma_a3[5] = -0.000278;
			return true;
		}
		if (f.is_open()) {
			fprintf(stdout, "Found configuration file! \n");
			f >> _cutParMap;
			//get SFcut parameters
			for (int iSector = 0; iSector < 6; iSector++) {
				char a0_SFmean[100], a1_SFmean[100], a2_SFmean[100], a3_SFmean[100];
				char a0_SFsigma[100], a1_SFsigma[100], a2_SFsigma[100], a3_SFsigma[100];
				sprintf(a0_SFmean, "%d_SFmean_a0", iSector+1);
				sprintf(a1_SFmean, "%d_SFmean_a1", iSector+1);
				sprintf(a2_SFmean, "%d_SFmean_a2", iSector+1);
				sprintf(a3_SFmean, "%d_SFmean_a3", iSector+1);
				sprintf(a0_SFsigma, "%d_SFsigma_a0", iSector+1);
				sprintf(a1_SFsigma, "%d_SFsigma_a1", iSector+1);
				sprintf(a2_SFsigma, "%d_SFsigma_a2", iSector+1);
				sprintf(a3_SFsigma, "%d_SFsigma_a3", iSector+1);
				sscanf(_cutParMap[a0_SFmean].c_str(),"%lf",&_mean_a0[iSector]);
				sscanf(_cutParMap[a1_SFmean].c_str(),"%lf",&_mean_a1[iSector]);
				sscanf(_cutParMap[a2_SFmean].c_str(),"%lf",&_mean_a2[iSector]);
				sscanf(_cutParMap[a3_SFmean].c_str(),"%lf",&_mean_a3[iSector]);
				sscanf(_cutParMap[a0_SFsigma].c_str(),"%lf",&_sigma_a0[iSector]);
				sscanf(_cutParMap[a1_SFsigma].c_str(),"%lf",&_sigma_a1[iSector]);
				sscanf(_cutParMap[a2_SFsigma].c_str(),"%lf",&_sigma_a2[iSector]);
				sscanf(_cutParMap[a3_SFsigma].c_str(),"%lf",&_sigma_a3[iSector]);
			}
			_nPbins = 15;
			returnVal =  true;
		} else {
			fprintf(stdout, "Not Found file with cut pars! \n");
			returnVal = false;
		}
		return returnVal;
	}
	static bool _passNphe(Int_t nPhe) {
		return ( nPhe > _cut_nPhe );
	}
	static bool _passCCtheta(Int_t sector, Int_t segment, Float_t thetaCC) {
		Int_t iSector = sector - 1;
		Float_t thetaLowFit = _a2_low[iSector]*(segment)*(segment) + _a1_low[iSector]*(segment) + _a0_low[iSector];
		Float_t thetaHighFit = _a2_high[iSector]*(segment)*(segment) + _a1_high[iSector]*(segment) + _a0_high[iSector];
		//printf("passCC: %f\t%f\t%f\n",thetaCC,thetaLowFit,thetaHighFit);
		return (thetaCC < thetaHighFit) && (thetaCC > thetaLowFit);
	}
	static bool _passSF(Int_t sector, Float_t p, Float_t sf) {
		bool returnVal = false;
		Int_t iSector = sector - 1;
		Float_t mean = _mean_a3[iSector]*p*p*p + _mean_a2[iSector]*p*p + _mean_a1[iSector]*p + _mean_a0[iSector];
		Float_t sigma = _sigma_a3[iSector]*p*p*p + _sigma_a2[iSector]*p*p + _sigma_a1[iSector]*p + _sigma_a0[iSector];
		Float_t SFhigh = mean + 3*sigma;
		Float_t SFlow = mean - 3*sigma;
		//printf("passSF: %f\t%f\t%f\n",sf,SFlow,SFhigh);
		return (sf >= SFlow && sf <= SFhigh);
	}
	static bool IsConfigLoaded = false;
	static bool IsElectron(Int_t sector, Float_t p, Float_t sf, Bool_t sim = kFALSE) {
		if ( !IsConfigLoaded) {
			_Get_SF_CutPars(sim);
			//_Get_CCtheta_CutPars();
			//_Get_nPhe_CutPars();
			IsConfigLoaded = true;
		}
		return ( _passSF(sector,p,sf) );
		//&& Threshold_pE(p) && _passNphe(nPhe) ); && _passCCtheta(sector, cc_segm, thetaCC));
	}
	static bool IsElectron(Int_t sector, Float_t p, Int_t cc_segm, Float_t sf, Float_t thetaCC, Int_t nPhe) {
		if ( !IsConfigLoaded ) {
			_Get_SF_CutPars();
			_Get_CCtheta_CutPars();
			_Get_nPhe_CutPars();
			IsConfigLoaded = true;
		}
		return ( _passSF(sector,p,sf) );
		//&& Threshold_pE(p) && _passNphe(nPhe) ); && _passCCtheta(sector, cc_segm, thetaCC));
	}
/* ************************************************************
 * END _Get_*_CutPar functions are from Arjun's HMeid class *
 * ************************************************************ */

	static bool IsBadScint(Int_t sector, Int_t paddle, Int_t run = 0) {
		bool badsc = false;
		switch(sector) {
		case 1:
			if ( paddle == 24 ) badsc = true;
			break;
		case 2:
			if ( paddle == 16 || paddle == 28
				  || (run>38449 && paddle==38) )
				badsc = true;
			break;
		case 3:
			if ( paddle == 2 || paddle == 11 || paddle == 24
				  || paddle == 27 || paddle == 28 || paddle == 40 )
				badsc = true;
			break;
		case 4:
			if ( paddle == 19 || paddle == 30 || paddle == 34
				  || (run>38450 && paddle==2) )
				badsc = true;
			break;
		case 5:
			if ( paddle == 2 || paddle == 18 || paddle == 20
				  || paddle == 40 || (run>38449 && paddle==34) )
				badsc = true;
			break;
		case 6:
			if ( paddle == 1 || paddle == 18 || paddle == 40 )
				badsc = true;
			break;
		}
		return badsc;
	}

}

#endif /* CUTS_H_ */
