#ifndef Q2BNG_H
#define Q2BNG_H

#include <TROOT.h>


/* 
[09-07-14]
This file contains relevant binning for making h8s:
	+ Coarse-W bins for which h8[W] are made.
	+ For each h8[W]
		+ Q2,W,M(M1 & M2),THETA,PHI,ALPHA
*/

//! 1. Set up Coarse-W binning

//! Use for coarse W bin ranges: [xmin,xmax)
//! To see study that I did to determine this binning, see:
//! obs/nb_h8_Coarse-W-bng
struct Wbng{
  Int_t nbins;
  Float_t xmin;
  Float_t xmax;
};
const Int_t NBINS_WCRS=6;//13;
const Wbng WCRSBIN[NBINS_WCRS]={
	{5,1.300,1.425},
	{6,1.425,1.575},
	{6,1.575,1.725},
	{5,1.725,1.850},
	{6,1.850,2.000},
	{5,2.000,2.125}
	/*{2.125,2.275},
	{2.275,2.425},
	{2.425,2.550},
	{2.550,2.700},
	{2.700,2.825},
	{2.825,2.950},
	{2.950,3.000}*/
};

//! 2. Q2 binning (see nb_h8_Q2-bng)
const Int_t NBINS_Q2=8; //! 0.5 GeV^2/bin
const Float_t XMIN_Q2=1.25;
const Float_t XMAX_Q2=5.25;


//! 3. W binning
//! + All that is needed is the binw:
//! 	+ xmin,xmax,nbins are determined by Coarse-W bins and the binw:
//! 		+ nbins=(CrsW.max-CrsW.min)/binw
//!   + Note that CrsW.max,CrsW.min are multiples of 0.025!
const Float_t BINW_W=0.025; //! 25 MeV/bin

//! 4. M binning (for details see nb_h8_Coarse-W-bng)
//! + All that is needed are the nbins:
//!		+ xmin,xmax,binw are determined by the Coarse-W bins 
const Int_t NBINS_M=14;

//! 5. THETA,PHI,ALPHA
const Int_t NBINS_THETA=10;
const Float_t XMIN_THETA=0;
const Float_t XMAX_THETA=180;

const Int_t NBINS_PHI=10;
const Float_t XMIN_PHI=0;
const Float_t XMAX_PHI=360;

const Int_t NBINS_ALPHA=10;
const Float_t XMIN_ALPHA=0;
const Float_t XMAX_ALPHA=360;

Int_t GetCrsWBinIdx(Float_t w){
	Int_t idx=9999;
	for(int i=0;i<NBINS_WCRS;i++){
		if(w>=WCRSBIN[i].xmin && w<WCRSBIN[i].xmax){
			idx=i;
			break;
		}
	}
	return idx;
}

#endif // Q2BNG_H