#ifndef Q2BNG_H
#define Q2BNG_H

#include <TROOT.h>


/* 
[09-07-14]
This file defines binning for Q2,W used for making h8s:
	+ Coarse W bins for which h8[W] are made.
	+ Q2 binning for h8[W]
	+ W binning for h8[W]
*/

//! 1. Set up coarse W binning
const Int_t kNCrsWBins=7;

//! Use for coarse W bin ranges: [xmin,xmax)
struct BinRange{
  Float_t xmin;
  Float_t xmax;
};

const BinRange kWCrsBin[kNCrsWBins]={
	{1.300,1.425},
	{1.425,1.575},
	{1.575,1.725},
	{1.725,1.850},
	{1.850,2.000},
	{2.000,2.125},
	{2.125,2.275},
	{2.275,2.425},
	{2.425,2,550},
	{2.550,2.700},
	{2.700,2.825},
	{2.825,2.950},
	{2.950,3.000}
};

//! 2. Set up Q2 binning (see nb_h8_Q2-bng)
const Int_t kNQ2Bins=8; //! 0.5 GeV^2/bin
const Float_t kQ2min=1.25
const Float_t kQ2max=5.25


//! W binning
//! + All that is needed is the binw, the limits are as per Coarse W bins
//! + nbins=(CrsW.max-CrsW.min)/binw
//!   + Note that CrsW.max,CrsW.min are multiples of 0.025!
const Int_t kWbinw=0.025 //! 25 MeV/bin

#endif // Q2BNG_H
