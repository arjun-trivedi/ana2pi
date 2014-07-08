#ifndef Q2BNG_H
#define Q2BNG_H

#include <TROOT.h>


/* 
This file defines the Coarse and Analysis binning for Q2W and should
be referenced by all the other code in the analysis that need to know
about the binning.
See ana2pi notes (tomboy) for details.
*/

//! 1. Set up coarse Q2W binning
const Int_t kQ2_NCrsBins=7;
const Int_t kW_NCrsBins=7;
const Int_t kQ2W_NCrsBins=49;

//! Use for coarse Q2,W,Q2W bin ranges
struct BinRange{
  Float_t xmin;
  Float_t xmax;
};

struct Q2WBinRange{
  Float_t q2min;
  Float_t q2max;
  Float_t wmin;
  Float_t wmax;
};

const BinRange kQ2_CrsBin[kQ2_NCrsBins]={
	{1.2,1.6},
	{1.6,2.0},
	{2.0,2.4},
	{2.4,3.0},
	{3.0,3.5},
	{3.5,4.2},
	{4.2,5.0}
};

//! Ensure that kW_CrsBin[kW_NCrsBins-1].xmax - kW_CrsBin[0].xmin is 
//! a multiple of 0.025 GeV, since these values are later used to 
//! to determine the number of analysis-bins i.e. for h8 binning
const BinRange kW_CrsBin[kW_NCrsBins]={
	{1.3,1.7},
	{1.7,2.0},
	{2.0,2.2},
	{2.2,2.4},
	{2.2,2.6},
	{2.2,2.8},
	{2.8,3.0} 
};

const Q2WBinRange kQ2W_CrsBin[kQ2W_NCrsBins]={
{kQ2_CrsBin[0].xmin,kQ2_CrsBin[0].xmax,kW_CrsBin[0].xmin,kW_CrsBin[0].xmax},
{kQ2_CrsBin[0].xmin,kQ2_CrsBin[0].xmax,kW_CrsBin[1].xmin,kW_CrsBin[1].xmax},
{kQ2_CrsBin[0].xmin,kQ2_CrsBin[0].xmax,kW_CrsBin[2].xmin,kW_CrsBin[2].xmax},
{kQ2_CrsBin[0].xmin,kQ2_CrsBin[0].xmax,kW_CrsBin[3].xmin,kW_CrsBin[3].xmax},
{kQ2_CrsBin[0].xmin,kQ2_CrsBin[0].xmax,kW_CrsBin[4].xmin,kW_CrsBin[4].xmax},
{kQ2_CrsBin[0].xmin,kQ2_CrsBin[0].xmax,kW_CrsBin[5].xmin,kW_CrsBin[5].xmax},
{kQ2_CrsBin[0].xmin,kQ2_CrsBin[0].xmax,kW_CrsBin[6].xmin,kW_CrsBin[6].xmax},


{kQ2_CrsBin[1].xmin,kQ2_CrsBin[1].xmax,kW_CrsBin[0].xmin,kW_CrsBin[0].xmax},
{kQ2_CrsBin[1].xmin,kQ2_CrsBin[1].xmax,kW_CrsBin[1].xmin,kW_CrsBin[1].xmax},
{kQ2_CrsBin[1].xmin,kQ2_CrsBin[1].xmax,kW_CrsBin[2].xmin,kW_CrsBin[2].xmax},
{kQ2_CrsBin[1].xmin,kQ2_CrsBin[1].xmax,kW_CrsBin[3].xmin,kW_CrsBin[3].xmax},
{kQ2_CrsBin[1].xmin,kQ2_CrsBin[1].xmax,kW_CrsBin[4].xmin,kW_CrsBin[4].xmax},
{kQ2_CrsBin[1].xmin,kQ2_CrsBin[1].xmax,kW_CrsBin[5].xmin,kW_CrsBin[5].xmax},
{kQ2_CrsBin[1].xmin,kQ2_CrsBin[1].xmax,kW_CrsBin[6].xmin,kW_CrsBin[6].xmax},

{kQ2_CrsBin[2].xmin,kQ2_CrsBin[2].xmax,kW_CrsBin[0].xmin,kW_CrsBin[0].xmax},
{kQ2_CrsBin[2].xmin,kQ2_CrsBin[2].xmax,kW_CrsBin[1].xmin,kW_CrsBin[1].xmax},
{kQ2_CrsBin[2].xmin,kQ2_CrsBin[2].xmax,kW_CrsBin[2].xmin,kW_CrsBin[2].xmax},
{kQ2_CrsBin[2].xmin,kQ2_CrsBin[2].xmax,kW_CrsBin[3].xmin,kW_CrsBin[3].xmax},
{kQ2_CrsBin[2].xmin,kQ2_CrsBin[2].xmax,kW_CrsBin[4].xmin,kW_CrsBin[4].xmax},
{kQ2_CrsBin[2].xmin,kQ2_CrsBin[2].xmax,kW_CrsBin[5].xmin,kW_CrsBin[5].xmax},
{kQ2_CrsBin[2].xmin,kQ2_CrsBin[2].xmax,kW_CrsBin[6].xmin,kW_CrsBin[6].xmax},

{kQ2_CrsBin[3].xmin,kQ2_CrsBin[3].xmax,kW_CrsBin[0].xmin,kW_CrsBin[0].xmax},
{kQ2_CrsBin[3].xmin,kQ2_CrsBin[3].xmax,kW_CrsBin[1].xmin,kW_CrsBin[1].xmax},
{kQ2_CrsBin[3].xmin,kQ2_CrsBin[3].xmax,kW_CrsBin[2].xmin,kW_CrsBin[2].xmax},
{kQ2_CrsBin[3].xmin,kQ2_CrsBin[3].xmax,kW_CrsBin[3].xmin,kW_CrsBin[3].xmax},
{kQ2_CrsBin[3].xmin,kQ2_CrsBin[3].xmax,kW_CrsBin[4].xmin,kW_CrsBin[4].xmax},
{kQ2_CrsBin[3].xmin,kQ2_CrsBin[3].xmax,kW_CrsBin[5].xmin,kW_CrsBin[5].xmax},
{kQ2_CrsBin[3].xmin,kQ2_CrsBin[3].xmax,kW_CrsBin[6].xmin,kW_CrsBin[6].xmax},

{kQ2_CrsBin[4].xmin,kQ2_CrsBin[4].xmax,kW_CrsBin[0].xmin,kW_CrsBin[0].xmax},
{kQ2_CrsBin[4].xmin,kQ2_CrsBin[4].xmax,kW_CrsBin[1].xmin,kW_CrsBin[1].xmax},
{kQ2_CrsBin[4].xmin,kQ2_CrsBin[4].xmax,kW_CrsBin[2].xmin,kW_CrsBin[2].xmax},
{kQ2_CrsBin[4].xmin,kQ2_CrsBin[4].xmax,kW_CrsBin[3].xmin,kW_CrsBin[3].xmax},
{kQ2_CrsBin[4].xmin,kQ2_CrsBin[4].xmax,kW_CrsBin[4].xmin,kW_CrsBin[4].xmax},
{kQ2_CrsBin[4].xmin,kQ2_CrsBin[4].xmax,kW_CrsBin[5].xmin,kW_CrsBin[5].xmax},
{kQ2_CrsBin[4].xmin,kQ2_CrsBin[4].xmax,kW_CrsBin[6].xmin,kW_CrsBin[6].xmax},

{kQ2_CrsBin[5].xmin,kQ2_CrsBin[5].xmax,kW_CrsBin[0].xmin,kW_CrsBin[0].xmax},
{kQ2_CrsBin[5].xmin,kQ2_CrsBin[5].xmax,kW_CrsBin[1].xmin,kW_CrsBin[1].xmax},
{kQ2_CrsBin[5].xmin,kQ2_CrsBin[5].xmax,kW_CrsBin[2].xmin,kW_CrsBin[2].xmax},
{kQ2_CrsBin[5].xmin,kQ2_CrsBin[5].xmax,kW_CrsBin[3].xmin,kW_CrsBin[3].xmax},
{kQ2_CrsBin[5].xmin,kQ2_CrsBin[5].xmax,kW_CrsBin[4].xmin,kW_CrsBin[4].xmax},
{kQ2_CrsBin[5].xmin,kQ2_CrsBin[5].xmax,kW_CrsBin[5].xmin,kW_CrsBin[5].xmax},
{kQ2_CrsBin[5].xmin,kQ2_CrsBin[5].xmax,kW_CrsBin[6].xmin,kW_CrsBin[6].xmax},

{kQ2_CrsBin[6].xmin,kQ2_CrsBin[6].xmax,kW_CrsBin[0].xmin,kW_CrsBin[0].xmax},
{kQ2_CrsBin[6].xmin,kQ2_CrsBin[6].xmax,kW_CrsBin[1].xmin,kW_CrsBin[1].xmax},
{kQ2_CrsBin[6].xmin,kQ2_CrsBin[6].xmax,kW_CrsBin[2].xmin,kW_CrsBin[2].xmax},
{kQ2_CrsBin[6].xmin,kQ2_CrsBin[6].xmax,kW_CrsBin[3].xmin,kW_CrsBin[3].xmax},
{kQ2_CrsBin[6].xmin,kQ2_CrsBin[6].xmax,kW_CrsBin[4].xmin,kW_CrsBin[4].xmax},
{kQ2_CrsBin[6].xmin,kQ2_CrsBin[6].xmax,kW_CrsBin[5].xmin,kW_CrsBin[5].xmax},
{kQ2_CrsBin[6].xmin,kQ2_CrsBin[6].xmax,kW_CrsBin[6].xmin,kW_CrsBin[6].xmax}
};

//! 2. Set up Ana Q2W binning
//! The following is used by DataAna::makeYields for h8 bng

//! Q2 variable binning
const Int_t kQ2_NAnaBins=kQ2_NCrsBins;
//const Float_t kQ2_AnaBins[kQ2_NAnaBins]={1.2,1.6,2.0,2.4,3.0,3.5,4.2,5.0};
const Float_t kQ2_AnaBins[kQ2_NAnaBins+1]={
	kQ2_CrsBin[0].xmin,
	kQ2_CrsBin[1].xmin,
	kQ2_CrsBin[2].xmin,
	kQ2_CrsBin[3].xmin,
	kQ2_CrsBin[4].xmin,
	kQ2_CrsBin[5].xmin,
	kQ2_CrsBin[6].xmin,
	kQ2_CrsBin[6].xmax
};
const Int_t kQ2_min=kQ2_CrsBin[0].xmin;
const Int_t kQ2_max=kQ2_CrsBin[kQ2_NCrsBins-1].xmax;

//~ W binning
const Int_t kW_min=kW_CrsBin[0].xmin;
const Int_t kW_max=kW_CrsBin[kW_NCrsBins-1].xmax;
const Float_t kW_binw=0.025;
const Int_t kW_NAnaBins=(kW_max-kW_min)/kW_binw;



//! W binning
Int_t kW_Nbins=80; //25 MeV/bin

#endif // Q2BNG_H
