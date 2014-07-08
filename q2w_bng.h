#ifndef Q2BNG_H
#define Q2BNG_H

#include <TROOT.h>

//! See ana2pi notes (tb) for details
const Int_t kNCrsQ2bins=7;
const Int_t kNCrsWbins=7;
const Int_t kNCrsQ2Wbins=49;

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

BinRange kQ2bin[kNCrsQ2bins]={
	{1.2,1.6},
	{1.6,2.0},
	{2.0,2.4},
	{2.4,3.0},
	{3.0,3.5},
	{3.5,4.2},
	{4.2,5.0}
};

BinRange kWbin[kNCrsWbins]={
	{1.3,1.7},
	{1.7,2.0},
	{2.0,2.2},
	{2.2,2.4},
	{2.2,2.6},
	{2.2,2.8},
	{2.8,3.0} 
};

//! Q2 variable binning
Int_t kQ2bin_nbins=7;
Float_t kQ2bin_bins[8]={1.2,1.6,2.0,2.4,3.0,3.5,4.2,5.0};

Q2WBinRange kQ2Wbin[kNCrsQ2Wbins]={
{kQ2bin[0].xmin,kQ2bin[0].xmax,kWbin[0].xmin,kWbin[0].xmax},
{kQ2bin[0].xmin,kQ2bin[0].xmax,kWbin[1].xmin,kWbin[1].xmax},
{kQ2bin[0].xmin,kQ2bin[0].xmax,kWbin[2].xmin,kWbin[2].xmax},
{kQ2bin[0].xmin,kQ2bin[0].xmax,kWbin[3].xmin,kWbin[3].xmax},
{kQ2bin[0].xmin,kQ2bin[0].xmax,kWbin[4].xmin,kWbin[4].xmax},
{kQ2bin[0].xmin,kQ2bin[0].xmax,kWbin[5].xmin,kWbin[5].xmax},
{kQ2bin[0].xmin,kQ2bin[0].xmax,kWbin[6].xmin,kWbin[6].xmax},


{kQ2bin[1].xmin,kQ2bin[1].xmax,kWbin[0].xmin,kWbin[0].xmax},
{kQ2bin[1].xmin,kQ2bin[1].xmax,kWbin[1].xmin,kWbin[1].xmax},
{kQ2bin[1].xmin,kQ2bin[1].xmax,kWbin[2].xmin,kWbin[2].xmax},
{kQ2bin[1].xmin,kQ2bin[1].xmax,kWbin[3].xmin,kWbin[3].xmax},
{kQ2bin[1].xmin,kQ2bin[1].xmax,kWbin[4].xmin,kWbin[4].xmax},
{kQ2bin[1].xmin,kQ2bin[1].xmax,kWbin[5].xmin,kWbin[5].xmax},
{kQ2bin[1].xmin,kQ2bin[1].xmax,kWbin[6].xmin,kWbin[6].xmax},

{kQ2bin[2].xmin,kQ2bin[2].xmax,kWbin[0].xmin,kWbin[0].xmax},
{kQ2bin[2].xmin,kQ2bin[2].xmax,kWbin[1].xmin,kWbin[1].xmax},
{kQ2bin[2].xmin,kQ2bin[2].xmax,kWbin[2].xmin,kWbin[2].xmax},
{kQ2bin[2].xmin,kQ2bin[2].xmax,kWbin[3].xmin,kWbin[3].xmax},
{kQ2bin[2].xmin,kQ2bin[2].xmax,kWbin[4].xmin,kWbin[4].xmax},
{kQ2bin[2].xmin,kQ2bin[2].xmax,kWbin[5].xmin,kWbin[5].xmax},
{kQ2bin[2].xmin,kQ2bin[2].xmax,kWbin[6].xmin,kWbin[6].xmax},

{kQ2bin[3].xmin,kQ2bin[3].xmax,kWbin[0].xmin,kWbin[0].xmax},
{kQ2bin[3].xmin,kQ2bin[3].xmax,kWbin[1].xmin,kWbin[1].xmax},
{kQ2bin[3].xmin,kQ2bin[3].xmax,kWbin[2].xmin,kWbin[2].xmax},
{kQ2bin[3].xmin,kQ2bin[3].xmax,kWbin[3].xmin,kWbin[3].xmax},
{kQ2bin[3].xmin,kQ2bin[3].xmax,kWbin[4].xmin,kWbin[4].xmax},
{kQ2bin[3].xmin,kQ2bin[3].xmax,kWbin[5].xmin,kWbin[5].xmax},
{kQ2bin[3].xmin,kQ2bin[3].xmax,kWbin[6].xmin,kWbin[6].xmax},

{kQ2bin[4].xmin,kQ2bin[4].xmax,kWbin[0].xmin,kWbin[0].xmax},
{kQ2bin[4].xmin,kQ2bin[4].xmax,kWbin[1].xmin,kWbin[1].xmax},
{kQ2bin[4].xmin,kQ2bin[4].xmax,kWbin[2].xmin,kWbin[2].xmax},
{kQ2bin[4].xmin,kQ2bin[4].xmax,kWbin[3].xmin,kWbin[3].xmax},
{kQ2bin[4].xmin,kQ2bin[4].xmax,kWbin[4].xmin,kWbin[4].xmax},
{kQ2bin[4].xmin,kQ2bin[4].xmax,kWbin[5].xmin,kWbin[5].xmax},
{kQ2bin[4].xmin,kQ2bin[4].xmax,kWbin[6].xmin,kWbin[6].xmax},

{kQ2bin[5].xmin,kQ2bin[5].xmax,kWbin[0].xmin,kWbin[0].xmax},
{kQ2bin[5].xmin,kQ2bin[5].xmax,kWbin[1].xmin,kWbin[1].xmax},
{kQ2bin[5].xmin,kQ2bin[5].xmax,kWbin[2].xmin,kWbin[2].xmax},
{kQ2bin[5].xmin,kQ2bin[5].xmax,kWbin[3].xmin,kWbin[3].xmax},
{kQ2bin[5].xmin,kQ2bin[5].xmax,kWbin[4].xmin,kWbin[4].xmax},
{kQ2bin[5].xmin,kQ2bin[5].xmax,kWbin[5].xmin,kWbin[5].xmax},
{kQ2bin[5].xmin,kQ2bin[5].xmax,kWbin[6].xmin,kWbin[6].xmax},

{kQ2bin[6].xmin,kQ2bin[6].xmax,kWbin[0].xmin,kWbin[0].xmax},
{kQ2bin[6].xmin,kQ2bin[6].xmax,kWbin[1].xmin,kWbin[1].xmax},
{kQ2bin[6].xmin,kQ2bin[6].xmax,kWbin[2].xmin,kWbin[2].xmax},
{kQ2bin[6].xmin,kQ2bin[6].xmax,kWbin[3].xmin,kWbin[3].xmax},
{kQ2bin[6].xmin,kQ2bin[6].xmax,kWbin[4].xmin,kWbin[4].xmax},
{kQ2bin[6].xmin,kQ2bin[6].xmax,kWbin[5].xmin,kWbin[5].xmax},
{kQ2bin[6].xmin,kQ2bin[6].xmax,kWbin[6].xmin,kWbin[6].xmax}
};

#endif // Q2BNG_H
