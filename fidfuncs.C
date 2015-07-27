#include "TF1.h"

Double_t getPaddleWeight(Int_t sector, Int_t paddle) {
  Bool_t pw[6][51];
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 51; j++) {
      pw[i][j] = 1;
    }
  }
  pw[0][24] = pw[0][42] = pw[0][46] = pw[0][47] = 0;
  pw[1][16] = pw[1][38] = pw[1][41] = pw[1][42] = pw[1][44] = pw[1][45] = pw[1][46] = pw[1][48] = 0;
  pw[2][2] = pw[2][11] = pw[2][24] = pw[2][25] = pw[2][27] = pw[2][28] = pw[2][40] = pw[2][41] = pw[2][42] = pw[2][43] = pw[2][44] = 0;
  pw[3][34] = pw[3][42] = pw[3][44] = pw[3][45] = pw[3][48] = 0;
  pw[4][2] = pw[4][18] = pw[4][34] = pw[4][36] = pw[4][40] = pw[4][42] = pw[4][43] = pw[4][44] = 0;
  pw[5][1] = pw[5][18] = pw[5][40] = pw[5][42] = pw[5][43] = pw[5][45] = pw[5][46] = pw[5][47] = 0;
  return (Double_t)(pw[sector-1][paddle]);
}

Double_t getPhiFid(Double_t p, Double_t theta, Int_t sector,
		Int_t tightness = 1, Int_t t0tightness = 0) {
	Int_t isector = sector-1;
	Double_t t1mt0[] = {8,8,8,8,8,8};
	Double_t ma[] = {1.5,1.5,1.5,1.5,1.5,1.5};
	Double_t mb[] = {0.1,0.1,0.1,0.1,0.1,0.1};
	Double_t t2[] = {70,70,70,70,70,70};
	Double_t t0parms[6][3] =
		{	{10.787299801658907, 22.45857319517163, 0.21749987294959106},
			{10.828070790230456, 22.12341012158817, 0.19010756666978365},
			{9.312911097075856, 30.973375501625757, 0.5039814192331635},
			{11.520929543525531, 19.709322349228394, 0.13831688498503214},
			{9.529958851709138, 29.55466910593996, 0.45512398922081015},
			{9.971998112496019, 28.44016526621845, 0.44311697599215194}
		};
	Double_t phi0parms[6][4] =
		{	{1.6839188190954342, 10.474164108316266, -3.0499977341743354, 3.432592304744964},
			{1.6762956840192975, 10.2529104749553, -3.449091515358539, 3.6491736492791502},
			{1.6999999999990663, 9.97998798303922, -2.5398914524479257, 1.928870428609786},
			{1.3000000000023624, 12.675184729042954, -5.858975050283604, 6.157806651701965},
			{1.6499998517208754, 11.32820223312787, -3.1583551082351176, 2.499451832577039},
			{1.3253604827603942, 11.464333746894148, -3.454330645184062, 3.0981926752851523}
		};
	//calculate min theta, parameter t0
	Double_t irat = 3375.0/2250.0;
	Double_t t0, t0a, t0b, t0c;
	t0a = t0parms[isector][0];
	t0b = t0parms[isector][1];
	t0c = t0parms[isector][2];
	t0 = t0a + t0b/((p+t0c)*irat);
	//calculate phi value at min theta, parameter phi0
	Double_t phi0, xl0l1, y0l0, ml0, ml1;
	xl0l1 = phi0parms[isector][0];
	y0l0 = phi0parms[isector][1];
	ml0 = phi0parms[isector][2];
	ml1 = phi0parms[isector][3];
	if (p <= xl0l1) phi0 = y0l0 + ml0*p;
	else phi0 = (y0l0 + ml0*p) + ml1*(p-xl0l1);
	//calculate phi boundary
	Double_t phib = 0;
	Double_t t1 = t0 + t1mt0[isector];
	if ( theta <= t0+t0tightness || theta >= t2[isector] ) phib = 0;
	else {
		Double_t phi1 = phi0 + ma[isector]*(t1-t0);
		if ( theta > t0 && theta < t1 ) phib = phi0 + ma[isector]*(theta-t0);
		else if ( theta >= t1 && theta < t2[isector] ) phib = phi1 + mb[isector]*(theta-t1);
	}
	phib -= tightness;
	phib = phib <= 0 ? 0 : phib;
	return phib;
}

Double_t dPhiFid(Double_t *x, Double_t *parms) {
	Double_t theta = x[0];
	Double_t p = parms[0];
	Int_t sector = parms[1];
	Double_t tightness = parms[2];
	Double_t t0tightness = parms[3];
	Int_t w = parms[4];
	Double_t phib = getPhiFid(p, theta, sector, tightness, t0tightness);
	return w*phib;
}

TF1* fPhiFid(Double_t p, Int_t sector, Int_t weight = 1,
		Double_t tightness = 0, Double_t t0tightness = 0) {
	TF1 *retFunc = new TF1("fphifid",dPhiFid,0,70,5);
	retFunc->SetParameters(p,sector,tightness,t0tightness,weight);
	return retFunc;
}

/****************************************************
[07-27-15] The following two functions
+ dPhiFid_mod
+ fPhiFid_mod
are modified versions for dPhiFid, fPhiFid that return
not the default phib=[-30,30], but a phib that is appropriate
for each sector.
*****************************************************/
Double_t dPhiFid_mod(Double_t *x, Double_t *parms) {
    Double_t theta = x[0];
    Double_t p = parms[0];
    Int_t sector = parms[1];
    Double_t tightness = parms[2];
    Double_t t0tightness = parms[3];
    Int_t w = parms[4];
    Double_t phib = getPhiFid(p, theta, sector, tightness, t0tightness);
	//! trivedia modified 
	float offst[]={0,60,120,180,240,300}; //offst to get -[30,30] to values appropriate per sector
	phib+=offst[sector-1];
	//
    return w*phib;
}

TF1* fPhiFid_mod(Double_t p, Int_t sector, Int_t weight = 1,
                Double_t tightness = 0, Double_t t0tightness = 0) {
	//! trivedia modified
	//TF1 *retFunc = new TF1("fphifid",dPhiFid,0,70,5);
    TF1 *retFunc = new TF1("fphifid_mod",dPhiFid_mod,0,70,5);
	//
    retFunc->SetParameters(p,sector,tightness,t0tightness,weight);
    return retFunc;
}

////////////////////////////////////////////////////////////////////////////
	
Bool_t inFid(Double_t p, Double_t theta, Double_t phi,
		Int_t sector, Double_t tightness = 1, Double_t t0tightness = 0) {
	Int_t isector = sector-1;
	Double_t phib = getPhiFid(p, theta, sector, tightness, t0tightness);
	//make sure phi is between -30 and 330 degrees
	phi = phi >= 330 ? phi-360 : (phi<-30 ? phi+360 : phi);
	//rotate to sector center = 0 degrees
	phi -= 60*isector;
	return (phi < phib && phi > -1*phib);  //symmetric boundary
	//return Fid::Instance()->InFid(p,theta,phi,sector,11,tightness);
}

