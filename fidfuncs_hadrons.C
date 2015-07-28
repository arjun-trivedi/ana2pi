#include "TF1.h"
using namespace TMath;

const Int_t ID_P = 2212;
const Int_t ID_PIP = 211;
const Int_t ID_PIM = -211;

//! + NPRT=2=proton and pip
//! + Fiducial cuts for pim are momentum dependent and therefore, different
const int NPRT=2; 
enum {P,PIP};
TString PRT_NAME[]={"p","pip"};

const int NDTYP=2;
enum {EXP,SIM};
TString DTYP_NAME[]={"exp","sim"};

const int NSCTR=6;

//! The following are taken from Evan's fid.hadrons.parms
float t0[NPRT][NDTYP]={
                       {8.12,9.38},
                       {8.12,6.88},
                       };

float Ah[NPRT][NDTYP][NSCTR]={
				{
					{23.925,22.143,24.602,23.923,23.265,23.343},
					{24.657,24.940,24.736,24.725,24.735,24.696}
				},
				{	{24.287,22.521,23.135,23.129,23.909,22.869},
					{24.536,24.372,24.542,24.504,24.490,24.297}
				}
			     };
float Bh[NPRT][NDTYP][NSCTR]={  
                                {
                                        {0.196,0.219,0.107,0.138,0.164,0.171},
                                        {0.112,0.097,0.104,0.106,0.104,0.105}
                                },
                                {       {0.167,0.212,0.131,0.157,0.138,0.174},
                                        {0.113,0.119,0.114,0.115,0.116,0.122}
                                }
                             };
float Ch[NPRT][NDTYP][NSCTR]={  
                                {
                                        {6.719,6.288,4.530,5.804,6.278,6.272},
                                        {4.423,2.995,3.344,3.777,3.614,3.634}
                                },
                                {       {7.128,7.319,6.634,6.907,6.027,6.836},
                                        {4.683,4.616,4.484,4.469,4.496,4.928}
                                }
                             };

float Al[NPRT][NDTYP][NSCTR]={
                                {
                                        {-23.355,-23.751,-23.264,-23.869,-22.186,-23.082},
                                        {-24.634,-24.748,-24.986,-24.713,-24.914,-24.706}
                                },
                                {       {-22.795,-23.142,-23.124,-23.510,-23.703,-22.760},
                                        {-24.690,-23.926,-23.704,-24.013,-24.151,-24.168}
                                }
                             };
float Bl[NPRT][NDTYP][NSCTR]={
                                {
                                        {0.143,0.128,0.158,0.140,0.232,0.162},
                                        {0.115,0.104,0.096,0.106,0.099,0.106}
                                },
                                {       {0.136,0.135,0.147,0.143,0.162,0.167},
                                        {0.107,0.124,0.130,0.124,0.121,0.118}
                                }
                             };
float Cl[NPRT][NDTYP][NSCTR]={
                                {
                                        {7.639,6.931,6.902,6.169,8.020,6.740},
                                        {4.782,3.574,2.995,3.863,3.150,3.727}
                                },
                                {       {7.695,7.511.7.084,7.050,7.956,7.525},
                                        {4.388,4.719,4.877,4.774,4.611,4.529}
                                }
                             };


Double_t dPhiFid_hdrn(Double_t *x, Double_t *parms) {
	Double_t theta = x[0];
	Double_t A = parms[0];
	Double_t B = parms[1];
	Double_t C = parms[2];
	Double_t phib = A*(1-Exp(-B*(theta-C)));
	return phib;
}

TF1* fPhiFid_hdrn_h(Int_t id, TString dtyp_name, Int_t sector, Int_t p) {	
	//! Get cut pars
	float _Ah=0;
	float _Bh=0;
	float _Ch=0;
	Int_t iprt=-1;
	Int_t idtyp=-1;

        if (id==ID_P){iprt=P;}
	else if (id==ID_PIP){iprt=PIP;}
	else if (id==ID_PIM){
		cout<<"Fiducial cuts for pi- not implemented"<<endl;
		return 0;
	}

	if (dtyp_name=="exp"){idtyp=EXP;}
	else if (dtyp_name=="sim"){idtyp=SIM;}
	else{printf("dtyp_name=%s not recognized",dtyp_name.Data());}

	_Ah=Ah[iprt][idtyp][sector-1];
	_Bh=Bh[iprt][idtyp][sector-1];
	_Ch=Ch[iprt][idtyp][sector-1];
	printf("Ah=%0.3f,Bh=%0.3f,Ch=%0.3f\n",_Ah,_Bh,_Ch);
	 
	TF1* retFunc = new TF1(TString::Format("fphifid_%s_%s_h",PRT_NAME[iprt],DTYP_NAME[idtyp]),dPhiFid_hdrn,0,70,3);
	retFunc->SetParameters(_Ah,_Bh,_Ch);
	return retFunc;
}

TF1* fPhiFid_hdrn_l(Int_t id, TString dtyp_name, Int_t sector, Int_t p) {
        //! Get cut pars
        float _Al=0;
        float _Bl=0;
        float _Cl=0;
        Int_t iprt=-1;
        Int_t idtyp=-1;

        if (id==ID_P){iprt=P;}
        else if (id==ID_PIP){iprt=PIP;}
        else if (id==ID_PIM){
                cout<<"Fiducial cuts for pi- not implemented"<<endl;
                return 0;
        }

        if (dtyp_name=="exp"){idtyp=EXP;}
        else if (dtyp_name=="sim"){idtyp=SIM;}
        else{printf("dtyp_name=%s not recognized",dtyp_name.Data());}

        _Al=Al[iprt][idtyp][sector-1];
        _Bl=Bl[iprt][idtyp][sector-1];
        _Cl=Cl[iprt][idtyp][sector-1];
        printf("Al=%0.3f,Bl=%0.3f,Cl=%0.3f\n",_Al,_Bl,_Cl);

        TF1* retFunc = new TF1(TString::Format("fphifid_%s_%s_l",PRT_NAME[iprt],DTYP_NAME[idtyp]),dPhiFid_hdrn,0,70,3);
        retFunc->SetParameters(_Al,_Bl,_Cl);
        return retFunc;
}

TLine* lt0(Int_t id, TString dtyp_name){
	//! Get t0
	float _t0=0;
	Int_t iprt=-1;
        Int_t idtyp=-1;

	if (id==ID_P){iprt=P;}
        else if (id==ID_PIP){iprt=PIP;}
        else if (id==ID_PIM){
                cout<<"Fiducial cuts for pi- not implemented"<<endl;
                return 0;
        }

        if (dtyp_name=="exp"){idtyp=EXP;}
        else if (dtyp_name=="sim"){idtyp=SIM;}
        else{printf("dtyp_name=%s not recognized",dtyp_name.Data());}
	
	_t0=t0[iprt][idtyp];
	TLine* retLine=new TLine(_t0,0,_t0,0);//y1 and y2 to be set by whatever uses lt0
	return retLine;
}

/****************************************************
[07-27-15] The following two functions
+ dPhiFid_hrdn_mod
+ fPhiFid_hrdn_h_mod
+ fPhiFid_hrdn_l_mod
are modified versions for dPhiFid, fPhiFid that return
not the default phib=[-30,30], but a phib that is appropriate
for each sector.
*****************************************************/
Double_t dPhiFid_hdrn_mod(Double_t *x, Double_t *parms) {
        Double_t theta = x[0];
        Double_t A = parms[0];
        Double_t B = parms[1];
        Double_t C = parms[2];
	//! trivedia modified
	Int_t sector=parms[3];
	//!
        Double_t phib = A*(1-Exp(-B*(theta-C)));
	//! trivedia modified 
        float offst[]={0,60,120,180,240,300}; //offst to get -[30,30] to values appropriate per sector
        phib+=offst[sector-1];
        //
        return phib;
}

TF1* fPhiFid_hdrn_h_mod(Int_t id, TString dtyp_name, Int_t sector, Int_t p) {
        //! Get cut pars
        float _Ah=0;
        float _Bh=0;
        float _Ch=0;
        Int_t iprt=-1;
        Int_t idtyp=-1;

        if (id==ID_P){iprt=P;}
        else if (id==ID_PIP){iprt=PIP;}
        else if (id==ID_PIM){
                cout<<"Fiducial cuts for pi- not implemented"<<endl;
                return 0;
        }

        if (dtyp_name=="exp"){idtyp=EXP;}
        else if (dtyp_name=="sim"){idtyp=SIM;}
        else{printf("dtyp_name=%s not recognized",dtyp_name.Data());}

        _Ah=Ah[iprt][idtyp][sector-1];
        _Bh=Bh[iprt][idtyp][sector-1];
        _Ch=Ch[iprt][idtyp][sector-1];
        printf("Ah=%0.3f,Bh=%0.3f,Ch=%0.3f\n",_Ah,_Bh,_Ch);

	//! trivedia modified
	//TF1* retFunc = new TF1(TString::Format("fphifid_%s_%s_h",PRT_NAME[iprt],DTYP_NAME[idtyp]),dPhiFid_hdrn_mod,0,70,3);
        //retFunc->SetParameters(_Ah,_Bh,_Ch);
        TF1* retFunc = new TF1(TString::Format("fphifid_%s_%s_h_mod",PRT_NAME[iprt],DTYP_NAME[idtyp]),dPhiFid_hdrn_mod,0,70,4);
        retFunc->SetParameters(_Ah,_Bh,_Ch,sector);
	//!
        return retFunc;
}

TF1* fPhiFid_hdrn_l_mod(Int_t id, TString dtyp_name, Int_t sector, Int_t p) {
        //! Get cut pars
        float _Al=0;
        float _Bl=0;
        float _Cl=0;
        Int_t iprt=-1;
        Int_t idtyp=-1;

        if (id==ID_P){iprt=P;}
        else if (id==ID_PIP){iprt=PIP;}
        else if (id==ID_PIM){
                cout<<"Fiducial cuts for pi- not implemented"<<endl;
                return 0;
        }

        if (dtyp_name=="exp"){idtyp=EXP;}
        else if (dtyp_name=="sim"){idtyp=SIM;}
        else{printf("dtyp_name=%s not recognized",dtyp_name.Data());}

        _Al=Al[iprt][idtyp][sector-1];
        _Bl=Bl[iprt][idtyp][sector-1];
        _Cl=Cl[iprt][idtyp][sector-1];
        printf("Al=%0.3f,Bl=%0.3f,Cl=%0.3f\n",_Al,_Bl,_Cl);

	//! trivedia modified
	//TF1* retFunc = new TF1(TString::Format("fphifid_%s_%s_l",PRT_NAME[iprt],DTYP_NAME[idtyp]),dPhiFid_hdrn_mod,0,70,3);
        //retFunc->SetParameters(_Al,_Bl,_Cl);
        TF1* retFunc = new TF1(TString::Format("fphifid_%s_%s_l_mod",PRT_NAME[iprt],DTYP_NAME[idtyp]),dPhiFid_hdrn_mod,0,70,4);
        retFunc->SetParameters(_Al,_Bl,_Cl,sector);
	//!
        return retFunc;
}

////////////////////////////////////////////////////////////////////////////

/*	
Bool_t inFid_hdrn(Double_t p, Double_t theta, Double_t phi,
		Int_t sector, Double_t tightness = 1, Double_t t0tightness = 0) {
	Int_t isector = sector-1;
	Double_t phib = getPhiFid(p, theta, sector, tightness, t0tightness);
	//make sure phi is between -30 and 330 degrees
	phi = phi >= 330 ? phi-360 : (phi<-30 ? phi+360 : phi);
	//rotate to sector center = 0 degrees
	phi -= 60*isector;
	return (phi < phib && phi > -1*phib);  //symmetric boundary
	//return Fid::Instance()->InFid(p,theta,phi,sector,11,tightness);
}*/

