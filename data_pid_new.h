#ifndef DATAPIDNEW_H
#define DATAPIDNEW_H

#include <TROOT.h>

class DataPidNew
{
public:
	DataPidNew();
	virtual ~DataPidNew();
	void Clear();
        Int_t h10IdxP;
        Int_t h10IdxPip;
        Int_t h10IdxPim;
        //! Following variables are added as
        //! Branches of Tree for study purpose
        //! Following are directly measured quantities for e-
        Float_t l_e;
        Float_t t_e;
        Float_t t_off;
        //! gpart
        Int_t gpart;
        //! gpart= npos+nneg+nzro
        Int_t npos;
        Int_t nneg;
        Int_t nzro;
        //! np,npip,npim
        Int_t np;
        Int_t npip;
        Int_t npim;
        //! Particles indexed by ntrk(+ve or -ve gpart particles)
        //! ntrk=Number of charged tracks
        static const Int_t kMaxTrack=20;
        Int_t ntrk;
        //! Following are directly measured quantities for ntrk
        Int_t q[kMaxTrack];
        Int_t dc[kMaxTrack];
        Int_t sc[kMaxTrack];
        Float_t p[kMaxTrack];
        Float_t l[kMaxTrack];
        Float_t t[kMaxTrack];
        Float_t b[kMaxTrack];
        Int_t id[kMaxTrack];
        Int_t sector[kMaxTrack];
        Int_t h10_idx[kMaxTrack];
        //! Following are quantities under particle assumption
        Float_t b_p[kMaxTrack];
        Float_t b_pip[kMaxTrack];
        Float_t b_pim[kMaxTrack];
        Float_t dt_p[kMaxTrack];
        Float_t dt_pip[kMaxTrack];
        Float_t dt_pim[kMaxTrack];
};

#endif 
