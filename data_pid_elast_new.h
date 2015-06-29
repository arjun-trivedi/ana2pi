#ifndef DATAPIDELASTNEW_H
#define DATAPIDELASTNEW_H

#include <TROOT.h>

class DataPidElastNew
{
public:
	DataPidElastNew();
	virtual ~DataPidElastNew();
	void Clear();
        Int_t h10IdxP;
        //! Following variables can be added as
        //! Branches of Tree for study purpose
        static const Int_t kMaxTrack=20;
        //! Following are directly measured quantities for e-
        Float_t l_e;
        Float_t t_e;
        Float_t t_off;
        //! ntrk=Number of charged tracks
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
