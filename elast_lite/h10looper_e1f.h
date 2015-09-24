//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep  6 09:24:36 2015 by ROOT version 5.33/03
// from TTree h10/CLAS
// found on file: /data/trivedia/e1f/elastic/sim/sim1.1/recon/1.root
//////////////////////////////////////////////////////////

#ifndef h10looper_e1f_h
#define h10looper_e1f_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h> //! Makes Info() work
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>

#include "constants.h"
using namespace E1F;

#include "mom_corr.cpp"
#include "fidfuncs.C"
#include "pid.h"

class h10looper_e1f {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //! To decode input h10type
   TString _expt;
   TString _dtyp;
   TString _seq;

   //! fout
   TFile* _fout;

   //! nentries_to_proc
   Long64_t _nentries_to_proc;

   //! cuts to be made in addition to 'dflt'
   //! eid: new cuts and corrections
   bool _use_cut_ECin_min; //!1
   bool _use_cut_ECfid; //!2
   bool _use_cut_zvtx;//!3
   bool _use_corr_sf_etot;//!4
   //! eid: to match Isupov
   bool _use_SChit;//!5
   bool _use_dc_stat;//!6
   //! Evans's EFID
   bool _use_ep_efid;//!7
   //! proton 
   bool _use_proton;//!8

   //! output objects
   //! Event level
   static const int NUM_EVT_STATS=4;
   enum {EVT_NULL, EVT_TRG, EVT_E, EVT_E_INFID, EVT_ELSTC};
   TH1D* _hevt;

   //! EID
   static const int NUM_EID_STATS=14;
   enum {EID_NULL, EID_TRG, EID_GPART0, EID_Q, 
         EID_HIT_DC, EID_HIT_CC, EID_HIT_SC, EID_HIT_EC, 
         EID_STAT, EID_DC_STAT, EID_P_MIN_ECTH, EID_ECIN_MIN, EID_EC_FID, EID_ZVTX, EID_SF};
   TH1D* _heid;
   //! for SF cut
   TF1** _sf_mean;
   TF1** _sf_min;
   TF1** _sf_max;
   //! for ECin min cut
   Float_t* _ECmin;
   //! for z-vtx cut
   Float_t* _z_vtx_min;
   Float_t* _z_vtx_max;

   //! EFID
   static const int NUM_EFID_STATS=2;
   enum {EFID_NULL, EFID_TOT, EFID_IN};
   TH1D* _hefid;

   //! mom corr
   TH2D* _hpcorr_dpVp;
   TH1D* _hpcorr_dcx;
   TH1D* _hpcorr_dcy;
   TH1D* _hpcorr_dcz;
   TH1D* _hpcorr_dp;
   //! for mom corr
   MomCorr_e1f* _pcorr;

   //! PID
   Pid* _pid_tool;
   static const int NUM_PID_STATS=2;
   enum {PID_NULL, PID_TOT, PID_PASS};
   TH1D* _hpid;
   
   //! delast
   TH1D* _hW;
   TH1D* _helast;
   static const int _NBINS=32;//64;
   static const float _THETA_MIN=14;
   static const float _THETA_MAX=46;
   TH1D** _hf;//!full sector
   TH1D** _hc;//!central sector
   //! elastic cut
   static const float _W_CUT_MIN=0.848;
   static const float _W_CUT_MAX=1.028;   

   //! ekin
   TLorentzVector _lvE0;
   TLorentzVector _lvP0;
   TLorentzVector* _lvE1;
   TLorentzVector* _lvQ;
   TLorentzVector* _lvW;
   float _Q2;
   float _W;
   float _theta;
   float _phi; //[-30,330]

   // Declaration of leaf types
   static const int _MAX_PARTS=100;
   Int_t           nprt;
   Int_t           pidpart[_MAX_PARTS];   //[nprt]
   Float_t         xpart[_MAX_PARTS];   //[nprt]
   Float_t         ypart[_MAX_PARTS];   //[nprt]
   Float_t         zpart[_MAX_PARTS];   //[nprt]
   Float_t         epart[_MAX_PARTS];   //[nprt]
   Float_t         pxpart[_MAX_PARTS];   //[nprt]
   Float_t         pypart[_MAX_PARTS];   //[nprt]
   Float_t         pzpart[_MAX_PARTS];   //[nprt]
   Float_t         qpart[_MAX_PARTS];   //[nprt]
   Int_t           flagspart[_MAX_PARTS];   //[nprt]
   UChar_t         npart;
   UChar_t         evstat;
   UInt_t          evntid;
   Char_t          evtype;
   Char_t          evntclas;
   Char_t          evthel;
   Int_t           evntclas2;
   Float_t         q_l;
   Float_t         t_l;
   Float_t         tr_time;
   Float_t         rf_time1;
   Float_t         rf_time2;
   Int_t           gpart;
   Short_t         id[_MAX_PARTS];   //[gpart]
   Char_t          stat[_MAX_PARTS];   //[gpart]
   UChar_t         dc[_MAX_PARTS];   //[gpart]
   UChar_t         cc[_MAX_PARTS];   //[gpart]
   UChar_t         sc[_MAX_PARTS];   //[gpart]
   UChar_t         ec[_MAX_PARTS];   //[gpart]
   UChar_t         lec[_MAX_PARTS];   //[gpart]
   Float_t         p[_MAX_PARTS];   //[gpart]
   Float_t         m[_MAX_PARTS];   //[gpart]
   Char_t          q[_MAX_PARTS];   //[gpart]
   Float_t         b[_MAX_PARTS];   //[gpart]
   Float_t         cx[_MAX_PARTS];   //[gpart]
   Float_t         cy[_MAX_PARTS];   //[gpart]
   Float_t         cz[_MAX_PARTS];   //[gpart]
   Float_t         vx[_MAX_PARTS];   //[gpart]
   Float_t         vy[_MAX_PARTS];   //[gpart]
   Float_t         vz[_MAX_PARTS];   //[gpart]
   Int_t           dc_part;
   UChar_t         dc_sect[_MAX_PARTS];   //[dc_part]
   UChar_t         dc_trk[_MAX_PARTS];   //[dc_part]
   Char_t          dc_stat[_MAX_PARTS];   //[dc_part]
   Float_t         dc_xsc[_MAX_PARTS];   //[dc_part]
   Float_t         dc_ysc[_MAX_PARTS];   //[dc_part]
   Float_t         dc_zsc[_MAX_PARTS];   //[dc_part]
   Float_t         dc_cxsc[_MAX_PARTS];   //[dc_part]
   Float_t         dc_cysc[_MAX_PARTS];   //[dc_part]
   Float_t         dc_czsc[_MAX_PARTS];   //[dc_part]
   Float_t         dc_xec[_MAX_PARTS];   //[dc_part]
   Float_t         dc_yec[_MAX_PARTS];   //[dc_part]
   Float_t         dc_zec[_MAX_PARTS];   //[dc_part]
   Float_t         dc_thcc[_MAX_PARTS];   //[dc_part]
   Float_t         dc_c2[_MAX_PARTS];   //[dc_part]
   Int_t           ec_part;
   UShort_t        ec_stat[_MAX_PARTS];   //[ec_part]
   UChar_t         ec_sect[_MAX_PARTS];   //[ec_part]
   Int_t           ec_whol[_MAX_PARTS];   //[ec_part]
   Int_t           ec_inst[_MAX_PARTS];   //[ec_part]
   Int_t           ec_oust[_MAX_PARTS];   //[ec_part]
   Float_t         etot[_MAX_PARTS];   //[ec_part]
   Float_t         ec_ei[_MAX_PARTS];   //[ec_part]
   Float_t         ec_eo[_MAX_PARTS];   //[ec_part]
   Float_t         ec_t[_MAX_PARTS];   //[ec_part]
   Float_t         ec_r[_MAX_PARTS];   //[ec_part]
   Float_t         ech_x[_MAX_PARTS];   //[ec_part]
   Float_t         ech_y[_MAX_PARTS];   //[ec_part]
   Float_t         ech_z[_MAX_PARTS];   //[ec_part]
   Float_t         ec_m2[_MAX_PARTS];   //[ec_part]
   Float_t         ec_m3[_MAX_PARTS];   //[ec_part]
   Float_t         ec_m4[_MAX_PARTS];   //[ec_part]
   Float_t         ec_c2[_MAX_PARTS];   //[ec_part]
   Int_t           sc_part;
   UChar_t         sc_sect[_MAX_PARTS];   //[sc_part]
   UChar_t         sc_hit[_MAX_PARTS];   //[sc_part]
   UChar_t         sc_pd[_MAX_PARTS];   //[sc_part]
   UChar_t         sc_stat[_MAX_PARTS];   //[sc_part]
   Float_t         edep[_MAX_PARTS];   //[sc_part]
   Float_t         sc_t[_MAX_PARTS];   //[sc_part]
   Float_t         sc_r[_MAX_PARTS];   //[sc_part]
   Float_t         sc_c2[_MAX_PARTS];   //[sc_part]
   Int_t           cc_part;
   UChar_t         cc_sect[_MAX_PARTS];   //[cc_part]
   UChar_t         cc_hit[_MAX_PARTS];   //[cc_part]
   Int_t           cc_segm[_MAX_PARTS];   //[cc_part]
   UShort_t        nphe[_MAX_PARTS];   //[cc_part]
   Float_t         cc_t[_MAX_PARTS];   //[cc_part]
   Float_t         cc_r[_MAX_PARTS];   //[cc_part]
   Float_t         cc_c2[_MAX_PARTS];   //[cc_part]
   Int_t           lac_part;
   Int_t           lec_sect[_MAX_PARTS];   //[lac_part]
   Int_t           lec_hit[_MAX_PARTS];   //[lac_part]
   Int_t           lec_stat[_MAX_PARTS];   //[lac_part]
   Float_t         lec_etot[_MAX_PARTS];   //[lac_part]
   Float_t         lec_t[_MAX_PARTS];   //[lac_part]
   Float_t         lec_r[_MAX_PARTS];   //[lac_part]
   Float_t         lec_x[_MAX_PARTS];   //[lac_part]
   Float_t         lec_y[_MAX_PARTS];   //[lac_part]
   Float_t         lec_z[_MAX_PARTS];   //[lac_part]
   Float_t         lec_c2[_MAX_PARTS];   //[lac_part]

   // List of branches
   TBranch        *b_nprt;   //!
   TBranch        *b_pidpart;   //!
   TBranch        *b_xpart;   //!
   TBranch        *b_ypart;   //!
   TBranch        *b_zpart;   //!
   TBranch        *b_epart;   //!
   TBranch        *b_pxpart;   //!
   TBranch        *b_pypart;   //!
   TBranch        *b_pzpart;   //!
   TBranch        *b_qpart;   //!
   TBranch        *b_flagspart;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_evstat;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_evtype;   //!
   TBranch        *b_evntclas;   //!
   TBranch        *b_evthel;   //!
   TBranch        *b_evntclas2;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_t_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_rf_time1;   //!
   TBranch        *b_rf_time2;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_id;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_ec;   //!
   TBranch        *b_lec;   //!
   TBranch        *b_p;   //!
   TBranch        *b_m;   //!
   TBranch        *b_q;   //!
   TBranch        *b_b;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_cy;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_dc_part;   //!
   TBranch        *b_dc_sect;   //!
   TBranch        *b_dc_trk;   //!
   TBranch        *b_dc_stat;   //!
   TBranch        *b_dc_xsc;   //!
   TBranch        *b_dc_ysc;   //!
   TBranch        *b_dc_zsc;   //!
   TBranch        *b_dc_cxsc;   //!
   TBranch        *b_dc_cysc;   //!
   TBranch        *b_dc_czsc;   //!
   TBranch        *b_dc_xec;   //!
   TBranch        *b_dc_yec;   //!
   TBranch        *b_dc_zec;   //!
   TBranch        *b_dc_thcc;   //!
   TBranch        *b_dc_c2;   //!
   TBranch        *b_ec_part;   //!
   TBranch        *b_ec_stat;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_ec_whol;   //!
   TBranch        *b_ec_inst;   //!
   TBranch        *b_ec_oust;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_ec_t;   //!
   TBranch        *b_ec_r;   //!
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_ec_m2;   //!
   TBranch        *b_ec_m3;   //!
   TBranch        *b_ec_m4;   //!
   TBranch        *b_ec_c2;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_hit;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_sc_c2;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_hit;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_cc_t;   //!
   TBranch        *b_cc_r;   //!
   TBranch        *b_cc_c2;   //!
   TBranch        *b_lac_part;   //!
   TBranch        *b_lec_sect;   //!
   TBranch        *b_lec_hit;   //!
   TBranch        *b_lec_stat;   //!
   TBranch        *b_lec_etot;   //!
   TBranch        *b_lec_t;   //!
   TBranch        *b_lec_r;   //!
   TBranch        *b_lec_x;   //!
   TBranch        *b_lec_y;   //!
   TBranch        *b_lec_z;   //!
   TBranch        *b_lec_c2;   //!

   h10looper_e1f(TString h10type,TChain* h10chain,TString fout_name, Long64_t nentries,
                 TString adtnl_cut_opt);
   virtual ~h10looper_e1f();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void setup_eid_cutpars(TString dtyp);
   void setup_adtnl_cut_opts(TString adtnl_cut_opt);

   bool evt_trigger_electron();
   void mom_corr_electron();
   bool electron_infid();
   void make_delast();

   void set_ekin();
   void reset_ekin();
   bool found_proton();
   
   int get_sector();
   void GetUVW(float xyz[3], float uvw[3]);
   bool pass_p_min_ECth();
   bool pass_ECin_min();
   bool pass_ECfid();
   bool pass_zvtx();
   bool pass_sf();
};

#endif

#ifdef h10looper_e1f_cxx
Int_t h10looper_e1f::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t h10looper_e1f::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void h10looper_e1f::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   if (_seq=="thrown"){
      fChain->SetBranchAddress("nprt", &nprt, &b_nprt);
      fChain->SetBranchAddress("pidpart", pidpart, &b_pidpart);
      fChain->SetBranchAddress("xpart", xpart, &b_xpart);
      fChain->SetBranchAddress("ypart", ypart, &b_ypart);
      fChain->SetBranchAddress("zpart", zpart, &b_zpart);
      fChain->SetBranchAddress("epart", epart, &b_epart);
      fChain->SetBranchAddress("pxpart", pxpart, &b_pxpart);
      fChain->SetBranchAddress("pypart", pypart, &b_pypart);
      fChain->SetBranchAddress("pzpart", pzpart, &b_pzpart);
      fChain->SetBranchAddress("qpart", qpart, &b_qpart);
      fChain->SetBranchAddress("flagspart", flagspart, &b_flagspart);
   }
   if (_seq=="recon"){
      fChain->SetBranchAddress("npart", &npart, &b_npart);
      fChain->SetBranchAddress("evstat", &evstat, &b_evstat);
      fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
      fChain->SetBranchAddress("evtype", &evtype, &b_evtype);
      fChain->SetBranchAddress("evntclas", &evntclas, &b_evntclas);
      fChain->SetBranchAddress("evthel", &evthel, &b_evthel);
      fChain->SetBranchAddress("evntclas2", &evntclas2, &b_evntclas2);
      fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
      fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
      fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
      fChain->SetBranchAddress("rf_time1", &rf_time1, &b_rf_time1);
      fChain->SetBranchAddress("rf_time2", &rf_time2, &b_rf_time2);
      fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
      fChain->SetBranchAddress("id", id, &b_id);
      fChain->SetBranchAddress("stat", stat, &b_stat);
      fChain->SetBranchAddress("dc", dc, &b_dc);
      fChain->SetBranchAddress("cc", cc, &b_cc);
      fChain->SetBranchAddress("sc", sc, &b_sc);
      fChain->SetBranchAddress("ec", ec, &b_ec);
      fChain->SetBranchAddress("lec", lec, &b_lec);
      fChain->SetBranchAddress("p", p, &b_p);
      fChain->SetBranchAddress("m", m, &b_m);
      fChain->SetBranchAddress("q", q, &b_q);
      fChain->SetBranchAddress("b", b, &b_b);
      fChain->SetBranchAddress("cx", cx, &b_cx);
      fChain->SetBranchAddress("cy", cy, &b_cy);
      fChain->SetBranchAddress("cz", cz, &b_cz);
      fChain->SetBranchAddress("vx", vx, &b_vx);
      fChain->SetBranchAddress("vy", vy, &b_vy);
      fChain->SetBranchAddress("vz", vz, &b_vz);
      fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
      fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
      fChain->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
      fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
      fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
      fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
      fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
      fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
      fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
      fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
      fChain->SetBranchAddress("dc_xec", dc_xec, &b_dc_xec);
      fChain->SetBranchAddress("dc_yec", dc_yec, &b_dc_yec);
      fChain->SetBranchAddress("dc_zec", dc_zec, &b_dc_zec);
      fChain->SetBranchAddress("dc_thcc", dc_thcc, &b_dc_thcc);
      fChain->SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
      fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
      fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
      fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
      fChain->SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
      fChain->SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
      fChain->SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
      fChain->SetBranchAddress("etot", etot, &b_etot);
      fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
      fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
      fChain->SetBranchAddress("ec_t", ec_t, &b_ec_t);
      fChain->SetBranchAddress("ec_r", ec_r, &b_ec_r);
      fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
      fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
      fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
      fChain->SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
      fChain->SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
      fChain->SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
      fChain->SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
      fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
      fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
      fChain->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
      fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
      fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
      fChain->SetBranchAddress("edep", edep, &b_edep);
      fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
      fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
      fChain->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
      fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
      fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
      fChain->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
      fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
      fChain->SetBranchAddress("nphe", nphe, &b_nphe);
      fChain->SetBranchAddress("cc_t", cc_t, &b_cc_t);
      fChain->SetBranchAddress("cc_r", cc_r, &b_cc_r);
      fChain->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
      fChain->SetBranchAddress("lac_part", &lac_part, &b_lac_part);
      fChain->SetBranchAddress("lec_sect", lec_sect, &b_lec_sect);
      fChain->SetBranchAddress("lec_hit", lec_hit, &b_lec_hit);
      fChain->SetBranchAddress("lec_stat", lec_stat, &b_lec_stat);
      fChain->SetBranchAddress("lec_etot", lec_etot, &b_lec_etot);
      fChain->SetBranchAddress("lec_t", lec_t, &b_lec_t);
      fChain->SetBranchAddress("lec_r", lec_r, &b_lec_r);
      fChain->SetBranchAddress("lec_x", lec_x, &b_lec_x);
      fChain->SetBranchAddress("lec_y", lec_y, &b_lec_y);
      fChain->SetBranchAddress("lec_z", lec_z, &b_lec_z);
      fChain->SetBranchAddress("lec_c2", lec_c2, &b_lec_c2);
   }
   Notify();
}

Bool_t h10looper_e1f::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void h10looper_e1f::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t h10looper_e1f::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef h10looper_e1f_cxx
