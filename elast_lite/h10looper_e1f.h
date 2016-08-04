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
//using namespace E1F; //! [11-12-15] now dynamically chosing namespace as per _expt at "appropriate places"

//! [10-11-15] 
//! + The following includes were moved to h10looper_e1f.cpp and
//!   replaced by Forward Declarations of Classes defined in the includes
//!   and used here.
//! + This was done to avoid "multiple declaration" errors from g++
//!   after implementing h10looper_2pi.h/.cpp that includes this file (since
//!   it derives from Clas h10looper_e1f) and therefore objects in these includes 
//!   were getting defined twice.
/*#include "mom_corr.cpp"
#include "fidfuncs.C"
#include "pid.h"*/
class MomCorr_e1f;
class Pid;

class h10looper_e1f {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   //! To decode input h10type
   TString _expt;
   TString _dtyp;
   TString _rctn;
   TString _seq;

   //! cutsncors: applicable only for Reconstructed Events
   bool _do_eid;
   bool _do_efid;
   bool _do_pcorr;
   bool _do_pid;
   bool _do_pfid;
   bool _do_evtsel_2pi;
   bool _do_evtsel_elast;
   bool _do_eff;
   bool _do_scpd;

   //! fout
   TFile* _fout;

   //! nentries_to_proc
   Long64_t _nentries_to_proc;

   //! additional options
   //! Numeric-coded options
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
   //! Evans's PFID
   bool _use_ep_pfid; //!8
   //! Variable binning in Q2 
   bool _use_Q2_var_binw_bng; //!9
   //! Topology t2(p,pip) [Default is to use t2'[p,pip,(pim)]]
   bool _use_t2; //!10
   //! MM2_cut_EI 
   bool _use_MM2_cut_EI; //! 11
   //! Reconcile
   //! + [12-06-15]
   //! + Till now the decision to Reconcile was based on 
   //!   expt and dtyp, and thus was done for expt:dtyp=e16:exp
   //! + However, after getting e16-sim-EI, I realize its needs to be Reconciled
   //! + Since I do not want to complicate the code by adding another possible 
   //!   value for dtyp other that exp and sim, I am going to use the following additional 
   //!   option when working with sim-EI or any other data that needs Reconciling.
   bool _do_reconcile; //! 12
   //! gpart cut for PID
   bool _use_gpart_pid; //! 13
   //! hit SC for PID
   bool _use_SChit_pid; //! 14
   //! stat for PID
   bool _use_stat_pid; //! 15
   //! Q2,W limits used for analysis
   bool _use_thesis_Q2W; //!16
   //! eff_scpd_at_mod //! 17
   bool _use_eff_scpd_at_mod;
   //! use_cut_ECfid_at_mod //! 18
   bool _use_cut_ECfid_at_mod;
   //! MM2_cut_SS //! conservatively chosen to [-0.16,0.16] => MM:[-0.4,0.4]
   bool _use_MM2_cut_SS; //! 19
   //! CC_cut_eff_lse/tgt //! 20/21
   bool _use_CC_cut_lse; //! 20
   bool _use_CC_cut_tgt; //! 21
   bool _use_CC_cut; //! this is set equal to _use_CC_cut_lse || _use_CC_cut_tgt
   //! The following options for using pmtC, which is normally ignored (no cut, wgt=1),
   //! will work only when option 20 or 21 are used:
   //! + av: cut,eff,wgt=av(pmtL,pmtR)
   //! + L: cut,eff,wgt=pmtL
   //! + R: cut,eff,wgt=pmtR
   bool _use_CC_cut_pmtC_av; //! 22
   bool _use_CC_cut_pmtC_L;  //! 23 
   bool _use_CC_cut_pmtC_R;  //! 24

   //! char-coded options
   bool _make_h10_skim_e;//! eid+efid; for Reco events
   /* [02-26-15]
   + Make skims for xsex=f(cutsncors) def.eq. Systematic Studies (SS)
   + Turn OFF all cuts that are used in SS (default values in parenthesis):
      1. ECfid:     (OFF)
      2. eff_scpd:  not use ':eff:scpd:' in proc-chain
      3. gpart_pid: (OFF)
      4. hitSC_pid: (ON) => adtnl_opt=':14:' used to turn OFF
      5. stat_pid:  (OFF) 
      6. MM2_cut:   Use wider-cut => adtnl_opt=':19:' for '_use_MM2_cut_SS'
   */
   bool _make_h10_skim_SS; //! 'eid:efid:pid:pfid:pcorr:evtsel_2pi:',':11:14:'
   //! [05-09-16]
   bool _make_h10_skim_eid;//! eid; for Reco events


   //! output objects
   //! + Only cuts-n-corrs objects common to 'elast' and '2pi':
   //!    + Event-level,EID,EFID,pcorr,PID 
   //!   are declared here.
   //!   (NOTE: Event-level and PID hists should be different for 'elast' and '2pi',
   //!    but they have been hacked to work for both)
   //!   
   //! + cuts-n-corrs objects that are different:
   //!    + PFID
   //!   are declared in h10looper_2pi.h
   
   //! Event level
   //! Following used for _rcnt=elast
   enum {EVT_NULL, EVT_TRG, EVT_E, EVT_E_INFID};
   static const int NUM_EVT_STATS_ELAST=4;
   enum {EVT_ELSTC=4};
   //! Following appended if _rctn=2pi
   static const int NUM_EVT_STATS_2PI=9;
   enum {EVT_P_PIP=4, EVT_P_PIP_INFID, EVT_INEFF, EVT_INSCPD,EVT_Q2W_KIN_PASS, EVT_2PI};
   //! common hist
   TH1D* _hevt;
   /*enum {EVT_NULL, EVT_TRG, EVT_E, EVT_E_INFID,EVT_ELSTC};
   static const int NUM_EVT_STATS=4;
   TH1D* _hevt;*/

   //! EID
   static const int NUM_EID_STATS=16;
   enum {EID_NULL, EID_TRG, EID_GPART0, EID_Q, 
         EID_HIT_DC, EID_HIT_CC, EID_HIT_SC, EID_HIT_EC, 
         EID_STAT, EID_DC_STAT, EID_P_MIN_ECTH, EID_ECIN_MIN, EID_EC_FID, EID_ZVTX, EID_NPHE, EID_SF, EID_E};
   TH1D* _heid;
   //! for p_min_ECth
   float _p_min_ECth;
   //! for ECin min cut
   Float_t* _ECmin;
   //! for ECfid cut
   float _Umin,_Umax;
   float _Vmin,_Vmax;
   float _Wmin,_Wmax;
   //! for SF cut
   TF1** _sf_mean;
   TF1** _sf_min;
   TF1** _sf_max;
   //! _hzvtxcorr[6][2] for zvtxcorr monitoring, before and after cut
   TH1F*** _hzvtxcorr;
   //! for z-vtx cut
   Float_t* _z_vtx_min;
   Float_t* _z_vtx_max;
   TH1F*** _hzvtxcut; //!_hzvtxcut[6][2]
   //! _hsf[6][2], SF for each sector, before and after cut
   TH2F*** _hsf;
   //! CC cut values[sct][sgm]
   float** _CC_cut_val;
   //! CC cut eff[sct][sgm][pmt]
   float*** _CC_cut_eff; 
   //! CC cut wgt[sct][sgm][pmt]
   float*** _CC_cut_wgt;
   //! _hnphe[6][18][2][2], hsf[sct][sgm][pmt][cut]
   TH1F***** _hnphe;
   
   //! EFID
   static const int NUM_EFID_STATS=2;
   enum {EFID_NULL, EFID_TOT, EFID_IN};
   TH1D* _hefid;
   //! _hefid_e[2], before and after cut
   TH2F** _hefid_e;
   
   //! pcorr objects[3] for e,p,pip
   TString* _pcorr_prtcl_names;
   TH2D** _hpcorr_dpVp;
   TH1D** _hpcorr_dcx;
   TH1D** _hpcorr_dcy;
   TH1D** _hpcorr_dcz;
   TH1D** _hpcorr_dp;
   //! for e1f:pcorr:e (currently e1f pcorr only for e)
   MomCorr_e1f* _pcorr_tool_e1f;

   //! PID
   Pid* _pid_tool;
   static const int NUM_PID_STATS=10;
   enum {PID_NULL, PID_TOT, PID_GPART_PASS, PID_Q_PASS, PID_HIT_DC, PID_HIT_SC, PID_STAT_GT0,
         PID_P_FOUND, PID_PIP_FOUND, PID_PIM_FOUND, PID_P_AND_PIP_FOUND};
   TH1D* _hpid;
   //! pid-cut prec and pstc hists (in each sector)
   //! + PID is done for each track by
   //!   1. Obtaining directly measured beta and mom.
   //!   2. Calculating dt vs. mom. based on particle hypothesis
   //!   3. Makeing cut on dt=0 for dt vs. mom. 
   //! *** prec *** 
   //! + beta vs. p for +tracks
   //! + dt vs. p for each +track under p and pip hypothesis
   //! + beta vs. p for -tracks
   //! + dt vs. p for each -track under pim hypothesis
   TH2D** _h_betaVp_pos;
   TH2D** _h_dtVp_pos_p;
   TH2D** _h_dtVp_pos_pip;
   TH2D** _h_betaVp_neg;
   TH2D** _h_dtVp_neg_pim;
   //! *** pstc ***
   //! + beta vs. p and dt vs. p for tracks after PID
   TH2D** _h_betaVp_p;
   TH2D** _h_dtVp_p;
   TH2D** _h_betaVp_pip;
   TH2D** _h_dtVp_pip;
   TH2D** _h_betaVp_pim;
   TH2D** _h_dtVp_pim;

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

   //! + Used for h10-skims
   //! + Currently: _make_h10_skim_e, [02-26-15] _make_h10_skim_SS
   TTree* _th10copy;

   //! ekin
   TLorentzVector _lvE0;
   TLorentzVector _lvP0;
   TLorentzVector _lvE1;
   TLorentzVector _lvQ;
   TLorentzVector _lvW;
   float _Q2;
   float _W;
   float _theta_e;
   float _phi_e; //[-30,330]
   float _p_e;
   int _sector_e,_sc_pd_e;
   int _h10idx_e;
   
   //! hadron kinematics
   //! Lab frame
   TLorentzVector _lvP;
   TLorentzVector _lvPip;
   TLorentzVector _lvPim;
   //! Varsets
   TLorentzVector _lvQCMS;
   TLorentzVector _lvP0CMS;
   TLorentzVector _lvPCMS;
   TLorentzVector _lvPipCMS;
   TLorentzVector _lvPimCMS;
   //! Lab frame (phi=[-30,330])
   float _p_p,_theta_p,_phi_p;
   float _p_pip,_theta_pip,_phi_pip;
   float _p_pim,_theta_pim,_phi_pim;
   //! sector,paddle
   int _sector_p,_sc_pd_p;
   int _sector_pip,_sc_pd_pip;
   int _sector_pim,_sc_pd_pim;
   //! Varsets kinematics (phi=[0,360])
   float _M_ppip,_M_ppim,_M_pippim;
   float _p_cms_p,_theta_cms_p,_phi_cms_p;
   float _p_cms_pip,_theta_cms_pip,_phi_cms_pip;
   float _p_cms_pim,_theta_cms_pim,_phi_cms_pim;
   float _alpha_1,_alpha_2,_alpha_3;
   //! h10idx
   int _h10idx_p,_h10idx_pip,_h10idx_pim;

   // Declaration of leaf types
   static const int _MAX_PARTS=100;
   //! PART Banks
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
   //! The following PART bank data was found in E16:exp!
   Int_t           Ipart10[_MAX_PARTS];   //[nprt]
   Float_t         Rpart11[_MAX_PARTS];   //[nprt]
   Float_t         Rpart12[_MAX_PARTS];   //[nprt]
   Int_t           Ipart13[_MAX_PARTS];   //[nprt]
   //! MCTK Banks
   Int_t mcnentr;
   UChar_t mcnpart;
   Int_t mcst[_MAX_PARTS]; //[mcnentr]
   Int_t mcid[_MAX_PARTS]; //[mcnentr]
   Int_t mcpid[_MAX_PARTS]; //[mcnentr]
   Float_t mctheta[_MAX_PARTS]; //[mcnentr]
   Float_t mcphi[_MAX_PARTS]; //[mcnentr]
   Float_t mcp[_MAX_PARTS]; //[mcnentr]
   Float_t mcm[_MAX_PARTS]; //[mcnentr]
   Float_t mcvx_x_el; //[mcnentr]
   Float_t mcvx_y_el; //[mcnentr]
   Float_t mcvx_z_el; //[mcnentr]
   //! SEB Banks
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
   //! Addition SEB Branches in E16:exp
   Float_t         rf_time;
   Int_t           l2bit;
   Int_t           l3bit;
   Int_t           hlsc;
   Int_t           intt;
   Int_t           st[_MAX_PARTS];   //[gpart]
   Int_t           tb_st[_MAX_PARTS];   //[dc_part]
   Float_t         dc_vx[_MAX_PARTS];   //[dc_part]
   Float_t         dc_vy[_MAX_PARTS];   //[dc_part]
   Float_t         dc_vz[_MAX_PARTS];   //[dc_part]
   Float_t         dc_vr[_MAX_PARTS];   //[dc_part]
   Float_t         tl1_cx[_MAX_PARTS];   //[dc_part]
   Float_t         tl1_cy[_MAX_PARTS];   //[dc_part]
   Float_t         tl1_cz[_MAX_PARTS];   //[dc_part]
   Float_t         tl1_x[_MAX_PARTS];   //[dc_part]
   Float_t         tl1_y[_MAX_PARTS];   //[dc_part]
   Float_t         tl1_z[_MAX_PARTS];   //[dc_part]
   Float_t         tl1_r[_MAX_PARTS];   //[dc_part]
   Float_t         lec_ein[_MAX_PARTS];   //[lac_part]
   Int_t           nschit;
   Int_t           scsect[_MAX_PARTS];   //[nschit]
   Int_t           schid[_MAX_PARTS];   //[nschit]
   Int_t           scpid[_MAX_PARTS];   //[nschit]
   Float_t         sct[_MAX_PARTS];   //[nschit]
   Float_t         sce[_MAX_PARTS];   //[nschit]
   Float_t         scx[_MAX_PARTS];   //[nschit]
   Float_t         scy[_MAX_PARTS];   //[nschit]
   Float_t         scz[_MAX_PARTS];   //[nschit]

   //for some Branches of SEB' Banks that are different
   struct tmp{
      Int_t id[_MAX_PARTS]; //[gpart]
      Int_t stat[_MAX_PARTS]; //[gpart]
      Int_t dc[_MAX_PARTS]; //[gpart]
      Int_t cc[_MAX_PARTS]; //[gpart]
      Int_t sc[_MAX_PARTS]; //[gpart]
      Int_t ec[_MAX_PARTS]; //[gpart]
      Int_t q[_MAX_PARTS]; //[gpart]
      Int_t dc_sect[_MAX_PARTS]; //[dc_part]
      Int_t dc_stat[_MAX_PARTS]; //[dc_part]
      Int_t ec_stat[_MAX_PARTS]; //[ec_part]
      Int_t ec_sect[_MAX_PARTS]; //[ec_part]
      Int_t sc_sect[_MAX_PARTS]; //[sc_part]
      Int_t sc_pd[_MAX_PARTS]; //[sc_part]
      Int_t sc_stat[_MAX_PARTS]; //[sc_part]
      Int_t cc_sect[_MAX_PARTS]; //[cc_part]
      Int_t nphe[_MAX_PARTS]; //[cc_part]
   } tmpVars;

   // List of branches
   //! PART Banks' Branches
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
   //! The following PART bank data was found in E16:exp!
   TBranch        *b_flagspart;   //!
   TBranch        *b_Ipart10;  //
   TBranch        *b_Rpart11;  //
   TBranch        *b_Rpart12;  //
   TBranch        *b_Ipart13;  //
   //! MCTK Bank's Branches 
   TBranch *b_mcnentr; //!
   TBranch *b_mcnpart; //!
   TBranch *b_mcst; //!
   TBranch *b_mcid; //!
   TBranch *b_mcpid; //!
   TBranch *b_mctheta; //!
   TBranch *b_mcphi; //!
   TBranch *b_mcp; //!
   TBranch *b_mcm; //!
   TBranch *b_mcvx_x_el; //!
   TBranch *b_mcvx_y_el; //!
   TBranch *b_mcvx_z_el; //!
   //! SEB Banks's Branches
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
   //! Addition SEB Branches in E16:exp
   TBranch        *b_rf_time;   //!
   TBranch        *b_l2bit;   //!
   TBranch        *b_l3bit;   //!
   TBranch        *b_hlsc;   //!
   TBranch        *b_intt;   //!
   TBranch        *b_st;   //!
   TBranch        *b_tb_st;   //!
   TBranch        *b_dc_vx;   //!
   TBranch        *b_dc_vy;   //!
   TBranch        *b_dc_vz;   //!
   TBranch        *b_dc_vr;   //!
   TBranch        *b_tl1_cx;   //!
   TBranch        *b_tl1_cy;   //!
   TBranch        *b_tl1_cz;   //!
   TBranch        *b_tl1_x;   //!
   TBranch        *b_tl1_y;   //!
   TBranch        *b_tl1_z;   //!
   TBranch        *b_tl1_r;   //!
   TBranch        *b_lec_ein;   //!
   TBranch        *b_nschit;   //!
   TBranch        *b_scsect;   //!
   TBranch        *b_schid;   //!
   TBranch        *b_scpid;   //!
   TBranch        *b_sct;   //!
   TBranch        *b_sce;   //!
   TBranch        *b_scx;   //!
   TBranch        *b_scy;   //!
   TBranch        *b_scz;   //!

   h10looper_e1f(TString h10type,TChain* h10chain,
                 TString cutsncors, 
                 TString fout_name, Long64_t nentries,
                 TString adtnl_opts);
   virtual ~h10looper_e1f();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void set_h10_SEB_BranchStatus(TTree* tr);
   void Reconcile();

   void setup_cutsncors(TString cutsncors);

   void setup_eid_cutpars(TString dtyp);
   void setup_eid_CC_cut_val(TString lvl);
   void setup_eid_CC_cut_eff(TString lvl);
   void setup_eid_CC_cut_wgt();
   void setup_adtnl_opts(TString adtnl_opts);

   bool pass_eid();
   bool pass_efid();
   void make_delast();

   bool pass_theta_vs_p(TString prtcl_name);
   bool is_scpd_bad(TString prtcl_name);

   void do_pcorr(); //![01-13-16] Works for {e1f,e16}*{2pi:top2',elast}
   void do_pcorr_helper(TString prtcl_name);
   void mom_corr_electron(); //! [01-13-16] Obsolete after do_pcorr 

   void do_zvtxcorr();

   void set_ekin();
   void reset_ekin();

   /*int found_proton();
   int found_pip();
   int found_pim();*/
   void fill_pid_prec_hists();
   void fill_pid_pstc_hists(int h10idx_p=-1, int h10idx_pip=-1);
   int found_hadron(TString hdrn_name);
   
   int get_sector(float phi); //! phi=[-30,330]
   void GetUVW(float xyz[3], float uvw[3]);
   bool pass_p_min_ECth();
   bool pass_ECin_min();
   bool pass_ECfid();
   bool pass_zvtx();
   bool pass_nphe();
   bool pass_sf();

   double get_CC_cut_wgt();

   //! Special functions to process sim-e16-EI
   //! + EI obtains "acceptance" as ST-pass-efid/ST-fullPS
   //!   and therefore, I need to implement code to do this step
   //! + Since this step is not a part of my process, I did not want
   //!   to disturb my code and therefore have implemented the following special
   //!   functions
   bool e16_ST_pass_efid(); 

   /*
   [02-28-16] hack-pcorr for '_make_h10_skim_SS'
   + The following functions are directly copied from 'do_pcorr()',
     'do_pcorr_helper()' and 'set_ekin()', respectively. However, the
     pcorr functions do not update the h10 variables and therefore directly return
     the momentum corrected Lorentz Vector that is then used by the set_ekin function.
   + This functionality was required for '_make_h10_skim_SS' because evtsel has
     to be done corrected momenta, however, the final h10 made has to contain
     uncorrected momenta. 
   + While the already existing functions could have been modified, or atleast the same names used but
     with different parameters, at this very late stage of my analysis, I decided to take a
     very conservative approach and keep these new functions completely separate, even with regard
     to their names so that their is no confusion at all.
   */
   void do_pcorr_but_not_update_h10(TLorentzVector &lvE,TLorentzVector &lvP,TLorentzVector &lvPip); //![01-13-16] Works for {e1f,e16}*{2pi:top2',elast}
   TLorentzVector do_pcorr_but_not_update_h10_helper(TString prtcl_name);
   void set_ekin_use_passed_lv(TLorentzVector lv);
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
      if(_rctn=="2pi"){
         fChain->SetBranchAddress("mcnentr", &mcnentr, &b_mcnentr);
         fChain->SetBranchAddress("mcnpart", &mcnpart, &b_mcnpart);
         fChain->SetBranchAddress("mcst", mcst, &b_mcst);
         fChain->SetBranchAddress("mcid", mcid, &b_mcid);
         fChain->SetBranchAddress("mcpid", mcpid, &b_mcpid);
         fChain->SetBranchAddress("mctheta", mctheta, &b_mctheta);
         fChain->SetBranchAddress("mcphi", mcphi, &b_mcphi);
         fChain->SetBranchAddress("mcp", mcp, &b_mcp);
         fChain->SetBranchAddress("mcm", mcm, &b_mcm);
         fChain->SetBranchAddress("mcvx_x_el", &mcvx_x_el, &b_mcvx_x_el);
         fChain->SetBranchAddress("mcvx_y_el", &mcvx_y_el, &b_mcvx_y_el);
         fChain->SetBranchAddress("mcvx_z_el", &mcvx_z_el, &b_mcvx_z_el);
      }else if (_rctn=="elast"){
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
   }
   if (_seq=="recon"){
      //! + First take care of Binding for some of the SEB Branches in h10 TTrees
      //!   that are generally can differ based on 'cooking'.
      //! + For my case, only some Branches for e16:exp are different.
      //! + The rest: e1f:exp,e1f:sim and e16:sim are same
      //! + For detailed notes see at-h10/formats/README
      if ( _do_reconcile || (_expt=="e16" && _dtyp=="exp") ){//! SEB' = Branches that differ for e16:exp
         std::cout<<"BINDING to tmpVars"<<std::endl;
         fChain->SetBranchAddress("id", tmpVars.id, &b_id);
         fChain->SetBranchAddress("stat", tmpVars.stat, &b_stat);
         fChain->SetBranchAddress("dc", tmpVars.dc, &b_dc);
         fChain->SetBranchAddress("cc", tmpVars.cc, &b_cc);
         fChain->SetBranchAddress("sc", tmpVars.sc, &b_sc);
         fChain->SetBranchAddress("ec", tmpVars.ec, &b_ec);
         fChain->SetBranchAddress("q", tmpVars.q, &b_q);
         fChain->SetBranchAddress("dc_sect", tmpVars.dc_sect, &b_dc_sect);
         fChain->SetBranchAddress("dc_stat", tmpVars.dc_stat, &b_dc_stat);
         fChain->SetBranchAddress("ec_stat", tmpVars.ec_stat, &b_ec_stat);
         fChain->SetBranchAddress("ec_sect", tmpVars.ec_sect, &b_ec_sect);
         fChain->SetBranchAddress("sc_sect", tmpVars.sc_sect, &b_sc_sect);
         fChain->SetBranchAddress("sc_pd", tmpVars.sc_pd, &b_sc_pd);
         fChain->SetBranchAddress("sc_stat", tmpVars.sc_stat, &b_sc_stat);
         fChain->SetBranchAddress("cc_sect", tmpVars.cc_sect, &b_cc_sect);
         fChain->SetBranchAddress("nphe", tmpVars.nphe, &b_nphe);
      }else{
         fChain->SetBranchAddress("id", id, &b_id);
         fChain->SetBranchAddress("stat", stat, &b_stat);
         fChain->SetBranchAddress("dc", dc, &b_dc);
         fChain->SetBranchAddress("cc", cc, &b_cc);
         fChain->SetBranchAddress("sc", sc, &b_sc);
         fChain->SetBranchAddress("ec", ec, &b_ec);
         fChain->SetBranchAddress("q", q, &b_q);
         fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
         fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
         fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
         fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
         fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
         fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
         fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
         fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
         fChain->SetBranchAddress("nphe", nphe, &b_nphe);         
      }
      //! Bind additional SEB branches and PART branches present in E16:exp
      if (_expt=="e16" && _dtyp=="exp"){
         //! additional SEB branches in E16:exp
         fChain->SetBranchAddress("rf_time", &rf_time, &b_rf_time);
         fChain->SetBranchAddress("l2bit", &l2bit, &b_l2bit);
         fChain->SetBranchAddress("l3bit", &l3bit, &b_l3bit);
         fChain->SetBranchAddress("hlsc", &hlsc, &b_hlsc);
         fChain->SetBranchAddress("intt", &intt, &b_intt);
         fChain->SetBranchAddress("st", st, &b_st);
         fChain->SetBranchAddress("tb_st", tb_st, &b_tb_st);
         fChain->SetBranchAddress("dc_vx", dc_vx, &b_dc_vx);
         fChain->SetBranchAddress("dc_vy", dc_vy, &b_dc_vy);
         fChain->SetBranchAddress("dc_vz", dc_vz, &b_dc_vz);
         fChain->SetBranchAddress("dc_vr", dc_vr, &b_dc_vr);
         fChain->SetBranchAddress("tl1_cx", tl1_cx, &b_tl1_cx);
         fChain->SetBranchAddress("tl1_cy", tl1_cy, &b_tl1_cy);
         fChain->SetBranchAddress("tl1_cz", tl1_cz, &b_tl1_cz);
         fChain->SetBranchAddress("tl1_x", tl1_x, &b_tl1_x);
         fChain->SetBranchAddress("tl1_y", tl1_y, &b_tl1_y);
         fChain->SetBranchAddress("tl1_z", tl1_z, &b_tl1_z);
         fChain->SetBranchAddress("tl1_r", tl1_r, &b_tl1_r);
         fChain->SetBranchAddress("lec_ein", lec_ein, &b_lec_ein);
         fChain->SetBranchAddress("nschit", &nschit, &b_nschit);
         fChain->SetBranchAddress("scsect", sc_sect, &b_sc_sect);
         fChain->SetBranchAddress("schid", schid, &b_schid);
         fChain->SetBranchAddress("scpid", scpid, &b_scpid);
         fChain->SetBranchAddress("sct", sct, &b_sct);
         fChain->SetBranchAddress("sce", sce, &b_sce);
         fChain->SetBranchAddress("scx", scx, &b_scx);
         fChain->SetBranchAddress("scy", scy, &b_scy);
         fChain->SetBranchAddress("scz", scz, &b_scz);
         //! PART banks found in E16:exp
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
         fChain->SetBranchAddress("Ipart10", Ipart10, &b_Ipart10); 
         fChain->SetBranchAddress("Rpart11", Rpart11, &b_Rpart11);
         fChain->SetBranchAddress("Rpart12", Rpart12, &b_Rpart12);
         fChain->SetBranchAddress("Ipart13", Ipart13, &b_Ipart13);
      }
      //! SEB branches present in all data (expt:dtyp:seq:rctn) EXCEPT when
      //! seq=="thrown"
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
      
      fChain->SetBranchAddress("lec", lec, &b_lec);
      fChain->SetBranchAddress("p", p, &b_p);
      fChain->SetBranchAddress("m", m, &b_m);
      
      fChain->SetBranchAddress("b", b, &b_b);
      fChain->SetBranchAddress("cx", cx, &b_cx);
      fChain->SetBranchAddress("cy", cy, &b_cy);
      fChain->SetBranchAddress("cz", cz, &b_cz);
      fChain->SetBranchAddress("vx", vx, &b_vx);
      fChain->SetBranchAddress("vy", vy, &b_vy);
      fChain->SetBranchAddress("vz", vz, &b_vz);
      fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
      
      fChain->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
      
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
      fChain->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
      
      fChain->SetBranchAddress("edep", edep, &b_edep);
      fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
      fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
      fChain->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
      fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
      
      fChain->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
      fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
      
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

void h10looper_e1f::Reconcile() {
   for (int i = 0; i < gpart; i++) {
      id[i] = tmpVars.id[i];
      stat[i] = tmpVars.stat[i];
      dc[i] = tmpVars.dc[i];
      cc[i] = tmpVars.cc[i];
      sc[i] = tmpVars.sc[i];
      ec[i] = tmpVars.ec[i];
      q[i] = tmpVars.q[i];
      /*std::cout<<"*** i="<<i<<"***"<<std::endl;
      std::cout<<"tmpvars.q="<<tmpVars.q[i]<<std::endl;
      std::cout<<"q="<<q[i]<<std::endl;
      if(q[i]==-1){
         std::cout<<"q[i]==-1"<<std::endl;
      }*/
   }
   for (int i = 0; i < dc_part; i++) {
      dc_sect[i] = tmpVars.dc_sect[i];
      dc_stat[i] = tmpVars.dc_stat[i];
   }
   for (int i = 0; i < ec_part; i++) {
      ec_stat[i] = tmpVars.ec_stat[i];
      ec_sect[i] = tmpVars.ec_sect[i];
   }
   for (int i = 0; i < sc_part; i++) {
      sc_stat[i] = tmpVars.sc_stat[i];
      sc_sect[i] = tmpVars.sc_sect[i];
      sc_pd[i] = tmpVars.sc_pd[i];
   }
   for (int i = 0; i < cc_part; i++) {
      cc_sect[i] = tmpVars.cc_sect[i];
      nphe[i] = tmpVars.nphe[i];
      /*std::cout<<"*** i="<<i<<"***"<<std::endl;
      std::cout<<"tmpVars.nphe="<<tmpVars.nphe[i]<<std::endl;
      std::cout<<"nphe="<<nphe[i]<<std::endl;*/
   }
}
#endif // #ifdef h10looper_e1f_cxx
