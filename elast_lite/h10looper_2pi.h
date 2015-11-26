//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep  6 09:24:36 2015 by ROOT version 5.33/03
// from TTree h10/CLAS
// found on file: /data/trivedia/e1f/elastic/sim/sim1.1/recon/1.root
//////////////////////////////////////////////////////////

#ifndef h10looper_2pi_h
#define h10looper_2pi_h

#include "h10looper_e1f.h"

#include <THnSparse.h>

class h10looper_2pi: public h10looper_e1f {
public :

   h10looper_2pi(TString h10type,TChain* h10chain,
                 TString cutsncors,
                 TString fout_name, Long64_t nentries,
                 TString adtnl_opts);
   virtual ~h10looper_2pi();
   virtual void     Loop();

   void reset_hkin();
   void set_hkin(int h10idx_p=-1, int h10idx_pip=-1); //! as per top2' logic

   bool pass_pid(int& h10idx_p, int& h10idx_pip);
   
   bool pass_pfid();
   bool proton_infid();
   bool pip_infid();

   void setup_d2pi();
   void fill_h8();

   //! + Following 2 functions taken directly from proc_d2pi.h
   //! + Modifications are purely for the sake of code readability
   float getPhi(TLorentzVector lv); //[0,360]; uses invTan
   float invTan(float y, float x); //angle in radians: [0, 2pi]
   float getAlpha(TVector3 uv_Gf,TVector3 uv_Gp,TVector3 uv_Bf,TVector3 uv_Bp);

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
   //! Varsets kinematics (phi=[0,360])
   float _M_ppip,_M_ppim,_M_pippim;
   float _p_cms_p,_theta_cms_p,_phi_cms_p;
   float _p_cms_pip,_theta_cms_pip,_phi_cms_pip;
   float _p_cms_pim,_theta_cms_pim,_phi_cms_pim;
   float _alpha_1,_alpha_2,_alpha_3;


   //! output objects

   //! PFID
   static const int NUM_PFID_STATS=4;
   enum {PFID_NULL, PFID_TOT, PFID_P_IN, PFID_PIP_IN, PFID_P_AND_PIP_IN};
   TH1D* _hpfid;
   //! _hpfid_p/pip[2], before and after cut
   TH2F** _hpfid_p;
   TH2F** _hpfid_pip;

   //! d2pi
   //! Q2-W
   TH2F* _hq2w_prec;
   TH2F* _hq2w_pstc;
   //! hmm2ppip,hmmpip
   TH1F* _hmm2_prec_fW;
   TH1F* _hmm_prec_fW;
   TH1F* _hmm2_pstc_fW;
   TH1F* _hmm_pstc_fW;
   //! hmm2ppip[NBINS_WCRS],hmmpip[NBINS_WCRS]
   TH1F** _hmm2_prec;
   TH1F** _hmm_prec;
   TH1F** _hmm2_pstc;
   TH1F** _hmm_pstc;
   //!h7[NBINS_WCRS][NVST]
   static const int NVST=3;
   THnSparse*** _h8;
   //! Q2,W cut values
   float _Q2_MIN,_Q2_MAX,_W_MIN,_W_MAX;
   //! MMcut value
   float _mm2ppip_l,_mm2ppip_h;
   
   
};

#endif


