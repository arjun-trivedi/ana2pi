//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Mar 25 06:48:22 2012 by ROOT version 5.33/01
// from TChain h10/
//////////////////////////////////////////////////////////

#ifndef SelSim_h
#define SelSim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include "data_omega.h"
#include "data_h10.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SelSim : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	DataOmega data;
	//tmp variables for data type compatibility
	Int_t id[DataH10::_MAX_PARTS]; //[gpart]
	Int_t stat[DataH10::_MAX_PARTS]; //[gpart]
	Int_t dc[DataH10::_MAX_PARTS]; //[gpart]
	Int_t cc[DataH10::_MAX_PARTS]; //[gpart]
	Int_t sc[DataH10::_MAX_PARTS]; //[gpart]
	Int_t ec[DataH10::_MAX_PARTS]; //[gpart]
	Int_t q[DataH10::_MAX_PARTS]; //[gpart]
	Int_t dc_sect[DataH10::_MAX_PARTS]; //[dc_part]
	Int_t dc_stat[DataH10::_MAX_PARTS]; //[dc_part]
	Int_t ec_stat[DataH10::_MAX_PARTS]; //[ec_part]
	Int_t ec_sect[DataH10::_MAX_PARTS]; //[ec_part]
	Int_t sc_sect[DataH10::_MAX_PARTS]; //[sc_part]
	Int_t sc_pd[DataH10::_MAX_PARTS]; //[sc_part]
	Int_t sc_stat[DataH10::_MAX_PARTS]; //[sc_part]
	Int_t cc_sect[DataH10::_MAX_PARTS]; //[cc_part]
	Int_t nphe[DataH10::_MAX_PARTS]; //[cc_part]

   // List of branches
   TBranch        *b_npart;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_t_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_id;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_ec;   //!
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
   TBranch        *b_dc_stat;   //!
   TBranch        *b_dc_xsc;   //!
   TBranch        *b_dc_ysc;   //!
   TBranch        *b_dc_zsc;   //!
   TBranch        *b_dc_cxsc;   //!
   TBranch        *b_dc_cysc;   //!
   TBranch        *b_dc_czsc;   //!
   TBranch        *b_ec_part;   //!
   TBranch        *b_ec_stat;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_vidmvrt;   //!
   TBranch        *b_ntrmvrt;   //!
   TBranch        *b_xmvrt;   //!
   TBranch        *b_ymvrt;   //!
   TBranch        *b_zmvrt;   //!
   TBranch        *b_ch2mvrt;   //!
   TBranch        *b_cxxmvrt;   //!
   TBranch        *b_cxymvrt;   //!
   TBranch        *b_cxzmvrt;   //!
   TBranch        *b_cyymvrt;   //!
   TBranch        *b_cyzmvrt;   //!
   TBranch        *b_stamvrt;   //!
   TBranch        *b_mcnentr;   //!
   TBranch        *b_mcnpart;   //!
   TBranch        *b_mcst;   //!
   TBranch        *b_mcid;   //!
   TBranch        *b_mcpid;   //!
   TBranch        *b_mctheta;   //!
   TBranch        *b_mcphi;   //!
   TBranch        *b_mcp;   //!
   TBranch        *b_mcm;   //!
   TBranch        *b_mcvx;   //!
   TBranch        *b_mcvy;   //!
   TBranch        *b_mcvz;   //!
   TBranch        *b_mctof;   //!
   TBranch        *b_nprt;   //!
   
   SelSim(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~SelSim() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) {
		Int_t retval = 0;
		if (fChain) {
			retval = fChain->GetTree()->GetEntry(entry, getall);
			Reconcile();
		  //return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
		 }
	 }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
	 
	 void Reconcile();

   ClassDef(SelSim,0);
};

#endif

#ifdef SelSim_cxx
void SelSim::Init(TTree *tree)
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
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("npart", &data.h10.npart, &b_npart);
   fChain->SetBranchAddress("q_l", &data.h10.q_l, &b_q_l);
   fChain->SetBranchAddress("t_l", &data.h10.t_l, &b_t_l);
   fChain->SetBranchAddress("tr_time", &data.h10.tr_time, &b_tr_time);
   fChain->SetBranchAddress("gpart", &data.h10.gpart, &b_gpart);
   fChain->SetBranchAddress("id", data.h10.id, &b_id);
   fChain->SetBranchAddress("stat", data.h10.stat, &b_stat);
   fChain->SetBranchAddress("dc", data.h10.dc, &b_dc);
   fChain->SetBranchAddress("cc", data.h10.cc, &b_cc);
   fChain->SetBranchAddress("sc", data.h10.sc, &b_sc);
   fChain->SetBranchAddress("ec", data.h10.ec, &b_ec);
   fChain->SetBranchAddress("p", data.h10.p, &b_p);
   fChain->SetBranchAddress("m", data.h10.m, &b_m);
   fChain->SetBranchAddress("q", data.h10.q, &b_q);
   fChain->SetBranchAddress("b", data.h10.b, &b_b);
   fChain->SetBranchAddress("cx", data.h10.cx, &b_cx);
   fChain->SetBranchAddress("cy", data.h10.cy, &b_cy);
   fChain->SetBranchAddress("cz", data.h10.cz, &b_cz);
   fChain->SetBranchAddress("vx", data.h10.vx, &b_vx);
   fChain->SetBranchAddress("vy", data.h10.vy, &b_vy);
   fChain->SetBranchAddress("vz", data.h10.vz, &b_vz);
   fChain->SetBranchAddress("dc_part", &data.h10.dc_part, &b_dc_part);
   fChain->SetBranchAddress("dc_sect", data.h10.dc_sect, &b_dc_sect);
   fChain->SetBranchAddress("dc_stat", data.h10.dc_stat, &b_dc_stat);
   fChain->SetBranchAddress("dc_xsc", data.h10.dc_xsc, &b_dc_xsc);
   fChain->SetBranchAddress("dc_ysc", data.h10.dc_ysc, &b_dc_ysc);
   fChain->SetBranchAddress("dc_zsc", data.h10.dc_zsc, &b_dc_zsc);
   fChain->SetBranchAddress("dc_cxsc", data.h10.dc_cxsc, &b_dc_cxsc);
   fChain->SetBranchAddress("dc_cysc", data.h10.dc_cysc, &b_dc_cysc);
   fChain->SetBranchAddress("dc_czsc", data.h10.dc_czsc, &b_dc_czsc);
   fChain->SetBranchAddress("dc_vx", data.h10.dc_vx, &b_dc_vx);
   fChain->SetBranchAddress("dc_vy", data.h10.dc_vy, &b_dc_vy);
   fChain->SetBranchAddress("dc_vz", data.h10.dc_vz, &b_dc_vz);
   fChain->SetBranchAddress("ec_part", &data.h10.ec_part, &b_ec_part);
   fChain->SetBranchAddress("ec_stat", data.h10.ec_stat, &b_ec_stat);
   fChain->SetBranchAddress("ec_sect", data.h10.ec_sect, &b_ec_sect);
   fChain->SetBranchAddress("etot", data.h10.etot, &b_etot);
   fChain->SetBranchAddress("ec_ei", data.h10.ec_ei, &b_ec_ei);
   fChain->SetBranchAddress("ec_eo", data.h10.ec_eo, &b_ec_eo);
   fChain->SetBranchAddress("sc_part", &data.h10.sc_part, &b_sc_part);
   fChain->SetBranchAddress("sc_sect", data.h10.sc_sect, &b_sc_sect);
   fChain->SetBranchAddress("sc_pd", data.h10.sc_pd, &b_sc_pd);
   fChain->SetBranchAddress("sc_stat", data.h10.sc_stat, &b_sc_stat);
   fChain->SetBranchAddress("sc_t", data.h10.sc_t, &b_sc_t);
   fChain->SetBranchAddress("sc_r", data.h10.sc_r, &b_sc_r);
   fChain->SetBranchAddress("cc_part", &data.h10.cc_part, &b_cc_part);
   fChain->SetBranchAddress("cc_sect", data.h10.cc_sect, &b_cc_sect);
   fChain->SetBranchAddress("cc_segm", data.h10.cc_segm, &b_cc_segm);
   fChain->SetBranchAddress("nphe", data.h10.nphe, &b_nphe);
   fChain->SetBranchAddress("vidmvrt", &data.h10.vidmvrt, &b_vidmvrt);
   fChain->SetBranchAddress("ntrmvrt", &data.h10.ntrmvrt, &b_ntrmvrt);
   fChain->SetBranchAddress("xmvrt", &data.h10.xmvrt, &b_xmvrt);
   fChain->SetBranchAddress("ymvrt", &data.h10.ymvrt, &b_ymvrt);
   fChain->SetBranchAddress("zmvrt", &data.h10.zmvrt, &b_zmvrt);
   fChain->SetBranchAddress("ch2mvrt", &data.h10.ch2mvrt, &b_ch2mvrt);
   fChain->SetBranchAddress("cxxmvrt", &data.h10.cxxmvrt, &b_cxxmvrt);
   fChain->SetBranchAddress("cxymvrt", &data.h10.cxymvrt, &b_cxymvrt);
   fChain->SetBranchAddress("cxzmvrt", &data.h10.cxzmvrt, &b_cxzmvrt);
   fChain->SetBranchAddress("cyymvrt", &data.h10.cyymvrt, &b_cyymvrt);
   fChain->SetBranchAddress("cyzmvrt", &data.h10.cyzmvrt, &b_cyzmvrt);
   fChain->SetBranchAddress("stamvrt", &data.h10.stamvrt, &b_stamvrt);
   fChain->SetBranchAddress("mcnentr", &data.h10.mcnentr, &b_mcnentr);
   fChain->SetBranchAddress("mcnpart", &data.h10.mcnpart, &b_mcnpart);
   fChain->SetBranchAddress("mcst", data.h10.mcst, &b_mcst);
   fChain->SetBranchAddress("mcid", data.h10.mcid, &b_mcid);
   fChain->SetBranchAddress("mcpid", data.h10.mcpid, &b_mcpid);
   fChain->SetBranchAddress("mctheta", data.h10.mctheta, &b_mctheta);
   fChain->SetBranchAddress("mcphi", data.h10.mcphi, &b_mcphi);
   fChain->SetBranchAddress("mcp", data.h10.mcp, &b_mcp);
   fChain->SetBranchAddress("mcm", data.h10.mcm, &b_mcm);
   fChain->SetBranchAddress("mcvx", data.h10.mcvx, &b_mcvx);
   fChain->SetBranchAddress("mcvy", data.h10.mcvy, &b_mcvy);
   fChain->SetBranchAddress("mcvz", data.h10.mcvz, &b_mcvz);
   fChain->SetBranchAddress("mctof", data.h10.mctof, &b_mctof);
}

Bool_t SelSim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef SelSim_cxx
