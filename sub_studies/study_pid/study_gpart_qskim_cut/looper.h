//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul 26 13:22:42 2015 by ROOT version 5.33/03
// from TTree t/TTree containing data for PID
// found on file: dpid_test.root
//////////////////////////////////////////////////////////

#ifndef looper_h
#define looper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class looper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         l_e;
   Float_t         t_e;
   Float_t         t_off;
   Int_t           gpart;
   Int_t           npos;
   Int_t           nneg;
   Int_t           nzro;
   Int_t           ntrk;
   Int_t           np;
   Int_t           npip;
   Int_t           npim;
   Int_t           dc[19];   //[ntrk]
   Int_t           sc[19];   //[ntrk]
   Int_t           q[19];   //[ntrk]
   Float_t         p[19];   //[ntrk]
   Float_t         l[19];   //[ntrk]
   Float_t         t[19];   //[ntrk]
   Float_t         b[19];   //[ntrk]
   Int_t           id[19];   //[ntrk]
   Int_t           sector[19];   //[ntrk]
   Int_t           h10_idx[19];   //[ntrk]
   //Float_t         b_p[19];   //[ntrk]
   Float_t         b_pip[19];   //[ntrk]
   Float_t         b_pim[19];   //[ntrk]
   Float_t         dt_p[19];   //[ntrk]
   Float_t         dt_pip[19];   //[ntrk]
   Float_t         dt_pim[19];   //[ntrk]

   // List of branches
   TBranch        *b_l_e;   //!
   TBranch        *b_t_e;   //!
   TBranch        *b_t_off;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_npos;   //!
   TBranch        *b_nneg;   //!
   TBranch        *b_nzro;   //!
   TBranch        *b_ntrk;   //!
   TBranch        *b_np;   //!
   TBranch        *b_npip;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_q;   //!
   TBranch        *b_p;   //!
   TBranch        *b_l;   //!
   TBranch        *b_t;   //!
   TBranch        *b_b;   //!
   TBranch        *b_id;   //!
   TBranch        *b_sector;   //!
   TBranch        *b_h10_idx;   //!
   //TBranch        *b_b_p;   //!
   TBranch        *b_b_pip;   //!
   TBranch        *b_b_pim;   //!
   TBranch        *b_dt_p;   //!
   TBranch        *b_dt_pip;   //!
   TBranch        *b_dt_pim;   //!

   looper(TTree *tree=0);
   virtual ~looper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef looper_cxx
looper::looper(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dpid_test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("dpid_test.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("dpid_test.root:/pid/cut");
      dir->GetObject("t",tree);

   }
   Init(tree);
}

looper::~looper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t looper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t looper::LoadTree(Long64_t entry)
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

void looper::Init(TTree *tree)
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

   fChain->SetBranchAddress("l_e", &l_e, &b_l_e);
   fChain->SetBranchAddress("t_e", &t_e, &b_t_e);
   fChain->SetBranchAddress("t_off", &t_off, &b_t_off);
   fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
   fChain->SetBranchAddress("npos", &npos, &b_npos);
   fChain->SetBranchAddress("nneg", &nneg, &b_nneg);
   fChain->SetBranchAddress("nzro", &nzro, &b_nzro);
   fChain->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("dc", dc, &b_dc);
   fChain->SetBranchAddress("sc", sc, &b_sc);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("l", l, &b_l);
   fChain->SetBranchAddress("t", t, &b_t);
   fChain->SetBranchAddress("b", b, &b_b);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("sector", sector, &b_sector);
   fChain->SetBranchAddress("h10_idx", h10_idx, &b_h10_idx);
   //fChain->SetBranchAddress("b_p", b_p, &b_b_p);
   fChain->SetBranchAddress("b_pip", b_pip, &b_b_pip);
   fChain->SetBranchAddress("b_pim", b_pim, &b_b_pim);
   fChain->SetBranchAddress("dt_p", dt_p, &b_dt_p);
   fChain->SetBranchAddress("dt_pip", dt_pip, &b_dt_pip);
   fChain->SetBranchAddress("dt_pim", dt_pim, &b_dt_pim);
   Notify();
}

Bool_t looper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void looper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t looper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef looper_cxx
