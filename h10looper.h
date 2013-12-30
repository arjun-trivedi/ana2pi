//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 26 13:31:00 2013 by ROOT version 5.33/03
// from TTree h10/CLAS
// found on file: qskim_38751.root
//////////////////////////////////////////////////////////

#ifndef h10looper_h
#define h10looper_h

#include "data_h10.h"
#include "ep_processor.h"
/*#include "proc_eid.h"
#include "proc_efid.h"
#include "proc_pid.h"
#include "proc_skim_q.h"
#include "proc_mom_cor.h"
#include "proc_top.h"
#include "proc_skim_q2w.h"
#include "proc_fill_skim.h"
#include "proc_copy_h10.h"*/


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class h10looper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   DataH10* dH10;
   //TString proc_list;
   EpProcessor* proc_chain;

   // Declaration of leaf types
   /*UChar_t         npart;
   UInt_t          evntid;
   Char_t          evthel;
   Float_t         q_l;
   Float_t         t_l;
   Float_t         tr_time;
   Int_t           gpart;
   Short_t         id[26];   //[gpart]
   Char_t          stat[26];   //[gpart]
   UChar_t         dc[26];   //[gpart]
   UChar_t         cc[26];   //[gpart]
   UChar_t         sc[26];   //[gpart]
   UChar_t         ec[26];   //[gpart]
   Float_t         p[26];   //[gpart]
   Float_t         m[26];   //[gpart]
   Char_t          q[26];   //[gpart]
   Float_t         b[26];   //[gpart]
   Float_t         cx[26];   //[gpart]
   Float_t         cy[26];   //[gpart]
   Float_t         cz[26];   //[gpart]
   Float_t         vx[26];   //[gpart]
   Float_t         vy[26];   //[gpart]
   Float_t         vz[26];   //[gpart]
   Int_t           dc_part;
   UChar_t         dc_sect[20];   //[dc_part]
   Char_t          dc_stat[20];   //[dc_part]
   Float_t         dc_xsc[20];   //[dc_part]
   Float_t         dc_ysc[20];   //[dc_part]
   Float_t         dc_zsc[20];   //[dc_part]
   Float_t         dc_cxsc[20];   //[dc_part]
   Float_t         dc_cysc[20];   //[dc_part]
   Float_t         dc_czsc[20];   //[dc_part]
   Int_t           ec_part;
   UShort_t        ec_stat[21];   //[ec_part]
   UChar_t         ec_sect[21];   //[ec_part]
   Float_t         etot[21];   //[ec_part]
   Float_t         ec_ei[21];   //[ec_part]
   Float_t         ec_eo[21];   //[ec_part]
   Float_t         ech_x[21];   //[ec_part]
   Float_t         ech_y[21];   //[ec_part]
   Float_t         ech_z[21];   //[ec_part]
   Int_t           sc_part;
   UChar_t         sc_sect[21];   //[sc_part]
   UChar_t         sc_pd[21];   //[sc_part]
   UChar_t         sc_stat[21];   //[sc_part]
   Float_t         sc_t[21];   //[sc_part]
   Float_t         sc_r[21];   //[sc_part]
   Int_t           cc_part;
   UChar_t         cc_sect[20];   //[cc_part]
   Int_t           cc_segm[20];   //[cc_part]
   UShort_t        nphe[20];   //[cc_part]

   // List of branches
   TBranch        *b_npart;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_evthel;   //!
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
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!*/

   h10looper(TTree *tree=0,DataH10* dataH10,EpProcessor* processor_chain);
   virtual ~h10looper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   //void SetupProcs();
};

#endif

#ifdef h10looper_cxx
h10looper::h10looper(TTree *tree, DataH10* dataH10, EpProcessor* processor_chain) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("qskim_38751.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("qskim_38751.root");
      }
      f->GetObject("h10",tree);

   }
   dH10 = dataH10;

   proc_chain = processor_chain;
   //SetupProcs();

   Init(tree);
}

h10looper::~h10looper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t h10looper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t h10looper::LoadTree(Long64_t entry)
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

void h10looper::Init(TTree *tree)
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


   dH10->Bind(fChain);
   /*fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
   fChain->SetBranchAddress("evthel", &evthel, &b_evthel);
   fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
   fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("stat", stat, &b_stat);
   fChain->SetBranchAddress("dc", dc, &b_dc);
   fChain->SetBranchAddress("cc", cc, &b_cc);
   fChain->SetBranchAddress("sc", sc, &b_sc);
   fChain->SetBranchAddress("ec", ec, &b_ec);
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
   fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
   fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
   fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
   fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
   fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
   fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
   fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
   fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
   fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
   fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
   fChain->SetBranchAddress("etot", etot, &b_etot);
   fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
   fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
   fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
   fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
   fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
   fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
   fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
   fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
   fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
   fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
   fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
   fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
   fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
   fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
   fChain->SetBranchAddress("nphe", nphe, &b_nphe);*/
   Notify();
}

Bool_t h10looper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void h10looper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t h10looper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

/*void h10looper::SetupProcs(){
   //Get list of processors
   TCollection *proc_list_tokens = proc_list.Tokenize(":");
   Info("SetupProcs", "List of processors received:\n");
   proc_list_tokens->Print();
         
   //instantiate topProc; by design, topProc "builds" a cascaded chain of Procs
   top_proc = new EpProcessor();
            
   if (proc_list_tokens->IsEmpty()) {
      Info("SetupProcs", "no processors specified; building default pipeline\n");
      EpProcessor *proc;
      proc = new ProcEid(fFileOut->mkdir("eid"), dAna, h10type->GetTitle());
      //procs.push_back(proc);
      top_proc->add(proc);
      //lastProcName= new TObjString("eid");
      //Info("SlaveBegin","last processor = %s", lastProcName->GetName());
   } else {
      TIter iter(proc_list_tokens);
      while(TObjString *obj_str = (TObjString*)iter.Next()) {
         EpProcessor *proc;
         TString str = obj_str->GetString();
         if (str.EqualTo("eid"))            proc = new ProcEid(mkdir("eid"), dAna, h10type->GetTitle());
         else if (str.EqualTo("eidmon"))     proc = new ProcEid(mkdir("eid"), dAna, h10type->GetTitle(),kTRUE);
         else if (str.EqualTo("eidmononly")) proc = new ProcEid(mkdir("eid"), dAna, h10type->GetTitle(),kTRUE,kTRUE);
         else if (str.EqualTo("efid"))       proc = new ProcEFid(mkdir("fid"), dAna, h10type->GetTitle());
         else if (str.EqualTo("efidmon"))    proc = new ProcEFid(mkdir("fid"), dAna, h10type->GetTitle(),kTRUE);
         else if (str.EqualTo("efidmononly"))proc = new ProcEFid(mkdir("fid"), dAna, h10type->GetTitle(),kTRUE,kTRUE);
         else if (str.EqualTo("qskim"))       proc = new ProcSkimQ(mkdir("qskim"), dAna, h10type->GetTitle());
         else if (str.EqualTo("mom"))      proc = new ProcMomCor(mkdir("mom"), dAna, h10type->GetTitle());
         else if (str.EqualTo("pid"))      proc = new ProcPid(mkdir("pid"), dAna, h10type->GetTitle());
         else if (str.EqualTo("pidmon"))     proc = new ProcPid(mkdir("pid"), dAna, h10type->GetTitle(),kTRUE);
         else if (str.EqualTo("pidmononly")) proc = new ProcPid(mkdir("pid"), dAna, h10type->GetTitle(),kTRUE,kTRUE);
         else if (str.EqualTo("top"))      proc = new ProcTop(mkdir("top"), dAna, h10type->GetTitle());
         else if (str.EqualTo("q2wskim")) proc = new ProcSkimQ2W(mkdir("q2wskim"), dAna, h10type->GetTitle());
         else if (str.EqualTo("fillskim"))   proc = new ProcFillSkim(mkdir("skim"), dAna, h10type->GetTitle());
         else if (str.EqualTo("copyh10")) proc = new ProcCopyH10(fFileOut, dAna, h10type->GetTitle());
         else {
            Info("Init","%s unrecognized processor\n",str.Data());
            continue;
         }
         //procs.push_back(proc);
         top_proc->add(proc);
         lastProcName = obj_str; //works because after last iteration, lastProcName is not updated
      }
      //Info("SlaveBegin","last processor = %s", lastProcName->GetName());
   }
}*/
#endif // #ifdef h10looper_cxx
