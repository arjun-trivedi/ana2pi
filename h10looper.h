//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 26 13:31:00 2013 by ROOT version 5.33/03
// from TTree h10/CLAS
// found on file: qskim_38751.root
//////////////////////////////////////////////////////////

#ifndef H10Looper_h
#define H10Looper_h

#include "data_h10.h"
#include "ep_processor.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class H10Looper {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   DataH10* dH10;
   DataAna* dAna;
   EpProcessor* proc_chain;

   Bool_t _use_q2w_elist;
   TFile* _f_q2w_el;
   TString _q2w;
   TEntryList* _el;

   H10Looper(TChain *tree=0,DataH10* dataH10=0,DataAna* dataAna=0,EpProcessor* processor_chain=0,
             Bool_t use_q2w_elist=kFALSE,TFile* f_q2w_el=NULL,TString q2w="");
   virtual ~H10Looper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop(Long64_t nentries);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef H10Looper_cxx
H10Looper::H10Looper(TChain *tree, DataH10* dataH10, DataAna* dataAna, EpProcessor* processor_chain,
                     Bool_t use_q2w_elist/*=kFALSE*/,TFile* f_q2w_el/*=NULL*/,TString q2w/*=""*/) : fChain(0) 
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
   dAna = dataAna;

   proc_chain = processor_chain;
   //SetupProcs();

   _use_q2w_elist=use_q2w_elist;
   _f_q2w_el=f_q2w_el;
   _q2w=q2w;
   _el=NULL;

   Init(tree);
}

H10Looper::~H10Looper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t H10Looper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t H10Looper::LoadTree(Long64_t entry)
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

void H10Looper::Init(TChain *tree)
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
   Notify();
}

Bool_t H10Looper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void H10Looper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t H10Looper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef H10Looper_cxx
