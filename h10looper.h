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

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class h10looper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   DataH10* dH10;
   EpProcessor* proc_chain;

   h10looper(TTree *tree=0,DataH10* dataH10,EpProcessor* processor_chain);
   virtual ~h10looper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
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

#endif // #ifdef h10looper_cxx
