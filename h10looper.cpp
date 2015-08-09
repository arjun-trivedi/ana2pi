#define H10Looper_cxx
#include "h10looper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TEntryList.h>

void H10Looper::Loop(Long64_t nentries)
{
   if (fChain == 0) return;

   //! Track memory usage
   MemInfo_t meminfo;
   gSystem->GetMemInfo(&meminfo);
   int mem_start=meminfo.fMemUsed+meminfo.fSwapUsed;


   Long64_t nentries_chain;// = fChain->GetEntries();
   if (_use_q2w_elist){
      _el=(TEntryList*)_f_q2w_el->Get(TString::Format("q2welist/%s",_q2w.Data()));
      nentries_chain = _el->GetN();
      fChain->SetEntryList(_el);
   } else {
      nentries_chain=fChain->GetEntries();
   }
   Long64_t nentries_to_proc=0;
   //nentries>nentries_chain?nentries_to_proc=nentries_chain:nentries_to_proc=nentries;//! old "logic", pre [08-09-15]
   nentries==0?nentries_to_proc=nentries_chain:nentries_to_proc=nentries;
   Info("H10Looper::Loop", "Number of entries in specified by user =  %llu",nentries);
   Info("H10Looper::Loop", "Number of entries in TChain =  %llu",nentries_chain);
   Info("H10Looper::Loop", "Number of entries to processess =  %llu",nentries_to_proc);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries_to_proc;jentry++) {
      // 1.Clear event data
      dH10->Clear();
      dAna->Clear();

      // 2.Load Branches into dH10 & set dH10::_ientry_h10chain
      Long64_t chainEntry=-1;
      if (_use_q2w_elist){
         Int_t treenum=0;
         Long64_t treeEntry = _el->GetEntryAndTree(jentry,treenum);
         chainEntry = treeEntry+fChain->GetTreeOffset()[treenum];
         // Info("H10Looper::Loop","listEntry=%lld, treeEntry=%lld, chainEntry=%lld, treenum=%d\n"
         //    ,jentry,treeEntry,chainEntry,treenum);
      }else{
         chainEntry=jentry;
      }
      Long64_t ientry = fChain->LoadTree(chainEntry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(chainEntry);   nbytes += nb;
      dH10->set_ientry_h10chain(chainEntry);
      
      // 2.1 If needed, Reconcile dH10
      if (dH10->rctn=="2pi_userana" || 
          dH10->rctn=="elast_userana" ||
          (dH10->expt=="e16" && dH10->dtyp=="exp")) {
         dH10->Reconcile();
      }
      
      if (jentry%100000==0) {
         Info("H10Looper::Loop", "Processing entry# %d(file=%s)",jentry,fChain->GetFile()->GetName());
         Info("H10Looper::Loop", "%5.2f%% entries processed",(float)jentry/nentries_to_proc*100);  
         gSystem->GetMemInfo(&meminfo);
         int mem=meminfo.fMemUsed+meminfo.fSwapUsed;
         Info("H10Looper::Loop", "Total memory used = %dMB",mem-mem_start);
      }
      //Reset proc_chain
      //3. Call proc_chain
      proc_chain->handle();
   }
   gSystem->GetMemInfo(&meminfo);
   int mem_end=meminfo.fMemUsed+meminfo.fSwapUsed;
   Info("H10Looper::Loop", "Total memory used = %dMB\n",mem_end-mem_start);
}
