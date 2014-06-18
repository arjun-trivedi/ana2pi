#define H10Looper_cxx
#include "h10looper.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

void H10Looper::Loop(Long64_t nentries)
{
   if (fChain == 0) return;

   //! Track memory usage
   MemInfo_t meminfo;
   gSystem->GetMemInfo(&meminfo);
   int mem_start=meminfo.fMemUsed+meminfo.fSwapUsed;


   //Long64_t nentries = fChain->GetEntriesFast();
   Int_t nentries_chain = fChain->GetEntries();
   Int_t nentries_to_proc=0;
   nentries>nentries_chain?nentries_to_proc=nentries_chain:nentries_to_proc=nentries;
   Info("H10Looper::Loop", "Number of entries in Chain =  %d\n",nentries_chain);
   Info("H10Looper::Loop", "Number of entries to processess =  %d\n",nentries_to_proc);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries_to_proc;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      // if (Cut(ientry) < 0) continue;

      // 1.Clear event data
      dH10->Clear();
      dAna->Clear();

      // 2.Load Branches into dH10
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // 2.1 If needed, Reconcile dH10
      if (dH10->rctn=="2pi_userana" || 
          dH10->rctn=="elast_userana" ||
          (dH10->expt=="e16" && dH10->dtyp=="exp")) {
         dH10->Reconcile();
      }
      
      if (jentry%100000==0) {
         Info("H10Looper::Loop", "Processing entry# %d\n",jentry);
         Info("H10Looper::Loop", "%5.2f%% entries processed\n",(float)jentry/nentries_to_proc*100);  
      }
      //Reset proc_chain
      //3. Call proc_chain
      proc_chain->handle();
   }
   gSystem->GetMemInfo(&meminfo);
   int mem_end=meminfo.fMemUsed+meminfo.fSwapUsed;
   Info("H10Looper::Loop", "Total memory used = %dMB\n",mem_end-mem_start);
}
